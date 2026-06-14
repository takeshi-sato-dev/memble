#!/usr/bin/env python3
"""
inject_posres.py

Inject staged, #ifdef-guarded [ position_restraints ] blocks into one or more
moleculetypes inside a GROMACS .itp, so a multi-stage equilibration can release
restraints gradually (the CHARMM-GUI Martini scheme: protein backbone and lipid
head beads restrained, force constants ramped down over the equilibration).

Position restraints are per-moleculetype with atom indices LOCAL to that
moleculetype, so this must edit each itp in place (operate on a LOCAL COPY of
any force-field itp, never the distributed original).

For each named moleculetype it finds the local indices of the requested bead
names in the [ atoms ] block and appends, at the end of that moleculetype:

    #ifdef POSRES_STEP1
    [ position_restraints ]
    <idx> 1 <fc> <fc> <fc>
    ...
    #endif
    #ifdef POSRES_STEP2
    ...

Each equilibration stage mdp then sets e.g.  define = -DPOSRES_STEP3  and only
that stage's block (with that stage's force constant) is active.

Usage (one moleculetype, beads BB, five stages):
  inject_posres.py --itp molecule_0.itp --mol molecule_0 --beads BB \\
      --stage POSRES_STEP1:1000 --stage POSRES_STEP2:500 \\
      --stage POSRES_STEP3:200 --stage POSRES_STEP4:100 \\
      --stage POSRES_STEP5:50

Repeat the call per itp / per moleculetype with the appropriate beads and the
SAME define names but each itp's own force constants.
"""

import argparse
import sys


def find_molecule_atom_indices(lines, mol, bead_names, resid_min=None, resid_max=None):
    """Return (insert_line_index, [local_atom_indices]) for moleculetype `mol`.

    insert_line_index is where to splice the restraint blocks (end of that
    moleculetype: just before the next [ moleculetype ] header, or EOF).
    Local atom indices are 1-based, as written in the itp [ atoms ] block.
    """
    sect = None
    in_target = False
    in_atoms = False
    local_idx = []
    pending_name = None
    insert_at = len(lines)
    seen_target = False

    i = 0
    while i < len(lines):
        raw = lines[i]
        s = raw.strip()
        # section header?
        if s.startswith("[") and s.endswith("]"):
            sect = s.strip("[] ").lower()
            if sect == "moleculetype":
                if in_target:
                    # next moleculetype starts -> end of our target block
                    insert_at = i
                    return insert_at, local_idx, seen_target
                # peek the moleculetype name on the next non-comment line
                j = i + 1
                while j < len(lines):
                    t = lines[j].strip()
                    if t and not t.startswith(";"):
                        pending_name = t.split()[0]
                        break
                    j += 1
                in_target = (pending_name == mol)
                if in_target:
                    seen_target = True
                in_atoms = False
            else:
                in_atoms = (sect == "atoms") and in_target
            i += 1
            continue

        if in_target and in_atoms and s and not s.startswith(";"):
            parts = s.split()
            # [ atoms ]: nr type resnr resid atomname cgnr charge [mass]
            if len(parts) >= 5:
                nr = parts[0]
                atomname = parts[4]
                in_range = True
                if resid_min is not None or resid_max is not None:
                    try:
                        rnr = int(parts[2])
                        if resid_min is not None and rnr < resid_min:
                            in_range = False
                        if resid_max is not None and rnr > resid_max:
                            in_range = False
                    except ValueError:
                        in_range = False
                if atomname in bead_names and in_range:
                    try:
                        local_idx.append(int(nr))
                    except ValueError:
                        pass
        i += 1

    return insert_at, local_idx, seen_target


def auto_head_bead(lines, mol):
    """Pick a restraint bead by priority among the molecule's atoms."""
    # re-scan to honor moleculetype boundaries correctly
    in_target = in_atoms = False
    names = []
    pending = False
    for raw in lines:
        t = raw.strip()
        if t.startswith("[") and t.endswith("]"):
            sect = t.strip("[] ").lower()
            if sect == "moleculetype":
                pending = True; in_target = False; in_atoms = False; continue
            in_atoms = (sect == "atoms") and in_target; continue
        if pending and t and not t.startswith(";"):
            in_target = (t.split()[0] == mol); pending = False; continue
        if in_target and in_atoms and t and not t.startswith(";"):
            p = t.split()
            if len(p) >= 5:
                names.append(p[4])
    for pref in ("PO4", "ROH", "GL1", "GL0", "NC3", "GM1", "AM1"):
        if pref in names:
            return pref
    return names[0] if names else None


def build_blocks(stages, indices):
    out = []
    for define, fc in stages:
        if fc <= 0:
            continue
        out.append("#ifdef %s\n" % define)
        out.append("[ position_restraints ]\n")
        out.append("; ai  funct  fx  fy  fz\n")
        for idx in indices:
            out.append("%6d 1 %d %d %d\n" % (idx, fc, fc, fc))
        out.append("#endif\n")
    out.append("\n")
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--itp", required=True, help="LOCAL COPY of the itp to edit in place")
    p.add_argument("--mol", required=True, help="moleculetype name to target")
    p.add_argument("--resid-min", type=int, default=None,
                   help="restrain only beads with itp resid >= this (TM core)")
    p.add_argument("--resid-max", type=int, default=None,
                   help="restrain only beads with itp resid <= this (TM core)")
    p.add_argument("--beads", required=True,
                   help="comma list of bead (atom) names to restrain, e.g. BB or PO4 or ROH")
    p.add_argument("--stage", action="append", required=True,
                   help="DEFINE:FC, repeatable, e.g. POSRES_STEP1:1000")
    args = p.parse_args()

    if args.beads.strip().lower() == "auto":
        with open(args.itp) as _fh:
            _lines = _fh.readlines()
        pick = auto_head_bead(_lines, args.mol)
        if not pick:
            sys.exit("ERROR: could not auto-pick a head bead for '%s'" % args.mol)
        bead_names = {pick}
        print("auto head bead for %s: %s" % (args.mol, pick))
    else:
        bead_names = set(b.strip() for b in args.beads.split(",") if b.strip())
    stages = []
    for st in args.stage:
        name, fc = st.split(":")
        stages.append((name.strip(), int(fc)))

    with open(args.itp) as fh:
        lines = fh.readlines()

    insert_at, indices, seen = find_molecule_atom_indices(
        lines, args.mol, bead_names, args.resid_min, args.resid_max)
    if not seen:
        sys.exit("ERROR: moleculetype '%s' not found in %s" % (args.mol, args.itp))
    if not indices:
        sys.exit("ERROR: no beads %s found in moleculetype '%s'"
                 % (sorted(bead_names), args.mol))

    blocks = build_blocks(stages, indices)
    new_lines = lines[:insert_at] + blocks + lines[insert_at:]
    with open(args.itp, "w") as fh:
        fh.writelines(new_lines)

    print("injected %d staged restraint block(s) on %d bead(s) of '%s' in %s"
          % (sum(1 for _, fc in stages if fc > 0), len(indices), args.mol, args.itp))


if __name__ == "__main__":
    main()
