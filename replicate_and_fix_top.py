#!/usr/bin/env python3
"""
replicate_and_fix_top.py

Two jobs, kept deliberately separate so each can be checked on its own:

  mode = replicate
      Read a single-chain coarse-grained peptide (PDB from martinize2 -x),
      make N translated copies on an XY grid, write one multi-chain PDB.
      The membrane-normal (z) is left untouched so every copy inserts into
      the same bilayer plane. Copies are placed far apart so insane builds
      one continuous bilayer around all of them.

  mode = fixtop
      Take the topol.top that insane wrote (which knows the lipid / water /
      ion counts but lists the protein as a single generic entry) and rebuild
      a clean topol.top with explicit Martini 3 includes and the correct
      protein copy count. Lipid / solvent / ion lines are copied verbatim
      from the insane top, so the counts stay exactly what insane placed.

This script does not call GROMACS, insane, or martinize2; it only moves
coordinates and rewrites a text topology. Run the two modes in order.
"""

import argparse
import sys


# ----------------------------------------------------------------------
# PDB parsing kept minimal and explicit (no external deps on the cluster).
# ----------------------------------------------------------------------
def read_pdb_atoms(path):
    """Return (header_lines_ignored, list_of_atom_record_strings)."""
    atoms = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atoms.append(line.rstrip("\n"))
    if not atoms:
        sys.exit("ERROR: no ATOM/HETATM records found in %s" % path)
    return atoms


def get_xyz(atom_line):
    # PDB fixed columns: x 31-38, y 39-46, z 47-54 (1-indexed)
    x = float(atom_line[30:38])
    y = float(atom_line[38:46])
    z = float(atom_line[46:54])
    return x, y, z


def set_xyz(atom_line, x, y, z):
    return "%s%8.3f%8.3f%8.3f%s" % (
        atom_line[:30], x, y, z, atom_line[54:]
    )


def set_chain(atom_line, chain_char):
    # chain ID is column 22 (1-indexed) -> index 21
    line = atom_line
    if len(line) < 22:
        line = line + " " * (22 - len(line))
    return line[:21] + chain_char + line[22:]


def centroid(atoms):
    n = len(atoms)
    sx = sy = sz = 0.0
    for a in atoms:
        x, y, z = get_xyz(a)
        sx += x
        sy += y
        sz += z
    return sx / n, sy / n, sz / n


def replicate(in_pdb, out_pdb, offsets_nm):
    """offsets_nm: list of (dx, dy) in nm for each copy. z is preserved."""
    atoms = read_pdb_atoms(in_pdb)
    cx, cy, _cz = centroid(atoms)  # center only in XY; keep z as-is
    chains = "ABCDEFGH"
    if len(offsets_nm) > len(chains):
        sys.exit("ERROR: more copies than chain IDs available")

    out = []
    serial = 0
    for i, (dx_nm, dy_nm) in enumerate(offsets_nm):
        dx = dx_nm * 10.0  # nm -> angstrom (PDB units)
        dy = dy_nm * 10.0
        for a in atoms:
            x, y, z = get_xyz(a)
            # recenter this copy on the XY origin, then shift to its grid slot
            nx = (x - cx) + dx
            ny = (y - cy) + dy
            serial += 1
            line = set_xyz(a, nx, ny, z)
            line = set_chain(line, chains[i])
            out.append(line)

    with open(out_pdb, "w") as fh:
        fh.write("TITLE     EGFR TM-JM x%d for insane (Martini 3, -GM3)\n"
                 % len(offsets_nm))
        for line in out:
            fh.write(line + "\n")
        fh.write("END\n")
    print("wrote %d copies (%d atoms) to %s"
          % (len(offsets_nm), serial, out_pdb))


# ----------------------------------------------------------------------
# Topology rebuild.
# ----------------------------------------------------------------------
def parse_insane_molecules(insane_top):
    """Return ordered list of (molname, count) from the [ molecules ] block.

    insane lists the protein it was given plus every lipid / solvent / ion
    it placed. We keep everything that is NOT a protein entry verbatim.
    """
    mols = []
    in_block = False
    with open(insane_top) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(";"):
                continue
            low = s.lower().replace(" ", "")
            if low.startswith("[molecules]"):
                in_block = True
                continue
            if in_block:
                if s.startswith("["):
                    break
                parts = s.split()
                if len(parts) >= 2:
                    mols.append((parts[0], parts[1]))
    if not mols:
        sys.exit("ERROR: no [ molecules ] block parsed from %s" % insane_top)
    return mols


def fixtop(insane_top, out_top, protein_itp, protein_name, n_protein,
           include_lines, protein_token):
    """Rebuild topol.top with explicit M3 includes and N protein copies.

    protein_token: the molecule name insane used for the protein in its
    [ molecules ] block (so we can drop / replace it). If insane was run
    WITHOUT a protein (-f omitted) this can be left empty and all parsed
    molecules are treated as non-protein.
    """
    mols = parse_insane_molecules(insane_top)

    non_protein = []
    for name, count in mols:
        if protein_token and name == protein_token:
            continue
        non_protein.append((name, count))

    with open(out_top, "w") as fh:
        fh.write("; topol.top rebuilt by replicate_and_fix_top.py\n")
        fh.write("; EGFR TM-JM x%d in CHOL:DIPC:DPSM, Martini 3, -GM3\n\n"
                 % n_protein)
        for inc in include_lines:
            fh.write('#include "%s"\n' % inc)
        fh.write('#include "%s"\n\n' % protein_itp)
        fh.write("[ system ]\n")
        fh.write("EGFR TM-JM x%d raft Martini3 -GM3\n\n" % n_protein)
        fh.write("[ molecules ]\n")
        fh.write("%-10s %d\n" % (protein_name, n_protein))
        for name, count in non_protein:
            fh.write("%-10s %s\n" % (name, count))
    print("wrote %s : %s x%d + %d non-protein species"
          % (out_top, protein_name, n_protein, len(non_protein)))


def compute_grid(n, box_x, box_y, margin):
    """Auto place n copies on a near-square grid inside (box_x, box_y) nm.

    Returns list of (dx, dy) in nm at cell centers, leaving `margin` nm from the
    box edges. Removes the need to hand-write --offsets for a given copy count.
    """
    import math
    ncols = int(math.ceil(math.sqrt(n)))
    nrows = int(math.ceil(n / ncols))
    usable_x = box_x - 2 * margin
    usable_y = box_y - 2 * margin
    offsets = []
    for r in range(nrows):
        for c in range(ncols):
            if len(offsets) >= n:
                break
            dx = margin + (usable_x * (c + 0.5) / ncols)
            dy = margin + (usable_y * (r + 0.5) / nrows)
            offsets.append((dx, dy))
    return offsets


# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="mode", required=True)

    pr = sub.add_parser("replicate")
    pr.add_argument("--in", dest="in_pdb", required=True,
                    help="single-chain CG peptide PDB (martinize2 -x output)")
    pr.add_argument("--out", dest="out_pdb", required=True,
                    help="output multi-chain PDB for insane -f")
    pr.add_argument("--offsets", default=None,
                    help="explicit semicolon list of dx,dy in nm "
                         "(overrides --grid), e.g. '10,10;30,10;10,30;30,30'")
    pr.add_argument("--grid", type=int, default=None,
                    help="auto-place this many copies on a near-square grid")
    pr.add_argument("--box", default=None,
                    help="box X,Y in nm for --grid auto placement, e.g. '40,40'")
    pr.add_argument("--margin", type=float, default=8.0,
                    help="nm kept clear from box edges for --grid (default 8)")

    pf = sub.add_parser("fixtop")
    pf.add_argument("--insane-top", required=True,
                    help="topol.top written by insane")
    pf.add_argument("--out", dest="out_top", required=True)
    pf.add_argument("--protein-itp", required=True,
                    help="protein itp from martinize2, e.g. molecule_0.itp")
    pf.add_argument("--protein-name", required=True,
                    help="[ moleculetype ] name inside the protein itp")
    pf.add_argument("--n-protein", type=int, required=True)
    pf.add_argument("--protein-token", default="",
                    help="protein name as it appears in the insane top "
                         "[ molecules ] block (to drop it); empty if none")
    pf.add_argument("--include", action="append", default=[],
                    help="M3 itp to #include; repeatable, in order")

    args = p.parse_args()

    if args.mode == "replicate":
        if args.offsets:
            offsets = []
            for chunk in args.offsets.split(";"):
                chunk = chunk.strip()
                if not chunk:
                    continue
                dx, dy = chunk.split(",")
                offsets.append((float(dx), float(dy)))
        elif args.grid:
            if not args.box:
                sys.exit("ERROR: --grid requires --box X,Y")
            bx, by = (float(v) for v in args.box.split(","))
            offsets = compute_grid(args.grid, bx, by, args.margin)
        else:
            sys.exit("ERROR: give --offsets or --grid with --box")
        if not offsets:
            sys.exit("ERROR: no offsets produced")
        replicate(args.in_pdb, args.out_pdb, offsets)

    elif args.mode == "fixtop":
        fixtop(args.insane_top, args.out_top, args.protein_itp,
               args.protein_name, args.n_protein, args.include,
               args.protein_token)


if __name__ == "__main__":
    main()
