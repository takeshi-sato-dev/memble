#!/usr/bin/env python3
"""Shift the protein in z by a small offset, for fine tuning how deep the
protein sits in the membrane.

The build centers the transmembrane core at the bilayer midplane. A researcher
may want the protein a little higher or lower than that, found by trying a few
values and rebuilding. This helper shifts only the protein beads (amino acid
residues) in z by --dz nanometers, leaving lipids, water, and ions in place.
The later declash step relaxes the lipids around the shifted protein, and the
position restraint reference is taken from the shifted coordinates, so the
offset is held through equilibration.

Default behavior with --dz 0 is a no-op copy, so the build can always call it.
"""
import argparse

# Amino acid residue names that mark a protein bead in a Martini structure.
AA = {
    "GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "SER",
    "THR", "CYS", "TYR", "ASN", "GLN", "ASP", "GLU", "LYS", "ARG", "HIS",
    "HSD", "HSE", "HSP", "HID", "HIE", "HIP",
}


def shift_gro(lines, dz):
    """Return new gro lines with protein bead z shifted by dz (nm).

    A gro atom line is fixed width: residue number in columns 1 to 5, residue
    name in 6 to 10, atom name in 11 to 15, atom number in 16 to 20, then x, y,
    z each in an 8.3 field. Only the z field of protein beads is changed, and
    the rest of every line is preserved exactly.
    """
    n = int(lines[1])
    head = lines[:2]
    body = lines[2:2 + n]
    box = lines[2 + n:]
    out_body = []
    moved = 0
    for ln in body:
        resn = ln[5:10].strip().upper()
        if resn in AA and len(ln) >= 44:
            z = float(ln[36:44]) + dz
            ln = ln[:36] + "%8.3f" % z + ln[44:]
            moved += 1
        out_body.append(ln)
    return head + out_body + box, moved


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--dz", type=float, required=True,
                    help="z shift in nm applied to protein beads (0 is a no-op)")
    ap.add_argument("--out", default=None,
                    help="output gro path; defaults to overwriting --gro")
    args = ap.parse_args()

    out = args.out or args.gro
    if args.dz == 0.0:
        if args.out and args.out != args.gro:
            with open(args.gro) as fh:
                data = fh.read()
            with open(out, "w") as fh:
                fh.write(data)
        print("shift_protein_z: dz=0, no shift applied")
        return

    with open(args.gro) as fh:
        lines = fh.read().splitlines()
    new_lines, moved = shift_gro(lines, args.dz)
    with open(out, "w") as fh:
        fh.write("\n".join(new_lines) + "\n")
    print("shift_protein_z: moved %d protein beads by dz=%.3f nm" % (moved, args.dz))


if __name__ == "__main__":
    main()
