# psf, crd and pdb for viewing
#
# Sourced by memble.sh. Not a script of its own: it reads and writes the
# variables memble.sh set up, in the order memble.sh sources it.

# ====================================================================
# 7. PSF / CRD / PDB (ParmEd)
# ====================================================================
# ParmEd writes the PSF, the CRD and a PDB for a viewer. The system runs without
# them, so a missing ParmEd is reported and the build carries on.
if ! "$PY" -c 'import parmed' >/dev/null 2>&1; then
  echo ">>> ParmEd is not installed, so system.psf, system.crd and system_parmed.pdb"
  echo "    were not written. The system itself is complete and runs without them."
  echo "    To get the viewer files:  pip install ParmEd"
else
if ! "$PY" - "${ALL[*]}" "$PROT_BLOCKS" "${PART_COUNTS[*]}" <<'PYEOF'
import sys, os, re, string, parmed as pmd
lipids = set(sys.argv[1].split())
prot_blocks = [int(x) for x in sys.argv[2].split()] if sys.argv[2].strip() else []
pcounts = [int(x) for x in sys.argv[3].split()] if len(sys.argv) > 3 and sys.argv[3].strip() else []
water  = {"W", "WF"}
ions   = {"NA", "CL", "ION", "NA+", "CL-"}
SEG_MEMB = "MEMB"; SEG_SOLV = "SOLV"; SEG_ION = "ION"  # CHARMM-GUI conventions

# ParmEd's GROMACS reader only supports 3-point vsite type 1, but Martini 3
# cholesterol (and some other lipids) use other virtual-site constructions.
# Those sections are irrelevant for PSF/segid generation, so build a flattened,
# vsite-stripped copy of the topology just for ParmEd. The real system.top
# (with vsites intact) is what GROMACS uses; this copy never touches the run.
def _inline(path, seen):
    ap = os.path.abspath(path); base = os.path.dirname(ap); out = []
    if ap in seen: return out
    seen.add(ap)
    with open(ap) as fh:
        for ln in fh:
            m = re.match(r'\s*#include\s+"([^"]+)"', ln)
            if m:
                inc = m.group(1); cand = None
                for d in ([''] if os.path.isabs(inc) else [os.getcwd(), base]):
                    p = inc if os.path.isabs(inc) else os.path.join(d, inc)
                    if os.path.exists(p): cand = p; break
                if cand: out += _inline(cand, seen); continue
            out.append(ln)
    return out
flat = _inline('system.top', set())
keep = []; skip = False
for ln in flat:
    st = ln.strip()
    if st.startswith('['):
        skip = st.strip('[] ').lower().startswith('virtual_sites')
    # Drop the data lines of virtual_sites sections (parmed cannot parse them),
    # but always keep preprocessor directives (#ifdef/#ifndef/#else/#endif/
    # #define) so their pairing stays balanced. A virtual_sites section can sit
    # between a #ifdef and its #endif; dropping the directive lines too would
    # leave an orphan #endif and break the parmed read.
    if skip and not st.startswith('#'):
        continue
    keep.append(ln)
open('system_parmed.top', 'w').write(''.join(keep))
top = pmd.load_file('system_parmed.top', xyz='system.gro', parametrize=False)
# protein residues, in order: ncopy TM-JM copies then each partner; one PRO chain per block
prot = [r for r in top.residues if r.name not in lipids
        and r.name not in water and r.name not in ions]
block_sizes = prot_blocks + pcounts
idx = 0; chain = 0
for b in block_sizes:
    for _ in range(b):
        if idx < len(prot):
            prot[idx].segid = "PRO" + string.ascii_uppercase[min(chain, 25)]; idx += 1
    chain += 1
for r in prot[idx:]:
    r.segid = "PRO" + string.ascii_uppercase[min(chain, 25)]
for r in top.residues:
    if r.name in lipids: r.segid = SEG_MEMB
    elif r.name in water: r.segid = SEG_SOLV
    elif r.name in ions: r.segid = SEG_ION
top.save('system.psf', overwrite=True); top.save('system_vmd.psf', vmd=True, overwrite=True)
top.save('system.crd', format='charmmcrd', overwrite=True); top.save('system_parmed.pdb', overwrite=True)
segs = sorted({r.segid for r in top.residues})
print('ParmEd wrote psf/crd/pdb; segids =', segs)
PYEOF
then :; else
  echo ">>> ParmEd could not write system.psf, system.crd and system_parmed.pdb."
  echo "    Those files are for a viewer. The system itself is complete and runs"
  echo "    without them, so the build carries on."
fi
fi
"$GMX" editconf -f system.gro -o system.pdb >/dev/null 2>&1 || true

