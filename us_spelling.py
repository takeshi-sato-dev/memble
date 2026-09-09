"""Rewrite British spellings as American ones, in place.

Usage:  python3 us_spelling.py README.md tests/test_pipeline.py
"""
import re, sys
MAP = [
 (r'\bminimis(ation|ations)\b', r'minimiz\1'),
 (r'\bminimise\b','minimize'), (r'\bminimises\b','minimizes'), (r'\bminimised\b','minimized'), (r'\bminimising\b','minimizing'),
 (r'\bcentre\b','center'), (r'\bcentres\b','centers'), (r'\bcentred\b','centered'), (r'\bcentring\b','centering'),
 (r'\bCentre\b','Center'), (r'\bCentred\b','Centered'), (r'\bCentring\b','Centering'),
 (r'\brecentring\b','recentering'), (r'\bRecentring\b','Recentering'),
 (r'\bbehaviour\b','behavior'), (r'\bbehaviours\b','behaviors'),
 (r'\bneighbour\b','neighbor'), (r'\bneighbours\b','neighbors'), (r'\bneighbouring\b','neighboring'),
 (r'\bneighbourhood\b','neighborhood'), (r'\bneighbourhoods\b','neighborhoods'),
 (r'\blicence\b','license'),
 (r'\bgrey\b','gray'), (r'\bGrey\b','Gray'),
 (r'\boptimiser\b','optimizer'), (r'\boptimisation\b','optimization'),
 (r'\bnormalise\b','normalize'), (r'\bnormalised\b','normalized'),
 (r'\banalyse\b','analyze'), (r'\banalysed\b','analyzed'), (r'\banalysing\b','analyzing'),
 (r'\bmodelling\b','modeling'), (r'\blabelled\b','labeled'), (r'\blabelling\b','labeling'),
 (r'\btowards\b','toward'), (r'\bwhilst\b','while'), (r'\bamongst\b','among'),
]
def convert(t):
    for a,b in MAP: t = re.sub(a,b,t)
    return t
for f in sys.argv[1:]:
    s=open(f).read(); n=convert(s)
    if n!=s: open(f,'w').write(n); print("changed", f)
    else: print("unchanged", f)
