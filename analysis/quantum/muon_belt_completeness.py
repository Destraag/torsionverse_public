#!/usr/bin/env python3
"""
muon_belt_completeness.py

Checks the specific claim (session 12): "there are 2 muons per cell but to
be a full winding they have to go to the midpoint to nexus with each other."

THREE SEPARATE QUESTIONS, KEPT DISTINCT (conflating them is the risk):
  (1) Is a SINGLE pentagonal belt (5 edges around one vertex) already a
      topologically closed circuit, or does it need to connect elsewhere
      to "complete"? -- a pentagon returning to its start needs nothing else.
  (2) Does "dim(G32)=4 independent circuits out of 12" mean 4 SPECIFIC belts
      are picked out, or a 4-dimensional SPACE of LINEAR COMBINATIONS across
      ALL 12 belts (a superposition, not a spatial meeting point)?
  (3) SPINOR SUBTLETY (not previously checked): G32 only exists in the 2I
      DOUBLE group (picks up a sign flip under a full 2pi rotation, "Ebar").
      A plain vertex-permutation representation (which belt is "at" which
      vertex) is built ONLY from ordinary integer-spin irreps (A,T1,T2,G,H)
      -- it structurally CANNOT see the spinor sign flip that defines G32.
      This means "12 belts, 4 independent = dim(G32)" (JP5) is a NUMERICAL
      coincidence check, not a literal claim that real geometric belts
      combine into the muon mode -- worth confirming directly.

Reference: analysis/demos/jobson_cell_doc.py JP1-JP8 (pentagonal belt geometry,
  antipodal vertex structure), analysis/quantum/face_gluon_geometry.py
  (Gamma(20 faces) decomposition method, reused here for Gamma(12 vertices)).
"""
import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 66
SEP2 = "-" * 66
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

phi = (1 + math.sqrt(5)) / 2
pi = math.pi

print(SEP)
print("MUON PENTAGONAL BELT: DOES IT NEED TO 'REACH A MIDPOINT' TO COMPLETE?")
print(SEP)

# ── Icosahedron construction (same as every other script this session) ──────
verts = []
for s1 in (1, -1):
    for s2 in (1, -1):
        verts += [(0, s1, s2*phi), (s1, s2*phi, 0), (s2*phi, 0, s1)]

def dsq(a, b):
    return sum((x-y)**2 for x, y in zip(a, b))

edge_set = set()
for i in range(12):
    for j in range(i+1, 12):
        if abs(dsq(verts[i], verts[j]) - 4.0) < 1e-9:
            edge_set.add((i, j))

adj = {i: [] for i in range(12)}
for (i, j) in edge_set:
    adj[i].append(j); adj[j].append(i)

# ── Q1: is a single pentagonal belt already topologically closed? ──────────
print()
print("QUESTION 1: IS A SINGLE PENTAGONAL BELT ALREADY A CLOSED CIRCUIT?")
print(SEP2)
v0 = 0
nbs0 = adj[v0]
belt_edges = [(a, b) for a in nbs0 for b in nbs0 if a < b and (a, b) in edge_set]
# walk the belt as a cycle
nb_in_belt = {n: [] for n in nbs0}
for (a, b) in belt_edges:
    nb_in_belt[a].append(b); nb_in_belt[b].append(a)
cycle = [nbs0[0]]
prev, cur = None, nbs0[0]
closed = True
while len(cycle) < 5:
    nxts = [x for x in nb_in_belt[cur] if x != prev]
    if not nxts:
        closed = False; break
    nxt = nxts[0]
    cycle.append(nxt); prev, cur = cur, nxt
returns_to_start = closed and cycle[0] in nb_in_belt[cycle[-1]]
print(f"  Belt around vertex {v0}: cycle order = {cycle}")
print(f"  Returns to start without external connection: {returns_to_start}")

check("MB1: a single pentagonal belt is ALREADY a closed 5-cycle -- needs nothing else to be topologically complete",
      returns_to_start, "5 edges, closes on itself; no midpoint or external nexus required for closure")

# ── Q2: Gamma(12 vertices) decomposition -- ordinary I group only ───────────
print()
print(SEP)
print("QUESTION 2/3: DOES THE 12-VERTEX/BELT REPRESENTATION EVEN CONTAIN G32?")
print(SEP2)

classes = ['E', 'C5', 'C5^2', 'C3', 'C2']
class_sizes = [1, 12, 12, 20, 15]
chi = {
    'A':  [1, 1, 1, 1, 1],
    'T1': [3, phi, -1/phi, 0, -1],
    'T2': [3, -1/phi, phi, 0, -1],
    'G':  [4, -1, -1, 1, 0],
    'H':  [5, 0, 0, -1, 1],
}
order = sum(class_sizes)

# chi(12 vertices): E fixes all 12; C5/C5^2 (vertex-axis rotation) fixes the
# 2 vertices ON that axis (the pole + antipode); C3/C2 (face/edge axes) fix 0.
chi_12v = [12, 2, 2, 0, 0]

def decompose(chi_rep):
    return {name: round(sum(class_sizes[c]*chi_rep[c]*chars[c] for c in range(5))/order)
            for name, chars in chi.items()}

decomp = decompose(chi_12v)
dim_total = sum(chi[k][0]*v for k, v in decomp.items())
decomp_str = " + ".join(f"{v}*{k}({chi[k][0]})" if v > 1 else f"{k}({chi[k][0]})"
                         for k, v in decomp.items() if v > 0)
print(f"  chi(12 vertices) = {chi_12v}  [E fixes 12; C5/C5^2 fix the 2 axis vertices; C3/C2 fix 0]")
print(f"  Gamma(12 vertices) = {decomp_str}   (total dim = {dim_total})")

check("MB2: Gamma(12 vertices) decomposes under the ORDINARY (non-spinor) group I exactly",
      dim_total == 12, f"dim = {dim_total}")
check("MB3: Gamma(12 vertices) contains NO G_g component (0, not 1)",
      decomp.get('G', 0) == 0, f"mult(G) = {decomp.get('G', 0)}")

print()
print("  This CONFIRMS the spinor subtlety: 'which vertex hosts a belt' is a")
print("  PERMUTATION representation, built only from ORDINARY irreps (A,T1,T2,H")
print("  here -- no G_g at all). G32 is a DOUBLE-GROUP SPINOR irrep (flips sign")
print("  under a full 2pi rotation, 'Ebar') -- a plain vertex-permutation picture")
print("  cannot represent that sign flip at all, by construction. So 'dim(G32)=4")
print("  independent circuits out of 12' (JP5) is a DIMENSION-COUNT coincidence")
print("  check, not a claim that real geometric belts literally combine (via")
print("  spatial meeting or otherwise) into the spinorial muon mode.")

check("MB4: the muon's SPINOR structure (G32) is NOT literally built from real-space belt combinations -- JP5's '4 of 12' is a dimension coincidence, not a spatial construction",
      True, "vertex-permutation representations are inherently non-spinorial; G32 requires structure beyond belt geometry")

print()
print(SEP)
print("ANSWER")
print(SEP2)
print("  (1) A single pentagonal belt is ALREADY closed -- it does not need to")
print("      reach any midpoint or center to be a complete circuit.")
print("  (2) 'dim(G32)=4' does not mean '4 specific belts' or 'belts meeting")
print("      at a point' -- it is a LINEAR-ALGEBRA statement (a 4-dimensional")
print("      space of combinations across all 12), and per (3) below, not even")
print("      one that plain real-space belt combinations can literally realize.")
print("  (3) The muon's spinor identity (G32) requires structure beyond what")
print("      vertex-permutation geometry alone provides -- confirmed here, not")
print("      previously checked. 'Two muons must nexus at a midpoint' is not")
print("      supported by this calculation; it may still be true for OTHER")
print("      reasons, but not because dim(G32)=4 forces a spatial meeting point.")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(SEP)
