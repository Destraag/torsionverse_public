"""
higgs_jamming_engagement_l18_selection_check.py
================================================
Tests the DYNAMICAL framing (author, 2026-09-22) for the T1g/A_g channel-count
cutoff found in higgs_ag_t1g_channel_count_v2.py: rather than picking l<=18
as one of three static cutoffs (12, 18, 20) that all happen to give ratio=2,
ask whether l=18 is SPECIFICALLY selected by the cell's approach to the
Maxwell-critical JAMMED state (n=18=dim(T_1g)*[dim(T_1g)+dim(T_2g)], the
same n=18 already independently validated as the gravity coupling exponent,
ih_lattice_phonon.py IP6-IP9) -- i.e. "Higgs mass = jammed state" (this
framework's SSB=JAMMING picture) should select l=18 as the physically
relevant stopping point, not 12 or 20.

METHOD (four honest, separately-checked pieces, not one combined claim):
  1. Re-ground n=18 directly (3V-E=6 Maxwell count, recomputed here, not just
     cited) and its 3(spatial)x6(soft-mode) structure [ih_lattice_phonon.py].
  2. Re-derive the A_g/T1g projection table (same character-projection method
     as higgs_ag_t1g_channel_count_v2.py) and cross-validate against its
     exact published l-lists -- this script does NOT introduce a new
     computation method here, only reuses/checks the established one.
  3. THE ACTUAL NEW TEST: is l=18 uniquely tied to Maxwell-rigidity, or is it
     just as arbitrary as 12 (=V, vertex count) and 20 (=F, face count) --
     both ALSO simple products of established I_h irrep dimensions? Full
     l=0..24 ratio landscape computed and cross-referenced against every
     dim(irrep_a)*dim(irrep_b) product, to see honestly whether 18 stands out
     for a MECHANISTIC reason (product involves the Maxwell-6 sum
     specifically) or merely a NUMERICAL one (any of several products land
     on an integer that happens to give ratio 2).
  4. GENUINE NEW CALCULATION: does the jamming transition, tracked via actual
     rigidity-matrix rank (reusing cell_self_stress_count.py's established
     method: rank(R) for the 12-vertex/30-edge network), engage GRADUALLY
     (which could motivate an intermediate l-cutoff) or as a hard all-or-
     nothing threshold at the full 30-edge count? Tested by removing each of
     the 30 edges in turn and recomputing rank -- deterministic, not sampled.

HONEST SCOPE, STATED UP FRONT (see also Section 5 verdict): even a fully
successful outcome here can only corroborate the FACTOR-OF-2 classification
between A_g (scalar) and T_1g (vector) channels -- doc_higgs.txt Section 3a's
linking-number theorem (torus-knot pi-rotation symmetry, exact, purely
topological) ALREADY derives this same factor of 2 independently, with no
l-cutoff needed at all. A ratio-only argument is structurally blind to
absolute magnitude: it cannot, even in the best case, derive WHY the base
coefficient is alpha/pi (rather than alpha/(2*pi) or any other value) --
that remains the actual open question, already honestly flagged in
doc_higgs.txt as a Standard Model import, not independently derived here.

Run: python analysis/higgs/higgs_jamming_engagement_l18_selection_check.py
"""

import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
import numpy as np

pi = math.pi
phi = (1 + math.sqrt(5)) / 2

SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "FAIL"
    results.append((name, s, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


print(SEP)
print("higgs_jamming_engagement_l18_selection_check.py")
print("Does the Maxwell-jamming transition specifically select l<=18,")
print("or is that cutoff just as arbitrary as l<=12 or l<=20?")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("SECTION 1: RE-GROUND n=18 DIRECTLY (not just cited)")
print(SEP2)
print()

verts_raw = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        verts_raw.append([0, s1, s2 * phi])
        verts_raw.append([s1, s2 * phi, 0])
        verts_raw.append([s2 * phi, 0, s1])
verts = np.array(verts_raw, dtype=float)
V = len(verts)

all_dists = sorted(np.linalg.norm(verts[i] - verts[j])
                    for i in range(V) for j in range(i + 1, V))
r_nn = all_dists[0]
tol = r_nn * 0.05
edges = [(i, j) for i in range(V) for j in range(i + 1, V)
         if abs(np.linalg.norm(verts[i] - verts[j]) - r_nn) < tol]
E = len(edges)

print(f"  V={V}, E={E}, 3V-E={3*V-E}")
check("J1 icosahedron: V=12, E=30, 3V-E=6 (Maxwell critical)",
      V == 12 and E == 30 and 3 * V - E == 6, f"V={V}, E={E}, 3V-E={3*V-E}")

dim = {'A_g': 1, 'T_1g': 3, 'T_2g': 3, 'G_g': 4, 'H_g': 5}
maxwell6 = dim['T_1g'] + dim['T_2g']
n_gravity = dim['T_1g'] * maxwell6
print(f"  Maxwell-6 = dim(T_1g)+dim(T_2g) = {maxwell6}")
print(f"  n = dim(T_1g) x Maxwell-6 = {dim['T_1g']} x {maxwell6} = {n_gravity}")
check("J2 n=18 = dim(T_1g) x [dim(T_1g)+dim(T_2g)] (ih_lattice_phonon.py IP8)",
      n_gravity == 18, f"3 x 6 = {n_gravity}")

# =============================================================================
print()
print(SEP2)
print("SECTION 2: RE-DERIVE A_g/T1g PROJECTION TABLE, CROSS-CHECK vs v2")
print(SEP2)
print()

CHI = {
    'A_g':  [1, 1, 1, 1, 1],
    'T_1g': [3, phi, 1 - phi, 0, -1],
    'T_2g': [3, 1 - phi, phi, 0, -1],
    'G_g':  [4, -1, -1, 1, 0],
    'H_g':  [5, 0, 0, -1, 1],
}
I_rotations = [(1, 0), (12, 2 * pi / 5), (12, 4 * pi / 5), (20, 2 * pi / 3), (15, pi)]
order_I = 60


def chi_Yl(l, theta):
    if abs(math.sin(theta / 2)) < 1e-12:
        return 2 * l + 1
    return math.sin((l + 0.5) * theta) / math.sin(theta / 2)


def n_irrep(l, irrep):
    total = sum(mult * c * chi_Yl(l, theta)
                for (mult, theta), c in zip(I_rotations, CHI[irrep]))
    val = total / order_I
    return int(round(val)) if abs(val - round(val)) < 1e-6 else 0


L_MAX = 24
table = {l: {irr: n_irrep(l, irr) for irr in CHI} for l in range(L_MAX + 1)}
Ag_at = [l for l in range(L_MAX + 1) if table[l]['A_g'] > 0]
T1g_at = [l for l in range(L_MAX + 1) if table[l]['T_1g'] > 0]

# Exact lists as published by higgs_ag_t1g_channel_count_v2.py (cross-check target)
Ag_at_v2 = [0, 6, 10, 12, 15, 16, 18, 20, 21, 22, 24]
T1g_at_v2 = [1, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

print(f"  A_g appears at l = {Ag_at}")
print(f"  T1g appears at l = {T1g_at}")
check("V1 A_g channel list matches higgs_ag_t1g_channel_count_v2.py exactly",
      Ag_at == Ag_at_v2, f"computed={Ag_at}")
check("V2 T1g channel list matches higgs_ag_t1g_channel_count_v2.py exactly",
      T1g_at == T1g_at_v2, f"computed={T1g_at}")

# =============================================================================
print()
print(SEP2)
print("SECTION 3: IS l=18 UNIQUELY TIED TO MAXWELL-RIGIDITY, OR JUST AS")
print("           ARBITRARY AS l=12 (=V) AND l=20 (=F)?")
print(SEP2)
print()

print(f"  {'l':>4}  {'n_Ag_cum':>9} {'n_T1g_cum':>10}  {'ratio':>7}")
ratio2_at = []
for l in range(1, L_MAX + 1):
    n_ag = sum(1 for x in Ag_at if x <= l)
    n_t1g = sum(1 for x in T1g_at if x <= l)
    ratio = n_t1g / n_ag if n_ag else float('nan')
    flag = "  <-- ratio=2" if n_ag and abs(ratio - 2.0) < 1e-9 else ""
    if flag:
        ratio2_at.append(l)
    print(f"  {l:>4}  {n_ag:>9} {n_t1g:>10}  {ratio:>7.4f}{flag}")

print()
print(f"  Full set of l giving EXACTLY ratio=2: {ratio2_at}")

# Candidate products of two established I_h irrep dimensions (excludes A_g=1,
# trivial), to see how many small integers <=24 are "just as simple" as 18.
irr_names = ['T_1g', 'T_2g', 'G_g', 'H_g']
products = {}
for a in irr_names:
    for b in irr_names:
        p = dim[a] * dim[b]
        if p <= L_MAX:
            products.setdefault(p, set()).add(tuple(sorted([a, b])))

print()
print("  Products dim(irrep_a)*dim(irrep_b) <= 24 (candidate 'simple' numbers):")
for p in sorted(products):
    print(f"    {p:>3} = {' , '.join(f'{a}*{b}' for a,b in sorted(products[p]))}")

# V=12 and F=20 are the icosahedron's own plain geometric invariants --
# unrelated to rigidity/jamming, just shape combinatorics (Euler's formula).
F = E - V + 2  # Euler's formula V-E+F=2 => F=E-V+2
print()
print(f"  Icosahedron F (faces, Euler's formula V-E+F=2): F = {F}")
check("S1 l=12 coincides with V (vertex count, plain geometry, NOT a rigidity count)",
      12 in ratio2_at and V == 12, f"V={V}")
check("S2 l=20 coincides with F (face count, plain geometry, NOT a rigidity count)",
      20 in ratio2_at and F == 20, f"F={F}")
check("S3 12 and 20 are ALSO simple products of two established irrep dims (3*4, 4*5)",
      12 in products and 20 in products,
      f"12={products.get(12)}, 20={products.get(20)}")
check("S4 18 is NOT a product of two DIFFERENT single irreps (no G_g/H_g factor pair gives 18)",
      18 not in products,
      f"18 absent from {sorted(products.keys())} -- 18 only arises as dim(T_1g)*[dim(T_1g)+dim(T_2g)]," 
      f" i.e. one irrep times the MAXWELL-RIGIDITY SUM, not a plain 2-irrep product")

print()
print("  HONEST READING: 12 and 20 are EQUALLY 'simple' numerically (V, F, and")
print("  2-irrep products) -- simplicity alone does not favor 18. What DOES")
print("  distinguish 18: only 18 is independently anchored to the MAXWELL-")
print("  RIGIDITY/JAMMING transition itself (validated separately as the")
print("  gravity exponent to 0.27%, ih_lattice_phonon.py IP7) -- 12 and 20 are")
print("  never used anywhere in this repo as a rigidity/engagement count, only")
print("  as static shape facts (V, F). The selection of l<=18 over l<=12/20 is")
print("  a MECHANISM-based argument (both concern the SAME jamming transition:")
print("  Higgs mass <-> Maxwell-critical jamming, per SSB=JAMMING), not a")
print("  numerical-coincidence one -- but it is still an ANALOGY/identification")
print("  (l-cutoff = engagement count), not an independent derivation of that")
print("  identification itself.")

# =============================================================================
print()
print(SEP2)
print("SECTION 4: DOES THE JAMMING TRANSITION ENGAGE GRADUALLY (l-cutoff-like)")
print("           OR AS AN ALL-OR-NOTHING THRESHOLD?")
print(SEP2)
print()


def rigidity_matrix(edge_list):
    R = np.zeros((len(edge_list), 3 * V))
    for k, (i, j) in enumerate(edge_list):
        rij = verts[j] - verts[i]
        rij_hat = rij / np.linalg.norm(rij)
        R[k, 3 * i:3 * i + 3] = rij_hat
        R[k, 3 * j:3 * j + 3] = -rij_hat
    return R


R_full = rigidity_matrix(edges)
rank_full = np.linalg.matrix_rank(R_full)
floppy_full = 3 * V - rank_full
print(f"  Full network (all {E} edges): rank(R)={rank_full}, floppy modes={floppy_full}")
check("R1 full network rank(R)=30, floppy=6 (matches cell_self_stress_count.py)",
      rank_full == 30 and floppy_full == 6, f"rank={rank_full}, floppy={floppy_full}")

# Deterministic single-edge-removal test (not sampled -- all 30 cases checked)
drop_ranks = []
for k in range(E):
    sub_edges = edges[:k] + edges[k + 1:]
    R_sub = rigidity_matrix(sub_edges)
    rank_sub = np.linalg.matrix_rank(R_sub)
    drop_ranks.append(rank_sub)

all_drop_to_29 = all(r == 29 for r in drop_ranks)
print(f"  Removing each single edge (all {E} cases): resulting rank(R) values = "
      f"{sorted(set(drop_ranks))}")
check("R2 EVERY single-edge removal drops rank to exactly 29 (floppy=7)",
      all_drop_to_29, f"ranks found: {sorted(set(drop_ranks))}")
check("R3 zero redundancy: no 29-edge subset reaches floppy=6 -- jamming is an "
      "all-or-nothing threshold at the FULL edge count, not a gradual ramp",
      all_drop_to_29,
      "every edge is independently load-bearing (isostatic, S=0 states of self-stress)")

print()
print("  CONSEQUENCE: there is no gradual, computed edge-engagement sequence")
print("  that approaches jamming smoothly -- the transition is a hard cliff at")
print("  m=30 (100% of constraints). This means Section 3's l<-->engagement-")
print("  count identification CANNOT be independently verified via an actual")
print("  intermediate dynamical sequence -- no such graded sequence exists in")
print("  the discrete constraint count. The l=18 selection stands ONLY on the")
print("  dimensional analogy (Section 1+3), not on a separately-computed")
print("  dynamical bridge.")

# =============================================================================
print()
print(SEP)
print("SECTION 5: FINAL HONEST VERDICT")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass}/{len(results)} checks PASS")
print()
print("  1. n=18's own structure (Section 1) and the A_g/T1g channel table")
print("     (Section 2) are both solid, cross-validated ground truth.")
print("  2. l=18 IS distinguishable from l=12/20 (Section 3) -- but by a")
print("     MECHANISM argument (same jamming transition as gravity's n=18),")
print("     not by numerical simplicity (12 and 20 are equally simple: V and")
print("     F, and equally valid 2-irrep products). This is genuine, if")
print("     partial, motivation -- not a coincidence, but also not a proof,")
print("     since it rests on identifying an l-cutoff with an engagement")
print("     count, an assumption this script does not independently derive.")
print("  3. Section 4 confirms no graded dynamical sequence exists to")
print("     independently verify that identification -- the real jamming")
print("     transition (edge/constraint count) is an all-or-nothing cliff,")
print("     not a ramp through intermediate l-like stages.")
print("  4. MOST IMPORTANT CAVEAT: even accepting l<=18 fully, this only")
print("     recovers the FACTOR-OF-2 (T1g channels = 2x A_g channels) --")
print("     doc_higgs.txt Section 3a's linking-number theorem (torus-knot")
print("     pi-rotation symmetry) ALREADY derives that same factor of 2")
print("     exactly, via a completely different, cleaner, purely topological")
print("     argument, with no l-cutoff needed at all. A ratio-only result")
print("     is structurally blind to absolute magnitude -- it cannot derive")
print("     WHY the base coefficient is alpha/pi rather than alpha/(2*pi) or")
print("     any other value. That remains the genuinely open question,")
print("     already honestly flagged in doc_higgs.txt as a Standard Model")
print("     import. This script's result -- even at its best -- is a second,")
print("     independent consistency check on the already-known factor of 2,")
print("     not progress toward deriving alpha/pi's own magnitude.")
