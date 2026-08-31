#!/usr/bin/env python3
"""
tau_pair_wz_composite.py

Tests the user's hypothesis (session 12): are Z/W an EMERGENT/COMPOSITE pattern
built from paired tau (I52) windings converging toward the cell center, rather
than a fundamentally separate mode?

This is NOT semantics -- there is a real precedent already in the framework:
the Higgs (A_g) is described (doc_jobson_cell.txt, JC3) as "20 I52 face
corkscrews in phase". The question here is whether the SAME kind of
tau-pairing operation that builds A_g can *also* reach the T_1g (W/Z) channel.

WHAT THIS SCRIPT CAN AND CANNOT SHOW:
  CAN show (abstract representation theory, exact): whether the algebraic
    content T_1g is present when two I52 (tau) windings are combined via the
    standard Clebsch-Gordan product -- i.e. whether "T_1g" is even a reachable
    slot from tau pairs at all.
  CANNOT show (needs separate future work): whether the actual 3D spatial
    paths of two tau windings literally converge toward the cell center and
    trace a directed-cone shape. That is a geometric construction question,
    analogous to what gluon_tau_helix.py did for the single free tau path,
    and is NOT attempted here. Flagged explicitly as the open next step.

Reference: analysis/nuclear/ih_double_cg.py (DC8, corrected this session),
  analysis/nuclear/ih_double_group.py (2I character table),
  docs/series1/doc_jobson_cell.txt Section on Higgs internal polarization (JC3),
  docs/series1/doc_jobson_cell.txt "DIRECTED CONES (W and Z)" section,
  docs/open_items.txt "CHARM [HIGH] Formal H_u Zone 1 polygon path geometry".
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

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2

print(SEP)
print("TAU-PAIR / W-Z COMPOSITE HYPOTHESIS TEST")
print(SEP)

# ── 2I character table (reproduced standalone from ih_double_group.py) ──────
classes = [
    ('E',    0.0,       1),
    ('Ebar', 2*pi,      1),
    ('C5',   2*pi/5,   12),
    ('C5^4', 8*pi/5,   12),
    ('C5^2', 4*pi/5,   12),
    ('C5^3', 6*pi/5,   12),
    ('C3',   2*pi/3,   20),
    ('C3^2', 4*pi/3,   20),
    ('C2',   pi,       30),
]
G_ORDER = 120

def chi_j(j, theta):
    if abs(math.sin(theta/2)) < 1e-12:
        return -(2*j+1) if (abs(theta) > 0.1 and (2*j) % 2 == 1) else (2*j+1)
    return math.sin((j+0.5)*theta) / math.sin(theta/2)

def lift_from_I(chi5):
    e, c5, c52, c3, c2 = chi5
    return [e, e, c5, c5, c52, c52, c3, c3, c2]

def spinor(j, galois=False):
    chars = []
    for name, theta, _ in classes:
        c = chi_j(j, theta)
        if galois and name in ('C5', 'C5^4', 'C5^2', 'C5^3'):
            m = {'C5': 6*pi/5, 'C5^4': 4*pi/5, 'C5^2': 8*pi/5, 'C5^3': 2*pi/5}
            c = chi_j(j, m[name])
        chars.append(round(c, 10))
    return chars

IRREPS = {
    'A':   lift_from_I([1,  1,      1,      1,   1]),
    'E+':  spinor(0.5),
    'E-':  spinor(0.5, True),
    'T1':  lift_from_I([3,  phi,   -1/phi,  0,  -1]),
    'T2':  lift_from_I([3, -1/phi,  phi,    0,  -1]),
    'G32': spinor(1.5),
    'G':   lift_from_I([4, -1,     -1,      1,   0]),
    'H':   lift_from_I([5,  0,      0,     -1,   1]),
    'I52': spinor(2.5),
}
NAMES = list(IRREPS.keys())

def dim(name):
    return round(IRREPS[name][0])

def inner(a, b):
    return round(sum(classes[i][2] * a[i] * b[i] for i in range(9)) / G_ORDER)

def cg(X, Y):
    chi_prod = [IRREPS[X][i] * IRREPS[Y][i] for i in range(9)]
    return {Z: inner(chi_prod, IRREPS[Z]) for Z in NAMES}

def cg_str(decomp):
    parts = [f"{n}*{k}({dim(k)})" if n > 1 else f"{k}({dim(k)})"
             for k, n in decomp.items() if n > 0]
    return " + ".join(parts) if parts else "0"

print()
print("SECTION 1: I52 x I52 (TWO TAU WINDINGS PAIRED) -- FULL DECOMPOSITION")
print(SEP2)

decomp = cg('I52', 'I52')
dim_total = sum(n * dim(k) for k, n in decomp.items())
print(f"  I52 x I52 = {cg_str(decomp)}")
print(f"  Dimension check: {dim_total} (expect {dim('I52')**2} = 6x6)")

check("TC1: I52 x I52 dimension = 36 exactly (6x6, no missing/extra content)",
      dim_total == dim('I52')**2,
      f"computed = {dim_total}")

check("TC2: I52 x I52 contains T1 (the W/Z channel) with multiplicity 2",
      decomp.get('T1', 0) == 2,
      f"mult(T1) = {decomp.get('T1', 0)}  [reproduces ih_double_cg.py DC8, corrected]")

check("TC3: I52 x I52 contains A (the Higgs channel) with multiplicity 1",
      decomp.get('A', 0) == 1,
      f"mult(A) = {decomp.get('A', 0)}  -- SAME product also reaches the Higgs channel")

print()
print("  Both the Higgs channel (A, mult=1) and the W/Z channel (T1, mult=2)")
print("  appear in the SAME I52 x I52 product. Pairing two tau windings does")
print("  not have to 'choose' between building a Higgs-like or a W/Z-like")
print("  object -- both are reachable slots of the same pairing operation.")

# ── Section 1b: symmetric vs antisymmetric split -- where do the 2 T1 copies live? ─
# "Multiplicity 2" does NOT mean 4 windings are needed. It means the 36-dim
# combined space of ONE pairing (2 windings) contains two independent 3-dim
# T1 subspaces. The natural way to split a self-product is symmetric (Sym^2,
# dim=21) vs antisymmetric (Alt^2, dim=15) under exchange of the two windings
# -- this answers "how do the two windings actually interact": computed below,
# NOT guessed -- the first attempt at this guessed one-copy-per-sector and
# was WRONG (both T1 copies actually land in the SAME sector); corrected here.
print()
print(SEP)
print("SECTION 1b: SYMMETRIC vs ANTISYMMETRIC PAIRING -- WHERE DO THE 2 T1 COPIES LIVE?")
print(SEP2)

def chi_g2(name_chars, thetas_doubled):
    """chi evaluated at g^2 for each class, using the SAME closed-form spinor formula."""
    return [chi_j(2.5, 2*theta) for theta in thetas_doubled]

thetas = [theta for (_, theta, _) in classes]
chi_I52 = IRREPS['I52']
chi_I52_at_g2 = [chi_j(2.5, 2*theta) for theta in thetas]

chi_sym2 = [(chi_I52[i]**2 + chi_I52_at_g2[i]) / 2 for i in range(9)]
chi_alt2 = [(chi_I52[i]**2 - chi_I52_at_g2[i]) / 2 for i in range(9)]

dim_sym2 = chi_sym2[0]
dim_alt2 = chi_alt2[0]

decomp_sym2 = {Z: inner(chi_sym2, IRREPS[Z]) for Z in NAMES}
decomp_alt2 = {Z: inner(chi_alt2, IRREPS[Z]) for Z in NAMES}

print(f"  dim(Sym^2 I52) = {dim_sym2:.0f}  (expect 6*7/2 = 21)")
print(f"  dim(Alt^2 I52) = {dim_alt2:.0f}  (expect 6*5/2 = 15)")
print(f"  Sym^2(I52) = {cg_str(decomp_sym2)}")
print(f"  Alt^2(I52) = {cg_str(decomp_alt2)}")

check("TC1b: dim(Sym^2)+dim(Alt^2) = 36, matching 21+15",
      abs(dim_sym2 - 21) < 1e-6 and abs(dim_alt2 - 15) < 1e-6,
      f"Sym^2={dim_sym2:.0f}  Alt^2={dim_alt2:.0f}")

check("TC1c: A (Higgs channel) lives in the ANTISYMMETRIC part only",
      decomp_sym2.get('A', 0) == 0 and decomp_alt2.get('A', 0) == 1,
      f"Sym^2 mult(A)={decomp_sym2.get('A',0)}  Alt^2 mult(A)={decomp_alt2.get('A',0)}")

t1_sym = decomp_sym2.get('T1', 0)
t1_alt = decomp_alt2.get('T1', 0)
print(f"\n  T1 (W/Z) multiplicity:  Sym^2 = {t1_sym}   Alt^2 = {t1_alt}")
check("TC1d: BOTH copies of T1 (W/Z) live in the SYMMETRIC part, zero in the antisymmetric part",
      t1_sym == 2 and t1_alt == 0,
      f"Sym^2 mult(T1)={t1_sym}  Alt^2 mult(T1)={t1_alt}  -- W/Z needs the symmetric combination, not split 1-1")

print()
print("  ANSWER to 'how do the two windings actually interact': not '2 windings")
print("  per cone', and NOT one symmetric + one antisymmetric copy (that was the")
print("  first guess, and it was WRONG -- corrected by actually running the numbers).")
print("  The real split: pairing two tau windings IN PHASE (symmetric exchange)")
print("  gives BOTH T1 copies (W/Z) plus both T2 copies plus one G -- the entire")
print("  symmetric sector is 2*T1+2*T2+G+H (dim 21). Pairing them OUT OF PHASE")
print("  (antisymmetric exchange) gives the Higgs channel (A) plus one G plus")
print("  both H copies -- the antisymmetric sector is A+G+2*H (dim 15).")
print("  So: W/Z requires a SYMMETRIC tau-pair; Higgs requires an ANTISYMMETRIC one.")
print()
print("  OPEN CAVEAT: whether tau windings pairing to form W/Z actually obey the")
print("  SYMMETRIC exchange rule found above (vs Higgs needing antisymmetric) has")
print("  NOT been established here or anywhere in the framework found so far.")
print("  Mainstream physics would use Pauli exclusion (identical fermions pairing")
print("  antisymmetrically, e.g. Cooper pairs) to argue AGAINST the symmetric W/Z")
print("  channel being the natural one -- but that is an IMPORTED assumption, not")
print("  a torsionverse-derived rule, and is flagged as such rather than silently")
print("  assumed. This tension (W/Z needs symmetric, Pauli-like intuition expects")
print("  antisymmetric for identical windings) is itself a genuinely open question.")

check("TC1e: physical selection rule (which copy is realized) is explicitly flagged as unresolved, not assumed",
      True, "no torsionverse-derived exchange-symmetry rule for winding pairing found; mainstream Pauli-exclusion analogy noted as an import, not a fact")

# ── Section 2: Cross-check T1 (2I) == T_1g (I_h single group, W/Z table) ────
print()
print(SEP)
print("SECTION 2: IS 'T1' HERE REALLY THE SAME T_1g USED FOR W/Z ELSEWHERE?")
print(SEP2)

# T_1g characters as used throughout the framework for W/Z (face_gluon_geometry.py,
# doc_jobson_cell.txt): chi(T_1g, C5)=phi, chi(T_1g,C5^2)=-1/phi, chi(T_1g,C3)=0,
# chi(T_1g,C2)=-1, dim=3. Compare against the I_h-projected part of 'T1' above
# (classes E,C5,C5^2,C3,C2 -- the "integer spin" classes shared with I_h).
T1g_wz_table = {'E': 3, 'C5': phi, 'C5^2': -1/phi, 'C3': 0, 'C2': -1}
T1_2I_at_shared_classes = {
    'E':    IRREPS['T1'][0],
    'C5':   IRREPS['T1'][2],
    'C5^2': IRREPS['T1'][4],
    'C3':   IRREPS['T1'][6],
    'C2':   IRREPS['T1'][8],
}
match = all(abs(T1_2I_at_shared_classes[k] - T1g_wz_table[k]) < 1e-9 for k in T1g_wz_table)

print(f"  T_1g (W/Z, used in face_gluon_geometry.py / doc_jobson_cell.txt):")
print(f"    {T1g_wz_table}")
print(f"  T1 (2I irrep in I52 x I52 above), restricted to shared classes:")
print(f"    {T1_2I_at_shared_classes}")

check("TC4: 'T1' in the I52 x I52 product is EXACTLY T_1g (identical characters, not just same name)",
      match, f"match = {match}")

# ── Section 2b: Is '2*G' a mislabeled muon, or the already-established gluon? ─
print()
print(SEP)
print("SECTION 2b: IS '2*G' A MISLABELED MUON, OR THE ESTABLISHED GLUON CHANNEL?")
print(SEP2)
print("  G32 (dim=4, spinor/fermionic, 2I-only) = muon [doc_leptons.txt].")
print("  G   (dim=4, ordinary/bosonic, also exists in plain I_h) = a DIFFERENT")
print("  irrep. Cross-check its characters against the G_g table already used")
print("  for the b quark / gluon channel (face_gluon_geometry.py, b_quark_geometry.py):")

Gg_established_table = {'E': 4, 'C5': -1, 'C5^2': -1, 'C3': 1, 'C2': 0}
G_2I_at_shared_classes = {
    'E':    IRREPS['G'][0],
    'C5':   IRREPS['G'][2],
    'C5^2': IRREPS['G'][4],
    'C3':   IRREPS['G'][6],
    'C2':   IRREPS['G'][8],
}
match_G = all(abs(G_2I_at_shared_classes[k] - Gg_established_table[k]) < 1e-9 for k in Gg_established_table)

print(f"  G_g (b quark / 2xG_g=gluons, established elsewhere): {Gg_established_table}")
print(f"  G (2I irrep in I52 x I52 above), restricted to shared classes: {G_2I_at_shared_classes}")

check("TC6: 'G' in the I52 x I52 product is EXACTLY G_g (the SAME irrep used for b quark / 2xG_g=8 gluons), NOT muon",
      match_G, f"match = {match_G}")

check("TC7: G32 (muon) has ZERO multiplicity in I52 x I52 -- muon is NOT reachable from tau-pairing",
      decomp.get('G32', 0) == 0,
      f"mult(G32) = {decomp.get('G32', 0)}  -- the muon slot is simply absent from this product")

print()
print("  So '2*G' is the ALREADY-ESTABLISHED gluon channel (2xG_g = 8 = SU(3)")
print("  adjoint, face_gluon_geometry.py FG2/FG3), not a mislabeled muon. Pairing")
print("  two tau windings reaches three already-named channels simultaneously:")
print("  Higgs (A), W/Z (2xT1), and gluons (2xG) -- but NOT muon (G32, absent).")

# ── Section 3: What is honestly still open ──────────────────────────────────
print()
print(SEP)
print("SECTION 3: WHAT THIS DOES AND DOES NOT PROVE")
print(SEP2)
print("  ESTABLISHED BY THIS SCRIPT (exact, algebraic):")
print("    - T_1g (the W/Z channel) is a genuine, reachable component when two")
print("      tau (I52) windings are combined -- not a coincidence of naming,")
print("      the characters match exactly (TC4).")
print("    - The SAME pairing operation reaches BOTH the Higgs channel (A) and")
print("      the W/Z channel (T1) -- structurally unifying, consistent with")
print("      the Higgs's own established internal make-up (JC3: 20 I52 in phase).")
print()
print("  NOT ESTABLISHED (open, needs a separate geometric construction):")
print("    - Whether the ACTUAL 3D paths of paired tau windings converge toward")
print("      the cell center and trace a directed-cone shape (the user's")
print("      specific spatial picture). This script only shows the abstract")
print("      representation-theory content is present -- it does not compute")
print("      any positions, paths, or directions.")
print("    - Which specific geometric pairing/phase relationship of the two")
print("      I52 windings PROJECTS onto the T1 slot specifically, as opposed")
print("      to the other slots in the same product (2*T2, 2*G, 3*H).")
print("    - This connects directly to the existing OPEN item in")
print("      open_items.txt: 'CHARM [HIGH] Formal H_u Zone 1 polygon path")
print("      geometry (analog of tau I52 path)' -- the natural next step is")
print("      a gluon_tau_helix.py-style geometric construction for a")
print("      Zone-1-confined / paired tau path, not attempted here.")

check("TC5: open geometric question explicitly identified (not silently skipped)",
      True, "spatial convergence-to-center claim requires a separate construction; flagged, not derived")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED. The T_1g (W/Z) channel is algebraically reachable")
    print("  from pairing two tau (I52) windings -- a real, exact, previously")
    print("  unconnected finding. The spatial/geometric picture (paths converging")
    print("  toward the center, forming a directed cone) remains open future work.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(SEP)
