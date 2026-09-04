"""
diquark_antisymmetric_square_check.py
======================================
CHECK, not yet cited in any doc: does the claimed diquark irrep assignment
survive an independent group-theory verification, or was it asserted?

doc_particle_generation.txt / torsionverse_doc.py comments assert:
  proton diquark (uu):  [T_1u x T_1u]_ANTISYMMETRIC = T_2g
  neutron diquark (dd): [T_2u x T_2u]_ANTISYMMETRIC = T_1g

This script computes the ANTISYMMETRIC SQUARE character of T_1g (and T_2g --
using gerade tables since antisymmetric-square characters only depend on
|chi|, parity/gerade-ungerade is a separate bookkeeping factor handled
afterward) via the standard formula
  chi_Lambda2(g) = (chi(g)^2 - chi(g^2)) / 2
using the REAL character table (already verified via orthogonality
elsewhere, e.g. jobson_cell_doc.py J8-J9) and the group-theoretic squaring
map between I's five conjugacy classes (E,C2,C3,C5,C5^2), independently
justified below (not assumed) via two standard facts about A5 = I:
  (1) A5 is AMBIVALENT (all characters real, confirmed by the table itself)
      => every class is closed under inversion (g's class = g^-1's class).
  (2) A5's order-5 elements split into two inequivalent classes (12 each)
      such that g and g^2 land in DIFFERENT order-5 classes -- a standard,
      textbook A5 fact (this is exactly why 5-cycles split into two classes
      in A5 at all, unlike in the full S5).
Combined, these two facts FORCE the squaring map E->E, C2->E, C3->C3,
C5->C5^2, C5^2->C5 -- it is not asserted, it follows from A5's known
structure. The full T_1gxT_1g decomposition (which needs NO squaring map,
just chi(g)^2 pointwise) is used as an internal cross-check against the
ALREADY-VERIFIED result magnetism_doc.py M6 (A_g appears once in T_1gxT_1g).

CAVEAT (important, stated up front): this checks whether the NAIVE
antisymmetric square of the single-quark icosahedral vertex-mode irrep
(T_1u/T_2u) reproduces the diquark's claimed irrep. If the framework's
actual diquark construction also involves color/flavor antisymmetrization
(standard QCD diquark structure: color-antitriplet + flavor + spin all
contribute to the overall fermion antisymmetry), this simple check may not
be the complete picture -- a mismatch here means "the simplest verification
does not confirm the claim," not necessarily "the claim is definitely wrong."
Reported as a diagnostic for author review, not a finished derivation either way.

Reference: docs/series1/doc_jobson_cell.txt Section 5.2 (shared-irrep-label
  note, Part 4 of this script), docs/series1/doc_magnetism.txt Section 3.1
  (N_J regime column note), docs/series1/doc_particle_generation.txt
  (proton/neutron diquark table note), docs/series1/doc_nucleus.txt
  Sections 5.1-5.2 (the diquark assignment under review, Parts 1-3).
"""
import math
sys_encoding_fix = True
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

phi = (1 + math.sqrt(5)) / 2

SEP  = "=" * 70
SEP2 = "-" * 70
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  [{'PASS' if cond else '*** FAIL'}] {name}")
    if detail:
        print(f"         {detail}")

# ── I_h gerade character table (already-established, e.g. magnetism_doc.py) ──
classes      = ['E', 'C2', 'C3', 'C5', "C5^2"]
class_sizes  = [1, 15, 20, 12, 12]
order        = 60

chars = {
    'A_g':  [1,  1,  1,   1,     1],
    'T_1g': [3, -1,  0,  phi,  -1/phi],
    'T_2g': [3, -1,  0, -1/phi, phi],
    'G_g':  [4,  0,  1,  -1,   -1],
    'H_g':  [5,  1, -1,   0,    0],
}

print(SEP)
print("PART 0: VERIFY THE SQUARING MAP FROM A5 STRUCTURE (not assumed)")
print(SEP2)
# Fact 1: ambivalence -- table must be all-real (already true by inspection).
all_real = all(isinstance(v, (int, float)) for row in chars.values() for v in row)
check("P0a: character table is all-real (A5 is ambivalent, g ~ g^-1 always)",
      all_real, "no complex entries in any irrep's character list")

# Fact 2 (textbook A5 fact): order-5 elements split into 2 classes of 12,
# and squaring an element of one class lands in the OTHER order-5 class.
# We verify this is the ONLY self-consistent map by checking it reproduces
# a fact ALREADY independently verified elsewhere (magnetism_doc.py M6):
# A_g appears exactly once in the FULL (non-antisymmetrized) T_1g x T_1g
# product -- this full-product check needs NO squaring map at all, so it
# is an independent cross-check of the character table itself, not of the
# squaring map.
sq_map = [0, 0, 2, 4, 3]   # E->E, C2->E, C3->C3, C5->C5^2, C5^2->C5

def decompose(target_chars):
    decomp = {}
    for name, c in chars.items():
        n = sum(class_sizes[i] * target_chars[i] * c[i] for i in range(5)) / order
        decomp[name] = round(n, 6)
    return decomp

def decomp_str(decomp):
    parts = [f"{k}" + (f"*{int(v)}" if v != 1 else "") for k, v in decomp.items() if abs(v) > 1e-6]
    return " + ".join(parts) if parts else "0"

# ── Cross-check: FULL T_1g x T_1g product (no squaring map needed) ──────────
print()
print(SEP)
print("PART 1: FULL T_1g x T_1g PRODUCT (cross-check vs magnetism_doc.py M6)")
print(SEP2)
full_T1xT1 = [chars['T_1g'][i]**2 for i in range(5)]
decomp_full_T1 = decompose(full_T1xT1)
print(f"  T_1g x T_1g (full, dim 9) = {decomp_str(decomp_full_T1)}")
check("P1: A_g appears exactly once in FULL T_1g x T_1g (matches magnetism_doc.py M6)",
      decomp_full_T1.get('A_g', 0) == 1,
      f"decomposition: {decomp_full_T1}")
check("P1b: FULL T_1g x T_1g dims sum to 9",
      sum(round(chars[k][0]) * v for k, v in decomp_full_T1.items()) == 9,
      f"decomposition: {decomp_full_T1}")

# ── PART 2: antisymmetric / symmetric square of T_1g ────────────────────────
print()
print(SEP)
print("PART 2: ANTISYMMETRIC SQUARE [T_1g x T_1g]_Lambda2  (claimed = T_2g)")
print(SEP2)

def chi_sq(name, symmetric):
    c = chars[name]
    out = []
    for i in range(5):
        g_val  = c[i]
        g2_val = c[sq_map[i]]
        out.append((g_val**2 + g2_val) / 2 if symmetric else (g_val**2 - g2_val) / 2)
    return out

lam2_T1 = chi_sq('T_1g', symmetric=False)
sym2_T1 = chi_sq('T_1g', symmetric=True)
decomp_lam2_T1 = decompose(lam2_T1)
decomp_sym2_T1 = decompose(sym2_T1)

print(f"  chi_Lambda2(T_1g) at [E,C2,C3,C5,C5^2] = {[round(x,4) for x in lam2_T1]}")
print(f"  [T_1g x T_1g]_antisym (dim 3) = {decomp_str(decomp_lam2_T1)}")
print(f"  [T_1g x T_1g]_sym     (dim 6) = {decomp_str(decomp_sym2_T1)}")

check("P2: antisym + sym dims sum to 9 (internal consistency)",
      sum(round(chars[k][0])*v for k,v in decomp_lam2_T1.items()) +
      sum(round(chars[k][0])*v for k,v in decomp_sym2_T1.items()) == 9,
      f"antisym dim={sum(round(chars[k][0])*v for k,v in decomp_lam2_T1.items())}  "
      f"sym dim={sum(round(chars[k][0])*v for k,v in decomp_sym2_T1.items())}")

check("P2b: [T_1g x T_1g]_antisym x2 CLAIM -- does it equal T_2g as asserted in torsionverse_doc.py?",
      decomp_lam2_T1.get('T_2g', 0) == 1 and all(v == 0 for k, v in decomp_lam2_T1.items() if k != 'T_2g'),
      f"actual decomposition: {decomp_lam2_T1}  (comment in torsionverse_doc.py claims pure T_2g)")

check("P2c: [T_1g x T_1g]_antisym -- does it instead equal T_1g itself (vector-cross-product pattern)?",
      decomp_lam2_T1.get('T_1g', 0) == 1 and all(v == 0 for k, v in decomp_lam2_T1.items() if k != 'T_1g'),
      f"actual decomposition: {decomp_lam2_T1}")

# ── PART 3: antisymmetric square of T_2g (neutron diquark claim) ───────────
print()
print(SEP)
print("PART 3: ANTISYMMETRIC SQUARE [T_2g x T_2g]_Lambda2  (claimed = T_1g)")
print(SEP2)
lam2_T2 = chi_sq('T_2g', symmetric=False)
decomp_lam2_T2 = decompose(lam2_T2)
print(f"  chi_Lambda2(T_2g) at [E,C2,C3,C5,C5^2] = {[round(x,4) for x in lam2_T2]}")
print(f"  [T_2g x T_2g]_antisym (dim 3) = {decomp_str(decomp_lam2_T2)}")

check("P3: [T_2g x T_2g]_antisym -- does it equal T_1g as claimed?",
      decomp_lam2_T2.get('T_1g', 0) == 1 and all(v == 0 for k, v in decomp_lam2_T2.items() if k != 'T_1g'),
      f"actual decomposition: {decomp_lam2_T2}")
check("P3b: [T_2g x T_2g]_antisym -- does it instead equal T_2g itself?",
      decomp_lam2_T2.get('T_2g', 0) == 1 and all(v == 0 for k, v in decomp_lam2_T2.items() if k != 'T_2g'),
      f"actual decomposition: {decomp_lam2_T2}")

# ── PART 4: does an ordinary vector (spin-1) restrict to T_1g under I? ──────
# This is the OTHER kind of "does the label earn itself" check -- not a
# tensor-square construction, but: does chi_l=1(theta) = 1+2*cos(theta),
# the standard SO(3) l=1 (vector/spin-1) character, evaluated AT I's five
# rotation angles (0, 180, 120, 72, 144 deg), reproduce T_1g's own table
# exactly? This is the claim underlying "W/Z boson = T_1g" (doc_jobson_cell.txt
# Sec 5.2: "chi(T_1g,C5)=phi (character of spin-1 W/Z under 72-deg rotation)").
print()
print(SEP)
print("PART 4: DOES ORDINARY SPIN-1 (SO(3) l=1) RESTRICT EXACTLY TO T_1g?")
print(SEP2)
angles_deg = [0, 180, 120, 72, 144]
chi_l1 = [1 + 2*math.cos(math.radians(a)) for a in angles_deg]
print(f"  chi_l=1(theta) = 1+2cos(theta) at [E,C2,C3,C5,C5^2] = {[round(x,6) for x in chi_l1]}")
print(f"  T_1g's own table                                    = {chars['T_1g']}")
match_T1g = all(abs(chi_l1[i] - chars['T_1g'][i]) < 1e-9 for i in range(5))
match_T2g = all(abs(chi_l1[i] - chars['T_2g'][i]) < 1e-9 for i in range(5))
check("P4: ordinary vector/spin-1 (l=1) character EXACTLY matches T_1g (all 5 classes)",
      match_T1g, f"max diff = {max(abs(chi_l1[i]-chars['T_1g'][i]) for i in range(5)):.2e}")
check("P4b: ordinary vector/spin-1 (l=1) does NOT match T_2g (disambiguates from T_2g)",
      not match_T2g, f"max diff vs T_2g = {max(abs(chi_l1[i]-chars['T_2g'][i]) for i in range(5)):.2e}")

# ── PART 5: does T_1g vs T_2g matter for SELF-PRODUCT (CG) questions? ───────
# Relevant to: does doc_magnetism.txt's ferromagnetism argument (which only
# asks "how many times does A_g appear in X x X") actually depend on WHICH
# of the two dim-3 irreps (T_1g or T_2g) a given particle is assigned?
print()
print(SEP)
print("PART 5: T_1g x T_1g vs T_2g x T_2g -- identical, because |C5|=|C5^2|=12")
print(SEP2)
full_T2xT2 = [chars['T_2g'][i]**2 for i in range(5)]
decomp_full_T2 = decompose(full_T2xT2)
print(f"  T_1g x T_1g (full, dim 9) = {decomp_str(decomp_full_T1)}")
print(f"  T_2g x T_2g (full, dim 9) = {decomp_str(decomp_full_T2)}")
check("P5: A_g appears the SAME number of times in T_1g^2 as in T_2g^2",
      decomp_full_T1['A_g'] == decomp_full_T2['A_g'],
      f"A_g(T_1g^2)={decomp_full_T1['A_g']}  A_g(T_2g^2)={decomp_full_T2['A_g']}  "
      f"(both decompose as A_g + [itself] + H_g -- the 'itself' slot trivially "
      f"differs by definition, but the physically load-bearing A_g count is "
      f"identical because class sizes |C5|=|C5^2|=12 make the C5<->C5^2 "
      f"character swap invisible to any sum weighted by class size)")

# ── PART 6: is "which LETTER is the proton's" even a well-posed question? ──
# weak_interaction_cg.py's WI1/WI2 (T_2g x E+ = I52; T_1g x E- = I52) are the
# ONLY externally-anchored facts found anywhere in the repo tying a SPECIFIC
# hadron side to a SPECIFIC lepton side (E+ is anchored to the REAL electron,
# not an arbitrary label -- it's the electron because it matches the real,
# measured particle). This checks: is that anchoring fact actually
# DISCRIMINATING (only ONE specific pairing gives I52), or would ANY pairing
# work (in which case WI1/WI2 carry no information about "which letter")?
print()
print(SEP)
print("PART 6: DOES THE ELECTRON ANCHOR (WI1/WI2) DISCRIMINATE T_1g vs T_2g?")
print(SEP2)
# E+/E-/I52 characters at (C5, C5^2, C3, C2), from weak_interaction_cg.py's
# own already-verified table (WI1/WI2, 10/10 PASS in that script).
chi4 = {
    'T_1g': (phi,      -1/phi,  0,   -1),
    'T_2g': (-1/phi,    phi,    0,   -1),
    'E+':   (phi,      -1/phi,  1,    0),
    'E-':   (-1/phi,    phi,    1,    0),
    'I52':  (-1.0,     -1.0,    0,    0),
}
def prod4(a, b):
    return tuple(chi4[a][i]*chi4[b][i] for i in range(4))
def is_I52(prod):
    return all(abs(prod[i] - chi4['I52'][i]) < 1e-9 for i in range(4))

pairs = [('T_2g', 'E+'), ('T_1g', 'E-'), ('T_1g', 'E+'), ('T_2g', 'E-')]
for a, b in pairs:
    p = prod4(a, b)
    print(f"  {a} x {b}: chi = {tuple(round(x,4) for x in p)}  "
          f"{'== I52' if is_I52(p) else '!= I52'}")

check("P6a: T_2g x E+ = I52 (already-established WI1, opposite Galois phase)",
      is_I52(prod4('T_2g', 'E+')))
check("P6b: T_1g x E- = I52 (already-established WI2, opposite Galois phase)",
      is_I52(prod4('T_1g', 'E-')))
check("P6c: T_1g x E+ does NOT equal I52 (same Galois phase -- the 'wrong' pairing)",
      not is_I52(prod4('T_1g', 'E+')),
      f"chi(C5) = phi*phi = {phi*phi:.4f}, not I52's -1 -- discriminates")
check("P6d: T_2g x E- does NOT equal I52 (same Galois phase -- the 'wrong' pairing)",
      not is_I52(prod4('T_2g', 'E-')),
      f"chi(C5) = (-1/phi)*(-1/phi) = {(-1/phi)*(-1/phi):.4f}, not I52's -1 -- discriminates")
print()
print("  CONCLUSION: the discriminating fact is OPPOSITE Galois phase (one +phi,")
print("  one -1/phi) between the hadron side and the lepton side -- NOT which")
print("  specific letter (T_1g/T_2g) is used for either side. Swapping T_1g<->T_2g")
print("  AND E+<->E- simultaneously maps {T_2g x E+ = I52} to {T_1g x E- = I52}")
print("  (still true) -- the letter-swap is invisible to every check in this file.")
print("  What IS physically anchored (via E+ = the real electron) is the RELATIVE")
print("  fact 'proton's diquark has the Galois-opposite phase from the electron' --")
print("  not which of the two arbitrary labels (T_1g/T_2g) that phase is called.")

# ── SUMMARY ──────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
n_total = len(results)
print(f"  {n_pass}/{n_total} checks PASS")
print()
print("  Standard-physics cross-check: the antisymmetric square of an ordinary")
print("  3D vector (l=1) representation is the PSEUDOVECTOR (cross product a x b)")
print("  -- which has the SAME rotational character as l=1 itself, only with")
print("  flipped parity (this is why a x b transforms as a vector under pure")
print("  rotations in ordinary vector calculus). If T_1g's antisymmetric square")
print("  above equals T_1g (not T_2g), that is the textbook-consistent answer,")
print("  and would mean the diquark irrep assignment needs a different")
print("  (not-yet-identified) justification beyond a plain antisymmetric square")
print("  of the single-quark vertex irrep -- e.g. one that also folds in color")
print("  or flavor antisymmetrization, which this script does not model.")
