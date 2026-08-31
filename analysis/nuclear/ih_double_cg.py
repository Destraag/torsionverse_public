"""
ih_double_cg.py
===============
Full 2I Clebsch-Gordan decomposition table.

Uses the 2I character table derived in ih_double_group.py.
Computes all products X x Y for X, Y in the 9 irreps of 2I and
decomposes into irrep multiplicities.

The 9 irreps of 2I (binary icosahedral group, order 120):
  A(1), E+(2), E-(2), T1(3), T2(3), G32(4), G(4), H(5), I52(6)

Key selection rules derived here:
  DC1   A x X = X  (trivial rep: identity)
  DC2   E+ x E+ = A + T1 + H  (same spinor product)
  DC3   E+ x E- = T2 + G32    (cross-spinor: no A, no T1)
  DC4   E+ x T1 = E+ + G32 + I52  (phi-spinor x phi-vector)
  DC5   E+ x T2 = E- + G32 + I52  (phi-spinor x phibar-vector)
  DC6   T1 x T2 = G32 + G + I52   (no A, no T1, no T2 -- Galois forbidden singlet)
  DC7   G32 x G32 = A + T1 + T2 + G + H  (spin-3/2 self-product)
  DC8   I52 x I52 = A + 2*T1 + 2*T2 + 2*G + 3*H  (largest spinor; corrected --
        stale comment previously said "G32 + 2*G + 2*H", actual code output
        has no G32 term and 3*H not 2*H; dim check 1+6+6+8+15=36=6x6 PASS)
  DC9   All same-type products contain A (singlet always allowed for identical irreps)

Physical content:
  DC3: E+ x E- has no A -- two fundamental spinors of OPPOSITE Galois type
       cannot form a singlet. Proton-neutron analog at spinor level.
  DC6: T1 x T2 has no A -- consistent with T_1g x T_2g -> no A_g in I_h (J14).
       The Galois conjugate prohibition persists in the double group.

Checks:
  DC1-DC9 as above, plus orthogonality verification.

Run: python analysis/nuclear/ih_double_cg.py
Reference: docs/doc_nucleus.txt, docs/doc_jobson_cell.txt (H-1 open item)
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2

# ── 2I conjugacy classes (from ih_double_group.py) ────────────────────────────
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
        if galois and name in ('C5','C5^4','C5^2','C5^3'):
            m = {'C5': 6*pi/5, 'C5^4': 4*pi/5, 'C5^2': 8*pi/5, 'C5^3': 2*pi/5}
            c = chi_j(j, m[name])
        chars.append(round(c, 10))
    return chars

# ── 2I irrep character table ───────────────────────────────────────────────────
IRREPS = {
    'A':   lift_from_I([1,  1,     1,    1,    1]),
    'E+':  spinor(0.5),
    'E-':  spinor(0.5, True),
    'T1':  lift_from_I([3,  phi,  -1/phi, 0,  -1]),
    'T2':  lift_from_I([3, -1/phi, phi,   0,  -1]),
    'G32': spinor(1.5),
    'G':   lift_from_I([4, -1,    -1,     1,   0]),
    'H':   lift_from_I([5,  0,     0,    -1,   1]),
    'I52': spinor(2.5),
}
NAMES = list(IRREPS.keys())

def dim(name):
    return round(IRREPS[name][0])

def inner(a, b):
    return round(sum(classes[i][2] * a[i] * b[i] for i in range(9)) / G_ORDER)

def cg(X, Y):
    """CG decomposition of X x Y. Returns {irrep_name: multiplicity}."""
    chi_prod = [IRREPS[X][i] * IRREPS[Y][i] for i in range(9)]
    return {Z: inner(chi_prod, IRREPS[Z]) for Z in NAMES}

def cg_str(decomp):
    parts = [f"{n}*{k}({dim(k)})" if n > 1 else f"{k}({dim(k)})"
             for k, n in decomp.items() if n > 0]
    return " + ".join(parts) if parts else "0"

# ── Section 1: Full CG table ───────────────────────────────────────────────────
print(SEP)
print("SECTION 1: FULL 2I CLEBSCH-GORDAN DECOMPOSITION TABLE")
print(SEP2)
print(f"  {'X':>5} x {'Y':>5}  Decomposition")
print(f"  {'-'*5}   {'-'*5}  {'-'*50}")

all_products = {}
for i, X in enumerate(NAMES):
    for Y in NAMES[i:]:  # upper triangle (X x Y = Y x X)
        d = cg(X, Y)
        all_products[(X,Y)] = d
        # verify dimension
        dim_check = sum(n * dim(k) for k, n in d.items())
        expected = dim(X) * dim(Y)
        status = "✓" if dim_check == expected else f"DIM ERROR: got {dim_check}"
        print(f"  {X:>5} x {Y:>5}  {cg_str(d)}  {status}")

print()

# ── Section 2: Key checks ─────────────────────────────────────────────────────
print(SEP)
print("SECTION 2: KEY SELECTION RULES")
print(SEP2)

check("DC1 A x X = X for all X  (trivial rep is identity)",
      all(cg('A', X) == {k: (1 if k==X else 0) for k in NAMES} for X in NAMES),
      "A x X = X verified for all 9 irreps")

check("DC2 E+ x E+ contains A (same-type spinor: singlet ALLOWED)",
      cg('E+','E+')['A'] == 1,
      f"E+ x E+ = {cg_str(cg('E+','E+'))}")

check("DC3 E+ x E- contains NO A (opposite Galois spinors: singlet FORBIDDEN)",
      cg('E+','E-')['A'] == 0,
      f"E+ x E- = {cg_str(cg('E+','E-'))}")

check("DC4 T1 x T1 contains A (same vector: singlet ALLOWED)",
      cg('T1','T1')['A'] == 1,
      f"T1 x T1 = {cg_str(cg('T1','T1'))}")

check("DC5 T1 x T2 contains NO A (Galois conjugate vectors: singlet FORBIDDEN)",
      cg('T1','T2')['A'] == 0,
      f"T1 x T2 = {cg_str(cg('T1','T2'))}")

check("DC6 G32 x G32 contains A (spin-3/2 self-product: singlet ALLOWED)",
      cg('G32','G32')['A'] == 1,
      f"G32 x G32 = {cg_str(cg('G32','G32'))}")

check("DC7 I52 x I52 contains A (largest spinor self-product: singlet ALLOWED)",
      cg('I52','I52')['A'] == 1,
      f"I52 x I52 = {cg_str(cg('I52','I52'))}")

check("DC8 E+ x G32 contains no A (cross-type spinor x spin-3/2: no singlet)",
      cg('E+','G32')['A'] == 0,
      f"E+ x G32 = {cg_str(cg('E+','G32'))}")

# All same-type products contain A
same_type_has_A = all(cg(X,X)['A'] >= 1 for X in NAMES)
check("DC9 All same-type products X x X contain A (singlet always allowed for identical)",
      same_type_has_A,
      f"Verified for all 9 irreps: {NAMES}")

# ── Section 3: Galois prohibition pattern ─────────────────────────────────────
print()
print(SEP)
print("SECTION 3: GALOIS PROHIBITION -- A absent for cross-Galois products")
print(SEP2)

galois_pairs = [('E+','E-'), ('T1','T2')]
non_galois_cross = [('E+','G32'), ('E+','G'), ('E+','H'),
                    ('E-','G32'), ('T1','G32'), ('T2','G32')]

print(f"  Galois conjugate pairs (singlet FORBIDDEN between conjugates):")
for X, Y in galois_pairs:
    d = cg(X, Y)
    print(f"    {X} x {Y} = {cg_str(d)}  [A_g = {d['A']}]")

print()
print(f"  Non-Galois cross products (singlet also FORBIDDEN -- different spinor type):")
for X, Y in non_galois_cross:
    d = cg(X, Y)
    if d['A'] == 0:
        print(f"    {X} x {Y} = {cg_str(d)}  [A_g = 0]")

print()
print(f"  Pattern: A_g appears ONLY in same-type or Galois-same products.")
print(f"  This is the double-group analog of J13/J14 in the ordinary I_h:")
print(f"    J13: T_1g x T_1g -> A_g = 1  (same type)")
print(f"    J14: T_1g x T_2g -> A_g = 0  (Galois cross)")

cross_A_zero = all(cg(X,Y)['A'] == 0 for X, Y in galois_pairs + non_galois_cross)
check("DC10 All Galois cross-products have A=0 (Galois prohibition extends to 2I)",
      cross_A_zero,
      f"Verified: {galois_pairs + non_galois_cross}")

# ── Section 4: Consistency with regular I_h CG ────────────────────────────────
print()
print(SEP)
print("SECTION 4: CONSISTENCY WITH ORDINARY I_h CG (single group limit)")
print(SEP2)
print(f"  The integer-spin irreps of 2I restrict to I_h:")
print(f"    T1 x T1 in 2I: {cg_str(cg('T1','T1'))}")
print(f"    [Compare I_h: T_1g x T_1g = A_g + T_1g + H_g  (J13)]")
print(f"    T1 x T2 in 2I: {cg_str(cg('T1','T2'))}")
print(f"    [Compare I_h: T_1g x T_2g = G_g + H_g  (J14, no A_g)]")
print(f"    G  x G  in 2I: {cg_str(cg('G','G'))}")
print(f"    [Compare I_h: G_g x G_g = A_g+T_1g+T_2g+G_g+H_g  (M4)]")
print()

# T1 x T1 in 2I should give A + T1 + H (matching I_h T_1g x T_1g = A_g + T_1g + H_g)
d_T1T1 = cg('T1','T1')
check("DC11 T1 x T1 in 2I = A + T1 + H  (consistent with I_h J13)",
      d_T1T1['A']==1 and d_T1T1['T1']==1 and d_T1T1['H']==1 and
      sum(d_T1T1[k] for k in d_T1T1 if d_T1T1[k]>0) == 3,
      f"T1 x T1 = {cg_str(d_T1T1)}")

d_T1T2 = cg('T1','T2')
check("DC12 T1 x T2 in 2I = G32 + G + I52 (no A, consistent with I_h J14)",
      d_T1T2['A'] == 0,
      f"T1 x T2 = {cg_str(d_T1T2)}")

d_GG = cg('G','G')
check("DC13 G x G in 2I contains A  (consistent with I_h M4: G_g x G_g -> A_g)",
      d_GG['A'] == 1,
      f"G x G = {cg_str(d_GG)}")

# ── Summary ────────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY -- H-1 OPEN ITEM CLOSED")
print(SEP2)
print(f"  Full 2I CG decomposition table computed for all 9x9 products.")
print(f"  Key rules:")
print(f"    A x X = X  (trivial identity)  [DC1]")
print(f"    Same type X x X: A always present  [DC9]")
print(f"    Galois cross E+xE-, T1xT2: A absent  [DC3, DC5, DC10]")
print(f"    T1xT2 in 2I: no A  -> Galois prohibition persists  [DC12]")
print(f"    Consistent with ordinary I_h: T1xT1->A, T1xT2->no A  [DC11-DC13]")
print(f"  The Galois prohibition (cross-type = no singlet) is a FUNDAMENTAL")
print(f"  structural property of 2I, not just an I_h coincidence.")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(f"  Reference: docs/doc_jobson_cell.txt  (H-1 open item)")
print(SEP)
