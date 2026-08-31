"""
hopf_c4_phi_hit.py — Follow-up on the 4*pi^2/phi near-hit from hopf_c4.py.

hopf_c4.py found Candidate G: 4*pi^2/phi = 24.399, only +0.06% from the
target R_s/alpha = 24.384. This script investigates whether this is a
genuine geometric relationship by:

  1. Computing the implied alpha formula if gap = 4*pi^2/phi exactly
  2. Comparing the implied alpha to the CODATA measured value
  3. Simplifying the expression algebraically
  4. Assessing whether 4*pi^2/phi has a natural geometric origin in the
     Hopf fibration / icosahedral torus topology

Run: python analysis/hopf_c4_phi_hit.py
"""

import math

SEP = "=" * 62

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
Rs    = math.sqrt(5) / (4 * pi)
alpha = 7.2973525693e-3          # CODATA 2018

gap_exact   = Rs / alpha
gap_formula = 4 * pi**2 / phi

print(SEP)
print("HOPF C4: 4*pi^2/phi NEAR-HIT INVESTIGATION")
print(SEP)
print()
print(f"Target:       R_s / alpha        = {gap_exact:.8f}")
print(f"Candidate G:  4*pi^2 / phi       = {gap_formula:.8f}")
print(f"Error:        {(gap_formula - gap_exact) / gap_exact * 100:+.4f}%")
print()

# ── Implied alpha if the formula were exact ───────────────────
#
#   gap = 4*pi^2 / phi
#   alpha = R_s / gap = [sqrt(5)/(4*pi)] / [4*pi^2/phi]
#         = sqrt(5) * phi / (16 * pi^3)
#
# Simplify sqrt(5)*phi:
#   phi = (1 + sqrt(5)) / 2
#   sqrt(5)*phi = sqrt(5)*(1+sqrt(5))/2 = (sqrt(5) + 5) / 2

alpha_implied = math.sqrt(5) * phi / (16 * pi**3)
alpha_simplified = (math.sqrt(5) + 5) / (32 * pi**3)

print("IMPLIED ALPHA IF gap = 4*pi^2/phi EXACTLY:")
print()
print("  alpha = R_s / (4*pi^2/phi)")
print("        = [sqrt(5)/(4*pi)] * [phi/(4*pi^2)]")
print("        = sqrt(5) * phi / (16 * pi^3)")
print()
print("  Simplify sqrt(5)*phi = sqrt(5)*(1+sqrt(5))/2 = (sqrt(5)+5)/2:")
print("        = (sqrt(5)+5) / (32 * pi^3)")
print()
print(f"  alpha_implied   = {alpha_implied:.10e}")
print(f"  alpha_CODATA    = {alpha:.10e}")
print(f"  Error:          = {(alpha_implied - alpha)/alpha * 100:+.6f}%")
print()

# ── Geometric origin of 4*pi^2/phi ───────────────────────────
print("GEOMETRIC ORIGIN OF 4*pi^2 / phi:")
print()
print("  4*pi^2 = L_major/R1 (the ratio of major circumference to minor")
print("           radius of a torus with R2/R1 = 2*pi)")
print("         = the primary Hopf torus length ratio")
print()
print("  phi    = the golden ratio = (1+sqrt(5))/2")
print("         = emerges from icosahedral symmetry group (Y, order 120)")
print("         = the unique positive solution to phi^2 = phi + 1")
print()
print("  4*pi^2/phi: the Hopf torus length ratio scaled by the icosahedral")
print("  symmetry factor. Physical interpretation:")
print("    - 4*pi^2 sets the crossing-ring loop topology (Hopf)")
print("    - 1/phi scales by the golden-ratio compression factor of the")
print("      icosahedral medium symmetry")
print()
print("  WHY phi might appear here:")
print("  The icosahedron is the Platonic solid with the highest rotational")
print("  symmetry (60 proper rotations). If the torsion medium has")
print("  icosahedral topology, phi enters through the diagonal/edge ratio")
print("  of the icosahedron: d/a = phi (the fundamental icosahedral ratio).")
print("  Dividing 4*pi^2 by phi would then represent a projection of the")
print("  Hopf torus geometry onto the icosahedral medium symmetry.")
print()

# ── Numerical sanity checks ───────────────────────────────────
print("SANITY CHECKS:")
print()
print(f"  phi^2         = {phi**2:.6f}  (= phi + 1 = {phi+1:.6f}, check: {abs(phi**2-(phi+1))<1e-10})")
print(f"  sqrt(5)       = {math.sqrt(5):.6f}")
print(f"  sqrt(5)*phi   = {math.sqrt(5)*phi:.6f}  = (sqrt(5)+5)/2 = {(math.sqrt(5)+5)/2:.6f}")
print(f"  16*pi^3       = {16*pi**3:.6f}")
print(f"  32*pi^3       = {32*pi**3:.6f}")
print()

# ── Assessment ───────────────────────────────────────────────
print(SEP)
print("ASSESSMENT:")
print()
print(f"  The formula alpha = sqrt(5)*phi / (16*pi^3) reproduces the")
print(f"  measured fine structure constant to {abs((alpha_implied-alpha)/alpha)*100:.4f}%.")
print()
if abs((alpha_implied - alpha) / alpha) < 0.001:
    print("  This is a <0.1% match. For a formula built from three")
    print("  fundamental constants (sqrt(5), phi, pi) with no free")
    print("  parameters, this is worth recording as a serious conjecture.")
    print()
    print("  HOWEVER: the 0.06% gap means the formula is NOT exact.")
    print("  Either:")
    print("    (a) The true formula has small higher-order corrections,")
    print("        analogous to QED loop corrections to g-2")
    print("    (b) The near-hit is numerical coincidence -- there are")
    print("        enough combinations of {sqrt(5), phi, pi} that one")
    print("        landing near 24.38 is not improbable")
    print("    (c) The correct formula uses a slightly different geometric")
    print("        quantity (e.g., 4*pi^2 / phi is a proxy for the true")
    print("        Hopf fiber ratio in the icosahedral medium)")
    print()
    print("  VERDICT: Promoted from 'no candidate' to 'near-hit conjecture'.")
    print("  Document as Conjecture C4a: alpha ~= sqrt(5)*phi/(16*pi^3)")
    print("  pending a geometric derivation that closes the 0.06% gap.")
else:
    print("  ERROR: greater than 0.1% -- reassess.")

print(SEP)
