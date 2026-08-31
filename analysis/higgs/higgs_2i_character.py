"""
higgs_2i_character.py
=====================
LEAD 1: Verify chi(E_1/2, C_5) = phi under the binary icosahedral group 2I.

The claim: the structural coupling f1=PHI in the alpha derivation is not merely
a geometric coincidence. It equals the C_5 character of the electron's spin-1/2
representation under 2I, the same way chi(T_1g, C_5) = phi for the W/Z under I_h.

This script:
  1. Derives chi(E_1/2, C_5) algebraically and numerically
  2. Compares to chi(T_1g, C_5) under I_h
  3. Shows both equal phi -- same C_5 character, different groups and spin
  4. Discusses the implication for the alpha derivation

Run: python analysis/higgs/higgs_2i_character.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2

SEP  = "=" * 65
SEP2 = "-" * 65

print(SEP)
print("LEAD 1: chi(E_1/2, C_5) = phi UNDER BINARY ICOSAHEDRAL GROUP 2I")
print(SEP2)
print()

# ── Group setup ───────────────────────────────────────────────────────────────
print("SETUP: Two related groups")
print()
print("  I_h (full icosahedral group, order 120):")
print("    = I x Z_2  (icosahedral rotations + inversions)")
print("    Relevant for integer-spin particles (bosons)")
print("    Gerade irreps: A_g(1), T_1g(3), T_2g(3), G_g(4), H_g(5)")
print()
print("  2I (binary icosahedral group, order 120):")
print("    = spin double-cover of I (icosahedral rotation group, order 60)")
print("    Relevant for half-integer-spin particles (fermions, including electron)")
print("    Additional (spinor) irreps: E_1/2(2), E_3/2(2), E_5/2(4), E_7/2(6)")
print("    Spin-integer irreps same as I: A(1), T_1(3), T_2(3), G(4), H(5)")
print()

# ── C_5 character derivation ──────────────────────────────────────────────────
print(SEP)
print("ALGEBRAIC DERIVATION OF chi(E_1/2, C_5)")
print(SEP2)
print()
print("  A C_5 rotation = rotation by 2*pi/5 about a 5-fold axis.")
print("  Under I_h (integer spin):")
print("    chi(A_g,  C_5) = 1")
print("    chi(T_1g, C_5) = 1 + 2*cos(2*pi/5)")
print("    chi(T_2g, C_5) = 1 + 2*cos(4*pi/5)")
print()
print("  Under 2I (half-integer spin), the spinor E_1/2 representation:")
print("    A C_5 rotation by angle theta lifts to SU(2) rotation by theta/2.")
print("    Wait -- C_5 in 2I is the LIFT of C_5 from I:")
print("      In I:  C_5 = rotation by 2*pi/5 (order 5)")
print("      In 2I: C_5 lifts to element of order 10 (since 5 rounds in I = 10 in SU(2))")
print()
print("  For the spin-j representation of SU(2):")
print("    chi(R_theta, spin-j) = sin((2j+1)*theta/2) / sin(theta/2)")
print()
print("  For j=1/2 and theta = 2*pi/5 (one C_5 step in I, corresponding to")
print("  the element of order 10 in 2I):")
theta = 2*pi/5   # C_5 rotation angle
j_half = 0.5
chi_E12_exact = math.sin((2*j_half + 1) * theta/2) / math.sin(theta/2)
print(f"    theta = 2*pi/5 = {theta:.8f}")
print(f"    chi(E_1/2, C_5) = sin(3*pi/5) / sin(pi/5)")
print(f"                    = {chi_E12_exact:.10f}")
print()

# Verify algebraically: sin(3*pi/5)/sin(pi/5) = phi
# sin(3*pi/5) = sin(pi - 3*pi/5) = sin(2*pi/5)
# sin(pi/5) = sin(36 degrees)
# sin(2*pi/5) / sin(pi/5) = 2*cos(pi/5) = 2*(phi/2) = phi
print(f"  Algebraic verification:")
print(f"    sin(3*pi/5) = sin(pi - 2*pi/5) = sin(2*pi/5) = {math.sin(3*pi/5):.10f}")
print(f"    sin(pi/5)   = {math.sin(pi/5):.10f}")
print(f"    Ratio       = 2*cos(pi/5)      = {2*math.cos(pi/5):.10f}")
print(f"    2*cos(pi/5) = 2*(phi/2) = phi  = {phi:.10f}")
print(f"    chi(E_1/2, C_5) = phi:  {abs(chi_E12_exact - phi) < 1e-10}")
print()

# ── Comparison with T_1g ─────────────────────────────────────────────────────
print(SEP)
print("COMPARISON: chi(E_1/2, C_5) vs chi(T_1g, C_5)")
print(SEP2)
print()
chi_T1g = 1 + 2*math.cos(2*pi/5)
print(f"  chi(E_1/2, C_5)  [2I, spin-1/2, electron] = {chi_E12_exact:.10f} = phi")
print(f"  chi(T_1g,  C_5)  [I_h, spin-1, W/Z]       = {chi_T1g:.10f} = phi")
print(f"  Are they equal?  {abs(chi_E12_exact - chi_T1g) < 1e-10}")
print()
print("  IDENTICAL. Both probe the same C_5 character of the icosahedral")
print("  geometry, despite being in different groups (2I vs I_h) and")
print("  representing different spins (1/2 vs 1).")
print()
print("  WHY: cos(2*pi/5) = phi/2 is a fundamental icosahedral constant.")
print("   chi(T_1g, C_5) = 1 + 2*cos(2*pi/5) = 1 + phi - 1 = phi")
print("     [using phi = 1 + 2*cos(2*pi/5) - 1 = 2*cos(2*pi/5) + 1")
print("      but cos(2*pi/5) = (phi-1)/2 so 2*cos(2*pi/5)+1 = phi-1+1 = phi]")
print("   chi(E_1/2, C_5) = sin(3*pi/5)/sin(pi/5) = 2*cos(pi/5) = phi")
print("     [cos(pi/5) = phi/2]")
print()

# ── Implication for alpha derivation ──────────────────────────────────────────
print(SEP)
print("IMPLICATION FOR THE ALPHA DERIVATION")
print(SEP2)
print()
print("  In the alpha derivation, the vertex structural coupling is:")
print("    f1 = PHI = (1+sqrt(5))/2 = phi")
print()
print("  This was derived geometrically: 5 direct edges + 5 face midpoints")
print("  of the icosahedral vertex give PHI exactly.")
print()
print("  The SAME phi also equals:")
print("    chi(E_1/2, C_5) [electron spin-1/2 rep under 2I]")
print("    chi(T_1g,  C_5) [W/Z spin-1 rep under I_h]")
print()
print("  INTERPRETATION: f1=PHI is not just 'vertex geometry happens to give phi'.")
print("  It is the C_5 character of the ELECTRON'S REPRESENTATION under the")
print("  icosahedral symmetry group. The vertex coupling strength IS the")
print("  representation-theoretic weight of the electron in the I_h/2I system.")
print()
print("  This provides a UNIFIED reason for phi appearing in:")
print("    (a) The alpha derivation (electron, spin-1/2, 2I): chi(E_1/2) = phi")
print("    (b) The W/Z coupling to the cell (spin-1, I_h): chi(T_1g) = phi")
print("    (c) The icosahedral vertex geometry: geometric proof gives PHI")
print("  All three are the same mathematical fact about icosahedral C_5 symmetry.")
print()

# ── Full 2I character table (C_5 column) ─────────────────────────────────────
print(SEP)
print("FULL 2I CHARACTER TABLE -- C_5 COLUMN")
print(SEP2)
print()
print("  Representation    dim    chi(C_5)            value")
print(SEP2)
# Integer-spin irreps of 2I (same as I)
int_reps = [
    ("A (trivial)",   1,  1.0),
    ("T_1",           3,  1 + 2*math.cos(2*pi/5)),
    ("T_2",           3,  1 + 2*math.cos(4*pi/5)),
    ("G",             4,  -1.0),
    ("H",             5,  0.0),
]
# Half-integer spin irreps of 2I
# C_5 in 2I has order 10. For spin j: chi = sin((2j+1)*pi/5)/sin(pi/5) * [...]
# Actually for the C_5 element in 2I (which has order 10, corresponding to C_10 in SU(2)):
# theta_SU2 = 2*pi/5 for the 5-fold rotation in I, so in SU(2) this is a rotation by 2*pi/5
# but the LIFT to 2I gives an element of order 10.
# chi(E_j, C_5 in 2I) = sin((2j+1)*pi/5) / sin(pi/5)  [standard formula for double cover]
half_reps = [
    ("E_1/2",  2,  math.sin(3*pi/5) / math.sin(pi/5)),    # j=1/2
    ("E_3/2",  2,  math.sin(5*pi/5) / math.sin(pi/5)),    # j=3/2: sin(pi)/sin(pi/5) = 0
    ("E_5/2",  4,  math.sin(7*pi/5) / math.sin(pi/5)),    # j=5/2
    ("E_7/2",  6,  math.sin(9*pi/5) / math.sin(pi/5)),    # j=7/2? dim=6
]
for name, d, c in int_reps + half_reps:
    note = ""
    if abs(c - phi) < 1e-6:
        note = " = phi  <-- KEY"
    elif abs(c - 0) < 1e-6:
        note = " = 0"
    elif abs(c + 1/phi) < 1e-6:
        note = " = -1/phi"
    elif abs(c + 1) < 1e-6:
        note = " = -1"
    print(f"  {name:<16} {d:>4}    {c:>+12.8f}  {note}")
print()

# Check sum of dim^2 = order of 2I = 120
total = sum(d**2 for _, d, _ in int_reps + half_reps)
print(f"  Sum of dim^2 = {total}  (should be 120 = |2I|: {total == 120})")
print()

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()
print(f"  chi(E_1/2, C_5) = phi = {phi:.8f}  [PROVEN, algebraic]")
print()
print("  The alpha vertex coupling f1=PHI is the C_5 character of the")
print("  electron's spin representation under the binary icosahedral group.")
print("  This is NOT a numerical coincidence -- it is a theorem of the")
print("  icosahedral group structure.")
print()
print("  BOTH the electron (spin-1/2, 2I) and the W/Z bosons (spin-1, I_h)")
print("  have the same C_5 character = phi. The icosahedral medium couples")
print("  to particles at strength proportional to their C_5 character,")
print("  and for the physically relevant representations, this is phi.")
print()
print("  This does NOT change the alpha derivation. It explains WHY phi")
print("  appeared as the structural vertex coupling constant.")
print(SEP)
