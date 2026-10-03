#!/usr/bin/env python3
"""
dipole_dipole_interaction_energy_check.py

DIRECT extension of analysis/nuclear/bar_magnet_pole_pressure_model.py:
that script established the SHAPE of a single two-pole pressure dipole's
own field/gradient (1/r^2 potential-like, 1/r^3 field-like). This script
asks the question actually needed by entanglement_phase_lock.py's Zone 3
inter-particle COUPLING ENERGY: when TWO such dipoles interact (not one
dipole's self-field observed from a point), does the resulting
INTERACTION ENERGY fall as 1/r^2 or 1/r^3?

MOTIVATION (notes/open_items/false_positive_scan_series1_resolved.txt
[DOC-ELECTRON-TIER3]; analysis/nuclear/bound_neutron_g_factor_distance_
dependence.py Section 8): two competing exponents were left as "neither
independently settled" -- n=3 (treating Zone 3 as a dipole FIELD, one
power steeper than its own potential) vs n=2 (bare reuse of an
intermediate internal CO-ROTATION VELOCITY falloff, v(r)~1/r^2, which is
a dimensionally different quantity -- m/s, not energy, with no stated
conversion). This script settles which is correct FOR THE ENERGY, by
computing it directly, rather than assuming either.

METHOD: build TWO independent two-pole pressure dipoles (dipole B fixed,
dipole A as a "test" dipole at center-to-center separation R), both using
the SAME already-established single-source Coulomb pressure formula
(doc_magnetism Section 1.2: P(r) = Q*K/(4*pi*r)) -- no new physics, no
fitted parameters. The interaction energy uses the SAME recipe already
used throughout this framework for electrostatic energy (U = charge *
potential-of-the-other-source, summed over every charge):

  U(R) = Q*P_B(z_A+) - Q*P_B(z_A-)

This is an EXACT closed form (no far-field approximation baked in) -- it
reduces algebraically to:

  U(R) = (Q^2*K)/(4*pi) * [ 2/R - 1/(R+ell) - 1/(R-ell) ]

Taylor-expanding this EXACT expression for R >> ell gives
U(R) ~ -(Q^2*ell^2*K)/(2*pi*R^3) + O(1/R^5) -- i.e. 1/R^3 is the
PREDICTED leading-order behavior. This script verifies that numerically
(not assumed) via log-log slope fitting, mirroring bar_magnet_pole_
pressure_model.py's own methodology, and additionally checks how close to
onset (R ~ ell, relevant to entanglement_phase_lock.py evaluating E_Z3
starting AT r=r_p itself, not deep in the far field) the asymptotic law
already holds.

Reference: analysis/nuclear/bar_magnet_pole_pressure_model.py;
  analysis/nuclear/bound_neutron_g_factor_distance_dependence.py Section 8;
  analysis/quantum/entanglement_phase_lock.py
"""
import math
import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP = "=" * 72
SEP2 = "-" * 72
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

pi = math.pi
eps_0 = 8.8541878128e-12
K = 1 / eps_0  # established bulk modulus = 1/eps_0 (doc_magnetism Section 1.1)

print(SEP)
print("DIPOLE-DIPOLE INTERACTION ENERGY -- DOES IT FALL AS 1/r^2 OR 1/r^3?")
print(SEP)

ell = 1.0       # each dipole's own internal pole separation (arbitrary units -- shape check)
Q_pole = 1.0    # pole strength (arbitrary units, same convention as bar_magnet_pole_pressure_model.py)

def pressure_at(point, Q, pos):
    """Single-source pressure, reusing the established Coulomb-pressure formula
    P(r) = Q*K/(4*pi*r) directly (doc_magnetism Section 1.2)."""
    r_vec = np.array(point) - np.array(pos)
    r = np.linalg.norm(r_vec)
    return Q * K / (4 * pi * r)

def dipole_B_pressure(point):
    """Dipole B, fixed: +Q at (0,0,+ell/2), -Q at (0,0,-ell/2)."""
    posN = np.array([0, 0, ell / 2])
    posS = np.array([0, 0, -ell / 2])
    return pressure_at(point, Q_pole, posN) + pressure_at(point, -Q_pole, posS)

def interaction_energy_axial(R):
    """Dipole A (same orientation as B) centered at separation R along the
    SAME axis as B -- U = Q*[P_B(A+) - P_B(A-)], the standard
    charge-times-potential-of-the-other-source recipe."""
    posA_plus  = np.array([0, 0, R + ell / 2])
    posA_minus = np.array([0, 0, R - ell / 2])
    return Q_pole * dipole_B_pressure(posA_plus) - Q_pole * dipole_B_pressure(posA_minus)

def interaction_energy_perp(R):
    """Dipole A (same orientation/moment direction as B -- internal
    separation still along z) with its CENTER offset from B's center along
    x by R (separation vector perpendicular to both dipole moments)."""
    posA_plus  = np.array([R, 0, ell / 2])
    posA_minus = np.array([R, 0, -ell / 2])
    return Q_pole * dipole_B_pressure(posA_plus) - Q_pole * dipole_B_pressure(posA_minus)

def exact_axial_formula(R):
    """Closed-form algebraic reduction of interaction_energy_axial(R), an
    independent cross-check of the direct numerical sum above."""
    return (Q_pole**2 * K / (4 * pi)) * (2 / R - 1 / (R + ell) - 1 / (R - ell))

def analytic_leading_order_axial(R):
    """Leading Taylor term for R >> ell (standard aligned-axial point-dipole
    formula, p = Q_pole*ell): U ~ -2*p^2*K/(4*pi*R^3)."""
    p = Q_pole * ell
    return -(2 * p**2 * K / (4 * pi)) / R**3

# =============================================================================
print()
print(SEP2)
print("PART 0: EXACT CLOSED-FORM REDUCTION MATCHES DIRECT 4-CHARGE SUM")
print(SEP2)
R_test = 50.0 * ell
direct = interaction_energy_axial(R_test)
closed_form = exact_axial_formula(R_test)
check("D0 direct numerical 4-charge sum matches the algebraic closed-form "
      "reduction (both are the SAME exact calculation, two independent code paths)",
      abs(direct - closed_form) < 1e-9 * abs(closed_form),
      f"direct={direct:.6e}, closed-form={closed_form:.6e}")

# =============================================================================
print()
print(SEP2)
print("PART 1: AXIAL (HEAD-TO-TAIL) INTERACTION ENERGY FALLOFF")
print(SEP2)
print("  Same far-field sampling convention as bar_magnet_pole_pressure_model.py")
print("  (R/ell = 10 to 320): fitting the log-log slope of |U(R)|.")

R_values = np.array([10.0, 20.0, 40.0, 80.0, 160.0, 320.0]) * ell
U_axial = np.array([interaction_energy_axial(R) for R in R_values])

log_R = np.log(R_values)
log_U_axial = np.log(np.abs(U_axial))
slope_axial, _ = np.polyfit(log_R, log_U_axial, 1)
print(f"  Fitted log-log slope (axial |U(R)| vs R): {slope_axial:.6f}")
print(f"  Expected: -3.0 (dipole-dipole interaction ENERGY, NOT -2.0)")
check("E1 axial dipole-dipole interaction energy falls as 1/R^3 (not 1/R^2), "
      "computed directly from the medium's own established single-source law",
      abs(slope_axial - (-3.0)) < 0.01, f"slope = {slope_axial:.4f}")

analytic_vals = np.array([analytic_leading_order_axial(R) for R in R_values])
rel_err_far = abs((U_axial[-1] - analytic_vals[-1]) / analytic_vals[-1])
check("E2 at the largest sampled separation (R=320*ell), the exact numerical "
      "value agrees with the standard textbook aligned-dipole formula "
      "(U ~ -2*p1*p2*K/(4*pi*R^3)) to within 0.1%",
      rel_err_far < 0.001, f"relative difference = {rel_err_far*100:.6f}%")

# =============================================================================
print()
print(SEP2)
print("PART 2: PERPENDICULAR (SIDE-BY-SIDE) INTERACTION ENERGY FALLOFF")
print(SEP2)
print("  Same parallel dipole orientation, but separation vector PERPENDICULAR")
print("  to the dipole axis -- checks the 1/R^3 power law is not an artifact")
print("  of the specific axial geometry.")

U_perp = np.array([interaction_energy_perp(R) for R in R_values])
log_U_perp = np.log(np.abs(U_perp))
slope_perp, _ = np.polyfit(log_R, log_U_perp, 1)
print(f"  Fitted log-log slope (perpendicular |U(R)| vs R): {slope_perp:.6f}")
check("E3 perpendicular-geometry dipole-dipole interaction energy ALSO falls "
      "as 1/R^3 (same power, different sign/coefficient than axial) -- the "
      "1/R^3 law is a property of the dipole-dipole interaction itself, not "
      "one specific orientation",
      abs(slope_perp - (-3.0)) < 0.01, f"slope = {slope_perp:.4f}")
check("E4 axial and perpendicular configurations have OPPOSITE sign (matching "
      "the standard dipole-dipole fact: aligned dipoles attract head-to-tail "
      "and repel side-by-side, or vice versa -- some definite sign flip, not "
      "the same sign)",
      (U_axial[-1] > 0) != (U_perp[-1] > 0),
      f"U_axial(320*ell)={U_axial[-1]:.4e}, U_perp(320*ell)={U_perp[-1]:.4e}")

# =============================================================================
print()
print(SEP2)
print("PART 3: HOW CLOSE TO r=r_p (R ~ ell, NOT FAR FIELD) DOES 1/R^3 HOLD?")
print(SEP2)
print("  entanglement_phase_lock.py evaluates E_Z3(r) STARTING AT r=r_p itself")
print("  (Zone 3 onset), not deep in the far field. If the dipole's own")
print("  internal scale ell is comparable to r_p, R/ell is not large at onset,")
print("  so Part 1's far-field fit is not automatically applicable there.")
print("  This reports where the LOCAL log-log slope actually settles near -3,")
print("  as an honest characterization rather than a hidden assumption.")

R_near = np.array([1.5, 2.0, 3.0, 5.0, 8.0, 15.0, 30.0, 60.0]) * ell
U_near = np.array([interaction_energy_axial(R) for R in R_near])
local_slopes = []
for i in range(len(R_near) - 1):
    s = (math.log(abs(U_near[i + 1])) - math.log(abs(U_near[i]))) / \
        (math.log(R_near[i + 1]) - math.log(R_near[i]))
    local_slopes.append(s)
    print(f"  R/ell in [{R_near[i]:5.1f}, {R_near[i+1]:5.1f}]: local slope = {s:.4f}")

check("E5 by R/ell=30-60 (still well short of Part 1's 10-320 far-field "
      "range), the local slope is already within 1% of the asymptotic -3.0",
      abs(local_slopes[-1] - (-3.0)) < 0.03, f"local slope at R/ell=30-60: {local_slopes[-1]:.4f}")

# =============================================================================
print()
print(SEP)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"RESULT: {len(results)} checks  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
    print()
    print("  CONCLUSION: the Zone 3 inter-particle COUPLING ENERGY exponent")
    print("  question is now DERIVED, not merely re-asserted. Computing the")
    print("  DIRECT interaction energy between two two-pole pressure dipoles --")
    print("  built from nothing but the already-established single-source")
    print("  Coulomb pressure law (doc_magnetism Section 1.2), no new postulates")
    print("  -- gives 1/R^3 in BOTH tested orientations, matches the standard")
    print("  textbook aligned-dipole coefficient to <0.001% in the far field,")
    print("  and the asymptotic power law is already accurate to ~1% by")
    print("  R~30-60*ell, not only in an extreme far-field limit. The n=2")
    print("  alternative (bare reuse of the v(r)~1/r^2 VELOCITY profile for an")
    print("  ENERGY) was never itself an energy calculation -- just a")
    print("  dimensionally mismatched reuse of a different quantity's exponent")
    print("  with no stated conversion. n=3 is the physically DERIVED choice")
    print("  for the Zone 3 coupling ENERGY specifically.")
else:
    print()
    print("  NOT ALL CHECKS PASSED -- see FAIL detail above; the exponent")
    print("  question should NOT be treated as settled by this script.")
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAIL: {name}")
            print(f"        {detail}")
print(SEP)
