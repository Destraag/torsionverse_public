#!/usr/bin/env python3
"""
bar_magnet_pole_pressure_model.py

FIRST piece of the population-imbalance/flow-geometry model for a line
(bar) magnet's N/S mechanism (mag_jamming_mechanism_investigation.txt).

SCOPE OF THIS SCRIPT (deliberately narrow): build the SIMPLEST possible
version of the mechanism -- two point pressure sources ("poles"), +Q_pole
at the N end and -Q_pole at the S end of a line magnet, reusing the
ALREADY-ESTABLISHED single-charge pressure formula directly (doc_magnetism
Section 1.2: P(r) = Q*K/(4*pi*r) = Q/(4*pi*eps_0*r), the same formula
already used for the Coulomb potential) -- and check whether this reduces,
in the far-field, to the CORRECT power-law behavior for a magnetic dipole
BEFORE attempting to derive Q_pole from any electron-population physics.

WHY THE POLE MODEL, NOT YET THE CURRENT-LOOP MODEL: mainstream E&M has TWO
equivalent external-field models for a bar magnet (Wikipedia, "Magnetic
dipole"): the pole model (two fictitious opposite magnetic charges) and
the Amperian current-loop model (real circulating current). Both give the
SAME external field at large distance. The pole model is the simpler one
to build first with already-established machinery (it is LITERALLY the
same math as two opposite electric charges, already fully derived here
via the Coulomb pressure Green's function) -- a natural first check before
attempting the harder, genuinely new current/population-imbalance version.

CRITICAL DIMENSIONAL CARE (the reason for PART 2 below): pressure here
plays the role of ELECTRIC POTENTIAL (Section 1.2's own "voltage IS
pressure"), not of the E or B FIELD. A potential-like dipole (two opposite
point sources a fixed distance apart) falls off as 1/r^2 in the far field,
NOT 1/r^3 -- the standard B-field 1/r^3 result (Wikipedia, "Magnetic
dipole") applies to a FIELD quantity (a derivative of potential), one
power of r steeper. This script checks BOTH: the raw pressure dipole
(expect 1/r^2) and its gradient (expect 1/r^3, the real target for
comparison to B) -- getting this distinction right, not conflating
pressure with field, is the whole point of Part 2.

NOT YET ATTEMPTED (honestly out of scope here): deriving Q_pole itself
from electron population/flux at the magnet's ends. This script only
checks that the SIMPLEST possible two-source pressure model has the
right mathematical shape to be a candidate at all, before investing in
that harder derivation.

Reference: notes/open_items/mag_jamming_mechanism_investigation.txt
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
K = 1 / eps_0  # established bulk modulus = 1/eps_0 (Section 1.1)

print(SEP)
print("BAR MAGNET TWO-POLE PRESSURE MODEL -- FAR-FIELD POWER-LAW CHECK")
print(SEP)

# ── Geometry: N pole at +ell/2, S pole at -ell/2, along z-axis ──────────────
ell = 1.0       # pole separation (arbitrary units -- this is a SHAPE check, not absolute)
Q_pole = 1.0    # pole "strength" (arbitrary units -- not yet derived from electron population)

def pressure_at(point, Q, pos):
    """Single-source pressure, reusing the established Coulomb-pressure formula
    P(r) = Q*K/(4*pi*r) directly (doc_magnetism Section 1.2)."""
    r_vec = np.array(point) - np.array(pos)
    r = np.linalg.norm(r_vec)
    return Q * K / (4 * pi * r)

def dipole_pressure(point):
    posN = np.array([0, 0, ell / 2])
    posS = np.array([0, 0, -ell / 2])
    return pressure_at(point, Q_pole, posN) + pressure_at(point, -Q_pole, posS)

# =============================================================================
print()
print("PART 1: RAW PRESSURE (POTENTIAL-LIKE) FALLOFF -- EXPECT 1/r^2, NOT 1/r^3")
print(SEP2)
print("  Sampling along the polar axis (theta=0, the 'axial' direction) at")
print("  increasing distance, fitting the log-log slope.")

r_values = np.array([10.0, 20.0, 40.0, 80.0, 160.0, 320.0]) * ell
axial_pressures = np.array([dipole_pressure([0, 0, r]) for r in r_values])

log_r = np.log(r_values)
log_p = np.log(np.abs(axial_pressures))
slope_pressure, intercept = np.polyfit(log_r, log_p, 1)
print(f"  r values (units of ell): {r_values/ell}")
print(f"  Fitted log-log slope (axial pressure vs r): {slope_pressure:.4f}")
print(f"  Expected: -2.0 (potential-like dipole, NOT -3.0)")
check("P1 raw pressure dipole falls off as 1/r^2 (potential-like), matching "
      "a genuine dipole of point sources, not yet the B-field target",
      abs(slope_pressure - (-2.0)) < 0.01, f"slope = {slope_pressure:.4f}")

# =============================================================================
print()
print("PART 2: GRADIENT OF PRESSURE (FIELD-LIKE) -- EXPECT 1/r^3, THE REAL B TARGET")
print(SEP2)
print("  B = curl(A) is a DERIVATIVE of the displacement field (Section 1.3),")
print("  not the pressure/potential itself -- so the correct comparison to the")
print("  standard B(r) ~ mu_0/(4pi) * [3*rhat*(m.rhat)-m]/r^3 result (Wikipedia,")
print("  'Magnetic dipole') is the GRADIENT of this pressure dipole, not the")
print("  raw pressure.")

def pressure_gradient_axial(r, h=1e-4):
    """Numerical d(pressure)/dz along the polar axis at distance r."""
    p_plus = dipole_pressure([0, 0, r + h])
    p_minus = dipole_pressure([0, 0, r - h])
    return (p_plus - p_minus) / (2 * h)

grad_values = np.array([pressure_gradient_axial(r) for r in r_values])
log_grad = np.log(np.abs(grad_values))
slope_grad, intercept_grad = np.polyfit(log_r, log_grad, 1)
print(f"  Fitted log-log slope (pressure GRADIENT vs r): {slope_grad:.4f}")
print(f"  Expected: -3.0 (matches the real B-field falloff)")
check("P2 gradient of the pressure dipole falls off as 1/r^3, matching the "
      "real external B-field falloff (Wikipedia dipole formula)",
      abs(slope_grad - (-3.0)) < 0.01, f"slope = {slope_grad:.4f}")

# =============================================================================
print()
print("PART 3: OFF-AXIS ANGULAR DEPENDENCE (theta=90deg, equatorial plane)")
print(SEP2)
print("  Standard dipole field has HALF the axial magnitude (same power law)")
print("  in the equatorial plane, with OPPOSITE sign for the gradient/field")
print("  component along the dipole axis -- check the shape, not just the r-scaling.")

r_equatorial = 100.0 * ell


def dipole_pressure_offaxis(r, theta):
    z = r * math.cos(theta)
    x = r * math.sin(theta)
    return dipole_pressure([x, 0, z])


p_axial = dipole_pressure_offaxis(r_equatorial, 0.0)
p_equatorial = dipole_pressure_offaxis(r_equatorial, pi / 2)
print(f"  Pressure at theta=0 (axial):  {p_axial:.6e}")
print(f"  Pressure at theta=90 (equatorial): {p_equatorial:.6e}")
print(f"  Ratio (should be ~0, equatorial pressure vanishes by symmetry "
      f"for a pure dipole -- pressure itself is a SCALAR odd function of z):")
check("P3 equatorial-plane pressure vanishes by symmetry (odd function of z, "
      "as expected for a scalar dipole potential)",
      abs(p_equatorial) < 1e-6 * abs(p_axial), f"|p_eq/p_axial| = {abs(p_equatorial/p_axial):.2e}")

# =============================================================================
print()
print(SEP)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"RESULT: {len(results)} checks  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
    print()
    print("  CONCLUSION: the simplest two-pole pressure model, built entirely from")
    print("  the ALREADY-ESTABLISHED Coulomb-pressure Green's function (no new")
    print("  physics), has the CORRECT mathematical shape to be a candidate for the")
    print("  external field: raw pressure falls as 1/r^2 (potential-like, as it")
    print("  must), and its gradient -- the field-like quantity that should compare")
    print("  to B -- falls as 1/r^3, matching the real external dipole field law.")
    print("  This is a NECESSARY, not sufficient, consistency check: it confirms")
    print("  the GEOMETRIC SETUP is right, not that Q_pole has any particular")
    print("  physical value yet. Deriving Q_pole from electron population/flux at")
    print("  the magnet's ends is the next, still-unattempted step.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAIL: {name}")
            print(f"        {detail}")
print(SEP)
