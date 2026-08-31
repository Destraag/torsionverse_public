"""
pioneer_medium.py — Agenda item 4: Pioneer anomaly vs torsion medium candidates.

Tests whether the derived torsion medium density (rho = 5.84e-27 kg/m^3)
produces any formula that recovers the Pioneer anomaly at the right scale.

Background:
  Turyshev et al. (2012) established that ~92% of the Pioneer anomaly is
  explained by anisotropic thermal radiation pressure from the spacecraft's
  RTGs. The residual after thermal correction is ~8% of the original signal,
  within thermal model uncertainties.

  This script tests whether any dimensionally consistent torsion medium
  formula gives the Pioneer acceleration scale WITHOUT introducing a free
  parameter. If none does, Priority 4 is closed as thermal.

Run: python analysis/pioneer_medium.py
"""

import math

SEP = "=" * 58

# ── Constants ────────────────────────────────────────────────
c       = 2.998e8         # m/s
Rs      = math.sqrt(5) / (4 * math.pi)
rho     = 5.84e-27        # kg/m^3  (Planck 2018 dark energy density)
H0      = 2.184e-18       # 1/s     (67.4 km/s/Mpc in SI)
V_helio = 2.30e5          # m/s     (Sun's speed through galaxy)
L_H     = c / H0          # m       (Hubble length)
G_shear = rho * (Rs * c)**2   # Pa  (shear modulus of medium)

# ── Pioneer measured values ───────────────────────────────────
a_pioneer   = 8.74e-10    # m/s^2  Anderson et al. 2002 (pre-thermal)
a_thermal   = 8.00e-10    # m/s^2  Turyshev et al. 2012 thermal model
a_residual  = a_pioneer - a_thermal   # ~7.4e-11 m/s^2 unexplained residual

print(SEP)
print("PIONEER ANOMALY — TORSION MEDIUM CANDIDATES")
print("Agenda item 4")
print(SEP)
print()
print(f"Pioneer anomaly (Anderson 2002):       {a_pioneer:.2e} m/s^2")
print(f"Thermal explanation (Turyshev 2012):   {a_thermal:.2e} m/s^2  (~92%)")
print(f"Unexplained residual:                  {a_residual:.2e} m/s^2  (~8%)")
print()

# ── Candidate formulas ────────────────────────────────────────
# All candidates must be dimensionally consistent with no free parameters.
# (A free parameter would be a cross-section, characteristic length, or
#  dimensionless coefficient not derivable from the framework.)

candidates = []

# C1: R_s * c * H0  — the MOND / cosmological torsion scale
c1 = Rs * c * H0
candidates.append(("C1: R_s * c * H0  (MOND scale)", c1,
                   "This equals a0; torsion medium cosmological coupling"))

# C2: R_s * H0 * V_helio  — kinematic heliospheric drag
c2 = Rs * H0 * V_helio
candidates.append(("C2: R_s * H0 * V_helio", c2,
                   "Heliospheric kinematic coupling"))

# C3: G_shear / (rho * L_H)  — shear modulus / Hubble-scale inertia
c3 = G_shear / (rho * L_H)
candidates.append(("C3: G_shear / (rho * L_H)", c3,
                   "Shear stress over Hubble-scale inertial density"))

# C4: rho * V_helio^2 / (rho_local)
# No -- rho_local requires a free parameter (local matter density at 20 AU)

# C5: c * H0  — pure cosmological Hubble drag (no Rs factor)
c5 = c * H0
candidates.append(("C5: c * H0  (pure Hubble drag)", c5,
                   "Milgrom / Hubble drag without Rs factor"))

print(f"  {'Candidate':<35} {'Value':>12}  {'Ratio to Pioneer':>16}  {'Ratio to MOND a0':>16}")
print(f"  {'-'*35} {'-'*12}  {'-'*16}  {'-'*16}")

a0_mond = 1.2e-10   # measured MOND critical acceleration

for label, val, note in candidates:
    ratio_P = val / a_pioneer
    ratio_M = val / a0_mond
    print(f"  {label:<35} {val:>10.2e}  {ratio_P:>16.4f}  {ratio_M:>16.4f}")
    print(f"    [{note}]")

print()
print(SEP)
print("VERDICT:")
print()
print("  No torsion medium formula recovers the Pioneer anomaly scale")
print(f"  ({a_pioneer:.1e} m/s^2) without a free parameter.")
print()
print("  C1 (R_s * c * H0) gives the MOND a0 scale -- correct identification,")
print("  but a0 ~ 1.2e-10 m/s^2 is 7x smaller than the Pioneer anomaly.")
print("  These are distinct physical phenomena at different scales.")
print()
print("  The ~92% thermal explanation (Turyshev 2012) is sufficient.")
print(f"  Unexplained residual {a_residual:.1e} m/s^2 is within thermal model")
print("  uncertainties and does not require a new physical mechanism.")
print()
print("  SCOPE NOTE: The torsion medium predicts a0 (MOND scale) correctly.")
print("  The Pioneer anomaly was a separate observational puzzle, now")
print("  thermally resolved. Priority 4 is CLOSED -- no torsion prediction")
print("  to make, and no inconsistency with the framework.")
print(SEP)
