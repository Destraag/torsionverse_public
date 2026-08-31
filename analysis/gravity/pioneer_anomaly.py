"""
pioneer_anomaly.py
Tests whether the torsion mass-loading mode can predict the Pioneer anomaly.

Torsion Universe Framework
Run: python analysis/pioneer_anomaly.py

Reads constants from analysis/constants.py.
See medium_chains.txt Section 14, Open Question Q1 for context.

Pioneer anomaly: a_P = 8.74e-10 m/s^2, roughly constant over 20-70 AU.
The Sun's spin-loading is only 1.8% saturated (flyby_anomaly.py),
so the Pioneer anomaly must be dominated by mass-loading, not spin-loading.

Sources:
  Turyshev et al. 2012, PRL 108:241101 (thermal recoil explanation)
  Anderson et al. 2002, Phys Rev D 65:082004 (original anomaly paper)
  Milgrom 1983, ApJ 270:365 (MOND)
"""

import math
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from constants import *

a_Pioneer = 8.74e-10   # m/s^2 measured (Anderson 2002)
distances_AU = [20, 30, 40, 50, 60, 70]

SEP2 = "=" * 62

print(SEP2)
print("PIONEER ANOMALY — MASS-LOADING PREDICTION CHECK")
print(SEP2)
print(f"  Measured: a_Pioneer = {a_Pioneer:.3e} m/s^2 (roughly constant 20-70 AU)")
print(f"           = {a_Pioneer/a0_m_s2:.2f} * a0")
print(f"  Key structural clue: approximately DISTANCE-INDEPENDENT over 50 AU range.")
print(f"  Any r-dependent formula is immediately suspect.")
print()

# ── F1: Simple mass-loading a = R_s * g_Newton ───────────────────────────────
print(SEP2)
print("  F1: a = R_s * g_Newton(r)  [R_s fraction of Newtonian gravity]")
print(SEP2)
print(f"  {'r (AU)':>8}  {'g_Newton':>12}  {'R_s*g_N':>12}  {'ratio to Pioneer':>18}")
for r_AU in distances_AU:
    r_m  = r_AU * AU_m
    g_N  = G * M_sun_kg / r_m**2
    F1   = Rs_exact * g_N
    print(f"  {r_AU:>8}  {g_N:>12.3e}  {F1:>12.3e}  {F1/a_Pioneer:>18.1f}x")
print(f"  VERDICT: F1 is thousands of times too large at 20 AU, falls steeply.")
print(f"           R-dependent. Does NOT explain Pioneer.")
print()

# ── F2: Geometric mean a = R_s * sqrt(g_Newton * c * H0) ────────────────────
print(SEP2)
print("  F2: a = R_s * sqrt(g_Newton * c*H0)  [geometric mean of local/cosmological]")
print(SEP2)
print(f"  {'r (AU)':>8}  {'g_Newton':>12}  {'F2':>12}  {'ratio':>18}  {'note'}")
for r_AU in distances_AU:
    r_m  = r_AU * AU_m
    g_N  = G * M_sun_kg / r_m**2
    F2   = Rs_exact * math.sqrt(g_N * cH0)
    note = "<-- closest" if abs(r_AU - 45) < 8 else ""
    print(f"  {r_AU:>8}  {g_N:>12.3e}  {F2:>12.3e}  {F2/a_Pioneer:>18.3f}x  {note}")
# Find crossover distance
r_cross_AU = None
for r_AU in range(1, 200):
    r_m = r_AU * AU_m
    g_N = G * M_sun_kg / r_m**2
    F2  = Rs_exact * math.sqrt(g_N * cH0)
    if abs(F2 - a_Pioneer) / a_Pioneer < 0.05:
        r_cross_AU = r_AU
        break
print(f"  F2 matches Pioneer to 5% at r ~ {r_cross_AU} AU" if r_cross_AU else
      "  F2 has no close-match distance within 200 AU")
print(f"  VERDICT: Closer than F1 but still r-dependent. Not a constant.")
print()

# ── F3: Standard MOND anomalous term ─────────────────────────────────────────
print(SEP2)
print("  F3: MOND anomalous term  a_MOND - a_Newton")
print("  (Pioneer at 20-70 AU is in the Newtonian regime: g_Newton >> a0)")
print(SEP2)
print(f"  {'r (AU)':>8}  {'g_Newton':>12}  {'g_N/a0':>10}  {'MOND correction':>18}  {'ratio':>8}")
for r_AU in distances_AU:
    r_m  = r_AU * AU_m
    g_N  = G * M_sun_kg / r_m**2
    x    = g_N / a0_m_s2
    mu   = x / math.sqrt(1 + x**2)        # standard MOND interpolation
    a_tot = g_N / mu
    corr = a_tot - g_N
    print(f"  {r_AU:>8}  {g_N:>12.3e}  {x:>10.1f}  {corr:>18.3e}  {corr/a_Pioneer:>8.5f}x")
print(f"  VERDICT: MOND correction is ~1e-4 of Pioneer at these distances.")
print(f"           Pioneer is in deep Newtonian regime. MOND cannot explain it.")
print()

# ── F4: Constant floor from Sun's kinematic torsion wake ─────────────────────
print(SEP2)
print("  F4: Sun's kinematic torsion wake through galaxy")
print("  The Sun orbits the galaxy at v_sun ~ 220 km/s.")
print("  The torsion medium's response to this motion may create a roughly")
print("  constant acceleration 'headwind' inside the heliosphere.")
print(SEP2)

v_sun_gal  = 220e3           # m/s  Sun's galactic orbital speed
r_helio_AU = 100             # AU   approximate heliosphere radius
r_helio_m  = r_helio_AU * AU_m

# Kinetic energy of the Sun relative to the galactic medium
E_kin_sun  = 0.5 * M_sun_kg * v_sun_gal**2    # J

# Torsion energy density in heliosphere if medium stores R_s fraction
V_helio    = (4/3) * pi * r_helio_m**3
rho_torsion = Rs_exact * E_kin_sun / V_helio   # J/m^3

print(f"  v_sun (galactic) = {v_sun_gal/1e3:.0f} km/s")
print(f"  E_kin(Sun)       = {E_kin_sun:.3e} J")
print(f"  Heliosphere vol  = {V_helio:.3e} m^3  (r = {r_helio_AU} AU)")
print(f"  Torsion energy density ~ R_s * E_kin / V = {rho_torsion:.3e} J/m^3")
print()
print(f"  To convert to acceleration we need the medium density (not yet defined).")
print(f"  Dimensional check: if medium density ~ rho_torsion / c^2:")
rho_medium_kgm3 = rho_torsion / c_m_s**2
print(f"    rho_medium ~ {rho_medium_kgm3:.3e} kg/m^3")
print(f"    For comparison: proton mass / (1 fm)^3 = {1.67e-27/(1e-15)**3:.3e} kg/m^3 (nuclear density)")
print(f"    Interplanetary medium density ~ 1e-20 to 1e-23 kg/m^3")
print()

# Descriptive n factor (from scale_check.py Part 5)
n_required = a_Pioneer / (Rs_exact * cH0)
print(f"  DESCRIPTIVE: a_Pioneer = {n_required:.2f} * R_s * c * H0  (from Part 5)")
print(f"  This says the local torsion energy density is {n_required:.2f}x cosmological baseline.")
print(f"  It is accurate but NOT a derivation of n from first principles.")
print()

# Can n be predicted from solar properties?
r_45AU = 45 * AU_m
g_N_45 = G * M_sun_kg / r_45AU**2
print(f"  Attempting to derive n = {n_required:.2f} from solar properties at 45 AU:")
candidates_n = [
    ("g_Newton(45AU) / a0",              g_N_45 / a0_m_s2),
    ("sqrt(g_Newton(45AU) / a0)",        math.sqrt(g_N_45 / a0_m_s2)),
    ("(g_Newton(45AU)/a0)^(1/3)",        (g_N_45 / a0_m_s2)**(1/3)),
    ("v_esc(45AU) / (R_s*c)",            math.sqrt(2*G*M_sun_kg/r_45AU) / (Rs_exact*c_m_s)),
    ("(v_sun/c) / R_s^2",               (v_sun_gal/c_m_s) / Rs_exact**2),
    ("2*pi / alpha",                     2*pi / alpha),
    ("R_s * (GM_sun/(r*c^2))^(-1)",      Rs_exact / (G*M_sun_kg/(r_45AU*c_m_s**2))),
]
print(f"  {'Expression':<36}  {'Value':>10}  {'Ratio to n_req':>16}")
print(f"  {'-'*36}  {'-'*10}  {'-'*16}")
for label, val in candidates_n:
    print(f"  {label:<36}  {val:>10.3f}  {val/n_required:>16.4f}")

print()

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP2)
print("PIONEER ANOMALY — SUMMARY")
print(SEP2)
print(f"  Measured: {a_Pioneer:.3e} m/s^2, constant 20-70 AU")
print()
print("  F1 (R_s * g_Newton):          FAIL -- too large, r-dependent")
print("  F2 (R_s * sqrt(g * cH0)):     FAIL -- r-dependent, not constant")
print("  F3 (MOND correction):          FAIL -- negligible at these distances")
print("  F4 (kinematic wake):           OPEN -- energy density set, accel needs medium density")
print()
print("  The Pioneer anomaly is an OPEN PROBLEM for the framework.")
print("  The n=7.22 description is accurate but underived.")
print()
print("  THREE PATHS TO RESOLUTION:")
print()
print("  PATH 1: Thermal recoil (conventional).")
print("    Turyshev et al. 2012 showed ~75% of Pioneer anomaly is thermal radiation")
print("    pressure from the spacecraft's RTGs. Whether the remaining ~25% is real")
print("    or systematic is contested. If fully thermal, the torsion framework")
print("    does not need to explain it.")
print()
print("  PATH 2: Torsion kinematic wake of the Sun.")
print("    Requires knowing the medium density. From the B meson radius prediction")
print(f"    (r_B = {Rs_exact*4.18/0.9:.3f} fm), the medium engages at the QCD string scale.")
print("    A medium density derivation from R_s and QCD parameters may be possible.")
print("    This is the blocking calculation.")
print()
print("  PATH 3: Modified interpolation function.")
print("    The torsion framework has not yet written the full MOND interpolation")
print("    function mu(x) using R_s. Standard MOND uses mu = x/sqrt(1+x^2).")
print("    A torsion-derived mu(x) that transitions at x = R_s (not x = 1) might")
print("    give a non-negligible correction at Pioneer distances.")
print(f"    At 45 AU: x = g_N/a0 = {g_N_45/a0_m_s2:.1f}. Transition at x = 1/R_s = {1/Rs_exact:.1f}.")
print(f"    Pioneer sits at x ~ 100-600 -- well above any transition at x = {1/Rs_exact:.0f}.")
print(f"    PATH 3 is unlikely unless the transition scale is much larger.")
print()
print("  RECOMMENDED ACTION: Determine whether the thermal explanation fully accounts")
print("  for Pioneer before investing in a torsion-based derivation of PATH 2.")
print("  If thermal is sufficient, Pioneer is not a test of the framework.")
print(SEP2)
