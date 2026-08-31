"""
flyby_anomaly.py
Torsion spin-loading prediction for the Earth flyby anomaly +
wrapper closure / rotation saturation threshold analysis.

Torsion Universe Framework
Run: python analysis/flyby_anomaly.py

Reads constants from analysis/constants.py.
See medium_chains.txt Sections 13-14 for full discussion.
See whitepaper.txt Appendix E for the R_s conjecture.

Sources:
  Anderson et al. 2008, Phys Rev Lett 100:091102
  Gravity Probe B: Everitt et al. 2011, PRL 106:221101
"""

import math
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from constants import *

# ============================================================
# PART A — ANDERSON FLYBY COEFFICIENT AND TORSION PREDICTION
# ============================================================

print(SEP)
print("FLYBY ANOMALY — PART A: ANDERSON COEFFICIENT vs TORSION PREDICTION")
print(SEP)

# Anderson 2008 empirical coefficient
K_Anderson = 2 * omega_E * R_E_m / c_m_s
print(f"  Anderson empirical K = 2*omega_E*R_E/c = {K_Anderson:.4e}")
print()

# Flyby data: Anderson et al. 2008, Table 1
# (name, dec_in deg, dec_out deg, dV_measured mm/s, V_inf km/s)
flybys = [
    ("Galileo I",   -12.52,  -34.26,   3.92,  8.949),
    ("NEAR",         20.76,  -72.07,  13.46,  6.851),
    ("Galileo II",  -34.26,   -4.87,  -4.60,  8.877),
    ("Rosetta I",     2.81,   34.29,   1.82,  3.863),
    ("Messenger",    31.44,  -31.92,   0.02,  4.056),
    ("Rosetta II",  -34.29,    2.81,   0.00,  8.396),
]

print(f"  {'Flyby':<14} {'dec_in':>7} {'dec_out':>8} {'dcos':>8} "
      f"{'dV_pred mm/s':>14} {'dV_meas mm/s':>14} {'ratio':>7}")
print(f"  {'-'*14} {'-'*7} {'-'*8} {'-'*8} {'-'*14} {'-'*14} {'-'*7}")

results = []
for name, dec_in, dec_out, dV_meas, V_inf_kms in flybys:
    V_inf = V_inf_kms * 1000
    c_in  = math.cos(math.radians(dec_in))
    c_out = math.cos(math.radians(dec_out))
    dcos  = c_in - c_out
    dV_pred = K_Anderson * V_inf * dcos * 1000   # mm/s
    ratio   = dV_pred / dV_meas if abs(dV_meas) > 0.1 else float('nan')
    results.append((name, dec_in, dec_out, dcos, dV_pred, dV_meas, ratio, V_inf))
    rs = f"{ratio:7.3f}" if not math.isnan(ratio) else "   (--)"
    print(f"  {name:<14} {dec_in:>7.2f} {dec_out:>8.2f} {dcos:>8.4f} "
          f"{dV_pred:>14.3f} {dV_meas:>14.3f} {rs:>7}")

rms = math.sqrt(sum((r[4]-r[5])**2 for r in results if abs(r[5])>0.1) /
                sum(1 for r in results if abs(r[5])>0.1))
print(f"\n  RMS residual (Anderson formula): {rms:.3f} mm/s")

# ── Framework prediction for K ────────────────────────────────────────────
print()
print(SEP)
print("  FRAMEWORK PREDICTION FOR K")
print(SEP)
print()
print("  Two loading modes:")
print("    Spin-loading: asymmetric (azimuthal drag, equatorial maximum)")
print("    Mass-loading: symmetric (all directions incl. poles), cancels on flyby")
print()
print("  GEOMETRIC NOTE: The torsion wrapper compresses poles too (candy-wrapper")
print("  twist convergence). The cos(dec) formula is Anderson's empirical fit from")
print("  a Lense-Thirring analogy. The FULL torsion geometry (polar convergent")
print("  compression + equatorial tangential drag) may produce a different geometric")
print("  factor, shifting c_torsion from R_s*c. See medium_chains.txt Sec 13.")
print()

K_theory_c = 2 * Rs_exact * omega_E * R_E_m / c_m_s
c_torsion_implied = 2 * Rs_exact * omega_E * R_E_m / K_Anderson
pct_off = (K_theory_c / K_Anderson - 1) * 100

print(f"  If c_torsion = c:")
print(f"    K_theory = 2*R_s*omega_E*R_E/c = {K_theory_c:.4e}")
print(f"    K_Anderson                     = {K_Anderson:.4e}")
print(f"    Deviation: {pct_off:+.2f}%")
print()
print(f"  Implied c_torsion for exact match:")
print(f"    c_torsion = 2*R_s*omega_E*R_E / K_Anderson")
print(f"              = {c_torsion_implied:.4e} m/s")
print(f"              = {c_torsion_implied/c_m_s:.6f} * c")
print(f"              = {c_torsion_implied/c_m_s / Rs_exact:.6f} * R_s * c")
print()
print(f"  R_s = {Rs_exact:.6f}")
print(f"  c_torsion / c = {c_torsion_implied/c_m_s:.6f}")
print(f"  Ratio c_torsion/(R_s*c) = {c_torsion_implied/(Rs_exact*c_m_s):.6f}")
print()
if abs(c_torsion_implied/c_m_s - Rs_exact) < 0.001:
    print("  KEY FINDING: c_torsion = R_s * c  (within 0.1%)")
    print("  R_s appears as the torsion propagation speed ratio -- 4th appearance.")
else:
    print(f"  c_torsion = {c_torsion_implied/c_m_s:.4f} * c")
    print(f"  This is {(c_torsion_implied/c_m_s/Rs_exact - 1)*100:+.2f}% from R_s.")

# ============================================================
# PART B — ROTATION SATURATION THRESHOLD
# ============================================================

print()
print(SEP)
print("FLYBY ANOMALY — PART B: WRAPPER CLOSURE / ROTATION SATURATION")
print(SEP)
print()
print(f"  Saturation condition:  v_rot / v_esc  >=  R_s = {Rs_exact:.4f}")
print(f"  Relaxation condition:  P_rot < tau_lower ~ 1e8 s (~1157 days)")
print(f"  Both must hold for full wrapper closure.")
print()

# (name, v_rot m/s, v_esc m/s, P_rot days, note)
bodies = [
    ("Sun",          1994.0,  617600.0,   25.40, "solar rotation"),
    ("Mercury",         3.0,    4250.0,   58.65, "slow rotator"),
    ("Venus",           1.8,   10360.0,  243.00, "very slow, P near tau_lower"),
    ("Earth",         465.1,   11186.0,    1.00, "Anderson flyby body"),
    ("Mars",          241.2,    5030.0,    1.03, ""),
    ("Jupiter",     12570.0,   59500.0,    0.41, "fast, fully saturated"),
    ("Saturn",       9870.0,   35500.0,    0.44, "fast, fully saturated"),
    ("Uranus",       2590.0,   21300.0,    0.72, "approaching threshold"),
    ("Neptune",      2680.0,   23500.0,    0.67, "approaching threshold"),
    ("Neutron star", 1.5e7,    1.8e8,      0.001,"ms pulsar estimate"),
]

tau_lower_s = 1e8
earth_sat   = 465.1 / 11186.0 / Rs_exact
K_max_full  = K_Anderson / earth_sat    # K if Earth were fully saturated

print(f"  {'Body':<14} {'v_rot/v_esc':>12} {'sat_frac':>10} {'P(days)':>9} "
      f"{'K_pred':>12}  {'Status'}")
print(f"  {'-'*14} {'-'*12} {'-'*10} {'-'*9} {'-'*12}  ------")

for name, v_rot, v_esc, p_rot, note in bodies:
    sat_frac = v_rot / v_esc / Rs_exact
    relax_ok = p_rot * 86400 < tau_lower_s
    k_pred   = K_max_full * min(sat_frac, 1.0) if relax_ok else 0.0
    if sat_frac >= 1.0 and relax_ok:
        status = "SATURATED"
    elif sat_frac >= 0.5:
        status = "NEAR"
    elif not relax_ok:
        status = "OPEN (slow)"
    else:
        status = "sub-sat"
    print(f"  {name:<14} {v_rot/v_esc:>12.5f} {sat_frac:>10.4f} {p_rot:>9.3f} "
          f"{k_pred:>12.4e}  {status}  {note}")

print()
print(f"  Earth fully-saturated K_max = {K_max_full:.4e}  "
      f"({K_max_full/K_Anderson:.1f}x Anderson K)")
print()

# Key predictions
print(SEP)
print("  KEY PREDICTIONS")
print(SEP)
jup_sat  = 12570.0 / 59500.0 / Rs_exact
sat_sat  = 9870.0  / 35500.0 / Rs_exact
mars_sat = 241.2   / 5030.0  / Rs_exact
print(f"  P1: Jupiter flyby  K ~ {K_max_full:.3e}  "
      f"({K_max_full/K_Anderson:.1f}x Earth Anderson K)  [sat_frac={jup_sat:.2f}]")
print(f"  P2: Saturn flyby   K ~ {K_max_full:.3e}  "
      f"({K_max_full/K_Anderson:.1f}x Earth Anderson K)  [sat_frac={sat_sat:.2f}]")
print(f"  P3: Mars flyby     K ~ {K_Anderson*mars_sat/earth_sat:.3e}  "
      f"({mars_sat/earth_sat:.2f}x Earth Anderson K)  [sat_frac={mars_sat:.2f}]")
print(f"  P4: Venus flyby    K ~ 0  (sat_frac=0.001, P=243d near tau_lower)")
print()
v_thresh_earth     = Rs_exact * 11186.0
P_thresh_earth_hrs = (2 * pi * R_E_m / v_thresh_earth) / 3600
print(f"  Earth saturation threshold: v_rot = {v_thresh_earth:.0f} m/s  "
      f"-> P_rot = {P_thresh_earth_hrs:.1f} hours")
print(f"  (Earth currently {v_thresh_earth/465.1:.1f}x too slow for full saturation)")
print()
print("  TESTABILITY: Cassini (Saturn) and Galileo (Jupiter) flyby datasets exist.")
print("  If K_Jupiter ~ 4x K_Anderson, this confirms saturation scaling prediction.")
print(SEP)

# ============================================================
# PART C — CROSS-BODY K SCALING: PREDICTION vs DATA AVAILABILITY
# ============================================================

print()
print(SEP)
print("FLYBY ANOMALY — PART C: CROSS-BODY K SCALING TEST")
print(SEP)
print()

# --- Framework prediction summary ---
sat_frac_E  = 465.1  / 11186.0 / Rs_exact    # 0.2337
sat_frac_J  = 12570.0 / 59500.0 / Rs_exact   # 1.187 -> capped at 1.0
sat_frac_S  = 9870.0  / 35500.0 / Rs_exact   # 1.563 -> capped at 1.0
K_max       = K_Anderson / sat_frac_E         # fully-saturated K
K_pred_J    = K_max                           # Jupiter: saturated
K_pred_S    = K_max                           # Saturn:  saturated

print("  FRAMEWORK PREDICTION (saturation scaling):")
print(f"    K_Earth(Anderson) = {K_Anderson:.4e}  [sat_frac = {sat_frac_E:.4f}]")
print(f"    K_Jupiter(pred)   = {K_pred_J:.4e}  [sat_frac = {sat_frac_J:.4f}, SATURATED]")
print(f"    K_Saturn(pred)    = {K_pred_S:.4e}  [sat_frac = {sat_frac_S:.4f}, SATURATED]")
print(f"    Ratio K_Jup/K_E   = {K_pred_J/K_Anderson:.2f}x")
print()

# --- Known Jupiter/Saturn flybys available for testing ---
# Trajectory data from published mission profiles and Jet Propulsion Lab records.
# V_inf values from JPL Horizons / published mission documentation.
# Dec_in/dec_out relative to each body's equatorial plane (Jupiter J2000 equator,
# Saturn equatorial plane) -- marked APPROX where not from primary analysis paper.
#
# GEOMETRY CAVEAT: The Anderson formula was derived empirically from Earth flybys
# at altitudes of ~960 km to ~3300 km (1.15-1.52 R_E).  For Jupiter flybys at
# 136 R_J (Cassini 2000), the same near-field formula almost certainly does not
# apply without a trajectory-integration correction.  The predicted dV below is
# therefore a NAIVE UPPER BOUND assuming the near-field formula extends to the
# actual flyby distance.  A proper test requires numerical integration of the
# torsion field along the actual hyperbolic trajectory.

print("  KNOWN JUPITER/SATURN FLYBYS (candidate tests):")
print()

# (mission, body, year, periapsis_Rbody, V_inf_km_s, dec_in_deg, dec_out_deg, dV_meas_note)
candidates = [
    # Cassini Jupiter flyby: Dec 30 2000, periapsis ~136 R_J (9.72e6 km)
    # V_inf from published trajectory: 10.41 km/s heliocentric Jupiter-relative
    # Dec_in/out relative to Jupiter equator: APPROX from approach geometry
    ("Cassini",    "Jupiter", 2000, 136.0,  10.41, -24.0,  26.0,
     "No published anomaly residual in Anderson 2008 format"),
    # Pioneer 11 Jupiter flyby: Dec 3 1974, periapsis ~1.6 R_J
    # V_inf ~11.5 km/s; dec estimates APPROX
    ("Pioneer 11", "Jupiter", 1974,   1.6,  11.50, -40.0,  60.0,
     "Thermal subtraction viable (Turyshev 2012 model); data quality TBD"),
    # Cassini Saturn flyby: NOT a Saturn gravity assist -- arrived in orbit 2004
    # Voyager 1 Saturn flyby: Nov 12 1980, periapsis ~3.1 R_S
    # V_inf ~5.9 km/s; dec APPROX
    ("Voyager 1",  "Saturn",  1980,   3.1,   5.90, -25.0,  65.0,
     "Pre-navigation-revolution precision; residuals not at mm/s level"),
    # Voyager 2 Saturn flyby: Aug 26 1981, periapsis ~2.7 R_S
    ("Voyager 2",  "Saturn",  1981,   2.7,   5.40, -20.0,  60.0,
     "Pre-navigation-revolution precision; residuals not at mm/s level"),
]

print(f"  {'Mission':<12} {'Body':<8} {'Yr':>4} {'r_peri/R':>9} {'V_inf':>7} "
      f"{'dec_in':>7} {'dec_out':>8} {'dV_pred mm/s':>14}")
print(f"  {'-'*12} {'-'*8} {'-'*4} {'-'*9} {'-'*7} {'-'*7} {'-'*8} {'-'*14}")

body_K = {"Jupiter": K_pred_J, "Saturn": K_pred_S}
for (mission, body, yr, r_peri, V_inf_kms, dec_in, dec_out, note) in candidates:
    V_inf   = V_inf_kms * 1000.0
    dcos    = math.cos(math.radians(dec_in)) - math.cos(math.radians(dec_out))
    K_body  = body_K[body]
    dV_pred = K_body * V_inf * dcos * 1000.0   # mm/s (naive, near-field formula)
    print(f"  {mission:<12} {body:<8} {yr:>4} {r_peri:>9.1f} {V_inf_kms:>7.2f} "
          f"{dec_in:>7.1f} {dec_out:>8.1f} {dV_pred:>14.1f}")

print()
print("  (dec angles relative to each body's equatorial plane; marked APPROX)")
print()

# --- Geometry correction factor ---
print("  GEOMETRY CORRECTION FOR DISTANT FLYBYS:")
print()
print("  The Anderson formula is derived from the near-field torsion gradient.")
print("  For a flyby at periapsis distance r_peri = n * R_body, the effective")
print("  K scales roughly as (R_body / r_peri)^2 relative to the near-surface")
print("  derivation.  This is a first-order correction only.")
print()

R_J_m = 6.991e7   # m
R_S_m = 5.823e7   # m

for (mission, body, yr, r_peri, V_inf_kms, dec_in, dec_out, note) in candidates:
    V_inf  = V_inf_kms * 1000.0
    dcos   = math.cos(math.radians(dec_in)) - math.cos(math.radians(dec_out))
    K_body = body_K[body]
    # Correction: Earth flybys were at ~1.15 R_E; scale to actual r_peri
    r_earth_ref = 1.15
    geom_factor = (r_earth_ref / r_peri) ** 2
    K_corrected = K_body * geom_factor
    dV_corr     = K_corrected * V_inf * dcos * 1000.0
    print(f"  {mission} {yr}: r_peri = {r_peri:.1f} R_{body[0]}  "
          f"geom_factor = ({r_earth_ref:.2f}/{r_peri:.1f})^2 = {geom_factor:.4e}  "
          f"dV_corrected = {dV_corr:.2f} mm/s")

print()

# --- Data availability verdict ---
print(SEP)
print("  DATA AVAILABILITY VERDICT")
print(SEP)
print()
print("  Cassini Jupiter (2000):")
print("    - BEST candidate: modern precision navigation, post-2000 Doppler quality")
print("    - V_inf and approach geometry reconstructible from SPICE kernels")
print("    - Published velocity residual: NOT in Anderson 2008 (Earth flybys only)")
print("    - Geometry-corrected predicted signal: ~0.01 mm/s (Cassini at 136 R_J)")
print("    - VERDICT: Signal too small at 136 R_J; geometry dilutes prediction")
print("      by factor ~1.4e-4 relative to near-surface formula.  Not a usable test.")
print()
print("  Pioneer 11 Jupiter (1974):")
print("    - periapsis 1.6 R_J -- NEAR-FIELD, geometry-corrected dV ~ 21 mm/s")
print("    - The Pioneer anomaly is thermal radiation pressure from the RTGs.")
print("      The RTG fuel (Pu-238, half-life 87.7 yr) decays slowly.  In 1974")
print("      the thermal anomaly was ~15% larger than the late-mission 8.74e-10 m/s^2,")
print("      but this is calculable from Turyshev et al. 2012 physics-based thermal model.")
print("    - Thermal contribution DURING the flyby arc (~1 hr near periapsis):")
print("        a_thermal * t_flyby ~ 9.5e-10 * 3600 ~ 3.4 mm/s contamination")
print("      This is much smaller than the predicted flyby signal of ~21 mm/s.")
print("      Signal-to-thermal ratio: ~6:1.  The thermal component IS separable")
print("      in principle using the Turyshev thermal model at the 1974 epoch.")
print("    - ACTUAL LIMITING FACTOR: 1974 Doppler tracking data quality.")
print("      DSN S-band tracking in 1974 had less precision than post-2000 missions.")
print("      Turyshev 2012 reprocessed the full dataset but the Jupiter flyby arc")
print("      was not the target.  Whether the 1974 residuals reach mm/s precision")
print("      after thermal subtraction requires checking against the raw archives.")
print("    - VERDICT: Potentially viable archival test.  Signal is large enough to")
print("      survive thermal subtraction.  Requires Turyshev thermal model at 1974")
print("      epoch + reprocessing of Jupiter flyby Doppler arc specifically.")
print()
print("  Voyager 1/2 Saturn (1980/1981):")
print("    - periapsis 3 R_S -- reasonable geometry")
print("    - Navigation precision of 1970s/80s: ~100 mm/s not ~1 mm/s")
print("    - VERDICT: Insufficient navigation precision for the predicted signal.")
print()
print("  CONCLUSION: No currently available dataset cleanly tests K_Jupiter/K_Earth")
print("  scaling.  Best existing candidate: Pioneer 11 Jupiter flyby (1974) --")
print("  near-field geometry (1.6 R_J), signal ~21 mm/s, thermal contamination ~3 mm/s,")
print("  separable in principle using Turyshev 2012 thermal model.  Requires archival")
print("  reprocessing of the 1974 Doppler arc with thermal subtraction as the target.")
print("  Ideal clean test: future Jupiter flyby at < 5 R_J with modern Doppler tracking.")
print(f"    K_Jupiter/K_Earth = {K_pred_J/K_Anderson:.2f}x  (saturation model)")
print()
print("  WHAT WOULD CONFIRM OR FALSIFY:")
print("    Confirm: dV residual at ~50-200 mm/s for a close Jupiter flyby")
print("    Falsify: dV residual < 5 mm/s for a close Jupiter flyby at < 5 R_J")
print("    (Mars flyby at K ~ 1.15x is feasible with upcoming Mars gravity assists)")
print(SEP)
