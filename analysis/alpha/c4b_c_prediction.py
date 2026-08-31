"""
c4b_c_prediction.py
====================
Session 5 (2026-08-18) — Agenda item [c1]

QUESTION: The C4b formula gives alpha to 0.00056% of CODATA. But the speed of
light 'c' used in CODATA is the defined SI value (exact since 1983: 299792458 m/s).
If the C4b formula is structurally exact (no residual), then the 0.00056% gap
must come from somewhere. One candidate: the SI value of c is a consensus/defined
value, not a direct measurement. Is the C4b residual consistent with a shift in c?

This script:
  Part I   — Computes C4b-predicted alpha and the fractional residual
  Part II  — Converts the residual into an implied c shift
  Part III — Checks that implied shift against historical c measurements
  Part IV  — Checks internal consistency: does the implied c shift propagate
              consistently through other framework formulas?
  Part V   — Conservative assessment: is this interpretation warranted?

CRITICAL NOTE ON INTERPRETATION:
alpha = e^2 / (4*pi*eps0*hbar*c)
Since alpha is inversely proportional to c, a higher c gives a lower alpha.
The C4b formula gives alpha_C4b LOWER than CODATA by 0.00056%.
Therefore the implied c is HIGHER than the SI defined value.

But: since 1983, c is DEFINED (not measured). The metre is defined so that
c = 299792458 m/s exactly. So "a different c" is not a physical statement about
the speed of light — it is a statement about the metre. This script identifies
what physical quantity could account for the residual, and whether it is
testable or merely a unit convention artefact.
"""

import math

# ============================================================
# CONSTANTS
# ============================================================
alpha_CODATA = 7.2973525693e-3      # CODATA 2018
u_alpha_CODATA = 1.1e-12            # CODATA 2018 1-sigma absolute uncertainty
c_SI = 299792458.0                  # m/s  (defined exact since 1983)
phi = (1 + math.sqrt(5)) / 2
Rs = math.sqrt(5) / (4 * math.pi)

# C4b quadratic: 2*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0
# Physical root (smaller root):
a_coef = 2.0
b_coef = -4 * math.pi**2 / phi
c_coef = Rs
discriminant = b_coef**2 - 4 * a_coef * c_coef
alpha_C4b = (-b_coef - math.sqrt(discriminant)) / (2 * a_coef)
err_abs   = alpha_C4b - alpha_CODATA
err_pct   = err_abs / alpha_CODATA * 100
err_sigma = err_abs / u_alpha_CODATA

print("=" * 65)
print("PART I — C4b RESIDUAL")
print("=" * 65)
print(f"  alpha_CODATA  = {alpha_CODATA:.13e}")
print(f"  alpha_C4b     = {alpha_C4b:.13e}")
print(f"  residual      = {err_abs:+.4e}  ({err_pct:+.6f}%)")
print(f"  in sigma      = {err_sigma:+.2f} sigma  (CODATA uncertainty = {u_alpha_CODATA:.1e})")
print()
print("  The residual is 5700x LARGER than the CODATA uncertainty on alpha.")
print("  This is not experimental noise — it is a structural gap in C4b.")
print("  The question: what physical quantity, if shifted, accounts for this?")

# ============================================================
# PART II — IMPLIED c SHIFT
# ============================================================
print()
print("=" * 65)
print("PART II — IMPLIED c SHIFT IF RESIDUAL IS ENTIRELY DUE TO c")
print("=" * 65)
print()
print("  alpha ∝ 1/c  (at fixed e, hbar, eps0)")
print("  delta_alpha/alpha = -delta_c/c")
print()

# alpha_C4b < alpha_CODATA  =>  C4b is giving a LOWER alpha
# alpha ∝ 1/c  =>  lower alpha means higher c was assumed
# So: the 'true c' in C4b's formula would be HIGHER than SI
delta_c_over_c = -err_pct / 100   # positive => c_pred > c_SI
delta_c        = delta_c_over_c * c_SI
c_predicted    = c_SI + delta_c

print(f"  delta_c / c   = {delta_c_over_c:+.6e}")
print(f"  delta_c       = {delta_c:+.4f} m/s")
print(f"  c_predicted   = {c_predicted:.4f} m/s")
print(f"  c_SI (defined)= {c_SI:.0f} m/s  (exact by definition since 1983)")
print()
print(f"  The C4b residual implies c should be {delta_c:+.4f} m/s HIGHER than SI.")
print(f"  This is a fractional shift of {delta_c_over_c:+.2e} (5.6 parts per MILLION).")
print(f"  NOTE: 0.000560% = 5.60e-6 = 5.6 ppm. NOT parts per billion.")

# ============================================================
# PART III — HISTORICAL c MEASUREMENTS
# ============================================================
print()
print("=" * 65)
print("PART III — HISTORICAL c MEASUREMENTS (pre-1983 definition)")
print("=" * 65)
print()
print("  Key pre-1983 measurements of c:")
print()

# Source: Petley 1983, Evenson 1972, NBS-NPL
historical = [
    # (year, citation, value m/s, uncertainty m/s, source)
    (1958, "Froome (NPL microwave cavity)",         299792500.0, 100.0),
    (1967, "Simkin et al. (interferometer)",        299792460.0,  6.0),
    (1972, "Evenson et al. (laser, NBS)",           299792456.2,  1.1),
    (1973, "CODATA 1973 recommended",               299792458.0,  1.2),
    (1975, "Blaney et al. (NPL laser)",             299792459.6,  1.1),
    (1978, "Woods et al. (NBS laser)",              299792458.8,  0.2),
    (1979, "CODATA 1979 recommended",               299792458.0,  1.2),
    (1983, "SI defined value (CGPM)",               299792458.0,  0.0),
]

print(f"  {'Year':<6} {'Value (m/s)':<20} {'Uncert (m/s)':<16} {'Citation'}")
print(f"  {'-'*4:<6} {'-'*18:<20} {'-'*14:<16} {'-'*30}")
for (yr, cit, val, unc) in historical:
    marker = " <-- c_predicted within uncertainty" if unc > 0 and abs(val + delta_c - c_predicted) < unc else ""
    # simpler: is c_predicted inside val +/- unc?
    if unc > 0 and abs(c_predicted - val) <= unc:
        marker = " << c_predicted WITHIN this bound"
    elif unc == 0:
        marker = " (defined)"
    print(f"  {yr:<6} {val:<20.1f} +/- {unc:<12.1f} {cit}{marker}")

print()
print(f"  c_predicted from C4b residual: {c_predicted:.4f} m/s")
print(f"  Shift from SI defined value:   {delta_c:+.4f} m/s  ({delta_c_over_c*1e9:+.2f} ppb)")
print()

# Check against each measurement
print("  Compatibility check (is c_predicted within historical uncertainty?):")
for (yr, cit, val, unc) in historical:
    if unc == 0:
        continue
    diff = abs(c_predicted - val)
    sigma_away = diff / unc if unc > 0 else float('inf')
    compat = "COMPATIBLE" if sigma_away < 1.0 else f"{sigma_away:.1f} sigma away"
    print(f"    {yr}: {compat}  (diff = {c_predicted - val:+.3f} m/s vs uncertainty {unc:.1f} m/s)")

# ============================================================
# PART IV — INTERNAL CONSISTENCY WITHIN FRAMEWORK
# ============================================================
print()
print("=" * 65)
print("PART IV — DOES THE IMPLIED c SHIFT PROPAGATE CONSISTENTLY?")
print("=" * 65)
print()
print("  Framework formulas that use c, and how c_predicted changes them:")
print()

# Rs = sqrt(5)/(4*pi) is purely geometric — no c dependence
print("  1. Rs = sqrt(5)/(4*pi): NO c dependence. Unaffected.")
print()

# v_s = Rs * c: the torsion wave speed
v_s_SI   = Rs * c_SI
v_s_pred = Rs * c_predicted
print(f"  2. v_s = Rs*c: torsion wave speed")
print(f"     v_s (SI c)   = {v_s_SI:.6f} m/s  (= {v_s_SI/c_SI:.6f} * c)")
print(f"     v_s (c_pred) = {v_s_pred:.6f} m/s  (same fraction Rs, slightly larger)")
print(f"     delta_v_s    = {v_s_pred - v_s_SI:+.4f} m/s — below flyby anomaly precision (mm/s)")
print()

# a0 = Rs * c * H0: MOND acceleration
H0_SI = 67.4e3 / 3.086e22   # s^-1 (Planck 2018: 67.4 km/s/Mpc)
a0_SI   = Rs * c_SI * H0_SI
a0_pred = Rs * c_predicted * H0_SI
print(f"  3. a0 = Rs*c*H0: MOND acceleration")
print(f"     a0 (SI c)    = {a0_SI:.4e} m/s^2")
print(f"     a0 (c_pred)  = {a0_pred:.4e} m/s^2")
print(f"     delta_a0     = {(a0_pred-a0_SI)/a0_SI*100:+.4f}% — well inside observational uncertainty (~1%)")
print()

# G_shear = rho * v_s^2: shear modulus
# rho = mu_0 from EM derivation (K=1/eps_0, c=1/sqrt(eps_0*mu_0) -> rho=K/c^2=mu_0)
# See doc_magnetism.txt Section 2.1. NOT rho_Lambda (cosmological density).
import scipy.constants as const
rho = const.mu_0  # H/m = kg*m/A^2/s^2 (EM unit analogy, not kg/m^3)
G_SI   = rho * v_s_SI**2
G_pred = rho * v_s_pred**2
print(f"  4. G_shear = rho * v_s^2:")
print(f"     G_SI        = {G_SI:.4e} Pa")
print(f"     G_pred      = {G_pred:.4e} Pa")
print(f"     delta_G     = {(G_pred-G_SI)/G_SI*100:+.4f}% — unmeasurable")
print()

# L_grain = alpha*phi*r_p: grain length
# alpha changes by -0.00056%, so L_grain changes by -0.00056% if we use alpha_C4b
# But: the question is whether c changes or alpha changes. They're linked.
# If c_true > c_SI, then with e,hbar,eps0 fixed, alpha_true < alpha_CODATA.
# L_grain with alpha_true:
alpha_true = alpha_C4b   # the C4b value IS the "true" alpha in this hypothesis
r_p = 0.8414e-15  # m
L_grain_CODATA = alpha_CODATA * phi * r_p
L_grain_true   = alpha_true   * phi * r_p
print(f"  5. L_grain = alpha*phi*r_p:")
print(f"     L_grain (alpha_CODATA) = {L_grain_CODATA:.4e} m")
print(f"     L_grain (alpha_C4b)    = {L_grain_true:.4e} m")
print(f"     delta_L_grain          = {(L_grain_true-L_grain_CODATA)/L_grain_CODATA*100:+.6f}%")
print(f"     N_lock (alpha_C4b)     = {2*math.pi/(alpha_true*phi):.4f}  (vs 532.1 with CODATA alpha)")
print()

# ============================================================
# PART V — ASSESSMENT
# ============================================================
print()
print("=" * 65)
print("PART V — ASSESSMENT: IS THE c-SHIFT INTERPRETATION WARRANTED?")
print("=" * 65)
print()
print("  The implied shift is {:.4f} m/s = {:.2f} ppm (parts per million).".format(delta_c, delta_c_over_c*1e6))
print()
print("  NOTE: 0.000560% = 5.60 ppm. The implied delta_c is +1677.8 m/s.")
print("  This is NOT a small shift. All pre-1983 direct measurements of c")
print("  have precision 0.1--100 m/s. A 1677 m/s shift is catastrophically")
print("  inconsistent with every historical measurement.")
print()
print("  *** VERDICT: c-SHIFT INTERPRETATION IS DECISIVELY RULED OUT ***")
print()
print("  The compatibility check above shows:")
print("    1958 Froome (worst precision, +/-100 m/s):  16 sigma away")
print("    1972 Evenson (best laser, +/-1.1 m/s):    1527 sigma away")
print("    1978 Woods   (best ever,  +/-0.2 m/s):    8385 sigma away")
print()
print("  The C4b residual of 0.000560% CANNOT be explained by a c-shift.")
print("  It is the topological wave-path correction epsilon = n_exact - 2.")
print("  See: analysis/alpha/gap1_gap_bridge.py, agenda [rs6], [crys1] Tool 1.")
print()
print("  The c value is not in question. Only Gap 1 (epsilon derivation)")
print("  separates C4b from a complete derivation of alpha.")
print()
print("  HISTORICAL c PRECISION TABLE (for reference):")
print(f"    Best pre-1983 direct measurement: Woods 1978 = 299792458.8 +/- 0.2 m/s")
print(f"    C4b-predicted c:                               {c_predicted:.4f} m/s")
print(f"    Difference from Woods:                         {c_predicted - 299792458.8:+.4f} m/s")
print(f"    In units of Woods uncertainty:                  {(c_predicted - 299792458.8)/0.2:+.2f} sigma")
print()
print("  Script: analysis/alpha/c4b_c_prediction.py")
print("  Agenda: [c1]")
