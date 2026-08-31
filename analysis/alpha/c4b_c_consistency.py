"""
c4b_c_consistency.py
=====================
Session 5 (2026-08-18) — Agenda item [c1] follow-up

QUESTION: c4b_c_prediction.py assumed alpha proportional to 1/c.
Is that correct? If we substitute c consistently through ALL places
it appears in the alpha expression, does the divergence reduce?

This script:
  Part I   — Correct SI algebra: what is alpha proportional to?
  Part II  — Sign correction: which direction does c_implied point?
  Part III — Does substituting c everywhere reduce the C4b gap?
  Part IV  — Fundamental argument: alpha is dimensionless (strongest)
  Part V   — Final verdict and connection to Gap 1
"""

import math

alpha_CODATA = 7.2973525693e-3       # CODATA 2018 (dimensionless, measured)
alpha_C4b    = 7.2973117300057e-3    # C4b quadratic root (dimensionless, predicted)
err_pct      = (alpha_C4b - alpha_CODATA) / alpha_CODATA * 100
c_SI         = 299792458.0           # m/s (defined)
hbar         = 1.054571817e-34       # J*s (defined post-2019 SI)
e_charge     = 1.602176634e-19       # C   (defined post-2019 SI)
mu_0         = 1.25663706212e-6      # N/A^2 (measured post-2019 SI)
eps0_SI      = 1.0 / (mu_0 * c_SI**2)  # F/m (derived from mu_0 and c)

print("=" * 65)
print("PART I — CORRECT SI ALGEBRA: alpha proportional to what?")
print("=" * 65)
print()
print("  alpha = e^2 / (4*pi*eps0*hbar*c)")
print()
print("  In SI, eps0 = 1/(mu_0*c^2). Substituting:")
print("  alpha = e^2 / (4*pi * (1/(mu_0*c^2)) * hbar * c)")
print("        = e^2 * mu_0 * c^2 / (4*pi * hbar * c)")
print("        = e^2 * mu_0 * c / (4*pi * hbar)")
print()
print("  Therefore: alpha is PROPORTIONAL TO c (not 1/c).")
print()
print("  NOTE: The original c4b_c_prediction.py had the sign wrong.")
print("  alpha_C4b < alpha_CODATA means c_implied < c_SI (lower, not higher).")
print()

# Recalculate c_implied with correct proportionality
delta_c_correct = err_pct / 100 * c_SI    # alpha prop to c => delta_alpha/alpha = delta_c/c
c_implied = c_SI + delta_c_correct
print(f"  err_pct           = {err_pct:+.6f}%  (C4b LOWER than CODATA)")
print(f"  delta_c           = {delta_c_correct:+.4f} m/s  (LOWER, not higher)")
print(f"  c_implied         = {c_implied:.4f} m/s")
print(f"  c_SI              = {c_SI:.0f} m/s")
print(f"  |shift|           = {abs(delta_c_correct):.1f} m/s  (same magnitude as before)")
print()
print("  The sign was wrong in the original script; the magnitude is the same.")
print("  1677.8 m/s shift, same historical incompatibility, opposite direction.")

print()
print("=" * 65)
print("PART II — WHAT IF c IS SUBSTITUTED CONSISTENTLY EVERYWHERE?")
print("=" * 65)
print()
print("  The question: if we use c_implied in ALL places in the alpha formula")
print("  (including eps0, which depends on c), does the C4b gap close?")
print()

# Compute alpha directly from fundamental constants at c_implied
# alpha = e^2 * mu_0 * c / (4*pi*hbar)
# This IS the correct formula (eps0 already eliminated via mu_0)
# Note: in post-2019 SI, e and hbar are DEFINED (exact). mu_0 is measured.
# So the only freedom is in c (or equivalently, mu_0).

for c_test, label in [
    (c_SI,         "c = c_SI (reference)"),
    (c_implied,    "c = c_implied (from C4b)"),
    (c_SI * 1.01,  "c = c_SI * 1.01 (1% higher)"),
    (c_SI * 0.99,  "c = c_SI * 0.99 (1% lower)"),
]:
    alpha_test = e_charge**2 * mu_0 * c_test / (4 * math.pi * hbar)
    diff_pct = (alpha_test - alpha_CODATA) / alpha_CODATA * 100
    print(f"  {label:<38} alpha = {alpha_test:.13e}  ({diff_pct:+.6f}%)")

print()
print("  The C4b formula gives alpha_C4b = 7.2973117...e-3.")
print("  The direct formula alpha = e^2*mu_0*c/(4pi*hbar) at c = c_implied")
print("  gives a DIFFERENT value because c_implied was derived from alpha.")
print("  Setting c_implied such that the formula reproduces alpha_C4b:")
print()
# Solve: alpha_C4b = e^2 * mu_0 * c_solve / (4*pi*hbar)
c_solve = alpha_C4b * 4 * math.pi * hbar / (e_charge**2 * mu_0)
print(f"  c needed to reproduce alpha_C4b = {c_solve:.4f} m/s")
print(f"  c_SI                             = {c_SI:.4f} m/s")
print(f"  difference                       = {c_solve - c_SI:+.4f} m/s")
print(f"  In sigma (Woods 1978, unc=0.2):  {(c_solve - c_SI)/0.2:+.1f} sigma")
print()
print("  Consistent substitution gives the SAME result: c must shift by ~-1677 m/s.")
print("  Substituting c everywhere does NOT reduce the divergence.")
print("  The magnitude is preserved; the sign is: c_implied is LOWER than c_SI.")

print()
print("=" * 65)
print("PART III — FUNDAMENTAL ARGUMENT (strongest)")
print("=" * 65)
print()
print("  Alpha is DIMENSIONLESS. Its value does not depend on unit choices.")
print()
print("  Independent measurements of alpha using different methods:")
print()
print("  1. Electron g-2 (Hanneke+2008, Harvard):")
print("       alpha^-1 = 137.035999084(51)")
print("       Method: QED calculation from anomalous magnetic moment.")
print("       Uses: only electron mass and magnetic moment ratio.")
print("       Sensitivity to c: NONE (pure dimensionless ratio).")
print()
print("  2. Atom interferometry / photon recoil (Parker+2018, Berkeley):")
print("       alpha^-1 = 137.035999046(27)")
print("       Method: measures h/m_Cs via photon recoil momentum.")
print("       Uses c: yes (to relate frequency to energy), but the")
print("       dimensionless alpha extracted is independent of c value.")
print()
print("  3. Quantum Hall resistance (von Klitzing constant):")
print("       R_K = h/e^2 = 25812.807... ohm  (measured, defined in 2019)")
print("       alpha = mu_0*c/(2*R_K)  -- uses c, but R_K is fixed by")
print("       QHE plateau which is a pure topological integer effect.")
print()
print("  All three methods agree on alpha to ~10^-10. If c were 1677 m/s")
print("  different, method 3 would give a DIFFERENT alpha from methods 1 and 2.")
print("  They agree. Therefore c is not wrong.")
print()

# Quantify: what would alpha_QHE be at c_solve?
R_K = 25812.807  # ohm (von Klitzing constant)
alpha_QHE_SI    = mu_0 * c_SI   / (2 * R_K)
alpha_QHE_solve = mu_0 * c_solve / (2 * R_K)
print(f"  alpha from QHE at c_SI:      {alpha_QHE_SI:.10e}  ({(alpha_QHE_SI - alpha_CODATA)/alpha_CODATA*100:+.6f}%)")
print(f"  alpha from QHE at c_implied: {alpha_QHE_solve:.10e}  ({(alpha_QHE_solve - alpha_CODATA)/alpha_CODATA*100:+.6f}%)")
print()
print("  If c were c_implied, QHE would give alpha_C4b — but g-2 and")
print("  atom interferometry would still give alpha_CODATA (no c dependence).")
print("  The three methods would DISAGREE by 1677.8 m/s worth of c — which")
print("  is 5.6 ppm. This disagreement has never been observed.")

print()
print("=" * 65)
print("PART IV — FINAL VERDICT")
print("=" * 65)
print()
print("  QUESTION: Does consistent c substitution reduce the C4b divergence?")
print("  ANSWER: NO. The magnitude is identical (~1677.8 m/s). Sign corrected.")
print()
print("  Is the c-shift interpretation ruled out?")
print("  ANSWER: YES, at two independent levels:")
print()
print("    Level 1 (historical):  |delta_c| = 1677.8 m/s is 8385 sigma from")
print("      the best pre-definition measurement (Woods 1978, +/-0.2 m/s).")
print("      Even the 1958 microwave cavity measurement (+/-100 m/s) is 16 sigma.")
print()
print("    Level 2 (fundamental): Alpha is dimensionless. Independent measurements")
print("      (g-2, atom interferometry, QHE) agree on alpha to ~10^-10 using")
print("      methods with different or no c dependence. A 5.6 ppm c shift would")
print("      cause a 5.6 ppm disagreement between QHE-derived alpha and g-2 alpha.")
print("      No such disagreement exists in the literature.")
print()
print("  CONCLUSION: The C4b residual of -0.000560% is a real physical gap.")
print("  It is Gap 1: the wave path correction epsilon = n_exact - 2 = 0.01869.")
print("  Until epsilon is derived topologically (WZW correlator, [crys1] Tool 1),")
print("  C4b is a numerically compelling conjecture, not a derivation.")
print()
print("  The c-shift line of inquiry is fully closed. No further investigation")
print("  warranted. Proceed to [crys1] holographic quasicrystal tools for Gap 1.")
print()
print("  Script: analysis/alpha/c4b_c_consistency.py")
print("  Agenda: [c1] follow-up (supersedes c4b_c_prediction.py Part V)")
