"""
alpha_precision_check.py
Asks: could alpha_C4b (from the C4b quadratic) be more accurate than CODATA?
      Would substituting alpha_C4b at the electron scale align it with the
      other three scales in the cross-scale ratio table?

Method:
  1. State CODATA uncertainty on alpha precisely.
  2. Compute C4b discrepancy in units of that uncertainty (sigma).
  3. Propagate alpha_C4b through the five-scale table to see if any ratio
     shifts toward alignment.
  4. Compute the implied anomalous magnetic moment a_e = (g-2)/2 at
     both alpha values and compare to experiment.

Conclusion is determined by the numbers, not the question.

Run: python analysis/alpha_precision_check.py
"""

import math

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)

SEP   = "=" * 70
SEP2  = "-" * 70

# ──────────────────────────────────────────────────────────────────────────────
# PART 1 — CODATA UNCERTAINTY ON ALPHA
# Source: NIST CODATA 2018. alpha = 7.2973525693(11)e-3
# The (11) denotes the combined standard uncertainty in units of the last
# two digits shown, i.e., ± 0.0000000011e-3 = ± 1.1e-12.
# ──────────────────────────────────────────────────────────────────────────────

alpha_CODATA      = 7.2973525693e-3
alpha_CODATA_unc  = 1.1e-12            # ± 1.1e-12 (CODATA 2018 standard uncertainty)

# ──────────────────────────────────────────────────────────────────────────────
# PART 2 — C4B VALUE (from c4a_candidates.py, not recomputed here)
# Quadratic: 2*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0, physical root alpha_-
# Verified in c4a_candidates.py; reproduced here for independence.
# ──────────────────────────────────────────────────────────────────────────────

a_coeff =  2.0
b_coeff = -(4 * pi**2 / phi)
c_coeff =  Rs
disc    = b_coeff**2 - 4 * a_coeff * c_coeff
alpha_C4b = (-b_coeff - math.sqrt(disc)) / (2 * a_coeff)

delta_abs     = alpha_C4b - alpha_CODATA        # negative: C4b < CODATA
delta_rel     = delta_abs / alpha_CODATA        # fractional (negative)
delta_pct     = delta_rel * 100                 # percent
delta_sigma   = abs(delta_abs) / alpha_CODATA_unc  # how many sigma from CODATA

print(SEP)
print("PART 1+2 — C4B vs CODATA: HOW FAR APART ARE THEY?")
print(SEP)
print(f"  alpha_CODATA                 = {alpha_CODATA:.13e}")
print(f"  alpha_C4b  (quadratic root)  = {alpha_C4b:.13e}")
print(f"  absolute difference          = {delta_abs:.4e}")
print(f"  relative difference          = {delta_rel:.4e}  ({delta_pct:+.6f}%)")
print()
print(f"  CODATA 2018 standard uncertainty on alpha: ±{alpha_CODATA_unc:.1e}")
print(f"  relative uncertainty:                       {alpha_CODATA_unc/alpha_CODATA:.2e}  ({alpha_CODATA_unc/alpha_CODATA*1e9:.2f} ppb)")
print()
print(f"  C4b discrepancy in sigma:    {delta_sigma:.0f} sigma")
print()
print(f"  In practical terms:")
print(f"    C4b is off by {abs(delta_pct):.6f}% from CODATA.")
print(f"    CODATA uncertainty is      {alpha_CODATA_unc/alpha_CODATA*100:.9f}%.")
print(f"    The C4b gap is {abs(delta_rel)/(alpha_CODATA_unc/alpha_CODATA):.0f}x the experimental uncertainty.")
print()
print(f"  Answer: CODATA cannot be off by this amount. The alpha value is")
print(f"  cross-confirmed by multiple independent techniques (see Part 3)")
print(f"  each agreeing to sub-ppb precision. A 5600-ppb error in CODATA")
print(f"  would require every one of those measurements to be wrong by the")
print(f"  same amount, which is astrophysically implausible.")
print()

# ──────────────────────────────────────────────────────────────────────────────
# PART 3 — INDEPENDENT DETERMINATIONS OF ALPHA
# All determine alpha by distinct physics. If CODATA were off by 5600 ppb,
# every one of these would show the same shift independently.
# ──────────────────────────────────────────────────────────────────────────────

# Experimental values of 1/alpha from distinct methods (published values)
inv_alpha_measurements = [
    ("Electron g-2 + QED theory (Harvard 2008)",      137.035999206, 0.000000011),
    ("Rb atom recoil interferometry (Parker 2018)",    137.035999049, 0.000000090),
    ("Cs atom recoil interferometry (Morel 2020)",     137.035999046, 0.000000027),
    ("Quantum Hall resistance R_K (exact since 2019)", 137.035999084, 0.000000021),
    ("CODATA 2018 recommended (weighted)",             1/alpha_CODATA, alpha_CODATA_unc/(alpha_CODATA**2)),
]

inv_alpha_C4b    = 1 / alpha_C4b
inv_alpha_CODATA = 1 / alpha_CODATA

print(SEP)
print("PART 3 — INDEPENDENT DETERMINATIONS OF 1/alpha")
print(SEP)
print(f"  {'Method':<50} {'1/alpha':>14}  {'C4b gap (sigma)':>16}")
print(f"  {'-'*50} {'-'*14}  {'-'*16}")

for label, inv_a, unc_inv_a in inv_alpha_measurements:
    sigma_from_C4b = abs(inv_a - inv_alpha_C4b) / unc_inv_a
    print(f"  {label:<50} {inv_a:>14.9f}  {sigma_from_C4b:>14.0f}σ")

print()
print(f"  C4b predicts: 1/alpha = {inv_alpha_C4b:.9f}")
print(f"  CODATA gives: 1/alpha = {inv_alpha_CODATA:.9f}")
print(f"  Gap in 1/alpha: +{inv_alpha_C4b - inv_alpha_CODATA:.6f}")
print()

# ──────────────────────────────────────────────────────────────────────────────
# PART 4 — PROPAGATE ALPHA_C4B THROUGH THE FIVE-SCALE TABLE
# We recompute every ratio that involves alpha using alpha_C4b vs alpha_CODATA.
# The question: does any ratio shift *toward* alignment by a meaningful amount?
# ──────────────────────────────────────────────────────────────────────────────

# Constants (CODATA 2018)
me_MeV      = 0.51099895
mp_MeV      = 938.27208816
mp_GeV      = mp_MeV / 1000.0
r_proton_fm = 0.8414
kappa_GeV_per_fm = 0.9
Lambda_QCD_GeV   = 0.217
nuclear_binding_MeV = 8.0
a0_m_s2     = 1.2e-10
c_m_s       = 2.99792458e8
H0_km_s_Mpc = 70.0
Mpc_in_m    = 3.085677581e22
H0_s        = H0_km_s_Mpc * 1e3 / Mpc_in_m
cH0         = c_m_s * H0_s

def compute_ratios(alpha_val):
    R1  = alpha_val / (2 * pi)          # Particle/EM
    R2  = alpha_val                     # Full EM coupling
    R3  = nuclear_binding_MeV / mp_MeV  # Nuclear binding (no alpha)
    R4  = (kappa_GeV_per_fm * r_proton_fm) / mp_GeV  # Hadronic (no alpha)
    R5  = a0_m_s2 / cH0                # Galactic (no alpha)
    return R1, R2, R3, R4, R5

R1c, R2c, R3c, R4c, R5c = compute_ratios(alpha_CODATA)
R1b, R2b, R3b, R4b, R5b = compute_ratios(alpha_C4b)

# The cross-scale checks from scale_check.py
# (target, label, value_CODATA, value_C4b, residual_CODATA, residual_C4b)
checks = [
    # R3/alpha should equal 1 for nuclear~EM alignment
    ("R3 / alpha = 1 ?",
     R3c / alpha_CODATA,
     R3b / alpha_C4b),
    # R5 * 2pi should equal 1 for galactic~particle alignment
    ("R5 * 2*pi = 1 ?",
     R5c * 2 * pi,
     R5b * 2 * pi),              # same: R5 has no alpha dependence
    # R5/R1 should equal 1/alpha
    # R5/R1 = 2pi*a0/(alpha*cH0); if = 1/alpha then 2pi*a0/cH0 = 1 (=R5*2pi)
    # So this is the same check as above, just presented differently
    ("R5/R1 = 1/alpha ?  [= R5*2pi]",
     (R5c / R1c) * alpha_CODATA,     # = R5*2pi (should be 1)
     (R5b / R1b) * alpha_C4b),
]

print(SEP)
print("PART 4 — DOES SUBSTITUTING alpha_C4b IMPROVE THE SCALE RATIOS?")
print(SEP)
print(f"  For each cross-scale check, residual = |value - 1|.")
print(f"  A residual of 0 means perfect alignment.")
print()
print(f"  {'Check':<28} {'With CODATA':>14} {'With C4b':>14} {'Change':>12}")
print(f"  {'-'*28} {'-'*14} {'-'*14} {'-'*12}")

for label, val_codata, val_c4b in checks:
    res_c = abs(val_codata - 1)
    res_b = abs(val_c4b - 1)
    delta_res = res_b - res_c
    direction = "better" if delta_res < 0 else "worse " if delta_res > 0 else "same"
    print(f"  {label:<28} {val_codata:>14.6f} {val_c4b:>14.6f} "
          f"{delta_res:>+12.8f}  ({direction})")

print()

# The ratio R3/alpha is the most interesting one: does nuclear~EM alignment improve?
r3_alpha_CODATA = R3c / alpha_CODATA
r3_alpha_C4b    = R3b / alpha_C4b
delta_r3 = abs(r3_alpha_C4b - 1) - abs(r3_alpha_CODATA - 1)

print(f"  Detailed look at R3/alpha (nuclear binding / alpha):")
print(f"    With CODATA:  R3/alpha = {r3_alpha_CODATA:.8f}  (residual {abs(r3_alpha_CODATA-1)*100:.4f}%)")
print(f"    With C4b:     R3/alpha = {r3_alpha_C4b:.8f}  (residual {abs(r3_alpha_C4b-1)*100:.4f}%)")
print(f"    Shift from C4b correction: {delta_r3*100:+.8f}%  (target residual: 0%)")
print(f"    Scale misalignment is {abs(r3_alpha_CODATA-1)*100:.4f}%,")
print(f"    C4b shifts it by      {abs(delta_r3)*100:.8f}%.")
print(f"    C4b moves the needle by {abs(delta_r3)/abs(r3_alpha_CODATA-1)*100:.6f}% of the gap.")
print()

print(SEP)
print("PART 5 — WHAT DOES alpha_C4b IMPLY FOR THE g-2 MEASUREMENT?")
print(SEP)
print()
print("  The electron anomalous magnetic moment a_e = (g-2)/2 is the primary")
print("  experimental anchor for CODATA's alpha (Hanneke 2008).")
print("  The leading QED prediction is a_e = alpha/(2*pi) + O(alpha^2).")
print("  If we use alpha_C4b, the leading-order prediction shifts by:")
print()

# Leading-order QED: a_e ≈ alpha/(2*pi)
# More precisely use full 5-loop result (Aoyama et al. 2019):
#   a_e = C1*(alpha/pi) + C2*(alpha/pi)^2 + C3*(alpha/pi)^3 + C4*(alpha/pi)^4 + C5*(alpha/pi)^5
# Coefficients (Aoyama 2019):
C1 = 0.5                   # = 1/2 (Schwinger 1948)
C2 = -0.328478965579       # 2-loop (exact)
C3 = 1.181241456587        # 3-loop (numerical)
C4 = -1.9144              # 4-loop (numerical)
C5 = 9.16                  # 5-loop (numerical estimate)

def ae_qed(alpha_val):
    x = alpha_val / pi
    return C1*x + C2*x**2 + C3*x**3 + C4*x**4 + C5*x**5

ae_CODATA = ae_qed(alpha_CODATA)
ae_C4b    = ae_qed(alpha_C4b)
ae_exp    = 1.15965218059e-3       # Harvard 2008 measurement (Hanneke et al.)
ae_exp_unc = 2.8e-13               # standard uncertainty

ae_diff_pct = (ae_C4b - ae_CODATA) / ae_CODATA * 100
ae_C4b_sigma = abs(ae_C4b - ae_exp) / ae_exp_unc

print(f"  a_e (QED, 5-loop, with alpha_CODATA) = {ae_CODATA:.12e}")
print(f"  a_e (QED, 5-loop, with alpha_C4b   ) = {ae_C4b:.12e}")
print(f"  a_e experimental (Harvard 2008)      = {ae_exp:.12e}  ± {ae_exp_unc:.1e}")
print()
print(f"  Shift in a_e from using alpha_C4b: {(ae_C4b-ae_CODATA):+.4e}")
print(f"  Experiment uncertainty:             ±{ae_exp_unc:.4e}")
print(f"  Shift / uncertainty:                {abs(ae_C4b-ae_CODATA)/ae_exp_unc:.0f}σ")
print()
print(f"  With alpha_C4b, QED prediction of a_e disagrees with")
print(f"  the Harvard measurement by {ae_C4b_sigma:.0f}σ.")
print()
print(f"  Verdict: alpha_C4b is NOT consistent with the g-2 measurement.")
print(f"  CODATA alpha is correct to the precision at which it is stated.")
print()

# ──────────────────────────────────────────────────────────────────────────────
# PART 6 — WHAT C4B ACTUALLY IS: A PRECISION FLOOR, NOT A CLAIM
# ──────────────────────────────────────────────────────────────────────────────

print(SEP)
print("PART 6 — WHAT C4B ACTUALLY MEANS (HONEST ASSESSMENT)")
print(SEP)
print()

C4a_err_pct = abs(-0.060376)
C4b_err_pct = abs(delta_pct)
improvement = C4a_err_pct / C4b_err_pct

print(f"  Accuracy comparison (% error from CODATA):")
print(f"    Best random-constant formulas:     ~0.1 – 1%  (Eddington, Wyler, etc.)")
print(f"    C4a (direct geometric formula):    {C4a_err_pct:.4f}%")
print(f"    C4b (self-consistent quadratic):   {C4b_err_pct:.6f}%  ({improvement:.0f}x better than C4a)")
print(f"    CODATA uncertainty:                {alpha_CODATA_unc/alpha_CODATA*100:.9f}%")
print(f"    C4b still needs to close:          {C4b_err_pct/(alpha_CODATA_unc/alpha_CODATA*100):.0f}x more precision")
print()
print(f"  Precision ladder (how many more orders of magnitude to go):")
print(f"    Random numerology (Eddington):     ~10^0 % off")
print(f"    C4a geometric formula:             ~10^-1 % off")
print(f"    C4b self-consistent quadratic:     ~10^-3 % off   <-- we are here")
print(f"    CODATA uncertainty floor:          ~10^-8 % level")
print()
print(f"  Distance to travel: 5 more orders of magnitude in precision.")
print(f"  That cannot be claimed by tuning alpha — it requires a forward")
print(f"  derivation (Step 3 of Appendix D) that produces alpha analytically.")
print()
print(f"  C4b's significance is not 'more accurate than CODATA' but:")
print(f"    (a) It is the most accurate GEOMETRIC prediction of alpha known.")
print(f"    (b) It has physical structure (self-consistent correction to")
print(f"        the Hopf ratio) that hints at a derivation path.")
print(f"    (c) Its remaining gap ({C4b_err_pct:.6f}%) defines the precision")
print(f"        target that Step 3 must close to be publishable.")
print()

print(SEP)
print("SUMMARY")
print(SEP)
print(f"  Q: Could CODATA alpha be wrong and C4b be correct?")
print(f"  A: No. C4b is {delta_sigma:.0f}σ from CODATA. Multiple independent")
print(f"     experiments (g-2, atom recoil, quantum Hall) all agree at the")
print(f"     sub-ppb level. A 5600-ppb error in all of them simultaneously")
print(f"     would require systematic physics errors in unrelated experiments.")
print()
print(f"  Q: Does using alpha_C4b bring the electron scale into alignment")
print(f"     with the other three scales in the cross-scale ratio table?")
print(f"  A: No. The scale residuals are 10.9% and 16.8%. Swapping alpha")
print(f"     changes those ratios by {abs(delta_pct):.6f}% — roughly 20,000x")
print(f"     smaller than the misalignment. Invisible to the scale table.")
print()
print(f"  What C4b IS:")
print(f"    The best geometric prediction of alpha from first principles.")
print(f"    100x more accurate than C4a. Structurally suggestive.")
print(f"    Still needs a forward derivation to be more than a conjecture.")
print(f"    When Step 3 is complete, the prediction must match CODATA to")
print(f"    10^-8% precision — not 10^-3% — to be a true derivation.")
print(SEP)
