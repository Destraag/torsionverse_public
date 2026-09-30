"""
weinberg_retraction_downstream_recompute.py
=============================================
After retracting the n*alpha^2*phi^2 Weinberg-angle two-loop term (see
notes/open_items/false_positive_scan_series1.txt [JC-WEINBERG-TUNED]),
recomputes the downstream quantities that used it: m_Z (via cos(theta_W))
and the SU(2) coupling g (via sin(theta_W)), using the honest one-loop
sin^2(theta_W) only.

Run: python analysis/higgs/weinberg_retraction_downstream_recompute.py
"""
import math

phi = (1 + 5 ** 0.5) / 2
alpha = 7.2973525693e-3
pi = math.pi
sqrt5 = math.sqrt(5)

SEP = "=" * 70

print(SEP)
print("DOWNSTREAM RECOMPUTE AFTER WEINBERG TWO-LOOP RETRACTION")
print(SEP)

# ── One-loop Weinberg angle (unchanged, honestly derived, GAP C) ────────────
cos_tw = math.sqrt(phi / sqrt5) * (1 + 5 * alpha)
sin2_tw = 1 - cos_tw ** 2
sin_tw = math.sqrt(sin2_tw)
print(f"\ncos(theta_W) = sqrt(phi/sqrt5)*(1+5*alpha) = {cos_tw:.10f}")
print(f"sin^2(theta_W) [one-loop only] = {sin2_tw:.10f}")
print(f"sin(theta_W)   [one-loop only] = {sin_tw:.10f}")

sin2_pdg = 0.222900
print(f"PDG sin^2(theta_W) = {sin2_pdg}")
print(f"Residual: {(sin2_tw - sin2_pdg)/sin2_pdg*100:+.4f}%  "
      f"({(sin2_tw - sin2_pdg)/3e-4:+.2f} sigma, using PDG unc ~3e-4)")

# ── m_Z from one-loop angle + PDG m_W ────────────────────────────────────────
mW_pdg = 80.3799
mZ_pdg = 91.1876
mZ_from_pdg_mW = mW_pdg / cos_tw
print(f"\nm_Z = m_W_PDG/cos(theta_W)_one-loop = {mZ_from_pdg_mW:.4f} GeV "
      f"(PDG {mZ_pdg}, gap {(mZ_from_pdg_mW-mZ_pdg)*1000:+.1f} MeV)")

# torsionverse-predicted m_W from Section 5a (g*v/2, using the same inputs
# doc_higgs.txt cites: m_W ~ 80.4 GeV conditional prediction)
mW_tv = 80.4  # as stated in doc_higgs.txt Section 5a / DOES NOT CLAIM
mZ_from_tv_mW = mW_tv / cos_tw
print(f"m_Z = m_W_torsionverse/cos(theta_W)_one-loop = {mZ_from_tv_mW:.4f} GeV "
      f"(gap {(mZ_from_tv_mW-mZ_pdg)*1000:+.1f} MeV)")

# ── g (SU(2) coupling) from one-loop sin(theta_W) ────────────────────────────
g_pred = math.sqrt(4 * pi * alpha) / sin_tw
v_ew = 246.2196
mW_from_g = g_pred * v_ew / 2
g_measured = 2 * mW_pdg / v_ew
print(f"\ng = sqrt(4*pi*alpha)/sin(theta_W)_one-loop = {g_pred:.6f}")
print(f"g (measured, 2*m_W_PDG/v) = {g_measured:.6f}")
print(f"Residual: {(g_pred-g_measured)/g_measured*100:+.4f}%")

# RG-run alpha_em(m_Z) variant, as doc_jobson_cell.txt Section 7 also quotes
alpha_mZ = 1/127.955  # standard RG-run value at m_Z scale, PDG
sin_tw_for_alt = sin_tw
g_pred_rg = math.sqrt(4 * pi * alpha_mZ) / sin_tw_for_alt
print(f"g [using RG-run alpha_em(m_Z)={alpha_mZ:.6f}] = {g_pred_rg:.6f}  "
      f"({(g_pred_rg-g_measured)/g_measured*100:+.4f}%)")

print(f"\n{SEP}")
print("SUMMARY (values to use in place of retracted two-loop-based figures)")
print(SEP)
print(f"  sin^2(theta_W) [one-loop]     = {sin2_tw:.6f}  ({(sin2_tw-sin2_pdg)/sin2_pdg*100:+.3f}%, -0.91 sigma)")
print(f"  m_Z [from PDG m_W]            = {mZ_from_pdg_mW:.4f} GeV  ({(mZ_from_pdg_mW-mZ_pdg)*1000:+.1f} MeV)")
print(f"  m_Z [from torsionverse m_W]   = {mZ_from_tv_mW:.4f} GeV  ({(mZ_from_tv_mW-mZ_pdg)*1000:+.1f} MeV)")
print(f"  g [from alpha(0), one-loop]    = {g_pred:.6f}  ({(g_pred-g_measured)/g_measured*100:+.2f}%)")
print(f"  g [from alpha(m_Z), one-loop]  = {g_pred_rg:.6f}  ({(g_pred_rg-g_measured)/g_measured*100:+.2f}%)")
