"""
higgs_h1_residual.py
====================
Investigates the remaining -0.15% gap in Gap H1:
  lambda_sub = (1-nu)/4 = 0.12909 vs lambda_SM = 0.12928  (-0.149%)

APPROACH:
  For bulk particles (electron): l=0 AND l=6 channels contribute -> L3 mean
  For sub-cell particles (Higgs): l=6 suppressed, only l=0 contributes
  The l=0 channel gives coupling f1 = PHI (golden ratio load path)
  A small l=0 vertex correction to lambda_sub may close the -0.15% gap.

LEADS:
  1. l=0 channel stiffness: lambda = lambda_sub * (1 + PHI * delta_k_l0)
  2. Alpha/pi correction propagated: lambda_sub * (1 + alpha/pi)?
  3. 1/N_J correction: lambda_sub * (1 + c * N_J)?

Run: python analysis/higgs/higgs_h1_residual.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

SEP  = "=" * 65
SEP2 = "-" * 65

Rs  = math.sqrt(5) / (4*math.pi)
nu  = (1 - 2*Rs**2) / (2*(1 - Rs**2))
L3  = (phi**3 + math.log(5)**3) / (phi**2 + math.log(5)**2)
n   = 2

lam_sub = (1 - nu) / 4   # 0.12909
lam_SM  = lam_SM          # 0.12928 from constants.py
gap_pct = (lam_SM / lam_sub - 1) * 100   # +0.149% needed

N_J_H = (hbar_c / (m_H_pdg22 * 1000)) / L_J  # 0.159 (sub-cell)

print(SEP)
print("GAP H1 RESIDUAL INVESTIGATION")
print(f"lambda_sub = {lam_sub:.8f}")
print(f"lambda_SM  = {lam_SM:.8f}")
print(f"Gap needed = {gap_pct:+.4f}%")
print(SEP)
print()

# ── Lead 1: l=0 channel stiffness increment ───────────────────────────────────
print("LEAD 1  l=0 channel vertex stiffness for sub-cell particles")
print(SEP2)
print()
print(f"  For bulk particles (N >> 1):")
print(f"    Both l=0 and l=6 channels active -> f_eff = L3(PHI,log5) = {L3:.6f}")
print(f"    Full correction: delta_n = L3 * delta_k -> n_exact = 2.01869")
print()
print(f"  For sub-cell particles (N < 1, Higgs N_J = {N_J_H:.4f}):")
print(f"    l=6 suppressed by (L_J/R_H)^6 = {(L_J/(hbar_c/(m_H_pdg22*1000)))**6:.2e}")
print(f"    Only l=0 channel active -> f_eff = f1 = PHI = {phi:.6f}")
print()
print(f"  IF lambda corrected by l=0 channel: lambda = lambda_sub * (1 + PHI * delta_k_l0)")
print(f"  delta_k_l0 needed = gap% / (PHI * 100) = {gap_pct/100:.6f} / {phi:.6f}")
delta_k_l0 = gap_pct / (phi * 100)
print(f"                    = {delta_k_l0:.8f}")
print()
delta_k_alpha = 0.01869 / L3
print(f"  Compare to alpha's delta_k = delta_n/L3 = {delta_k_alpha:.6f}")
print(f"  Ratio delta_k_l0 / delta_k_alpha = {delta_k_l0/delta_k_alpha:.6f}")
print(f"  = {delta_k_l0/delta_k_alpha:.4f}  -- is this alpha, Rs, or related?")
print(f"  alpha = {alpha:.6f}  (close? {abs(delta_k_l0/delta_k_alpha - alpha)/alpha*100:.1f}% off)")
print(f"  Rs    = {Rs:.6f}  (close? {abs(delta_k_l0/delta_k_alpha - Rs)/Rs*100:.1f}% off)")
print(f"  alpha/pi = {alpha/math.pi:.6f}")
print()

# ── Lead 2: alpha/pi correction to lambda ─────────────────────────────────────
print("LEAD 2  Does alpha/pi propagate from m_H correction to lambda?")
print(SEP2)
print()
print(f"  m_H = E_cell*(1+alpha/pi):  the scalar QED correction adds alpha/pi to m_H")
print(f"  Does the SAME correction apply to lambda_sub?")
print(f"  alpha/pi = {alpha/math.pi:.8f} = {alpha/math.pi*100:.4f}%")
lam_corrected_api = lam_sub * (1 + alpha/math.pi)
print(f"  lambda_sub * (1 + alpha/pi) = {lam_corrected_api:.8f}")
print(f"  vs lambda_SM = {lam_SM:.8f}")
print(f"  Deviation: {(lam_corrected_api/lam_SM-1)*100:+.4f}%  -- OVERCORRECTS by {abs((lam_corrected_api/lam_SM-1)*100):.3f}%")
print(f"  alpha/pi = {alpha/math.pi*100:.4f}% vs gap needed = {gap_pct:.4f}%")
print(f"  alpha/pi is {alpha/math.pi*100/gap_pct:.2f}x the needed gap -- too large by {alpha/math.pi*100/gap_pct:.2f}x")
print()

# ── Lead 3: 1/N_J correction ──────────────────────────────────────────────
print("LEAD 3  1/N_J finite-size correction for sub-cell particles")
print(SEP2)
print()
print(f"  N_J_Higgs = {N_J_H:.4f}  (very sub-cell)")
print(f"  Finite-size expansion: lambda = lambda_sub * (1 + c * N_J + ...)")
print(f"  For gap = {gap_pct:.4f}%: c = gap% / (N_J * 100) = {gap_pct/(N_J_H*100):.6f}")
c_Ngrain = gap_pct / (N_J_H * 100)
print(f"  c = {c_Ngrain:.6f}")
print(f"  Is c a known quantity?")
print(f"  Rs = {Rs:.6f}  (c/Rs = {c_Ngrain/Rs:.4f})")
print(f"  alpha = {alpha:.6f}  (c/alpha = {c_Ngrain/alpha:.4f})")
print(f"  Rs^2  = {Rs**2:.6f}  (c/Rs^2 = {c_Ngrain/Rs**2:.4f})")
print(f"  alpha/pi = {alpha/math.pi:.6f}  (close? {abs(c_Ngrain-alpha/math.pi)/(alpha/math.pi)*100:.1f}% off)")
print()

# ── Lead 4: exact zero-correction check ───────────────────────────────────────
print("LEAD 4  Is the gap within PDG measurement uncertainty?")
print(SEP2)
print()
print(f"  lambda_SM comes from m_H = {m_H_pdg22} +/- {m_H_pdg_unc} GeV and v_EW = {v_EW} GeV")
print(f"  delta_lambda/lambda = 2*delta_m_H/m_H = 2*{m_H_pdg_unc}/{m_H_pdg22} = {2*m_H_pdg_unc/m_H_pdg22*100:.3f}%")
delta_lam_meas = 2 * m_H_pdg_unc / m_H_pdg22 * 100
print(f"  Measurement uncertainty in lambda: {delta_lam_meas:.3f}%")
print(f"  Gap in lambda: {gap_pct:.4f}%")
print(f"  Gap / measurement unc = {gap_pct/delta_lam_meas:.2f}x")
print()
print(f"  The gap ({gap_pct:.3f}%) is {gap_pct/delta_lam_meas:.1f}x the measurement uncertainty.")
print(f"  It is NOT within measurement uncertainty -- a real correction is needed.")
print()

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print(f"  Gap needed: +{gap_pct:.4f}% to lambda_sub = {lam_sub:.6f}")
print()
print(f"  Lead 1 (l=0 channel): delta_k_l0 = {delta_k_l0:.6f}")
print(f"    ratio to alpha's delta_k: {delta_k_l0/delta_k_alpha:.4f}")
print(f"    Physical: l=0 vertex correction for sub-cell particles")
print(f"    STATUS: plausible, delta_k_l0 derivation needed")
print()
print(f"  Lead 2 (alpha/pi propagation): overcorrects by {alpha/math.pi*100/gap_pct:.1f}x")
print(f"    STATUS: NOT the right correction -- too large")
print()
print(f"  Lead 3 (1/N_J): c = {c_Ngrain:.6f}, not obviously clean")
print(f"    STATUS: possible but coefficient not yet derived")
print()
print(f"  Lead 4 (measurement): gap is {gap_pct/delta_lam_meas:.1f}x larger than uncertainty")
print(f"    STATUS: real correction needed, not measurement noise")
print()
print(f"  RECOMMENDED NEXT STEP: investigate l=0 vertex stiffness increment")
print(f"  for sub-cell particles. The l=0 channel IS active for Higgs (not suppressed).")
print(f"  delta_k_l0 = {delta_k_l0:.6f} needs a geometric/topological derivation.")
print(SEP)
