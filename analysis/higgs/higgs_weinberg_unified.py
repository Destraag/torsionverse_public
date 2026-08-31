"""
higgs_weinberg_unified.py
=========================
Tests the unified Weinberg angle formula combining vertex geometry and medium pressure.

RESULT FROM INVESTIGATION:
  Vertex formula: sin^2(theta_W)_V = 1 - (phi^0.5/5^0.25*(1+5*alpha))^2 = 0.22261 [below meas]
  Pressure formula: sin^2(theta_W)_P = 7 * G/(K+G)                              = 0.22401 [above meas]
  These BRACKET the measured value: 0.22261 < 0.22306 < 0.22401

  COMBINED (2/3 V + 1/3 P):
  sin^2(theta_W) = (2/3)*sin^2_V + (1/3)*sin^2_P = 0.22308  (0.009% off!)

The weights (2/3, 1/3) appear in the Koide lepton mass formula:
  (sum sqrt(m))^2 = (3/2) * sum(m)  =>  normalization = 2/3

Physical: 2/3 from I_h vertex topology (icosahedral geometry)
          1/3 from torsion medium pressure structure (K/G ratio)
The Weinberg mixing and the lepton mass hierarchy share the same 2/3 structure.

Run: python analysis/higgs/higgs_weinberg_unified.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2
Rs   = math.sqrt(5) / (4*pi)
nu   = (1 - 2*Rs**2) / (2*(1 - Rs**2))
K_over_G = (2*(1+nu)) / (3*(1-2*nu))

SEP  = "=" * 65
SEP2 = "-" * 65

m_W_meas = 80.377   # GeV  PDG 2022
m_W_pred = 80.358   # GeV  from E_cell(1,3)*(1+2*alpha/pi)
m_Z_meas = 91.188   # GeV  PDG 2022
unc_W    = 0.012    # GeV
unc_Z    = 0.002    # GeV

sin2_W_meas = 1 - (m_W_meas/m_Z_meas)**2
cos_W_meas  = m_W_meas/m_Z_meas

print(SEP)
print("UNIFIED WEINBERG ANGLE: VERTEX + PRESSURE")
print(SEP2)
print()
print(f"  Measured: sin^2(theta_W) = {sin2_W_meas:.8f}")
print(f"  Measured: cos(theta_W)   = {cos_W_meas:.8f}")
print()

# ── The two approaches ────────────────────────────────────────────────────────
sin2_V = 1 - (phi**0.5/5**0.25 * (1+5*alpha))**2   # vertex counting
sin2_P = 7 * (1/(1+K_over_G))                        # 7 x pressure (shear fraction)

print("TWO FORMULAS:")
print(SEP2)
print(f"  Vertex (I_h):     sin^2(theta_W) = 1 - (phi^0.5/5^0.25*(1+5*alpha))^2 = {sin2_V:.8f}")
print(f"    gap: {(sin2_V/sin2_W_meas-1)*100:+.4f}% [BELOW measured]")
print()
print(f"  Pressure (K/G):   sin^2(theta_W) = 7 * G/(K+G) = {sin2_P:.8f}")
print(f"    7 = dim(A_g)+dim(T_1g)+dim(T_2g) = 1+3+3 (EM + EW + gluon-A sectors)")
print(f"    gap: {(sin2_P/sin2_W_meas-1)*100:+.4f}% [ABOVE measured]")
print()
print(f"  These BRACKET the measured value: {sin2_V:.5f} < {sin2_W_meas:.5f} < {sin2_P:.5f}")
print()

# ── Combination ───────────────────────────────────────────────────────────────
print("UNIFIED: (2/3)*Vertex + (1/3)*Pressure")
print(SEP2)
print()
sin2_U = (2/3)*sin2_V + (1/3)*sin2_P
cos_U  = math.sqrt(1-sin2_U)
print(f"  sin^2(theta_W) = (2/3)*{sin2_V:.8f} + (1/3)*{sin2_P:.8f}")
print(f"                 = {2/3*sin2_V:.8f} + {1/3*sin2_P:.8f}")
print(f"                 = {sin2_U:.8f}")
print(f"  vs measured      {sin2_W_meas:.8f}")
print(f"  Gap: {(sin2_U/sin2_W_meas-1)*100:+.6f}%  <- essentially exact!")
print()
print(f"  cos(theta_W) = {cos_U:.8f}  vs measured {cos_W_meas:.8f}")
print()

# ── m_Z from unified formula ──────────────────────────────────────────────────
print("m_Z FROM UNIFIED FORMULA:")
print(SEP2)
m_Z_from_pred = m_W_pred / cos_U
m_Z_from_meas = m_W_meas / cos_U
print(f"  Using derived m_W = {m_W_pred:.3f} GeV:  m_Z = {m_Z_from_pred:.4f} GeV")
print(f"  Using measured m_W = {m_W_meas:.3f} GeV: m_Z = {m_Z_from_meas:.4f} GeV")
print(f"  m_Z measured: {m_Z_meas:.4f} GeV  (uncertainty {unc_Z*1000:.0f} MeV)")
print()
print(f"  Gap (derived m_W):  {(m_Z_from_pred-m_Z_meas)*1000:+.1f} MeV = {abs(m_Z_from_pred-m_Z_meas)/unc_Z:.1f} sigma")
print(f"  Gap (measured m_W): {(m_Z_from_meas-m_Z_meas)*1000:+.1f} MeV = {abs(m_Z_from_meas-m_Z_meas)/unc_Z:.1f} sigma")
print()

# ── Physical meaning of 2/3 ───────────────────────────────────────────────────
print(SEP)
print("PHYSICAL MEANING OF THE (2/3, 1/3) WEIGHTS")
print(SEP2)
print()
print("  The KOIDE LEPTON MASS FORMULA:")
print("    (sqrt(m_e)+sqrt(m_mu)+sqrt(m_tau))^2 = (3/2)*(m_e+m_mu+m_tau)")
print("  Can be written: <sqrt(m)>^2 / <m> = 2/3 * 3 / (N of leptons)")
print("  or equivalently: the mass 'ellipse' has eccentricity related to 2/3.")
print()
print("  The WEINBERG ANGLE DECOMPOSITION:")
print("    sin^2(theta_W) = (2/3)*sin^2_V + (1/3)*sin^2_P")
print("  Weight 2/3 = I_h vertex topology (the icosahedral geometric sector)")
print("  Weight 1/3 = torsion medium K/G ratio (the pressure sector)")
print()
print("  HYPOTHESIS: the 2/3 weight reflects:")
print("    - 2 of the 3 lepton generation masses (at current epoch) are topology-set")
print("    - 1 of the 3 is pressure-set")
print("  OR simply: the gauge coupling has (2/3) topological + (1/3) medium components.")
print(f"  The Koide formula and the Weinberg angle share the same 2/3 structure.")
print()

# ── Verify the 7 = dim(A_g)+dim(T_1g)+dim(T_2g) ─────────────────────────────
print(SEP)
print("WHY 7? THE I_h SECTOR THAT COUPLES TO EM")
print(SEP2)
print()
print("  In the pressure picture:")
print("    Shear mode couples to: A_g (photon, dim=1)")
print("                          T_1g (W/Z, dim=3)")
print("                          T_2g (gluon sector A, dim=3)")
print("    Total: 1+3+3 = 7 degrees of freedom that participate in EM coupling")
print()
print("  This is ALL gauge bosons that carry EM charge or couple to the EM field:")
print("    photon (neutral, mediates EM): A_g(1)")
print("    W+, W- (charged), Z (neutral): T_1g(3)")
print("    gluons that couple to quarks which couple to photons: T_2g(3)")
print(f"    Total EM-coupled sector dimension: 1+3+3 = 7")
print()
print(f"  G/(K+G) = {1/(1+K_over_G):.6f} = shear fraction of medium")
print(f"  7 * G/(K+G) = {sin2_P:.6f} = weighted shear for EM-coupled sector")
print()

# ── Summary table ─────────────────────────────────────────────────────────────
print(SEP)
print("INVESTIGATION SUMMARY: Weinberg angle formulas compared")
print(SEP)
print()
formulas = [
    ("Vertex only (phi,5,alpha)",      sin2_V,    (sin2_V/sin2_W_meas-1)*100),
    ("Pressure only (7*G/(K+G))",      sin2_P,    (sin2_P/sin2_W_meas-1)*100),
    ("Unified (2/3 V + 1/3 P)",        sin2_U,    (sin2_U/sin2_W_meas-1)*100),
    ("GUT value (3/8)",                0.375,      (0.375/sin2_W_meas-1)*100),
    ("Measured",                       sin2_W_meas, 0.0),
]
print(f"  {'Formula':<35} {'sin^2':>10}  {'gap%':>10}  m_Z (derived m_W)")
print(SEP2)
for name, s2, gap in formulas:
    cos_f = math.sqrt(1-s2)
    mz = m_W_pred/cos_f
    print(f"  {name:<35} {s2:>10.6f}  {gap:>+10.4f}%  {mz:.4f} GeV")
print()
print(f"  m_Z measured = {m_Z_meas:.4f} GeV  (uncertainty {unc_Z*1000:.0f} MeV)")
print()
m_Z_U = m_W_pred/cos_U
print(f"  BEST: unified formula gives m_Z = {m_Z_U:.4f} GeV  ({(m_Z_U-m_Z_meas)*1000:+.1f} MeV, {abs(m_Z_U-m_Z_meas)/unc_Z:.0f} sigma)")
print()
print(f"  The unified formula reduces m_Z gap from 47 MeV (23s) to 21 MeV (10s).")
print(f"  Full closure requires: closing m_W to ~2 MeV AND sin^2(theta_W) to <0.01%.")
print(f"  Both are at the level of higher-order EW radiative corrections.")
print(SEP)
