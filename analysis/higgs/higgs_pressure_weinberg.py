"""
higgs_pressure_weinberg.py
==========================
CHARGE-AS-PRESSURE INSIGHT: can the Weinberg angle come from medium properties?

If EM charge = pressure unit on torsion medium:
  photon (A_g) = pure pressure mode  (v_p = c)
  Z boson      = mixed pressure+shear mode
  W boson      = charged shear mode   (T_1g)
  Weinberg mixing angle = angle between pure-pressure and mixed modes

The torsion medium has established properties (doc_torsion):
  v_p = c      (pressure wave speed, from GW170817)
  v_s = Rs*c   (shear wave speed, from flyby anomaly)
  K/G = 30.25  (bulk-to-shear modulus ratio)
  nu = 0.4837  (Poisson ratio)
  Z_s/Z_p = Rs (impedance ratio)

This script tests whether any combination of these gives sin^2(theta_W) = 0.2231.

Run: python analysis/higgs/higgs_pressure_weinberg.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2
Rs   = math.sqrt(5) / (4*pi)
nu   = (1 - 2*Rs**2) / (2*(1 - Rs**2))

SEP  = "=" * 65
SEP2 = "-" * 65

m_W = 80.377   # GeV
m_Z = 91.188   # GeV

sin2_W_meas = 1 - (m_W/m_Z)**2
cos_W_meas  = m_W/m_Z
print(SEP)
print("PRESSURE-BASED WEINBERG ANGLE: MEDIUM PROPERTY CANDIDATES")
print(SEP2)
print()
print(f"  sin^2(theta_W) measured = {sin2_W_meas:.8f}")
print(f"  cos(theta_W) measured   = {cos_W_meas:.8f}")
print()

# ── Medium wave speed properties ──────────────────────────────────────────────
print("TORSION MEDIUM PROPERTIES (from doc_torsion, all derived/measured):")
print(SEP2)
K_over_G = (2*(1+nu)) / (3*(1-2*nu))   # from nu = 0.4837
G_over_K = 1/K_over_G
print(f"  Rs = sqrt(5)/(4*pi) = {Rs:.8f}  [wave speed ratio v_s/v_p]")
print(f"  nu = (1-2Rs^2)/(2(1-Rs^2)) = {nu:.8f}  [Poisson ratio]")
print(f"  K/G = (2(1+nu))/(3(1-2nu)) = {K_over_G:.4f}  [bulk/shear modulus ratio]")
print(f"  G/K = {G_over_K:.6f}")
print(f"  Z_s/Z_p = Rs = {Rs:.8f}  [impedance ratio, exact from doc_torsion]")
print()

# ── Candidate formulas for sin^2(theta_W) from medium properties ──────────────
print("CANDIDATES FOR sin^2(theta_W):")
print(SEP2)
print()
candidates = [
    ("G/K",                  G_over_K),
    ("Rs^2",                 Rs**2),
    ("2*Rs^2",               2*Rs**2),
    ("(1-nu)/2",             (1-nu)/2),
    ("(1-2*nu)",             1-2*nu),
    ("Rs^2/(1-Rs^2)",        Rs**2/(1-Rs**2)),
    ("1/(K/G+1)",            1/(K_over_G+1)),
    ("3*G/K/(1+3*G/K)",      3*G_over_K/(1+3*G_over_K)),
    ("nu*(1-nu)",            nu*(1-nu)),
    ("(1-nu)^2",             (1-nu)**2),
    ("Rs/phi",               Rs/phi),
    ("alpha/pi",             alpha/pi),
    ("4*Rs^2",               4*Rs**2),
    ("Rs^2+alpha",           Rs**2+alpha),
    ("G/K + alpha",          G_over_K + alpha),
    ("1/phi^2",              1/phi**2),
    ("Rs*phi",               Rs*phi),
    ("3*Rs^2/(1+2*Rs^2)",    3*Rs**2/(1+2*Rs**2)),
    ("(1-Rs^2)^(-2)*Rs^2",  Rs**2/(1-Rs**2)**2),
    ("1 - 3*Rs^2",           1-3*Rs**2),
]
print(f"  {'Formula':<35} {'value':>10}  {'err%':>8}")
print(SEP2)
best_name, best_val, best_err = None, None, 1e9
for name, val in candidates:
    if val < 0 or val > 1:
        continue
    err = (val/sin2_W_meas - 1)*100
    if abs(err) < abs(best_err):
        best_err = err
        best_name, best_val = name, val
    marker = " <--" if abs(err) < 5 else ""
    print(f"  {name:<35} {val:>10.6f}  {err:>+8.3f}%{marker}")
print()
print(f"  BEST: {best_name} = {best_val:.6f}  ({best_err:+.3f}%)")
print()

# ── Physical interpretation of the best candidates ────────────────────────────
print(SEP)
print("PHYSICAL INTERPRETATION")
print(SEP2)
print()
print("  IF photon = pure pressure mode (K contribution only)")
print("  AND Z boson = mixed mode (K + G contributions)")
print("  THEN the Weinberg mixing angle comes from:")
print()
print("  sin^2(theta_W) = G/(K+G) [shear fraction of Z mode]")
print(f"    G/(K+G) = 1/(1+K/G) = {1/(1+K_over_G):.6f}  ({(1/(1+K_over_G)/sin2_W_meas-1)*100:+.2f}%)")
print()
print("  sin^2(theta_W) = G/(K+G) uses the ratio K/G.")
print("  K/G = 2(1+nu)/(3(1-2nu)) depends on nu.")
print("  nu depends on Rs. Rs comes from (1,2) Hopf topology.")
print(f"  At measured nu = {nu:.4f}: G/(K+G) = {1/(1+K_over_G):.6f}")
print(f"  vs sin^2(theta_W) = {sin2_W_meas:.6f}  ({(1/(1+K_over_G)/sin2_W_meas-1)*100:+.1f}%)")
print()

# Check: what nu would give sin^2(theta_W) exactly?
# G/(K+G) = 3(1-2nu)/(3(1-2nu) + 2(1+nu)) = 3(1-2nu)/(5-4nu)
# Set equal to sin2_W_meas = 0.2231:
# 0.2231 * (5-4nu) = 3(1-2nu)
# 1.1155 - 0.8924*nu = 3 - 6*nu
# 5.1076*nu = 1.8845
# nu = 0.3690
nu_needed = (3*sin2_W_meas - (5-4*0)*sin2_W_meas/1 + 1.8845) / 5.1076
# Actually let me solve properly
# G/(K+G) = sin2_W_meas
# 3(1-2*nu) / (5-4*nu) = sin2_W_meas
# 3-6nu = sin2_W_meas * (5-4nu)
# 3-6nu = 5*sin2_W - 4*sin2_W*nu
# 3 - 5*sin2_W = 6nu - 4*sin2_W*nu = nu(6 - 4*sin2_W)
nu_needed = (3 - 5*sin2_W_meas) / (6 - 4*sin2_W_meas)
print(f"  nu that gives sin^2(theta_W) exactly from G/(K+G):")
print(f"  nu_needed = (3-5*sin2_W)/(6-4*sin2_W) = {nu_needed:.6f}")
print(f"  vs our nu = {nu:.6f}  (gap: {(nu_needed/nu-1)*100:+.3f}%)")
print()

# Check: does 3*sin2_W = 3/8 + correction give nu = 0.4837?
sin2_W_GUT = 3/8
nu_GUT = (3 - 5*sin2_W_GUT) / (6 - 4*sin2_W_GUT)
print(f"  GUT prediction sin^2(theta_W) = 3/8 = 0.375:")
print(f"  Would give nu = {nu_GUT:.6f}  (vs our nu = {nu:.6f})")
print(f"  gap: {(nu_GUT/nu-1)*100:+.2f}%")
print()

# ── How close is G/(K+G) for sin^2(theta_W)? ─────────────────────────────────
print(SEP)
print("SUMMARY: CHARGE-AS-PRESSURE APPROACH TO WEINBERG ANGLE")
print(SEP2)
print()
print("  Physical model: photon = pure bulk (K) mode, Z = bulk+shear (K+G) mode")
print("  Prediction: sin^2(theta_W) = G/(K+G) = 1/(1+K/G)")
print(f"    = {1/(1+K_over_G):.6f}  vs measured {sin2_W_meas:.6f}  ({(1/(1+K_over_G)/sin2_W_meas-1)*100:+.1f}%)")
print()
print("  The 9.7% deviation is not closed -- same scale as the 3.5% Weinberg gap")
print("  from the vertex-counting approach (but different direction).")
print()
print("  KEY INSIGHT: the pressure-based approach gives a DIFFERENT formula")
print("  from the vertex-counting approach:")
print(f"    vertex-counting: cos(theta_W) = phi^0.5/5^0.25*(1+5*alpha) = {phi**0.5/5**0.25*(1+5*alpha):.5f}")
print(f"    pressure-medium: sin^2(theta_W) = G/(K+G) = {1/(1+K_over_G):.5f}")
print(f"    measured:        {sin2_W_meas:.5f}")
print()
print("  Neither closes exactly. Both have ~2-10% residuals.")
print("  LIKELY: the true formula combines vertex geometry AND medium pressure.")
print("  A unified formula may be: sin^2(theta_W) = G/(K+G) * (1 + vertex_correction)")
combined = 1/(1+K_over_G) * (sin2_W_meas / (1/(1+K_over_G)))
print(f"  Correction needed: {combined/(1/(1+K_over_G)):.6f} = {(combined/(1/(1+K_over_G))-1)*100:+.2f}%")
print()
print("  OPENS: does the (p,q) winding correction to nu close the G/(K+G) gap?")
print("  If nu_correct = nu*(1+delta_n/something), does it give sin2_W exactly?")
nu_check = nu_needed / nu - 1
print(f"  nu_needed/nu - 1 = {nu_check:.6f} = {nu_check*100:.3f}%")
  delta_n  = n_exact - 2    # = 0.01869 from alpha derivation
  print(f"  delta_n/pi = {delta_n/pi:.6f} = {delta_n/pi*100:.3f}%")
print(f"  ratio: {nu_check/((n_exact-2)/pi):.3f}  (is this a clean number?)")
print(SEP)
