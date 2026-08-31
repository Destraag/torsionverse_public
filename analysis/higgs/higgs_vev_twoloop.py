"""
higgs_vev_twoloop.py
====================
Investigates the two-loop vev lead (R9): does adding alpha^2*phi^2 to the
spin correction close the vev gap, and is there a third term that closes
the remaining 0.66 MeV residual?

FINDING (R9, 2026-08-20):
  m_H = E_cell*(1 + alpha/pi + alpha^2*phi^2)  =>  v = 246.219 GeV  (gap -0.66 MeV)
  vs baseline:  m_H = E_cell*(1 + alpha/pi)     =>  v = 246.185 GeV  (gap -35 MeV)

QUESTION: what closes the remaining 0.66 MeV?
  - Is it alpha^3*phi^4 (continuing the pattern)?
  - Is it the full geometric series sum alpha^2*phi^2/(1-alpha*phi^2)?
  - Is it a known correction from elsewhere in the framework?

Run: python analysis/higgs/higgs_vev_twoloop.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import alpha, E_cell_GeV, phi, hbar_c, L_J

pi    = math.pi
Rs    = math.sqrt(5) / (4*pi)
nu    = (1 - 2*Rs**2) / (2*(1 - Rs**2))
lam   = (1 - nu) / 4
v_EW  = 246.2196   # GeV  from G_F (CODATA-2018 exact: 246.21965...)

n_exact = 2.01869   # from alpha derivation Gap 1
delta_n = n_exact - 2

SEP  = "=" * 70
SEP2 = "-" * 70

def vev(mH):     return mH / math.sqrt(2*lam)
def gap_MeV(v):  return (v - v_EW) * 1000
def sigma_v(v):  return gap_MeV(v) / (110 / math.sqrt(2*lam))  # prop from m_H unc

print(SEP)
print("TWO-LOOP VEV INVESTIGATION  --  doc_higgs R9")
print(SEP)
print(f"  E_cell = {E_cell_GeV:.9f} GeV")
print(f"  lambda = (1-nu)/4 = {lam:.10f}")
print(f"  v_EW   = {v_EW:.6f} GeV  (from G_F, CODATA-2018)")
print(f"  phi    = {phi:.10f}")
print(f"  phi^2  = {phi**2:.10f}  (= phi+1 exactly)")
print(f"  alpha*phi^2 = {alpha*phi**2:.8f}")
print()

print(SEP)
print("STEP 1  Baseline and one-loop")
print(SEP2)
mH_0 = E_cell_GeV * (1 + alpha/pi)
v0   = vev(mH_0)
print(f"  Baseline (1+a/pi):          m_H={mH_0:.6f}  v={v0:.6f}  gap={gap_MeV(v0):+.2f} MeV")

print()
print(SEP)
print("STEP 2  Two-loop: adding alpha^2*phi^2")
print(SEP2)
c2    = alpha**2 * phi**2
mH_2  = E_cell_GeV * (1 + alpha/pi + c2)
v2    = vev(mH_2)
print(f"  + alpha^2*phi^2 = {c2:.8e}")
print(f"  m_H = {mH_2:.9f} GeV")
print(f"  v   = {v2:.9f} GeV")
print(f"  gap = {gap_MeV(v2):+.4f} MeV   ({gap_MeV(v2)/gap_MeV(v0)*100:.2f}% of baseline)")
print()

print(SEP)
print("STEP 3  Third term: alpha^3*phi^4 (pattern continuation)")
print(SEP2)
c3     = alpha**3 * phi**4
mH_3   = E_cell_GeV * (1 + alpha/pi + c2 + c3)
v3     = vev(mH_3)
print(f"  phi^4 = {phi**4:.8f}  (= 3*phi+2 = {3*phi+2:.8f})")
print(f"  alpha^3*phi^4 = {c3:.8e}")
print(f"  Ratio alpha^3*phi^4 / alpha^2*phi^2 = alpha*phi^2 = {alpha*phi**2:.6f}")
print(f"  m_H = {mH_3:.9f} GeV")
print(f"  v   = {v3:.9f} GeV")
print(f"  gap = {gap_MeV(v3):+.4f} MeV")
print()

print(SEP)
print("STEP 4  Geometric series: alpha^2*phi^2 * 1/(1 - alpha*phi^2)")
print(SEP2)
ratio = alpha * phi**2   # common ratio of the geometric series
geo_sum = c2 / (1 - ratio)
mH_geo = E_cell_GeV * (1 + alpha/pi + geo_sum)
v_geo  = vev(mH_geo)
print(f"  Common ratio r = alpha*phi^2 = {ratio:.8f}")
print(f"  Sum = alpha^2*phi^2 / (1-alpha*phi^2) = {geo_sum:.8e}")
print(f"  = alpha^2*phi^2 + alpha^3*phi^4 + alpha^4*phi^6 + ...")
print(f"  m_H = {mH_geo:.9f} GeV")
print(f"  v   = {v_geo:.9f} GeV")
print(f"  gap = {gap_MeV(v_geo):+.4f} MeV")
print()

print(SEP)
print("STEP 5  n_exact correction COMBINED with alpha^2*phi^2")
print(SEP2)
# n_exact shifts the alpha/pi coefficient to (n_exact/2)*alpha/pi
mH_ne_c2 = E_cell_GeV * (1 + (n_exact/2)*alpha/pi + c2)
v_ne_c2  = vev(mH_ne_c2)
print(f"  n_exact/2 * alpha/pi + alpha^2*phi^2")
print(f"  m_H = {mH_ne_c2:.9f} GeV")
print(f"  v   = {v_ne_c2:.9f} GeV")
print(f"  gap = {gap_MeV(v_ne_c2):+.4f} MeV")
print()

print(SEP)
print("STEP 6  Claim 8 pattern: does alpha^2*phi^2 / alpha/pi = alpha*pi*phi^2?")
print(SEP2)
# Claim 8: 1 + alpha + alpha^2*phi  (ratio alpha*phi between terms)
# m_H series: alpha/pi + alpha^2*phi^2 + ?
# If ratio between m_H terms is also alpha*phi^2:
ratio_m8  = alpha * phi          # Claim 8 ratio
ratio_mH  = c2 / (alpha/pi)      # actual ratio of term2/term1 in m_H series
print(f"  Claim 8 term ratio: alpha*phi    = {ratio_m8:.8f}")
print(f"  m_H term2/term1 ratio:           = {ratio_mH:.8f}")
print(f"  = alpha*phi^2*pi = alpha*pi*phi^2 = {alpha*pi*phi**2:.8f}")
print(f"  Pattern: m_H ratio = Claim-8 ratio * pi * phi = {ratio_m8*pi*phi:.8f}")
print()

print(SEP)
print("STEP 7  What delta closes gap exactly (target = 0)?")
print(SEP2)
v_target = v_EW
mH_target = v_target * math.sqrt(2*lam)
delta_total = mH_target/E_cell_GeV - 1
delta_remaining = delta_total - (alpha/pi + c2)
print(f"  m_H needed for exact closure:  {mH_target:.9f} GeV")
print(f"  Correction needed (total):     {delta_total:.10f}")
print(f"  Already have (a/pi + a^2*phi^2):{alpha/pi + c2:.10f}")
print(f"  Remaining gap correction:       {delta_remaining:.4e}")
print()
print(f"  Compare to:")
print(f"    alpha^3*phi^4     = {alpha**3*phi**4:.4e}  ratio={delta_remaining/(alpha**3*phi**4):.4f}")
print(f"    alpha^3*phi^3     = {alpha**3*phi**3:.4e}  ratio={delta_remaining/(alpha**3*phi**3):.4f}")
print(f"    alpha^3*phi^2     = {alpha**3*phi**2:.4e}  ratio={delta_remaining/(alpha**3*phi**2):.4f}")
print(f"    alpha^2*delta_n/pi= {alpha**2*delta_n/pi:.4e}  ratio={delta_remaining/(alpha**2*delta_n/pi):.4f}")
print(f"    alpha*delta_n/(pi^2)={alpha*delta_n/pi**2:.4e}  ratio={delta_remaining/(alpha*delta_n/pi**2):.4f}")
print(f"    delta_n*alpha^2/pi= {delta_n*alpha**2/pi:.4e}  ratio={delta_remaining/(delta_n*alpha**2/pi):.4f}")
print(f"    Rs*alpha^2        = {Rs*alpha**2:.4e}  ratio={delta_remaining/(Rs*alpha**2):.4f}")
print(f"    alpha^2/(2*pi)    = {alpha**2/(2*pi):.4e}  ratio={delta_remaining/(alpha**2/(2*pi)):.4f}")
print()

print(SEP)
print("SUMMARY")
print(SEP2)
results = [
    ("Baseline  1+a/pi",                       mH_0,     v0),
    ("+ a^2*phi^2",                             mH_2,     v2),
    ("+ a^2*phi^2 + a^3*phi^4",                 mH_3,     v3),
    ("+ geo sum a^2*phi^2/(1-a*phi^2)",         mH_geo,   v_geo),
    ("n_exact + a^2*phi^2",                     mH_ne_c2, v_ne_c2),
]
for label, mH, v in results:
    print(f"  {label:<40} v={v:.6f}  gap={gap_MeV(v):+7.2f} MeV")
print()
print(f"  v_EW (G_F) = {v_EW:.6f} GeV")
print(f"  Free parameters in all above: 0")
