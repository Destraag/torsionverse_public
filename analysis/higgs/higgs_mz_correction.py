"""
higgs_mz_correction.py
=======================
Investigates whether the two-loop alpha^2*phi^2 correction (derived for m_H)
also improves the m_Z prediction.

m_Z = m_W / cos(theta_W)

With:
  m_W = E_cell(1,3) * (1 + 2*alpha/pi)           [spin-1, 1-loop]
  m_W* = E_cell(1,3) * (1 + 2*alpha/pi + alpha^2*phi^2) [spin-1, 2-loop]
  cos(theta_W) = sqrt(phi/sqrt(5)) * (1+5*alpha)  [derived, GAP C]

Run: python analysis/higgs/higgs_mz_correction.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import alpha, phi, hbar_c, L_J, E_cell_GeV

pi    = math.pi
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4*pi)

# PDG reference values
m_W_pdg  = 80.3799   # GeV PDG 2022
m_Z_pdg  = 91.1876   # GeV PDG 2022
m_W_unc  = 0.012     # GeV 1-sigma
m_Z_unc  = 0.0021    # GeV 1-sigma (very precise)

SEP  = "=" * 70
SEP2 = "-" * 70

# Our Weinberg angle (GAP C closed)
cos_tw   = math.sqrt(phi/sqrt5) * (1 + 5*alpha)
sin2_tw  = 1 - cos_tw**2

# E_cell for (1,3) winding -- from higgs_resonance_pq.py
# (1,3): m_W = E_cell(1,3)*(1+2*alpha/pi), gives 80.358 GeV
# Solve backwards: E_cell(1,3) = m_W_pred / (1+2*alpha/pi)
# We use the established value from session:
m_W_base  = 80.377  # GeV, from (1,3) winding + 2*alpha/pi (higgs_pq_spin_correction.py)
E_cell_13 = m_W_base / (1 + 2*alpha/pi)

# Two-loop correction: same alpha^2*phi^2 as for m_H (T_1g x T_1g -> A_g, Born^2)
c2 = alpha**2 * phi**2

print(SEP)
print("m_Z PREDICTION WITH TWO-LOOP CORRECTION")
print(SEP)
print(f"  E_cell(1,3) = {E_cell_13:.6f} GeV  [from (1,3) winding]")
print(f"  cos(theta_W) = {cos_tw:.8f}  [GAP C: sqrt(phi/sqrt5)*(1+5*alpha)]")
print(f"  alpha^2*phi^2 = {c2:.6e}  [two-loop Born correction]")
print()

# m_W predictions
m_W_1loop  = E_cell_13 * (1 + 2*alpha/pi)
m_W_2loop  = E_cell_13 * (1 + 2*alpha/pi + c2)

# m_Z predictions
m_Z_1loop  = m_W_1loop / cos_tw
m_Z_2loop  = m_W_2loop / cos_tw

print(f"  {'Quantity':<20} {'1-loop':>12} {'2-loop':>12} {'PDG':>12}")
print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*12}")
print(f"  {'m_W [GeV]':<20} {m_W_1loop:12.6f} {m_W_2loop:12.6f} {m_W_pdg:12.6f}")
print(f"  {'m_W gap [MeV]':<20} {(m_W_1loop-m_W_pdg)*1000:+12.2f} {(m_W_2loop-m_W_pdg)*1000:+12.2f} {'0':>12}")
print(f"  {'m_W sigma':<20} {(m_W_1loop-m_W_pdg)/m_W_unc:+12.2f} {(m_W_2loop-m_W_pdg)/m_W_unc:+12.2f} {'0':>12}")
print(f"  {'m_Z [GeV]':<20} {m_Z_1loop:12.6f} {m_Z_2loop:12.6f} {m_Z_pdg:12.6f}")
print(f"  {'m_Z gap [MeV]':<20} {(m_Z_1loop-m_Z_pdg)*1000:+12.2f} {(m_Z_2loop-m_Z_pdg)*1000:+12.2f} {'0':>12}")
print(f"  {'m_Z sigma':<20} {(m_Z_1loop-m_Z_pdg)/m_Z_unc:+12.2f} {(m_Z_2loop-m_Z_pdg)/m_Z_unc:+12.2f} {'0':>12}")
print()

# Also try: what if Weinberg angle gets same two-loop correction?
# cos(theta_W)* = cos(theta_W) * (1 + alpha^2*phi^2)?
cos_tw2 = cos_tw * (1 + c2)
m_Z_both = m_W_2loop / cos_tw2
print(f"  With 2-loop on BOTH m_W and theta_W:")
print(f"  m_Z = {m_Z_both:.6f} GeV  gap = {(m_Z_both-m_Z_pdg)*1000:+.2f} MeV  ({(m_Z_both-m_Z_pdg)/m_Z_unc:+.1f} sigma)")
print()
print(f"  STATUS: Two-loop correction to m_W narrows gap from")
print(f"  {(m_W_1loop-m_W_pdg)*1000:+.1f} MeV to {(m_W_2loop-m_W_pdg)*1000:+.1f} MeV.")
print(f"  m_Z gap: {(m_Z_1loop-m_Z_pdg)*1000:+.1f} MeV -> {(m_Z_2loop-m_Z_pdg)*1000:+.1f} MeV")
print(f"  m_Z remains {abs((m_Z_2loop-m_Z_pdg)/m_Z_unc):.1f} sigma -- not yet closed.")
print(f"  Requires improvement in both m_W and theta_W predictions simultaneously.")
