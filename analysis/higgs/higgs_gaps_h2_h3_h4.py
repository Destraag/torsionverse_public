"""
higgs_gaps_h2_h3_h4.py
=======================
Quick assessment of Gaps H2, H3, H4 before writing dedicated scripts.

H2: Branching ratios from icosahedral symmetry
H3: Yukawa coupling derivation from cell geometry
H4: Higgs width from jamming relaxation time

Run: python analysis/higgs/higgs_gaps_h2_h3_h4.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

SEP  = "=" * 65
SEP2 = "-" * 65

pi   = math.pi
phi  = (1+math.sqrt(5))/2
Rs   = math.sqrt(5)/(4*pi)
n    = 2

# ── GAP H3: Yukawa coupling from N_J ──────────────────────────────────────
print(SEP)
print("GAP H3  Yukawa coupling vs N_J")
print(SEP)
print()
print("  Yukawa: g_f = sqrt(2)*m_f/v  [SM definition, proportional to mass]")
print("  N_J: N_f = hbar_c/(m_f*L_J)  [inversely proportional to mass]")
print("  => g_f = sqrt(2)*hbar_c/(v * N_J * L_J)  [equivalent]")
print()

v = v_EW  # 246.22 GeV
lam_sub = (1-(1-2*Rs**2)/(2*(1-Rs**2)))/4

print(f"  With v = {v} GeV, L_J = {L_J:.6f} fm:")
print(f"  g_f = sqrt(2)*hbar_c / (v*L_J * N_J)")
print(f"      = sqrt(2)*{hbar_c:.4f} / ({v*1000:.1f}*{L_J:.6f} * N_J)")
coeff = math.sqrt(2)*hbar_c / (v*1000*L_J)
print(f"      = {coeff:.6f} / N_J")
print()

particles_yukawa = [
    ("top quark", 173000, 0.114),   # g_t = sqrt(2)*173000/246220 ≈ 0.995 (heavy!)
    ("b quark",   4180,   0.024),
    ("tau",       1777,   0.010),
    ("charm",     1280,   0.007),
    ("muon",       106,   0.0006),
]
print(f"  {'Particle':<12}  {'N_J':>8}  {'g_pred':>8}  {'g_SM':>8}  {'match':>8}")
print(f"  {'-'*12}  {'-'*8}  {'-'*8}  {'-'*8}  {'-'*8}")
for name, mass_MeV, g_SM in particles_yukawa:
    N = hbar_c / (mass_MeV * L_J)
    g_pred = coeff / N
    match_pct = (g_pred/g_SM - 1)*100
    print(f"  {name:<12}  {N:>8.2f}  {g_pred:>8.5f}  {g_SM:>8.5f}  {match_pct:>+7.3f}%")

print()
print("  VERDICT: g_f = coeff/N_J is EXACT by construction (circular).")
print("  g_f ∝ 1/N_J is just a restatement of mass ∝ 1/N_J.")
print("  H3 requires an INDEPENDENT derivation of why particles have")
print("  these masses from cell geometry -- that is a bigger problem.")
print("  H3 is NOT closeable from current framework alone.")
print()

# ── GAP H4: Higgs width ────────────────────────────────────────────────────────
print(SEP)
print("GAP H4  Higgs total width Gamma_H = 4.07 +/- 0.17 MeV")
print(SEP)
print()
print("  Naive: tau_relax = L_J/v_s = 1.86e-25 s vs tau_H = 1.6e-22 s (860x off)")
print()
print("  NEW CANDIDATE: Gamma_H = alpha^2 * m_H / phi")
print()

Gamma_pred = alpha**2 * m_H_pdg22 * 1000 / phi  # MeV
Gamma_SM   = 4.07    # MeV PDG 2022
Gamma_unc  = 0.17    # MeV 1-sigma

print(f"  alpha^2 = {alpha**2:.6e}")
print(f"  m_H = {m_H_pdg22*1000:.1f} MeV")
print(f"  phi = {phi:.6f}")
print(f"  alpha^2 * m_H / phi = {Gamma_pred:.4f} MeV")
print(f"  vs Gamma_H (PDG)    = {Gamma_SM:.4f} +/- {Gamma_unc:.2f} MeV")
print(f"  Deviation           = {(Gamma_pred-Gamma_SM)/Gamma_SM*100:+.2f}%")
print(f"  Significance        = {(Gamma_pred-Gamma_SM)/Gamma_unc:.2f} sigma")
print()
print(f"  Also try m_H_pred = {m_H_pred*1000:.1f} MeV:")
Gamma_pred2 = alpha**2 * m_H_pred * 1000 / phi
print(f"  alpha^2 * m_H_pred / phi = {Gamma_pred2:.4f} MeV  (dev {(Gamma_pred2-Gamma_SM)/Gamma_SM*100:+.2f}%)")
print()

# Physical argument
print(f"  Physical argument:")
print(f"  - Higgs decays as a scalar (spin-0) -> two EM vertices -> alpha^2")
print(f"  - Decay rate proportional to available energy -> m_H")
print(f"  - Icosahedral suppression by phi (inverse of inflation factor)")
print(f"  - Result: Gamma_H = alpha^2 * m_H / phi")
print()

# Check if this matches the bb contribution specifically
# H->bb dominates: Gamma(H->bb) = 3*g_b^2/(8*pi) * m_H * sqrt(1-4*m_b^2/m_H^2)
m_b = 4.180e3  # MeV
m_H_MeV = m_H_pdg22 * 1000
g_b = math.sqrt(2)*m_b/v_EW/1000  # dimensionless
beta_b = math.sqrt(1 - 4*m_b**2/m_H_MeV**2)
Gamma_bb = 3 * g_b**2 / (8*pi) * m_H_MeV * beta_b
print(f"  SM H->bb partial width: {Gamma_bb:.3f} MeV (tree level)")
print(f"  As fraction of total:   {Gamma_bb/Gamma_SM*100:.1f}% (SM says 58.1%)")
print(f"  (Discrepancy from QCD corrections ~20% to b Yukawa)")
print()
print(f"  VERDICT: Gamma_H = alpha^2*m_H/phi = {Gamma_pred:.3f} MeV matches PDG")
print(f"  to {abs(Gamma_pred-Gamma_SM)/Gamma_unc:.1f} sigma. H4 is a STRONG CANDIDATE.")
print(f"  Physical derivation needed: why phi suppression for Higgs decay?")
print()

# ── GAP H2: Branching ratios ───────────────────────────────────────────────────
print(SEP)
print("GAP H2  Branching ratios from icosahedral symmetry")
print(SEP)
print()
print("  I_h representations: A_g(1), T_1g(3), T_2g(3), G_g(4), H_g(5)")
print()
print("  Decay mode analysis:")
print(f"  H->bb:        3 colors x 1 flavor = 3 final states  (T_1g or T_2g?)")
print(f"  H->WW*:       2 W bosons (W+, W-)  = 2 final states  (not clean)")
print(f"  H->tautau:    1 lepton pair        = 1 final state   (A_g?)")
print(f"  H->ZZ*:       1 Z pair             = 1 final state   (A_g?)")
print(f"  H->cc:        3 colors x 1 flavor  = 3 final states  (T_1g or T_2g?)")
print()
print(f"  Ratio check: BR(bb)/BR(tautau) = {0.581/0.063:.1f}")
print(f"  Color ratio: 3/1 = 3  (explains only 3x, not {0.581/0.063:.1f}x)")
print(f"  Mass ratio:  (m_b/m_tau)^2 = {(4180/1777)**2:.1f}  (explains {(4180/1777)**2:.1f}x)")
print(f"  Color * mass: 3 * {(4180/1777)**2:.1f} = {3*(4180/1777)**2:.1f}  vs actual {0.581/0.063:.1f}x")
print()
print(f"  The SM branching ratios are primarily from Yukawa^2 x color,")
print(f"  not from icosahedral symmetry directly.")
print(f"  H2 requires group-theoretic decomposition of SM couplings")
print(f"  into I_h representations -- non-trivial calculation.")
print(f"  VERDICT: H2 NOT closeable in current session. Future work.")
print()

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print(f"  H2 (branching ratios):   OPEN -- requires group theory")
print(f"  H3 (Yukawa derivation):  NOT closeable -- circular in current framework")
print(f"  H4 (width):              STRONG CANDIDATE -- Gamma_H = alpha^2*m_H/phi")
print(f"                           = {Gamma_pred:.3f} MeV vs {Gamma_SM:.2f} +/- {Gamma_unc:.2f} MeV")
print(f"                           ({(Gamma_pred-Gamma_SM)/Gamma_unc:.1f} sigma)")
print(f"  H1 (lambda residual):    PINNED -- -0.15% gap, revisit after H4")
print(SEP)
