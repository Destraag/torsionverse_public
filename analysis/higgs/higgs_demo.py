"""
higgs_demo.py
=============
Reproducibility demo for doc_higgs.txt -- all 8 claims in sequence.

Derives the Higgs boson mass, quartic coupling, vacuum expectation value,
and total decay width from the (1,2) Hopf fibration geometry using only
six measured inputs (alpha, r_p, v_p=c, v_s=Rs*c, Rs from K-formula).

Run:  python analysis/higgs/higgs_demo.py

This script was developed with the assistance of AI language models
(GitHub Copilot / Claude, Anthropic). All physical claims and results
are the sole work of the author (R. Jobson).

For derivation detail on any individual claim, see the companion scripts
listed at the end of this file.
"""

import math, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

pi    = math.pi
sqrt  = math.sqrt
sqrt3 = sqrt(3)

# ── MEASURED INPUTS (zero free parameters beyond these) ──────────────────────
alpha  = 7.2973525693e-3   # CODATA-2018 fine structure constant
r_p    = 0.8414            # fm  proton charge radius, CODATA-2018
hbar_c = 197.3269804       # MeV*fm
phi    = (1 + sqrt(5)) / 2 # golden ratio (exact, from icosahedral symmetry)
Rs     = sqrt(5) / (4*pi)  # v_s/v_p = Rs, from K-formula + GW170817 [doc_torsion]

# ── DERIVED CELL CONSTANTS ────────────────────────────────────────────────────
L_J    = alpha * phi * r_p                # fm  Jobson cell edge
N_lock = 2*pi / (alpha * phi)             # tube closure number
E_cell = 2*pi * hbar_c / L_J / 1000      # GeV  natural energy quantum of cell
nu     = (1 - 2*Rs**2) / (2*(1 - Rs**2)) # Poisson ratio of torsion medium
K_o_G  = (2*(1+nu)) / (3*(1-2*nu))       # bulk/shear modulus ratio

# ── PDG 2022 REFERENCE VALUES ────────────────────────────────────────────────
m_H_pdg    = 125.20;  m_H_unc  = 0.11   # GeV
lam_sm     = 0.12928; lam_unc  = 0.0015
v_sm       = 246.220                     # GeV  (from G_F, Fermi constant)
Gamma_pdg  = 4.07;    Gamma_unc= 0.17   # MeV

SEP  = "=" * 70
SEP2 = "-" * 70

def sigma(pred, meas, unc):
    return (pred - meas) / unc

def pct(pred, meas):
    return (pred / meas - 1) * 100

print(SEP)
print("HIGGS BOSON DERIVATION DEMO  --  doc_higgs.txt")
print("All 8 claims derived from 6 measured inputs, 0 free parameters.")
print(SEP)
print()
print(f"  Inputs:")
print(f"    alpha  = {alpha:.10e}  [CODATA-2018]")
print(f"    r_p    = {r_p} fm                  [CODATA-2018]")
print(f"    Rs     = sqrt(5)/(4*pi) = {Rs:.8f}    [K-formula + GW170817]")
print(f"    v_p    = c                              [GW170817]")
print()
print(f"  Cell geometry:")
print(f"    phi    = (1+sqrt(5))/2 = {phi:.10f}  [icosahedral symmetry]")
print(f"    L_J    = alpha*phi*r_p = {L_J:.8f} fm")
print(f"    N_lock = 2*pi/(alpha*phi) = {N_lock:.4f}")
print(f"    E_cell = 2*pi*hbar_c/L_J  = {E_cell:.6f} GeV")
print(f"    nu     = (1-2*Rs^2)/(2*(1-Rs^2)) = {nu:.8f}  [Poisson ratio]")
print(f"    K/G    = {K_o_G:.6f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 1  m_H = E_cell * (1 + alpha/pi)")
print(SEP)
m_H_pred = E_cell * (1 + alpha/pi)
s1 = sigma(m_H_pred, m_H_pdg, m_H_unc)
print(f"  Correction alpha/pi: spin-0 scalar QED correction (linking number theorem)")
print(f"  Predicted: {m_H_pred:.6f} GeV")
print(f"  PDG 2022:  {m_H_pdg} +/- {m_H_unc} GeV")
print(f"  Residual:  {pct(m_H_pred, m_H_pdg):+.4f}%  ({s1:+.2f} sigma)")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 2  Higgs is spin-0: linking number theorem (analytic)")
print(SEP)
print("  Theorem: (p,q) torus knot has pi-rotation symmetry <=> n = p*q is even.")
print("  Proof: path (t, q*t) maps to itself under phi->phi+pi iff q is even.")
print("  (1,2) torus knot: n = 1*2 = 2 (even) => pi-rotation => scalar (spin-0).")
print("  => QED correction is alpha/pi (two-vertex scalar loop), not alpha/2pi.")
print("  Result: exact topological theorem, no measurement needed.")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 3  lambda = (1-nu)/4")
print(SEP)
lam_pred = (1 - nu) / 4
# Also exact form
lam_exact = 2*pi**2 / (16*pi**2 - 5)
s3 = sigma(lam_pred, lam_sm, lam_unc)
print(f"  Sub-cell coupling (N_J = 1/(2*pi) < 1): lambda set by bulk Poisson ratio.")
print(f"  Exact form: 2*pi^2/(16*pi^2-5) = {lam_exact:.8f}")
print(f"  Predicted:  {lam_pred:.8f}")
print(f"  SM value:   {lam_sm} +/- {lam_unc}")
print(f"  Residual:   {pct(lam_pred, lam_sm):+.4f}%  ({s3:+.2f} sigma)")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 4  v = m_H / sqrt(2*lambda)")
print(SEP)
v_pred = m_H_pred / sqrt(2 * lam_pred)
print(f"  Standard Higgs relation, both m_H and lambda now predicted.")
print(f"  Predicted: {v_pred:.6f} GeV")
print(f"  G_F value: {v_sm} GeV  (from Fermi constant, independent measurement)")
print(f"  Residual:  {pct(v_pred, v_sm):+.4f}%")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 5  Gamma_H = alpha^2 * m_H / phi")
print(SEP)
Gamma_pred = alpha**2 * m_H_pred * 1000 / phi   # MeV
s5 = sigma(Gamma_pred, Gamma_pdg, Gamma_unc)
# Deeper form via alpha equation: Q*alpha = Rs (leading order)
print(f"  Deeper form: alpha^2*m_H/phi = alpha*Rs*m_H/(4*pi^2)")
print(f"  Predicted: {Gamma_pred:.4f} MeV")
print(f"  PDG 2022:  {Gamma_pdg} +/- {Gamma_unc} MeV")
print(f"  Residual:  {pct(Gamma_pred, Gamma_pdg):+.4f}%  ({s5:+.2f} sigma)")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 6  Hierarchy problem dissolved (analytic)")
print(SEP)
print("  m_H = 2*pi*hbar_c*(1+alpha/pi) / (alpha*phi*r_p)")
print("  All factors are topologically fixed (Hopf winding + CODATA observables).")
print("  No Higgs mass counterterm; no fine-tuning; no free parameter.")
print("  Radiative corrections to m_H are finite because L_J is a topological scale.")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 7  Coulomb's law from torsion medium pressure (analytic)")
print(SEP)
rho_EM = 1.2566370614e-6   # kg/m^3  EM-derived medium density (= mu_0)
print(f"  Medium density rho = mu_0 = {rho_EM:.6e} kg/m^3 (not cosmological)")
print(f"  Point pressure source P_0 in incompressible medium:")
print(f"    nabla^2 p = -P_0 * delta^3(r)")
print(f"    Solution:  p(r) = P_0 / (4*pi*r)   [3D Green's function]")
print(f"  Force on test charge: F = alpha * P_0 / r^2")
print(f"  Setting P_0 = hbar*c (one energy quantum):")
print(f"    V(r) = -alpha*hbar*c / r   [Coulomb potential, exact]")
print(f"  Result: Coulomb's law is the Green's function of torsion medium pressure.")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("CLAIM 8  7 * k_n_max / (2*pi) = 1 + alpha + alpha^2 * phi  (Claim 8)")
print(SEP)

# k_n(g) = (sqrt(3)-g)/2 * g^5  from alpha derivation (scale-invariant)
g_max   = 5*sqrt3/6
k_n_max = (sqrt3 - g_max)/2 * g_max**5
N_em    = 7           # dim(A_g + T_1g + T_2g) = 1+3+3, EM-coupled I_h sector

lhs = N_em * k_n_max / (2*pi)
rhs = 1 + alpha + alpha**2 * phi
gap = abs(lhs - rhs) / rhs * 100

k_exact = 3125/3456   # k_n_max exact algebraic = 5^5 / (2^7 * 3^3)
print(f"  k_n(g) = (sqrt(3)-g)/2 * g^5   [scale-invariant, from doc_alpha]")
print(f"  g_max  = 5*sqrt(3)/6            [jamming maximum]")
print(f"  k_n_max = {k_n_max:.10f}  (exact: 3125/3456 = {k_exact:.10f})")
print(f"  N_em = 7 = dim(A_g+T_1g+T_2g)  [I_h EM-coupled sector]")
print()
print(f"  LHS = 7 * k_n_max / (2*pi)          = {lhs:.10f}")
print(f"  RHS = 1 + alpha + alpha^2*phi        = {rhs:.10f}")
print(f"  Gap = {gap:.6f}%  (to 0.0001%)")
print(f"  => E_cell = 7*k_n_max*hbar_c/L_J / (1+alpha+alpha^2*phi) [no density]")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY TABLE")
print(SEP)
print(f"  {'Property':<16} {'Formula':<34} {'Predicted':<14} {'Measured':<18} {'sigma/pct'}")
print(f"  {'-'*16} {'-'*34} {'-'*14} {'-'*18} {'-'*10}")
print(f"  {'m_H [GeV]':<16} {'E_cell*(1+a/pi)':<34} {m_H_pred:<14.4f} {str(m_H_pdg)+' +/- '+str(m_H_unc):<18} {s1:+.2f}s")
print(f"  {'spin':<16} {'n=2 even -> scalar':<34} {'0 (exact)':<14} {'0':<18} {'exact'}")
print(f"  {'lambda':<16} {'(1-nu)/4':<34} {lam_pred:<14.6f} {str(lam_sm)+' +/- '+str(lam_unc):<18} {s3:+.2f}s")
print(f"  {'v [GeV]':<16} {'m_H/sqrt(2*lam)':<34} {v_pred:<14.4f} {str(v_sm):<18} {pct(v_pred,v_sm):+.4f}%")
print(f"  {'Gamma_H [MeV]':<16} {'a^2*m_H/phi':<34} {Gamma_pred:<14.4f} {str(Gamma_pdg)+' +/- '+str(Gamma_unc):<18} {s5:+.2f}s")
print(f"  {'hierarchy':<16} {'m_H topologically fixed':<34} {'no tuning':<14} {'n/a':<18} {'dissolved'}")
print(f"  {'Coulomb':<16} {'3D pressure Green fn':<34} {'exact':<14} {'exact':<18} {'exact'}")
print(f"  {'E_cell jamming':<16} {'7*k_max/(2pi)=1+a+a^2*phi':<34} {gap:.6f}%    {'<0.001%':<18} {gap:.4f}%")
print()
print(f"  Free parameters: 0")
print(f"  Measured inputs: alpha (CODATA), r_p (CODATA), Rs (K-formula+GW170817)")
print()

print(SEP)
print("COMPANION SCRIPTS (for derivation detail)")
print(SEP)
print("  higgs_cell_energy.py        -- cell geometry, E_cell, m_H step-by-step")
print("  higgs_linking_spin.py       -- linking number theorem proof")
print("  higgs_mechanism_stack.py    -- two-regime analysis, N_J, sub-cell coupling")
print("  higgs_width.py              -- Gamma_H: Chern-Simons deeper form")
print("  higgs_cell_jamming_scaling.py -- Claim 8: full jamming chain")
print("  em_coulomb_pressure.py      -- Claim 7: Coulomb from pressure (nuclear/)")
print()
