"""
constants.py  --  Higgs boson and Jobson cell constants for Doc Higgs investigation.

Two parallel sets: Higgs (SM, measured) and Jobson cell (derived from geometry).
Purpose: check whether Higgs properties derive from cell geometry.
"""

import math

pi = math.pi

# ── TOPOLOGY (from (1,2) Hopf fibration -- all derived) ──────────────────────
phi    = (1 + math.sqrt(5)) / 2          # golden ratio = 1.6180...
alpha  = 7.2973525693e-3                 # fine structure constant (CODATA)
r_p    = 0.8414e-15                      # m  proton charge radius (CODATA)
hbar_c = 197.3269804                     # MeV*fm
hbar_c_Jm = 3.16153e-26                  # J*m

# ── JOBSON CELL (derived from topology + CODATA) ──────────────────────────────
L_J       = alpha * phi * (r_p * 1e15)   # fm  cell edge length = 0.00993 fm
L_J_m     = alpha * phi * r_p            # m   = 9.93e-18 m
N_lock    = 2 * pi / (alpha * phi)       # tube closure number = 532.1
E_cell_J  = 2 * pi * hbar_c_Jm / L_J_m  # J   cell energy (retained for reference)
E_cell_GeV = 2 * pi * hbar_c / L_J / 1000  # GeV; use MeV*fm path to avoid eV rounding

# Scalar QED correction for spin-0 particle (standard, not a fit)
alpha_pi  = alpha / pi                   # 0.002323 = 2x Schwinger correction

# Predicted Higgs mass
m_H_pred  = E_cell_GeV * (1 + alpha_pi)  # GeV = 125.089 GeV

# ── HIGGS BOSON (PDG 2022, all measured) ─────────────────────────────────────
m_H_pdg22   = 125.20    # GeV  PDG 2022 combined
m_H_pdg_unc = 0.11      # GeV  1-sigma uncertainty
m_H_old     = 125.09    # GeV  older PDG combined (pre-2022)
v_EW        = 246.22    # GeV  electroweak vacuum expectation value
lam_SM      = m_H_pdg22**2 / (2 * v_EW**2)  # Higgs quartic self-coupling = 0.1293
Gamma_H     = 4.1e-3    # GeV  Higgs total decay width
tau_H       = 6.582e-25 / Gamma_H  # s   Higgs lifetime ~ 1.6e-22 s

# ── HIGGS QUARTIC COUPLING LEAD ───────────────────────────────────────────────
# Candidate: lambda = phi/(4*pi) -- to be investigated
lam_phi4pi  = phi / (4 * pi)             # 0.12877 (-0.40% from SM)
v_pred_phi  = m_H_pred / math.sqrt(phi / (2 * pi))  # GeV: if lam = phi/(4pi)

# ── HIGGS DECAY BRANCHING RATIOS (PDG 2022) ───────────────────────────────────
# Format: (final state, mass_MeV, branching_ratio, N_J)
# N_J = R_compton/L_J = (hbar_c/mass)/L_J
def N_J(mass_MeV):
    if mass_MeV <= 0: return None
    return (hbar_c / mass_MeV) / L_J

higgs_decays = [
    # name,      mass_MeV,  BR,    notes
    ("bb",         4180,   0.581, "b quark -- BOUNDARY REGIME"),
    ("WW*",       80400,   0.214, "W boson -- SUB-CELL"),
    ("gg",            0,   0.082, "gluons -- massless, loop-induced"),
    ("tautau",     1777,   0.063, "tau lepton -- near-bulk"),
    ("ZZ*",       91200,   0.026, "Z boson -- sub-cell"),
    ("cc",         1280,   0.029, "charm quark -- near-bulk"),
    ("gammagamma",    0,   0.002, "photons -- massless, loop-induced"),
    ("Zgamma",        0,  0.0015, "Z+photon -- loop-induced"),
    ("mumu",        106,   0.0002,"muon -- bulk"),
]

# ── ICOSAHEDRAL DECOMPOSITIONS (I_h) ─────────────────────────────────────────
# Irreducible representations and multiplicities (from character table)
# A_g(1), T_1g(3), T_2g(3), G_g(4), H_g(5) + ungerade partners
Ih_reps = {
    "A_g":  1,   # totally symmetric -- already used in alpha derivation
    "T_1g": 3,   # 3 components
    "T_2g": 3,   # 3 components
    "G_g":  4,   # 4 components
    "H_g":  5,   # 5 components -- 5-fold, related to f2=log5
    "A_u":  1,
    "T_1u": 3,
    "T_2u": 3,
    "G_u":  4,
    "H_u":  5,
}
# Note: 1+3+3+4+5 = 16 per parity sector; total dim = 32 per l-shell
# The H_g (5-fold) dominates -- 5 b quarks in proton sea? coincidence?

# ── KNOWN PARTICLE N_J TABLE ────────────────────────────────────────────────────
particle_masses = {
    "top quark":    173000,
    "Higgs boson":  125090,
    "Z boson":       91200,
    "W boson":       80400,
    "b quark":        4180,
    "tau lepton":     1777,
    "c quark":        1280,
    "s quark":         100,
    "proton":          938,
    "muon":            106,
    "pion":            140,
    "electron":        0.511,
}
