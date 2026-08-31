"""
constants.py
Shared physical constants for torsion universe analysis scripts.

All values: CODATA 2018 / PDG 2022 unless noted.
Import with: from constants import *

IMPORTANT: Rs_exact = sqrt(5)/(4*pi) is the conjectured medium saturation
constant. See medium_chains.txt and whitepaper.txt Appendix E for derivation
and evidence. It is NOT a measured constant -- it is the best-fit exact form.
"""

import math

pi      = math.pi
phi     = (1 + math.sqrt(5)) / 2   # golden ratio

# ── Medium saturation constant (conjectured) ────────────────────────────────
Rs_exact = math.sqrt(5) / (4 * pi)   # = 0.17794

# ── Dimensionless ────────────────────────────────────────────────────────────
alpha           = 7.2973525693e-3      # fine structure constant

# ── Speeds / cosmological ────────────────────────────────────────────────────
c_m_s           = 2.99792458e8        # speed of light (m/s)
H0_km_s_Mpc     = 70.0                # Hubble constant (km/s/Mpc)
Mpc_in_m        = 3.085677581e22      # 1 Megaparsec in metres
H0_s            = H0_km_s_Mpc * 1e3 / Mpc_in_m   # H0 in s^-1
cH0             = c_m_s * H0_s        # c * H0 in m/s^2

# ── Particle / EM ───────────────────────────────────────────────────────────
me_MeV          = 0.51099895          # electron rest mass (MeV)

# ── Nuclear / hadronic ───────────────────────────────────────────────────────
mp_MeV          = 938.27208816        # proton rest mass (MeV)
mp_GeV          = mp_MeV / 1000.0
r_proton_fm     = 0.8414              # proton charge radius (fm, CODATA 2018)
kappa_GeV_per_fm = 0.9                # QCD string tension (GeV/fm)
Lambda_QCD_GeV  = 0.217               # QCD confinement scale (GeV, MS-bar)
nuclear_binding_MeV = 8.0             # avg binding energy per nucleon (MeV)
sigma_piN_MeV   = 45.0                # pion-nucleon sigma term (MeV)
m_b_GeV         = 4.18                # bottom quark mass (GeV, PDG 2022)

# ── Galactic / MOND ──────────────────────────────────────────────────────────
a0_m_s2         = 1.2e-10             # MOND critical acceleration (m/s^2)

# ── Gravitational ───────────────────────────────────────────────────────────
G               = 6.67430e-11         # gravitational constant (m^3 kg^-1 s^-2)
M_sun_kg        = 1.98892e30          # solar mass (kg)
M_E_kg          = 5.9722e24           # Earth mass (kg)
R_E_m           = 6.371e6             # Earth radius (m)
J_E             = 5.861e33            # Earth angular momentum (kg m^2 s^-1)
omega_E         = 7.2921150e-5        # Earth rotation rate (rad/s)
AU_m            = 1.495978707e11      # 1 AU in metres

# ── Derived cluster B ratios (for cross-checks) ──────────────────────────────
R_nuclear   = nuclear_binding_MeV / sigma_piN_MeV
R_hadronic  = (kappa_GeV_per_fm * r_proton_fm) / m_b_GeV
R_galactic  = a0_m_s2 / cH0
R_s_mean    = (R_nuclear + R_hadronic + R_galactic) / 3

SEP = "=" * 62

if __name__ == "__main__":
    print("Torsion Universe — Shared Constants")
    print(SEP)
    print(f"  Rs_exact  = sqrt(5)/(4*pi) = {Rs_exact:.6f}")
    print(f"  alpha     = {alpha:.10f}")
    print(f"  R_nuclear = {R_nuclear:.6f}  (8 MeV / sigma_piN)")
    print(f"  R_hadronic= {R_hadronic:.6f}  (kappa*rp / m_b)")
    print(f"  R_galactic= {R_galactic:.6f}  (a0 / cH0)")
    print(f"  R_s_mean  = {R_s_mean:.6f}  (cluster B centre)")
    print(SEP)
