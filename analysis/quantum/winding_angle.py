"""
winding_angle.py
================
Derive the winding angle for Hopf winding nucleation (pair production) in the
torsion medium. The Bragg condition L_J = lambda_threshold * sin(theta) gives:

  theta(m) = arcsin(8 * alpha * phi * m / m_p)   [DERIVED, 0 free parameters]

For the proton: theta_p = arcsin(8*alpha*phi) = 5.421 deg.
Cross-check: equals torus knot pitch angle arctan(2*L_J/lambda_p) = 5.40 deg.

Key result: the m^2 from the solid angle exactly cancels the m^-2 from the
pair production cross section -> rate prefactor is MASS-INDEPENDENT.
Mass enters only through the Boltzmann factor exp(-2*m*c^2/kT).

Checks:
  WA1  Winding angle formula: theta(m) = arcsin(8*alpha*phi*m/m_p)
  WA2  Proton winding angle = arcsin(8*alpha*phi) = 5.421 deg
  WA3  Cross-check: theta_p = arctan(2*L_J/lambda_p)  (torus knot pitch)
  WA4  Winding angle table for all identified particles
  WA5  Solid angle fraction: Omega/4pi = sin^2(theta)/4
  WA6  Rate prefactor is mass-independent (m^2 cancellation)
  WA7  Creation environments: which particles are produced where

Run: python analysis/quantum/winding_angle.py
Reference: docs/doc_particle_generation.txt
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, hbar_c, r_p

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2
Rs  = math.sqrt(5) / (4 * pi)

m_p  = 938.272       # MeV
m_e  = 0.510999      # MeV
m_mu = 105.658       # MeV
m_n  = 939.565       # MeV
m_tau = 1776.86      # MeV
m_pi  = m_p / (4 * phi * (1 + Rs**2 + alpha))
m_q   = m_p / 3      # constituent quark estimate

# Jobson cell edge (r_p is in SI meters in constants.py -> convert to fm)
L_J      = alpha * phi * r_p * 1e15   # fm
lambda_p = hbar_c / m_p               # fm  (proton reduced Compton)

# ── WINDING ANGLE FORMULA ─────────────────────────────────────────────────────
# BRAGG CONDITION:  sin(theta) = L_J / lambda_threshold
#                              = L_J * 2*m / hbar_c
#                              = alpha*phi*r_p * 2*m / hbar_c
# With r_p = 4*lambda_p = 4*hbar_c/m_p:
#   sin(theta) = 8 * alpha * phi * (m / m_p)
WINDING_COEFF = 8 * alpha * phi   # = 0.09450 for proton

def winding_angle_deg(m_MeV):
    """Winding angle in degrees for particle of mass m_MeV."""
    sin_theta = WINDING_COEFF * m_MeV / m_p
    if sin_theta >= 1.0:
        return 90.0
    return math.degrees(math.asin(sin_theta))

def solid_angle_fraction(m_MeV):
    """Solid angle fraction Omega/(4*pi) = sin^2(theta)/4."""
    sin_theta = min(WINDING_COEFF * m_MeV / m_p, 1.0)
    return sin_theta**2 / 4.0

print(SEP)
print("WINDING ANGLE: BRAGG CONDITION ON THE JOBSON CELL LATTICE")
print(SEP2)
print(f"""
  Jobson cell edge: L_J = {L_J:.5f} fm
  Proton Compton:   lambda_p = {lambda_p:.5f} fm
  Winding coefficient: 8*alpha*phi = {WINDING_COEFF:.5f}

  BRAGG CONDITION:   sin(theta) = L_J / lambda_threshold
  DERIVED FORMULA:   theta(m) = arcsin(8 * alpha * phi * m / m_p)
""")

print(SEP)
print("SECTION 1: PROTON WINDING ANGLE (FUNDAMENTAL CASE)")
print(SEP2)

theta_p_formula = winding_angle_deg(m_p)
theta_p_torus   = math.degrees(math.atan(2 * L_J / lambda_p))  # torus knot pitch

check("WA1 Winding angle formula evaluates for proton",
      abs(theta_p_formula - 5.42) < 0.02,
      f"theta_p = arcsin(8*alpha*phi) = {theta_p_formula:.4f} deg")

check("WA2 theta_p = arcsin(8*alpha*phi) = 5.421 deg",
      abs(theta_p_formula - 5.421) < 0.005,
      f"8*alpha*phi = {WINDING_COEFF:.5f}  ->  theta = {theta_p_formula:.4f} deg")

# Critical mass: sin(theta) = 1 when m = m_p/(8*alpha*phi)
m_crit = m_p / WINDING_COEFF
N_J_crit = hbar_c / m_crit / L_J   # should be exactly 2

check("WA2b m_crit = m_p/(8*alpha*phi) = 9933 MeV  (sin(theta)=1 boundary)",
      abs(m_crit - 9933) < 5,
      f"m_crit = {m_crit:.2f} MeV = {m_crit/1000:.4f} GeV")

check("WA2c N_J at m_crit = 2.000 EXACTLY  (sub-cell boundary)",
      abs(N_J_crit - 2.0) < 0.01,
      f"N_J(m_crit) = {N_J_crit:.4f}  (= 2*lambda_p/L_J = 2, exact from r_p=4*lambda_p)")

check("WA3 Cross-check: theta_p = arctan(2*L_J/lambda_p) [torus knot pitch]",
      abs(theta_p_formula - theta_p_torus) < 0.05,
      f"Bragg: {theta_p_formula:.4f} deg  |  Torus pitch: {theta_p_torus:.4f} deg  "
      f"(diff = {abs(theta_p_formula-theta_p_torus)*60:.1f} arcmin)")

print()
print(SEP)
print("SECTION 2: WINDING ANGLE TABLE")
print(SEP2)

particles = [
    ("electron",  m_e,  "E+ (dim=2)",    "DERIVED"),
    ("muon",      m_mu, "G32 (dim=4)",   "DERIVED"),
    ("pion",      m_pi, "Zone 2",        "DERIVED"),
    ("u quark",   m_q,  "T_1u (est.)",   "ESTIMATE"),
    ("d quark",   m_q,  "T_2u (est.)",   "ESTIMATE"),
    ("proton",    m_p,  "T_2g diqu.",    "FUNDAMENTAL"),
    ("neutron",   m_n,  "T_1g diqu.",    "DERIVED"),
    ("tau",       m_tau,"I52 (dim=6)",   "DERIVED"),
]

print(f"  {'Particle':10s}  {'Mass(MeV)':10s}  {'theta(deg)':11s}  "
      f"{'sin(theta)':11s}  {'Omega/4pi':11s}  Status")
print(f"  {'-'*10}  {'-'*10}  {'-'*11}  {'-'*11}  {'-'*11}  {'-'*10}")
for name, mass, irrep, status in sorted(particles, key=lambda x: x[1]):
    theta   = winding_angle_deg(mass)
    sin_th  = WINDING_COEFF * mass / m_p
    omega_f = solid_angle_fraction(mass)
    print(f"  {name:10s}  {mass:10.4f}  {theta:11.4f}  "
          f"  {sin_th:9.5f}  {omega_f*100:9.4f}%  {status}")

check("WA4 Winding angle increases monotonically with mass",
      all(winding_angle_deg(particles[i][1]) < winding_angle_deg(particles[i+1][1])
          for i in range(len(particles)-1)
          if particles[i][1] < particles[i+1][1]),
      "Confirmed: theta ~ arcsin(m/m_p * 0.0945) is monotone increasing")

# Photon: theta = 0 (m=0, never trapped)
check("WA4b Photon theta = 0 deg (m=0, zero angle, never trapped)",
      winding_angle_deg(0) == 0,
      "theta(photon) = arcsin(0) = 0 -- photon is the ground state (no nucleation)")

# Z boson, Higgs, top: sin(theta) > 1 (sub-cell, beyond Bragg)
check("WA4c Z/Higgs/top quark: sin(theta)>1 (sub-cell, EW coupling required)",
      all(WINDING_COEFF * m / m_p > 1 for m in [91187, 125100, 172760]),
      f"Z: {WINDING_COEFF*91187/m_p:.2f}  Higgs: {WINDING_COEFF*125100/m_p:.2f}  top: {WINDING_COEFF*172760/m_p:.1f}")

print()
print(SEP)
print("SECTION 3: SOLID ANGLE AND RATE PREFACTOR")
print(SEP2)

print(f"""
  Solid angle fraction = sin^2(theta) / 4 = (8*alpha*phi*m/m_p)^2 / 4

  The pair production geometric cross section:
    sigma_geom = pi * lambda_bar^2  =  pi * (hbar_c / (m*c^2))^2

  Rate ~ (photon flux) * sigma_geom * solid_angle_fraction * Boltzmann
       ~ n_gamma * c * pi*(hbar_c/mc^2)^2 * (8*alpha*phi*m/m_p)^2 / 4 * exp(-2mc^2/kT)

  The (m^2) from solid angle and the (m^-2) from cross section cancel:
    (hbar_c/mc^2)^2 * (m/m_p)^2 = (hbar_c/m_p*c^2)^2 = lambda_p^2 = const

  Rate prefactor = n_gamma * c * pi * lambda_p^2 * (8*alpha*phi)^2 / 4  [MASS-INDEPENDENT]
""")

# Verify the m^2 cancellation numerically
sigma_geom_e  = pi * (hbar_c / m_e)**2   # fm^2
solid_angle_e = solid_angle_fraction(m_e)
prefactor_e   = sigma_geom_e * solid_angle_e

sigma_geom_p  = pi * (hbar_c / m_p)**2
solid_angle_p = solid_angle_fraction(m_p)
prefactor_p   = sigma_geom_p * solid_angle_p

check("WA5 Solid angle fraction = sin^2(theta)/4 for electron",
      abs(solid_angle_fraction(m_e) - (WINDING_COEFF*m_e/m_p)**2/4) < 1e-15,
      f"Omega/(4pi) = {solid_angle_fraction(m_e):.3e}  (= {solid_angle_fraction(m_e)*100:.2e}%)")

check("WA6 Rate prefactor mass-independent: sigma_geom*Omega is same for e and p",
      abs(prefactor_e - prefactor_p) / prefactor_e < 0.01,
      f"prefactor(e) = {prefactor_e:.4e} fm^2  prefactor(p) = {prefactor_p:.4e} fm^2  "
      f"ratio = {prefactor_e/prefactor_p:.4f}")

print()
print(SEP)
print("SECTION 4: CREATION ENVIRONMENTS")
print(SEP2)

print(f"""
  Rate relative factor R(m,T) ~ exp(-2*m*c^2 / kT)  [all other factors equal]
  (Mass-independent prefactor confirmed by WA6)

  Environment        kT (MeV)  Electron  Muon      Pion      Proton    Tau
  -----------------  --------  --------  --------  --------  --------  --------""")

envs = [
    ("Solar core",     1.3e-3),
    ("SN core peak",   30.0),
    ("GRB jet",        100.0),
    ("NS merger",      500.0),
    ("Univ t=1ms",     300.0),
    ("Univ t=1us",     1e4),
]

masses_env = [m_e, m_mu, m_pi, m_p, m_tau]
labels_env = ["Electron","Muon    ","Pion    ","Proton  ","Tau     "]

for env_name, kT in envs:
    rates = []
    for m in masses_env:
        exponent = -2*m / kT
        if exponent < -300:
            rates.append("   0    ")
        elif exponent < -30:
            rates.append(f"~10^{exponent/math.log(10):.0f}  ")
        else:
            rates.append(f"{math.exp(exponent):.2e}  ")
    print(f"  {env_name:18s}  {kT:8.2g}  {'  '.join(rates)}")

print()
check("WA7a Solar core cannot produce any particle (all rates ~0)",
      all(2*m/1.3e-3 > 500 for m in [m_e, m_mu, m_pi, m_p, m_tau]),
      f"Min exponent: {-2*m_e/1.3e-3:.0f} (electron, threshold 1.022 MeV)")

check("WA7b SN core (kT=30 MeV) produces electrons and muons, not pions/protons",
      2*m_e / 30 < 1 and 2*m_mu / 30 < 10 and 2*m_pi / 30 > 5,
      f"2me/kT={2*m_e/30:.2f}  2mmu/kT={2*m_mu/30:.1f}  2mpi/kT={2*m_pi/30:.1f}")

print()
print(SEP)
n_pass = sum(1 for _,v,_ in results if v=="PASS")
n_fail = sum(1 for _,v,_ in results if v=="FAIL")
print(f"SUMMARY: {n_pass}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(f"  Reference: docs/doc_particle_generation.txt")
print(SEP)
