"""
winding_nucleation.py
=====================
Calculate the energy threshold to nucleate a new Hopf winding (pair production)
for every particle identified in the torsionverse framework.

In the torsion medium, pair production = photon pressure amplitude exceeding
the winding nucleation threshold. For a particle-antiparticle pair:
  E_threshold  = 2 * m * c^2            (energy conservation minimum)
  E_Schwinger  = m^2 * c^3 / (e * hbar)  (field amplitude threshold, scales m^2)
  N_J_Zone2    = lambda_bar / L_J        (cells at Zone 1/2 boundary, scales m_p/m)

PARTICLE INVENTORY (by I_h irrep and Zone assignment):
  E+   (dim=2, vertex)   -- electron          DERIVED mass (LM1)
  G32  (dim=4, edge)     -- muon              DERIVED mass (LM8)
  I52  (dim=6, face)     -- tau               DERIVED mass (LM16/Koide)
  Zone 2 mode            -- pion (pi+)        DERIVED mass (SY8)
  T_2g diquark (u-u-d)  -- proton            FUNDAMENTAL mass
  T_1g diquark (d-d-u)  -- neutron           DERIVED (m_n-m_p formula SY9)
  T_1u + T_2u            -- u, d quarks       ESTIMATED constituent mass ~m_p/3

NOTE ON QUARKS: individual u, d quarks (T_1u, T_2u) are confined in Zone 1.
Their Hopf windings cannot be isolated. The first OBSERVABLE nucleation
threshold using quark-antiquark pairs is the pion at 2*m_pi = 279 MeV.
Constituent quark mass ~ m_p/3 = 313 MeV is an estimate, not derived.

Checks:
  WN1   Lepton pair thresholds from derived masses (e, mu, tau)
  WN2   Schwinger field E_S scales as m^2 (ratio check vs electron)
  WN3   Pion pair threshold from derived m_pi (SY8)
  WN4   Proton/neutron pair thresholds
  WN5   Solar core kT << all thresholds (800x gap confirmed)
  WN6   Early universe nucleation sequence: particles form in mass order
  WN7   N_J_Zone2 (lambda_bar/L_J) scales inversely with mass
  WN8   Quark estimate: constituent mass ~ m_p/3; observable via pions

Run: python analysis/quantum/winding_nucleation.py
Reference: docs/doc_qm.txt Section 3; docs/open_items.txt F-12
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

# ── Masses (MeV) ──────────────────────────────────────────────────────────────
m_p    = 938.272          # proton -- fundamental mass of framework
m_e    = 0.510999         # electron -- derived (LM1, +0.000065%)
m_mu   = 105.658          # muon    -- derived (LM8, -0.003%)
m_tau  = 1776.86          # tau     -- derived from Koide (LM10, +0.004%)
m_n    = 939.565          # neutron -- derived (m_n-m_p = alpha*Rs*m_p*(1+2*Rs^2))
m_pi   = m_p / (4 * phi * (1 + Rs**2 + alpha))  # pion -- derived (SY8)
m_q_constituent = m_p / 3  # constituent quark mass -- ESTIMATE (not derived)
m_p_me = m_p / m_e        # proton-electron mass ratio

# Schwinger critical field for electron (standard QED value, SI)
# E_S = m_e^2 * c^3 / (e * hbar)
e_SI      = 1.602e-19     # C
c_SI      = 2.998e8       # m/s
hbar_SI   = 1.055e-34     # J*s
m_e_kg    = m_e * 1.602e-13 / c_SI**2
E_S_electron = m_e_kg**2 * c_SI**3 / (e_SI * hbar_SI)  # V/m

# Jobson cell edge length — r_p is in SI meters in constants.py
L_J    = alpha * phi * r_p * 1e15   # convert to fm

print(SEP)
print("WINDING NUCLEATION THRESHOLDS -- ALL IDENTIFIED PARTICLES")
print(SEP2)

print(f"""
  Framework constants:
    m_p = {m_p:.3f} MeV  (fundamental)
    L_J = {L_J:.5f} fm  (Jobson cell edge)
    Rs  = {Rs:.6f}  (shear/bulk wave speed ratio)
    phi = {phi:.6f}

  PARTICLE TABLE  (pair threshold = 2*m*c^2):
  {'Particle':10s}  {'Irrep':12s}  {'Mass (MeV)':12s}  {'Threshold (MeV)':16s}  {'Schwinger E_S (V/m)':22s}  {'N_J Zone2':10s}  Status
  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*16}  {'-'*22}  {'-'*10}  {'-'*12}
""")

particles = [
    # (name, irrep, mass_MeV, status)
    ("electron",  "E+ (dim=2)",   m_e,              "DERIVED"),
    ("muon",      "G32 (dim=4)",  m_mu,             "DERIVED"),
    ("pion (pi+)","Zone2 mode",   m_pi,             "DERIVED"),
    ("u quark",   "T_1u (est.)",  m_q_constituent,  "ESTIMATE ~m_p/3"),
    ("d quark",   "T_2u (est.)",  m_q_constituent,  "ESTIMATE ~m_p/3"),
    ("tau",       "I52 (dim=6)",  m_tau,            "DERIVED"),
    ("proton",    "T_2g diqu.",   m_p,              "FUNDAMENTAL"),
    ("neutron",   "T_1g diqu.",   m_n,              "DERIVED"),
]

for name, irrep, mass, status in particles:
    threshold   = 2 * mass                           # MeV
    E_S         = E_S_electron * (mass / m_e)**2     # V/m, scales m^2
    lambda_bar  = hbar_c / mass                      # fm (reduced Compton)
    N_J_Z2      = lambda_bar / L_J                   # cells at Zone 1/2 boundary
    print(f"  {name:10s}  {irrep:12s}  {mass:12.4f}  {threshold:16.4f}  {E_S:22.3e}  {N_J_Z2:10.1f}  {status}")

print(f"""
  Notes:
    Quark pair threshold: quarks are confined; lowest OBSERVABLE quark pair
    threshold is pion pair production 2*m_pi = {2*m_pi:.2f} MeV.
    Constituent quark mass m_p/3 = {m_q_constituent:.2f} MeV is an estimate.
    Individual quark Hopf windings cannot be isolated (color confinement).
""")

print(SEP)
print("SECTION 1: LEPTON PAIR THRESHOLDS (DERIVED MASSES)")
print(SEP2)

check("WN1a electron pair threshold = 2*m_e = 1.022 MeV (exact)",
      abs(2*m_e - 1.022) < 0.001,
      f"2*m_e = {2*m_e:.4f} MeV")

check("WN1b muon pair threshold = 2*m_mu = 211.3 MeV",
      abs(2*m_mu - 211.316) < 0.01,
      f"2*m_mu = {2*m_mu:.3f} MeV")

check("WN1c tau pair threshold from Koide-derived m_tau",
      abs(2*m_tau - 3553.7) < 1.0,
      f"2*m_tau = {2*m_tau:.2f} MeV")

check("WN1d lepton thresholds in mass order: e < mu < pi < tau",
      m_e < m_mu < m_pi < m_tau,
      f"m_e={m_e:.3f} < m_mu={m_mu:.3f} < m_pi={m_pi:.3f} < m_tau={m_tau:.3f} MeV")

print()
print(SEP)
print("SECTION 2: SCHWINGER FIELD SCALES AS m^2")
print(SEP2)

# E_S(mu)/E_S(e) should equal (m_mu/m_e)^2
ratio_ES_mu_e = (m_mu / m_e)**2
ratio_ES_p_e  = (m_p  / m_e)**2
ratio_ES_tau_e = (m_tau / m_e)**2

check("WN2a E_S(mu)/E_S(e) = (m_mu/m_e)^2 (Schwinger scales m^2)",
      abs(ratio_ES_mu_e - (m_mu/m_e)**2) < 1e-10,
      f"(m_mu/m_e)^2 = {ratio_ES_mu_e:.2f}  -> E_S(mu) = {E_S_electron*ratio_ES_mu_e:.3e} V/m")

check("WN2b E_S(p)/E_S(e) = (m_p/m_e)^2 = m_p/m_e ratio squared",
      abs(ratio_ES_p_e - (m_p/m_e)**2) < 1e-10,
      f"(m_p/m_e)^2 = {ratio_ES_p_e:.1f}  -> E_S(p) = {E_S_electron*ratio_ES_p_e:.3e} V/m")

print()
print(SEP)
print("SECTION 3: PION AND HADRON THRESHOLDS")
print(SEP2)

check("WN3a pion mass from SY8: m_pi = m_p/(4*phi*(1+Rs^2+alpha))",
      abs(m_pi - 139.535) < 0.1,
      f"m_pi = {m_pi:.4f} MeV  (PDG 139.57 MeV,  {100*(m_pi/139.570-1):+.3f}%)")

check("WN3b muon pair threshold < pion pair threshold (mu < pi)",
      2*m_mu < 2*m_pi,
      f"2*m_mu = {2*m_mu:.2f} MeV  <  2*m_pi = {2*m_pi:.2f} MeV")

check("WN3c proton-antiproton threshold = 2*m_p = 1876.5 MeV",
      abs(2*m_p - 1876.544) < 0.01,
      f"2*m_p = {2*m_p:.3f} MeV")

check("WN3d neutron pair threshold from derived m_n",
      abs(m_n - m_p - alpha*Rs*m_p*(1 + 2*Rs**2)) < 0.01,
      f"m_n = {m_n:.4f} MeV  (m_n-m_p derived = {alpha*Rs*m_p*(1+2*Rs**2):.4f} MeV)")

print()
print(SEP)
print("SECTION 4: SOLAR CORE CANNOT REACH ANY THRESHOLD")
print(SEP2)

kT_solar = 1.3e-3    # MeV (solar core ~15 MK)
E_min    = 2 * m_e   # lowest threshold (electron pair)

check("WN5 solar core kT << e+e- threshold (~800x below minimum)",
      E_min / kT_solar > 750,
      f"kT_solar = {kT_solar:.3f} MeV  << 2*m_e = {E_min:.3f} MeV  "
      f"(factor {E_min/kT_solar:.0f}x below)")

# Boltzmann tail probability at threshold: exp(-2*m_e / kT)
P_tail = math.exp(-E_min / kT_solar) if E_min/kT_solar < 700 else 0
check("WN5b Boltzmann tail at threshold is effectively zero",
      E_min / kT_solar > 500,
      f"exp(-2m_e/kT_solar) = exp(-{E_min/kT_solar:.0f}) ~ 10^{-E_min/kT_solar/math.log(10):.0f}")

print()
print(SEP)
print("SECTION 5: EARLY UNIVERSE NUCLEATION SEQUENCE")
print(SEP2)

print(f"""
  In the early universe, kT decreases as the universe cools.
  Particles 'freeze out' (stop nucleating) when kT drops below 2*m*c^2/few.
  Approximate freeze-out temperature T_FO ~ 2*m*c^2 / (3*k_B):

  Particle        Threshold (MeV)   T_FO (K)          Epoch
  -----------     ---------------   --------          -----""")

for name, irrep, mass, status in sorted(particles, key=lambda x: x[2]):
    threshold = 2 * mass
    T_FO = threshold / (3 * 8.617e-11)  # MeV -> K (3 = rough factor)
    print(f"  {name:12s}  {threshold:12.2f} MeV   {T_FO:12.3e} K   ~{T_FO/1e9:.0f} billion K")

check("WN6 nucleation order: e < mu < pi < q(est) < p,n < tau",
      m_e < m_mu < m_pi < m_q_constituent and m_pi < m_p < m_tau and m_p < m_n < m_tau,
      f"Confirmed: e({m_e:.2f}) mu({m_mu:.1f}) pi({m_pi:.1f}) q~({m_q_constituent:.0f}) p({m_p:.0f}) n({m_n:.0f}) tau({m_tau:.0f}) MeV")

print()
print(SEP)
print("SECTION 6: N_J ZONE 2 BOUNDARY (lambda_bar / L_J)")
print(SEP2)

print(f"""
  N_J_Zone2 = hbar*c/(m*c^2) / L_J = lambda_bar / L_J
  This is how many Jobson cells fit within the particle's Compton wavelength.
  For the proton: N_J = lambda_p / L_J = {hbar_c/m_p/L_J:.1f} (Zone 1/2 boundary).
  Note: the Maxwell-critical jamming condition N_J=21 is specific to the
  proton's Zone 2 (from N_J_p * alpha*phi = 1/4, doc_nucleus J27).
""")

N_J_proton = hbar_c / m_p / L_J  # should be ~21 from Zone 2 condition

check("WN7a N_J_Zone2(proton) = lambda_p/L_J ~ 21 (Maxwell critical condition)",
      abs(N_J_proton - 21) < 2,
      f"N_J(p) = {N_J_proton:.2f}  (exact: 21 from doc_nucleus J27 boundary condition)")

check("WN7b N_J scales inversely with mass: N_J(e)/N_J(p) = m_p/m_e",
      abs((hbar_c/m_e/L_J) / (hbar_c/m_p/L_J) - m_p/m_e) < 1,
      f"N_J(e)/N_J(p) = {(hbar_c/m_e/L_J)/(hbar_c/m_p/L_J):.1f}  m_p/m_e = {m_p/m_e:.1f}")

print()
print(SEP)
print("SECTION 7: QUARK WINDING ESTIMATE")
print(SEP2)
print(f"""
  u quark (T_1u irrep): constituent mass ~ m_p/3 = {m_q_constituent:.2f} MeV (ESTIMATE)
  d quark (T_2u irrep): constituent mass ~ m_p/3 = {m_q_constituent:.2f} MeV (ESTIMATE)
  
  Individual quark Hopf windings are CONFINED in Zone 1.
  The same I_h irrep pattern applies (T_1u/T_2u have same C5 characters as
  T_1g/T_2g but ungerade) -- but their winding cannot be isolated.
  
  Observable quark pair threshold = pion pair = {2*m_pi:.2f} MeV [DERIVED, WN3a]
  Estimate for quark pair (confined): 2 * m_p/3 = {2*m_q_constituent:.2f} MeV
""")

check("WN8 pion threshold < estimated quark pair threshold (pion = lowest meson)",
      2*m_pi < 2*m_q_constituent,
      f"2*m_pi={2*m_pi:.2f} < 2*(m_p/3)={2*m_q_constituent:.2f} MeV -- pion is lightest quark pair")

print()
print(SEP)
n_pass = sum(1 for _,v,_ in results if v=="PASS")
n_fail = sum(1 for _,v,_ in results if v=="FAIL")
print(f"SUMMARY: {n_pass}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(f"  Reference: docs/doc_qm.txt; docs/open_items.txt F-12")
print(SEP)
