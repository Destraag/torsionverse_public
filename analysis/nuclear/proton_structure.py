"""
proton_structure.py
===================
Derives the geometric and pressure structure of the proton in the torsion
medium using all established framework constants.

QUESTION ADDRESSED (from nuclear_pressure.txt):
  What does the proton look like in the torsion medium?
  How does it maintain the negative pressure (Coulomb well)?
  What is the precise role of Jobson cell jamming at the boundary?

ESTABLISHED INPUTS (all from prior scripts, no new assumptions):
  - L_J = alpha*phi*r_p  [Jobson cell edge, constants.py]
  - E_cell = 2*pi*hbar_c/L_J  [cell binding energy, constants.py]
  - N_J = hbar_c / (m * L_J)  [jamming number, rs_v3.py]
  - K = 1/eps_0  [bulk modulus, doc_magnetism M2]
  - rho = mu_0  [medium density, doc_magnetism M1]
  - V(r) = alpha*hbar_c/r  [Coulomb potential, doc_higgs C7]

Run: python analysis/nuclear/proton_structure.py
Reference: docs/nuclear_pressure.txt
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, phi, hbar_c, r_p, E_cell_GeV, L_J as L_J_fm
# constants.py: r_p in metres, L_J in fm, hbar_c in MeV*fm

pi      = math.pi
m_p     = 938.272      # MeV  proton mass (CODATA)
m_e     = 0.51100      # MeV  electron mass (CODATA)
m_b     = 4180.0       # MeV  b quark mass
eps_0   = 8.8542e-12   # F/m
mu_0    = 4*pi*1e-7    # kg/m^3  = rho_medium
c_si    = 2.99792e8    # m/s
hbar_c_Jm = 3.16153e-26 # J·m

r_p_fm  = r_p * 1e15  # convert metres → fm  (r_p = 0.8414e-15 m → 0.8414 fm)

SEP  = "=" * 65
SEP2 = "-" * 65

# ── Derived constants ─────────────────────────────────────────────────────────
# L_J_fm imported directly from constants (already in fm)
E_cell_MeV = E_cell_GeV * 1000    # MeV

# ── N_J formula: hbar_c / (m * L_J)  [from rs_v3.py] ────────────────────────
N_J_p = hbar_c / (m_p * L_J_fm)   # proton
N_J_e = hbar_c / (m_e * L_J_fm)   # electron
N_J_b = hbar_c / (m_b * L_J_fm)   # b quark

results = []
def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

# ── SECTION 1: N_J scale verification ────────────────────────────────────────
print(SEP)
print("SECTION 1: N_J SCALE HIERARCHY  (N_J = hbar_c / (m * L_J))")
print(SEP2)
print(f"  L_J  = alpha*phi*r_p = {L_J_fm:.5f} fm")
print(f"  E_cell               = {E_cell_MeV:.1f} MeV = {E_cell_GeV:.2f} GeV")
print()
print(f"  {'Particle':<12} {'m (MeV)':>12} {'N_J':>10}  Regime")
print(f"  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*20}")
for name, m, nj_known, label in [
    ("b quark",   m_b,  4.75,  "boundary (sub-grain)"),
    ("proton",    m_p,  21.0,  "boundary (Maxwell critical)"),
    ("electron",  m_e,  38870, "deep bulk"),
]:
    nj = hbar_c / (m * L_J_fm)
    print(f"  {name:<12}  {m:>12.1f}  {nj:>10.1f}  {label}")
print()

check("PS1 N_J(proton) = hbar_c/(m_p*L_J) ~ 21  [boundary, Maxwell critical]",
      abs(N_J_p - 21.0)/21.0 < 0.02,
      f"N_J_p = {N_J_p:.2f}  (target 21.0)")
check("PS2 N_J(electron) = hbar_c/(m_e*L_J) ~ 38870  [deep bulk]",
      abs(N_J_e - 38870)/38870 < 0.01,
      f"N_J_e = {N_J_e:.0f}  (target 38870)")
check("PS3 N_J(b quark) = hbar_c/(m_b*L_J) ~ 4.75  [boundary]",
      abs(N_J_b - 4.75)/4.75 < 0.05,
      f"N_J_b = {N_J_b:.3f}  (target 4.75)")

# ── SECTION 2: Proton geometry in medium ─────────────────────────────────────
print()
print(SEP)
print("SECTION 2: PROTON GEOMETRIC STRUCTURE IN THE TORSION MEDIUM")
print(SEP2)

lambda_p_fm   = hbar_c / m_p               # reduced Compton wavelength (fm)
r_p_over_lambda = r_p_fm / lambda_p_fm
V_proton_fm3  = (4/3) * pi * r_p_fm**3   # fm^3
V_cell_fm3    = L_J_fm**3                 # fm^3
N_cells_disp  = V_proton_fm3 / V_cell_fm3

lambda_p_cells = lambda_p_fm / L_J_fm    # = N_J_p
r_p_cells      = r_p_fm / L_J_fm

print(f"  Proton charge radius:       r_p     = {r_p_fm:.4f} fm = {r_p_cells:.1f} L_J cells")
print(f"  Proton Compton wavelength:  lambda_p = {lambda_p_fm:.4f} fm = {lambda_p_cells:.1f} L_J cells")
print(f"  Ratio r_p / lambda_p        = {r_p_over_lambda:.3f}")
print()
print(f"  SCALE LAYERS (outward from proton centre):")
print(f"    r < lambda_p = {lambda_p_fm:.4f} fm  (N_J < 1):  SUB-CELL -- medium absent, quark confinement zone")
print(f"    r ~ lambda_p = {lambda_p_fm:.4f} fm  (N_J = 21): BOUNDARY -- Maxwell jamming (3V-E=6)")
print(f"    lambda_p < r < r_p                TRANSITION -- cells lock in I_h chirality")
print(f"    r > r_p      = {r_p_fm:.4f} fm  (N_J > {r_p_cells:.0f}): BULK -- cells free, carry Coulomb 1/r field")
print()
print(f"  Excluded volume:            V_p    = {V_proton_fm3:.4f} fm^3 = {V_proton_fm3*1e-45:.3e} m^3")
print(f"  Jobson cell volume:         V_cell = {V_cell_fm3:.4e} fm^3")
print(f"  Displaced cells:            N_disp = {N_cells_disp:.3e}")
print()

check("PS4 r_p / lambda_p = 4.0 (proton charge radius = 4 Compton wavelengths)",
      abs(r_p_over_lambda - 4.0)/4.0 < 0.02,
      f"r_p/lambda_p = {r_p_over_lambda:.4f}  (target 4.000)")

# ── SECTION 3: Pressure field at the proton boundary ─────────────────────────
print()
print(SEP)
print("SECTION 3: COULOMB PRESSURE WELL AT THE PROTON SURFACE")
print(SEP2)

V_at_rp   = alpha * hbar_c / r_p_fm           # MeV  Coulomb potential at r = r_p
a_0_fm    = hbar_c / (m_e * alpha)             # fm  Bohr radius = 52917 fm
V_at_a0   = alpha * hbar_c / a_0_fm            # MeV  potential at Bohr radius (= 2 * 13.6 eV)

# Pressure = epsilon_0 * E^2 / 2 at r = r_p
e_charge  = 1.602e-19   # C
r_p_m_si  = r_p_fm * 1e-15  # m
E_field_rp = alpha * hbar_c_Jm / (e_charge * r_p_m_si**2)   # V/m
P_maxwell  = 0.5 * eps_0 * E_field_rp**2                 # Pa  Maxwell stress

# QCD bag constant for comparison (standard value ~60 MeV/fm^3)
B_qcd_MeV_fm3 = 60.0     # MeV/fm^3  (standard MIT bag model)
B_qcd_Pa = B_qcd_MeV_fm3 * 1.602e-13 / (1e-15)**3  # Pa

print(f"  Coulomb potential at r = r_p:  V(r_p) = {V_at_rp:.4f} MeV")
print(f"  Coulomb potential at r = a_0:  V(a_0) = {V_at_a0*1e6:.4f} eV  (= 13.6 eV expected)")
print(f"  Electric field at r = r_p:     E(r_p) = {E_field_rp:.4e} V/m")
print(f"  Maxwell stress pressure:       P(r_p) = {P_maxwell:.4e} Pa")
print(f"  MIT bag constant (QCD):        B_QCD  = {B_qcd_Pa:.4e} Pa")
print(f"  P(r_p) / B_QCD = {P_maxwell/B_qcd_Pa:.2e}")
print()
print(f"  NOTE: P(r_p) >> B_QCD. The Coulomb surface pressure at r_p is the")
print(f"  electric-field Maxwell stress, much stronger than the bag constant.")
print(f"  The bag constant is the net INTERIOR pressure, a softer confinement.")
print(f"  In the torsion medium: interior pressure = 0 (no cells inside lambda_p),")
print(f"  bag constant = E_cell/V_cell = {E_cell_MeV/V_cell_fm3:.3e} MeV/fm^3")
print()

check("PS5 V(a_0) = alpha*hbar_c/a_0 = 27.21 eV = 2x13.6 eV  [virial theorem: E = V/2]",
      abs(V_at_a0*1e6 - 27.211)/27.211 < 0.01,
      f"V(a_0) = {V_at_a0*1e6:.4f} eV  (a_0 = {a_0_fm:.1f} fm,  expected 27.211 eV)")

# ── SECTION 4: Mechanism of sustained negative pressure ──────────────────────
print()
print(SEP)
print("SECTION 4: HOW THE PROTON MAINTAINS NEGATIVE PRESSURE")
print(SEP2)
print("""
  CHAIN (each step established):

  1. EXCLUDED VOLUME [established, exclusion gravity section]:
     The proton occupies r < r_p with zero medium inside.
     V_excluded = {V:.4f} fm^3.

  2. MAXWELL JAMMING AT r ~ lambda_p [N_J = 21, established]:
     The cells at the proton's Compton wavelength scale are exactly at
     3V-E=6 (Maxwell critical). They CANNOT flow inward to fill the void.
     They are locked -- not over-constrained (can't move at all), not
     under-constrained (floppy), but marginally stable.

  3. BULK CELLS TRANSMIT THE RESTORING FORCE [K = 1/eps_0, proven]:
     The medium's bulk modulus K = 1/eps_0 creates an inward restoring
     force on every shell of cells at r > lambda_p. They push inward
     but are stopped by the jammed boundary cells.

  4. NET RESULT -- SUSTAINED 1/r PRESSURE GRADIENT [C7, proven]:
     The jammed boundary transmits the force as a static pressure field.
     V(r) = -alpha*hbar*c/r  [the Coulomb well, Green's function of K]
     This field persists indefinitely with no energy input because the
     jammed state at N_J = 21 is the lowest-energy configuration.

  5. MAINTENANCE ENERGY = ZERO:
     The proton's rest mass m_p*c^2 = {mp:.1f} MeV is the INITIAL cost of
     creating the excluded volume -- the elastic strain energy stored in
     the displaced cells. Once created, the jammed boundary holds the
     configuration. No ongoing energy input is needed.

  6. THE (1,2) HOPF WINDING [argued, from doc_alpha]:
     The inward vs outward direction (charge sign) is set by the winding
     topology. (1,2) winding = inward = positive charge = proton.
     (2,1) mirror = outward = negative charge = electron.
     This is argued but not yet formally derived in this framework.
""".format(V=V_proton_fm3, mp=m_p))

# ── SECTION 5: Electron contrast ─────────────────────────────────────────────
print(SEP)
print("SECTION 5: ELECTRON CONTRAST -- 12-BOUNCE STANDING WAVE")
print(SEP2)
print(f"""
  The electron is NOT an excluded volume. N_J_e = {N_J_e:.0f} >> 1 means the
  cells flow FREELY through the electron's wave packet. No jamming, no
  frozen boundary.

  Instead (from doc_alpha): the electron is a STANDING WAVE that bounces
  12 times per cycle against the I_h icosahedral cell vertices. Each of
  the 12 vertex contacts transfers outward momentum to the surrounding
  cells. Time-averaged over one cycle: net outward pressure (positive
  pressure source).

  This is the opposite mechanism:
    Proton: EXCLUDED VOLUME + JAMMED BOUNDARY → static inward well
    Electron: STANDING WAVE BOUNCING × 12 → dynamic outward pressure

  Both produce a 1/r pressure field (Coulomb's law holds for both)
  because both are described by the same Green's function V(r) ~ 1/r.
  The COUPLING alpha is the same for both because it comes from the
  (1,2)/(2,1) Hopf topology -- same magnitude, opposite direction.

  Electron cavity size = Bohr radius a_0 = {hbar_c/(m_e*alpha):.0f} fm = {hbar_c/(m_e*alpha)/L_J_fm:.0f} L_J cells
  This is the resonant cavity where the 12-bounce standing wave fits.
""")

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
passed = sum(1 for _,s,_ in results if s=="PASS")
failed = sum(1 for _,s,_ in results if s=="FAIL")
print(f"  Total checks: {len(results)}   PASS: {passed}   FAIL: {failed}")
if failed == 0:
    print()
    print("  ALL CHECKS PASSED.")
    print("  Reference: docs/nuclear_pressure.txt")
    print("  Key numbers for the proton structure section:")
    print(f"    N_J_p = {N_J_p:.2f}  (boundary, Maxwell critical)")
    print(f"    lambda_p = {lambda_p_fm:.4f} fm  (Compton wavelength = 21 cells)")
    print(f"    r_p = {r_p_fm:.4f} fm  (charge radius = {r_p_cells:.0f} cells)")
    print(f"    r_p/lambda_p = {r_p_over_lambda:.3f}  (charge radius = 4x Compton)")
    print(f"    V_excluded = {V_proton_fm3:.4f} fm^3  ({N_cells_disp:.2e} cells displaced)")
    print(f"    E_cell = {E_cell_MeV:.1f} MeV  (binding energy at boundary)")
    print(f"    E_cell/m_p = {E_cell_MeV/m_p:.2f}  (cell energy / proton mass)")
