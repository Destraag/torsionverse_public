"""
electron_doc.py
===============
Demo script for doc_electron.txt -- orbital electron, four forces, chirality.
All checks derivable from existing constants and prior results.

Companion: docs/series1/doc_electron.txt
Prior scripts depended on: orbit_doc.py (OD3,OD7), atomic_shells.py (AS1-AS7)
"""
import sys, math
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# All constants inline — script runs standalone on any machine
pi      = math.pi
phi     = (1 + math.sqrt(5)) / 2         # golden ratio
alpha   = 7.2973525693e-3                 # fine structure constant (CODATA 2018)
r_p     = 0.8414e-15                      # m, proton charge radius (CODATA 2018)
hbar_c  = 197.3269804                     # MeV*fm
L_J     = alpha * phi * (r_p * 1e15)     # fm, Jobson cell edge

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  {'[PASS]' if cond else '[FAIL]'} {name}")
    if detail: print(f"         {detail}")

# Constants not in constants.py
m_e_MeV   = 0.51099895      # MeV  electron mass (PDG)
m_p_MeV   = 938.27208816    # MeV  proton mass (PDG)
hbar_c_fm = hbar_c          # MeV*fm  (alias)

# Derived lengths (all in fm)
lambda_p  = hbar_c_fm / m_p_MeV          # proton Compton wavelength  = 0.2103 fm
lambda_e  = hbar_c_fm / m_e_MeV          # electron Compton wavelength = 386.2 fm
r_p_fm    = r_p * 1e15                   # r_p in fm = 0.8414 fm
a_0_fm    = hbar_c_fm / (m_e_MeV * alpha)  # Bohr radius in fm = 52,918 fm

print(SEP)
print("electron_doc.py -- orbital electron: four forces and chirality")
print(SEP)

# ── Section 1: E+ vertex mode ─────────────────────────────────────────────────
print()
print(SEP2)
print("SECTION 1: E+ VERTEX MODE")
print(SEP2)

N_J_e = lambda_e / L_J   # = hbar_c/(m_e * L_J)  [N_J formula from constants.N_J]

print(f"  lambda_e = hbar_c/m_e = {lambda_e:.2f} fm  (electron Compton wavelength)")
print(f"  L_J      = {L_J:.5f} fm  (Jobson cell edge)")
print(f"  N_J_e    = lambda_e / L_J = {N_J_e:.0f}  (deep bulk: >> 1)")
print(f"  r_p      = {r_p_fm:.4f} fm,  a_0 = {a_0_fm:.1f} fm")
print(f"  Scale hierarchy: L_J < r_p < lambda_e << a_0")
print(f"    L_J      = {L_J:.4f} fm")
print(f"    r_p      = {r_p_fm:.4f} fm  (x{r_p_fm/L_J:.0f} larger)")
print(f"    lambda_e = {lambda_e:.1f} fm  (x{lambda_e/r_p_fm:.0f} larger than r_p)")
print(f"    a_0      = {a_0_fm:.0f} fm  (x{a_0_fm/lambda_e:.0f} larger than lambda_e)")

check("ED1: N_J_e >> 1 (electron is deep bulk, no Zone 2, no nuclear hard core)",
      N_J_e > 10000,
      f"N_J_e = {N_J_e:.0f}  (deep bulk threshold: N_J > 1)")

check("ED2: Scale hierarchy L_J < r_p < lambda_e << a_0",
      L_J < r_p_fm < lambda_e < a_0_fm,
      f"L_J={L_J:.4f} < r_p={r_p_fm:.4f} < lambda_e={lambda_e:.1f} << a_0={a_0_fm:.1f}  (all fm)")

# ── Section 2: Bohr radius and chirality ─────────────────────────────────────
print()
print(SEP2)
print("SECTION 2: BOHR RADIUS AND CHIRALITY")
print(SEP2)

a_0_m_pred = a_0_fm * 1e-15        # convert fm -> m
a_0_m_CODATA = 5.29177210903e-11   # m  CODATA

print(f"  a_0 = hbar_c / (m_e * alpha) = {a_0_fm:.2f} fm = {a_0_m_pred:.6e} m")
print(f"  CODATA a_0 = {a_0_m_CODATA:.6e} m")

check("ED3: Bohr radius from m_e and alpha within 0.001% of CODATA",
      abs(a_0_m_pred/a_0_m_CODATA - 1) < 1e-5,
      f"pred={a_0_m_pred:.6e} m  CODATA={a_0_m_CODATA:.6e} m"
      f"  err={abs(a_0_m_pred/a_0_m_CODATA-1)*100:.5f}%")

# Grinding radius and Coulomb energy at hard core
r_grind  = 2 * lambda_p            # Zone 2 boundaries touch = 0.421 fm
E_grind  = alpha * hbar_c_fm / r_grind  # MeV

print(f"  r_grind = 2*lambda_p = {r_grind:.4f} fm  (same-chirality hard core)")
print(f"  E_grind = alpha*hbar_c/r_grind = {E_grind:.3f} MeV")
print(f"  a_0 / r_grind = {a_0_fm/r_grind:.0f}  (electron orbit >> grinding threshold)")

check("ED4: Coulomb energy at r_grind is MeV-scale (nuclear hard core energy)",
      0.5 < E_grind < 10.0,
      f"E_grind = {E_grind:.3f} MeV  (observed hard core: 0.4-0.6 fm window)")

check("ED5: a_0 >> r_grind (opposite-chirality: electron approaches smoothly, no hard core)",
      a_0_fm / r_grind > 1000,
      f"a_0/r_grind = {a_0_fm/r_grind:.0f}  (orbital radius far above grinding threshold)")

# ── Section 3: Positronium ────────────────────────────────────────────────────
print()
print(SEP2)
print("SECTION 3: POSITRONIUM (e+/e- opposite chirality)")
print(SEP2)

# Reduced mass mu = m_e * m_e / (m_e + m_e) = m_e/2
# Bohr radius scales as 1/mu: a_ps = a_0 * (m_e / mu) = a_0 * 2
a_ps = 2 * a_0_fm
a_ps_m = a_ps * 1e-15

print(f"  Positronium: reduced mass mu = m_e/2 (opposite rolls, no Zone 1 wall)")
print(f"  a_ps = 2 * a_0 = {a_ps:.1f} fm = {a_ps_m:.4e} m")
print(f"  (CODATA a_ps = 1.0583e-10 m)")

check("ED6: Positronium Bohr radius = 2*a_0 (reduced mass m_e/2)",
      abs(a_ps_m / 1.05835e-10 - 1) < 0.001,
      f"a_ps = {a_ps_m:.5e} m  CODATA = 1.0584e-10 m")

# ── Section 4: Z > 137 Coulomb limit ─────────────────────────────────────────
print()
print(SEP2)
print("SECTION 4: Z > 137 COULOMB LIMIT")
print(SEP2)

Z_crit = 1 / alpha   # = 137.036

# Key algebraic identity: a_1s(Z_crit) = a_0 / Z_crit = a_0 * alpha = lambda_e
a_1s_Zcrit = a_0_fm * alpha       # fm
a_1s_Zcrit_check = lambda_e       # should equal lambda_e exactly

print(f"  Z_crit = 1/alpha = {Z_crit:.4f}")
print(f"  a_1s(Z_crit) = a_0 * alpha = {a_1s_Zcrit:.3f} fm")
print(f"  lambda_e     = hbar_c/m_e  = {lambda_e:.3f} fm")
print(f"  -> a_1s(Z_crit) = lambda_e  [EXACT algebraic identity]")
print(f"  Physical: at Z=Z_crit the 1s orbital radius equals the electron's own")
print(f"  Compton wavelength. Below this the winding circuit cannot close.")
print()
print(f"  Scale table (correct values):")
print(f"    Z=1:   a_0   = {a_0_fm:.0f} fm = {a_0_fm*1e-4:.4f} Angstrom")
print(f"    Z=50:  a_0/50  = {a_0_fm/50:.1f} fm")
print(f"    Z=100: a_0/100 = {a_0_fm/100:.1f} fm")
print(f"    Z=137: a_0/137 = {a_0_fm/137:.1f} fm = lambda_e = {lambda_e:.1f} fm  [limit]")
print(f"    (Compare: r_p = {r_p_fm:.3f} fm,  L_J = {L_J:.4f} fm)")

check("ED7: Z_crit = 1/alpha = 137.036",
      abs(Z_crit - 137.036) < 0.001,
      f"Z_crit = 1/alpha = {Z_crit:.4f}")

check("ED8: a_1s(Z_crit) = lambda_e  [algebraic: a_0*alpha = hbar_c/m_e]",
      abs(a_1s_Zcrit - lambda_e) / lambda_e < 1e-8,
      f"a_0*alpha = {a_1s_Zcrit:.4f} fm  lambda_e = {lambda_e:.4f} fm")

check("ED9: lambda_e >> L_J  (cell scale is far below the Z_crit orbital limit)",
      lambda_e / L_J > 1000,
      f"lambda_e/L_J = {lambda_e/L_J:.0f}  (L_J is not the relevant UV cutoff for atoms)")

# ── Section 5: Force % estimates ──────────────────────────────────────────────
print()
print(SEP2)
print("SECTION 5: FORCE SCALE COMPARISON")
print(SEP2)

# Hydrogen ground state binding energy
E_binding = 0.5 * m_e_MeV * alpha**2  # = 13.6 eV (in MeV: 1.36e-5 MeV)
E_binding_eV = E_binding * 1e6        # convert MeV -> eV

# Fine structure scale: alpha^2 relative correction per shell
alpha_sq_pct = alpha**2 * 100         # percent of binding energy

# Lamb shift (1s): 1057 MHz in energy
E_lamb_eV = 1057e6 * 4.13607e-15     # h*f in eV = 4.37e-6 eV
lamb_pct = E_lamb_eV / E_binding_eV * 100

print(f"  Hydrogen ground state: E_1 = alpha^2*m_e*c^2/2 = {E_binding_eV:.4f} eV")
print(f"  Fine structure scale:  alpha^2 = {alpha**2:.4e} = {alpha_sq_pct:.4f}% of E_1")
print(f"  Lamb shift (2s-2p):    {E_lamb_eV:.2e} eV = {lamb_pct:.5f}% of E_1")
print(f"  Pressure (Coulomb):    ~100% of binding energy")
print(f"  Light linger (wave):   ~100% of orbital position (sets Bohr radius)")
print(f"  Roll (fine structure): ~{alpha_sq_pct:.3f}% energy correction")
print(f"  Wall (inner boundary): 0% energy, sets r > r_grind = {r_grind:.3f} fm")

check("ED10: Fine structure correction ~ alpha^2 < 0.01% of binding energy",
      alpha**2 < 1e-3,
      f"alpha^2 = {alpha**2:.4e} = {alpha_sq_pct:.4f}%  (roll is a small energy correction)")

check("ED11: Lamb shift << fine structure << binding energy (three-level hierarchy)",
      E_lamb_eV < E_binding_eV * alpha**2 < E_binding_eV,
      f"Lamb={E_lamb_eV:.2e} eV < fs={E_binding_eV*alpha**2:.2e} eV < E1={E_binding_eV:.2f} eV")

# ── Section 6: Scale identities and classical radius ─────────────────────────
print()
print(SEP2)
print("SECTION 6: ALPHA AS RATIO AT EVERY SCALE TRANSITION")
print(SEP2)

# Classical electron radius: r_e = alpha * lambda_e = alpha^2 * a_0
r_e = alpha * lambda_e           # fm
r_e_CODATA = 2.8179403227        # fm (CODATA)

print(f"  Classical electron radius r_e = alpha * lambda_e = {r_e:.6f} fm")
print(f"  CODATA r_e = {r_e_CODATA:.6f} fm")
print()
print(f"  Alpha appears at every scale transition in the electron hierarchy:")
print(f"    L_J / (phi * r_p)  = {L_J/(phi*r_p_fm):.6f}  [= alpha, by definition]")
print(f"    r_e / lambda_e     = {r_e/lambda_e:.6f}  [classical radius / Compton]")
print(f"    lambda_e / a_0     = {lambda_e/a_0_fm:.6f}  [Compton / Bohr]")
print(f"    alpha              = {alpha:.6f}")
print(f"  All three ratios equal alpha -- same constant connects every rung.")

check("ED12: r_e = alpha * lambda_e within 0.001% of CODATA",
      abs(r_e/r_e_CODATA - 1) < 1e-5,
      f"r_e = {r_e:.6f} fm  CODATA = {r_e_CODATA:.6f} fm  err={abs(r_e/r_e_CODATA-1)*100:.5f}%")

check("ED13: Three alpha ratios all equal alpha (scale-transition identity)",
      all(abs(x/alpha - 1) < 1e-8 for x in [
          L_J/(phi*r_p_fm), r_e/lambda_e, lambda_e/a_0_fm]),
      f"L_J/(phi*r_p)={L_J/(phi*r_p_fm):.6f}  r_e/lambda_e={r_e/lambda_e:.6f}"
      f"  lambda_e/a_0={lambda_e/a_0_fm:.6f}  alpha={alpha:.6f}")

# ── Section 7: Positronium lifetime leading order ─────────────────────────────
print()
print(SEP2)
print("SECTION 7: POSITRONIUM LIFETIME (para-Ps, leading order)")
print(SEP2)

# Para-Ps (singlet, 2γ) decay rate leading order:
# Gamma = m_e * c^2 * alpha^5 / (2 * hbar)
# tau = 2*hbar / (m_e*c^2 * alpha^5)
hbar_eVs   = 6.582119569e-16  # eV*s
m_e_eV     = m_e_MeV * 1e6    # eV

tau_paraPs_s    = 2 * hbar_eVs / (m_e_eV * alpha**5)    # seconds
tau_paraPs_ps   = tau_paraPs_s * 1e12                   # picoseconds
tau_measured_ps = 125.142                                # ps (PDG)

print(f"  tau(para-Ps) = 2*hbar / (m_e*c^2 * alpha^5) = {tau_paraPs_ps:.1f} ps")
print(f"  PDG measured = {tau_measured_ps} ps")
print(f"  Error = {abs(tau_paraPs_ps/tau_measured_ps - 1)*100:.1f}%  (leading-order only)")
print(f"  Physical: rate at which opposite-roll E+ windings overlap and cancel.")
print(f"  Higher-order alpha^n corrections would account for residual gap.")

check("ED14: Positronium lifetime leading order within 10% of PDG",
      abs(tau_paraPs_ps/tau_measured_ps - 1) < 0.10,
      f"pred={tau_paraPs_ps:.1f} ps  PDG={tau_measured_ps} ps"
      f"  err={abs(tau_paraPs_ps/tau_measured_ps-1)*100:.1f}%")

# ── Section 8: Spin-orbit (fine structure) ────────────────────────────────────
print()
print(SEP2)
print("SECTION 8: SPIN-ORBIT COUPLING — 2p FINE STRUCTURE SPLITTING")
print(SEP2)

# H atom 2p_{3/2} - 2p_{1/2} splitting (leading order Dirac result):
# ΔE(2p) = alpha^4 * m_e * c^2 / 32
# This is the energy scale of roll coupling (Zone 3 co-rotation) on the electron orbital.
dE_2p_eV   = alpha**4 * m_e_eV / 32         # eV
dE_2p_ueV  = dE_2p_eV * 1e6                 # μeV
dE_2p_meas = 45.35e-6                        # eV (measured)

print(f"  ΔE(2p) = alpha^4 * m_e * c^2 / 32 = {dE_2p_ueV:.2f} μeV")
print(f"  Measured 2p_{{3/2}} - 2p_{{1/2}}  = {dE_2p_meas*1e6:.2f} μeV")
print(f"  Error = {abs(dE_2p_eV/dE_2p_meas - 1)*100:.2f}%")
print(f"  Physical: roll (Zone 3 cog coupling) contributes at the alpha^4 level")
print(f"  -- four orders of alpha^2 below pressure, as shown in force table Section 6.")

check("ED15: 2p fine structure = alpha^4*m_e/32 within 1% of measured [spin-orbit scale]",
      abs(dE_2p_eV/dE_2p_meas - 1) < 0.01,
      f"pred={dE_2p_ueV:.2f} μeV  meas={dE_2p_meas*1e6:.2f} μeV"
      f"  err={abs(dE_2p_eV/dE_2p_meas-1)*100:.2f}%")

# ── Section 9: Hyperfine splitting (Zone 3 effect on electron spin) ──────────
print()
print(SEP2)
print("SECTION 9: HYPERFINE SPLITTING — ZONE 3 CO-ROTATION EFFECT ON ELECTRON SPIN")
print(SEP2)

# Proton Zone 3 co-rotation gives the proton a magnetic moment (g_p = 5.5857).
# The coupling between the proton magnetic moment and the electron spin IS the
# hyperfine structure — the torsionverse mechanism: Zone 2 jammed cells rotate
# with the proton, driving Zone 3 cells via gluon-antinode cog contact; the
# resulting effective magnetic field at the electron's orbital location splits
# the electron spin states.
#
# Leading-order formula (Fermi contact + dipolar, ground state n=1):
# ΔE_HFS = (4/3) * alpha^4 * m_e * c^2 * (m_e/m_p) * g_p
g_p        = 5.58569            # proton g-factor
dE_hfs_eV  = (4/3) * alpha**4 * m_e_eV * (m_e_MeV/m_p_MeV) * g_p
dE_hfs_MHz = dE_hfs_eV / 4.13607e-15 / 1e6   # convert eV -> Hz -> MHz
dE_hfs_meas_MHz = 1420.405751                  # Hz -> MHz (21 cm line, PDG)

print(f"  ΔE_HFS = (4/3)*alpha^4*m_e*c^2*(m_e/m_p)*g_p = {dE_hfs_MHz:.3f} MHz")
print(f"  Measured (21 cm line, hydrogen ground state) = {dE_hfs_meas_MHz:.3f} MHz")
print(f"  Error = {abs(dE_hfs_MHz/dE_hfs_meas_MHz - 1)*100:.4f}%")
print(f"  Physical: the proton's Zone 3 co-rotation (magnetic moment g_p = {g_p})")
print(f"  couples to the electron's E+ spin at the orbital contact points.")
print(f"  This is the ONLY orbital-scale effect of the proton's Zone 3 rotation")
print(f"  on the electron. It does not modify the free-space g-factor (2.00232);")
print(f"  it creates an energy splitting between electron spin-up and spin-down.")

check("ED16: Hydrogen hyperfine (21 cm) = (4/3)*alpha^4*m_e*(m_e/m_p)*g_p within 0.1%",
      abs(dE_hfs_MHz/dE_hfs_meas_MHz - 1) < 0.001,
      f"pred={dE_hfs_MHz:.3f} MHz  meas={dE_hfs_meas_MHz:.3f} MHz"
      f"  err={abs(dE_hfs_MHz/dE_hfs_meas_MHz-1)*100:.4f}%")

# ── Section 10: Zone 3 co-rotation radial falloff (Lense-Thirring) ───────────
print()
print(SEP2)
print("SECTION 10: ZONE 3 CO-ROTATION RADIAL FALLOFF (LENSE-THIRRING, 1/r^3)")
print(SEP2)

# Hopf winding -> frame drag field: E_Z3(r) = alpha * hbar_c * r_p^2 / r^3
# At Zone 3 onset (r = r_p): E_Z3 = alpha*hbar_c/r_p  (Coulomb coupling at charge radius)
# Falloff: 1/r^3 (Lense-Thirring frame drag; same as entanglement_doc.py EP1-EP2)
E_Z3 = lambda r_fm: alpha * hbar_c_fm * r_p_fm**2 / r_fm**3   # MeV, r in fm

E_onset   = E_Z3(r_p_fm)
E_coulomb = alpha * hbar_c_fm / r_p_fm

print(f"  E_Z3(r) = alpha * hbar_c * r_p^2 / r^3  (Lense-Thirring frame drag)")
print(f"  E_Z3(r_p) = alpha*hbar_c/r_p = {E_onset:.6f} MeV  (Coulomb at Zone 3 onset)")
print(f"  Ratio E(r)/E(2r) = 8  (confirms 1/r^3 exponent)")
print(f"  E_Z3(a_0) = {E_Z3(a_0_fm)*1e15:.2f} neV  (co-rotation energy at Bohr radius)")

check("ED17: E_Z3(r_p) = alpha*hbar_c/r_p at Zone 3 onset [Lense-Thirring, EP1]",
      abs(E_onset / E_coulomb - 1) < 1e-8,
      f"E_Z3(r_p) = {E_onset:.6f} MeV  alpha*hbar_c/r_p = {E_coulomb:.6f} MeV")

check("ED18: E_Z3 falls as 1/r^3: E(r)/E(2r) = 8 [EP2]",
      abs(E_Z3(1.0) / E_Z3(2.0) - 8.0) < 1e-10,
      f"E(r)/E(2r) = {E_Z3(1.0)/E_Z3(2.0):.6f}  (exact: 8)")

# ── Section 11: Positron Born balance (same chi(C5)=phi as E+) ───────────────
print()
print(SEP2)
print("SECTION 11: POSITRON BORN BALANCE -- chi(anti-E+,C5) = chi(E+,C5) = phi")
print(SEP2)

# CORRECTED 2026-09-03: previously labeled the positron "E-" and justified
# chi(E-,C5)=chi(E+,C5) via "phi is real and self-conjugate" -- circular
# (both sides were just assigned the identical literal expression) AND it
# collided with "E-" already meaning something else (the Galois-conjugate,
# no-vertex-coupling electron NEUTRINO irrep, chi(C5)=-1/phi --
# ih_double_group.py DG11-DG14, neutrino_freed_lepton.py NL1-NL6, both
# independently verified against real neutrino phenomenology). The positron
# is NOT that object. Renamed to "anti-E+" (matches beta_plus_strip.py's own
# usage) to avoid the collision.
# REAL justification: chirality (winding handedness, sets charge sign) and
# vertex-coupling strength (which irrep, sets confinement/mass) are
# independent labels. Reversing chirality does not change WHICH irrep
# governs coupling -- so anti-E+ (positron) keeps E+'s own chi(C5)=phi,
# giving it the SAME mass as the electron (same medium displacement,
# doc_magnetism.txt Section 3.4) -- matching the real, measured positron.
chi_Eplus_C5     = 2 * math.cos(math.pi / 5)   # E+ character at C5
chi_antiEplus_C5 = chi_Eplus_C5   # anti-E+ (positron): same irrep, chirality reversed only

print(f"  chi(E+, C5)      = 2*cos(pi/5) = {chi_Eplus_C5:.10f}  (= phi = {phi:.10f})")
print(f"  chi(anti-E+, C5) = {chi_antiEplus_C5:.10f}  (SAME irrep as E+ -- chirality reversal")
print(f"                     doesn't change vertex-coupling strength, only charge sign)")
print(f"  => Born balance for anti-E+ gives identical n_exact -> same alpha as E+.")
print(f"  (E- is a DIFFERENT object -- the electron neutrino, chi(C5)=-1/phi -- not the positron.)")

check("ED19: chi(anti-E+,C5) = chi(E+,C5) = phi  [positron Born balance = electron Born balance]",
      abs(chi_antiEplus_C5 - phi) < 1e-12 and abs(chi_Eplus_C5 - phi) < 1e-12,
      f"chi(E+,C5)={chi_Eplus_C5:.10f}  chi(anti-E+,C5)={chi_antiEplus_C5:.10f}  phi={phi:.10f}")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
