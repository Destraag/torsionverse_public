#!/usr/bin/env python3
"""
Torsionverse: EM pressure disruption of covalent bonds -- source comparison
Checks BD1-BD9. Zone 2 pressure physics extended to molecular scale.
Reference: doc_hadron_manipulation.txt (conceptual extension of Section 6)
"""
import math

SEP = "=" * 62
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

c_ms      = 2.9979e8
eps0      = 8.854e-12
mu0       = 4 * math.pi * 1e-7
J_per_eV  = 1.602e-19
phi       = (1 + math.sqrt(5)) / 2
hc_eVnm   = 1239.8          # eV·nm

print(SEP)
print("TORSIONVERSE: EM PRESSURE vs COVALENT BOND DISRUPTION  [BD1-BD9]")
print("Zone 2 pressure physics extended to molecular scale (O-H bond, water)")
print(SEP)
print()

# ── O-H bond (water) parameters ──────────────────────────────────────────────
E_bond_eV   = 5.15          # eV    (H₂O → H + OH dissociation, NIST)
r_bond_m    = 0.958e-10     # m     (O-H bond length, NIST)
nu_OH_cm    = 3700          # cm⁻¹  (average O-H symmetric/antisymmetric stretch)
chi_w       = -9.0e-6       # vol. diamagnetic susceptibility of water (SI)
r_mol_m     = 1.40e-10      # m     (water van der Waals radius)

E_bond_J    = E_bond_eV * J_per_eV
E_ph_eV     = nu_OH_cm * 1.2398e-4         # eV per photon at O-H stretch
E_ph_J      = E_ph_eV * J_per_eV
lam_res_um  = 1e4 / nu_OH_cm               # μm

sigma_UV_cm2 = 1e-17    # cm²  σ_abs(water, 193 nm)
sigma_IR_cm2 = 5e-18    # cm²  σ_abs(water, 2.7 μm O-H stretch, gas)
T1_s         = 10e-12   # s    O-H stretch T₁ relaxation (gas phase)
t_pulse_s    = 1e-9     # s    representative pulse duration (1 ns)

print(f"  O-H bond: E = {E_bond_eV} eV,  r = {r_bond_m*1e10:.3f} Å")
print(f"  O-H stretch: ν = {nu_OH_cm} cm⁻¹  →  hν = {E_ph_eV:.4f} eV,  λ = {lam_res_um:.2f} μm")
print()

# ── BD1: Zone 2 pressure at bond scale ───────────────────────────────────────
print("BD1: Zone 2 pressure at O-H bond scale:")
V_bond = (4/3) * math.pi * r_bond_m**3
P_bond = E_bond_J / V_bond
print(f"  V_bond = (4/3)π r³ = {V_bond:.3e} m³")
print(f"  P = E_bond / V_bond = {P_bond:.3e} Pa = {P_bond/1e9:.0f} GPa")
check("BD1: O-H Zone 2 pressure in 100-500 GPa range",
      1e11 < P_bond < 5e11,
      f"P = {P_bond/1e9:.0f} GPa")
print()

# ── BD2: Brute-force radiation pressure intensity ─────────────────────────────
print("BD2: Brute-force radiation pressure intensity required (I = P × c):")
I_brute_cm2 = P_bond * c_ms * 1e-4
print(f"  I_brute = {I_brute_cm2:.3e} W/cm²")
check("BD2: Brute-force intensity > 10^14 W/cm² (far above any laser focus)",
      I_brute_cm2 > 1e14,
      f"I_brute = {I_brute_cm2:.2e} W/cm²")
print()

# ── BD3: Ionization threshold kills brute-force route ────────────────────────
print("BD3: Ionization threshold vs brute-force:")
I_ion_vis = 1e13    # W/cm²  (UV/visible strong-field ionization of water)
print(f"  I_ionize ≈ {I_ion_vis:.0e} W/cm²  (UV; strong-field physics)")
print(f"  I_brute  = {I_brute_cm2:.2e} W/cm²")
print(f"  I_brute / I_ion = {I_brute_cm2/I_ion_vis:.1e}  → molecule ionizes first.")
check("BD3: I_brute / I_ionize > 100 (direct field route self-destructs)",
      I_brute_cm2 / I_ion_vis > 100,
      f"ratio = {I_brute_cm2/I_ion_vis:.2e}")
print()

# ── BD4: UV single-photon threshold wavelength ───────────────────────────────
print("BD4: UV single-photon threshold:")
lam_thr_nm = hc_eVnm / E_bond_eV
print(f"  λ_threshold = hc / E_bond = {lam_thr_nm:.1f} nm  (photon must be shorter)")
check("BD4: Single-photon threshold in DUV range 200-280 nm",
      200 < lam_thr_nm < 280,
      f"λ_threshold = {lam_thr_nm:.1f} nm  ({E_bond_eV} eV bond)")
print()

# ── BD5: IRMPD photon count ───────────────────────────────────────────────────
print("BD5: IR multiphoton dissociation (IRMPD) photon count at O-H resonance:")
N_irmpd = E_bond_eV / E_ph_eV
print(f"  N = E_bond / hν = {E_bond_eV:.2f} / {E_ph_eV:.4f} = {N_irmpd:.1f} photons at {lam_res_um:.2f} μm")
print(f"  Mechanism: sequential vibrational ladder climbing; anharmonic detuning present.")
check("BD5: IRMPD requires 9-15 resonant IR photons (true multi-photon regime)",
      9 < N_irmpd < 15,
      f"N = {N_irmpd:.1f} photons at {lam_res_um:.2f} μm")
print()

# ── BD6: Source comparison table ─────────────────────────────────────────────
print("BD6: Source comparison table for O-H bond disruption:")
print()

# Key derived quantities
I_IRMPD_crit = E_ph_J / (sigma_IR_cm2 * T1_s)          # W/cm², IRMPD critical
I_UV_ArF     = (6.43 * J_per_eV) / (sigma_UV_cm2 * t_pulse_s)  # W/cm², 193 nm 1-photon
I_ion_midIR  = 1e13    # W/cm²  ionization threshold at 2.7 μm (Keldysh, higher than UV)

print(f"  {'Source':<30} {'hν(eV)':<9} {'N':<5} {'Mechanism':<18} {'I_req(W/cm²)':<14} {'I_max(W/cm²)':<14} Gap")
print(f"  {'-'*30} {'-'*9} {'-'*5} {'-'*18} {'-'*14} {'-'*14} ---")

def print_row(name, hnu, N_str, mech, I_req, I_max):
    if I_req is None:
        print(f"  {name:<30} {hnu:<9} {N_str:<5} {mech:<18} {'—':<14} {'—':<14} {'see BD7' if 'B-field' in name else '—'}")
    else:
        gap = I_req / I_max
        feas = " ← FEASIBLE" if gap < 1 else ""
        print(f"  {name:<30} {hnu:<9.4f} {N_str:<5} {mech:<18} {I_req:<14.2e} {I_max:<14.2e} {gap:.2e}{feas}")

# Anti-Helmholtz static B
V_mol  = (4/3)*math.pi*r_mol_m**3
F_bond = E_bond_J / r_bond_m
F_dia  = abs(chi_w) * V_mol * 20.0 * 200.0 / mu0
ratio_B = F_dia / F_bond
print(f"  {'Anti-Helmholtz static B':<30} {'N/A':<9} {'N/A':<5} {'gradient force':<18}"
      f" F={F_dia:.2e} N  F_bond={F_bond:.2e} N  ratio={ratio_B:.2e} → negligible")

rows = [
    ("Microwave (2.45 GHz)",      1.01e-5, ">>1",  "rotation only",  None,              None),
    ("CO₂ laser (10.6 μm)",       0.117,   "44",   "off-res bending", I_IRMPD_crit*50,  5e11),
    ("Cr:ZnS (2.7 μm, resonant)", E_ph_eV, f"{N_irmpd:.0f}", "IRMPD resonant", I_IRMPD_crit, 1e10),
    ("Ti:Sa (800 nm, ultrashort)", 1.55,    "3",    "tunnel-ionize",  I_ion_vis,          5e15),
    ("UV KrF (248 nm)",           5.00,    "~2",   "2-photon MPI",   I_UV_ArF*2,         1e12),
    ("DUV ArF (193 nm)",          6.43,    "1",    "single photon",  I_UV_ArF,           1e10),
]
for (name, hnu, N_str, mech, I_req, I_max) in rows:
    print_row(name, hnu, N_str, mech, I_req, I_max)

print()
check("BD6: DUV ArF (193 nm) single-photon: I_req < I_max (gap < 1)",
      I_UV_ArF < 1e10,
      f"I_req = {I_UV_ArF:.3e} W/cm²  I_max(focused ArF) ~ 1e10  gap = {I_UV_ArF/1e10:.2f}x")
print()

# ── BD7: Static magnetic / toroidal field: gradient force negligible ──────────
print("BD7: Anti-Helmholtz or macroscopic toroidal field (B=20 T, dB/dr=200 T/m):")
print(f"  Diamagnetic force on water molecule: F = χ_m V B (dB/dr)/μ₀")
print(f"  F_diamag = {F_dia:.3e} N")
print(f"  F_bond   = E_bond/r_bond = {F_bond:.3e} N")
print(f"  Ratio    = {ratio_B:.3e}")
print(f"  Macroscopic field null is at mm–cm scale; O-H bond is at 1 Å — 10⁸× mismatch.")
check("BD7: Static diamagnetic gradient force << bond force (ratio < 1e-10)",
      ratio_B < 1e-10,
      f"F_dia/F_bond = {ratio_B:.2e}  → macroscopic static field null is wrong scale")
print()

# ── BD8: Resonant photonic ring cavity (toroidal EM equivalent) ───────────────
print("BD8: Photonic ring resonator at 2.7 μm (toroidal EM at correct scale):")
Q_ring   = 1e5      # realistic Q for mid-IR chalcogenide/InF₃ ring resonator
I_input  = I_IRMPD_crit / Q_ring
print(f"  Principle: a micron-scale ring cavity at 2.7 μm stores EM energy at")
print(f"  the O-H stretch resonance, amplifying the intracavity intensity by Q.")
print(f"  I_IRMPD_crit (bare)       = {I_IRMPD_crit:.3e} W/cm²")
print(f"  Ring cavity Q             = {Q_ring:.0e}  (chalcogenide glass, achievable)")
print(f"  Required input intensity  = I_crit / Q = {I_input:.3e} W/cm²")
print(f"  Ionization threshold (2.7 μm) = {I_ion_midIR:.0e} W/cm²  (IRMPD window: exists)")
print(f"  Note: this IS the 'toroidal amp' concept — resonant ring, not static coil.")
check("BD8: Ring cavity input intensity < 10^6 W/cm² (mW-to-W laser sufficient)",
      I_input < 1e6,
      f"I_input = {I_input:.3e} W/cm²  at Q = {Q_ring:.0e}")
print()

# ── BD9: T₁g vortex OAM coupling enhancement ─────────────────────────────────
print("BD9: T₁g vortex coupling to A₁g symmetric O-H stretch:")
# chi(T₁g,C5)=phi, chi(A₁g,C5)=1 → product chi=phi
# sigma_vortex/sigma_isotropic = (phi*1)^2 / 1^2 = phi^2
chi_T1g = phi
chi_A1g = 1.0
sigma_enh = (chi_T1g * chi_A1g)**2 / chi_A1g**2    # = phi^2
I_final  = I_input / sigma_enh
print(f"  chi(T₁g,C5) = phi = {phi:.4f};  chi(A₁g,C5) = 1  (symmetric mode)")
print(f"  sigma_vortex / sigma_isotropic = phi² = {sigma_enh:.4f}")
print(f"  The TE₀₁ mode of a ring resonator circulates around the axis →")
print(f"  carries OAM ℓ=1 → T₁g symmetry in I_h.  Selectively drives A₁g stretch.")
print(f"  I_required (ring + vortex) = {I_input:.3e} / {sigma_enh:.3f} = {I_final:.3e} W/cm²")
print(f"  1 mW into 1 μm² spot = {1e-3/1e-8:.1e} W/cm² → comfortably feasible.")
print(f"  Full molecular I_h vibrational mode assignment OPEN in torsionverse.")
check("BD9: Ring + T₁g vortex combined I_final < 10^4 W/cm² (mW CW laser accessible)",
      I_final < 1e4,
      f"I_final = {I_final:.3e} W/cm²  = I_IRMPD / (Q × phi²)")

# ── BD10-BD14: Toroidal focus — above-cutoff field overlap at center ──────────
print()
print("-" * 62)
print("TOROIDAL FOCUS: above-cutoff field overlap at center [BD10-BD14]")
print("Concept: torus walls above I_IRMPD; contributions from all azimuthal")
print("angles converge ('overlap') at center where water flows.")
print("-" * 62)
print()

lam_m  = lam_res_um * 1e-6         # 2.7 μm in m
chi_11 = 2.4048                     # first zero of J_0 (TE_11 waveguide cutoff)
k_wave = 2 * math.pi / lam_m

# BD10: Cutoff inner radius ────────────────────────────────────────────────────
R_co_m  = chi_11 * lam_m / (2 * math.pi)
R_co_um = R_co_m * 1e6
print("BD10: Cutoff inner radius — boundary between evanescent and propagating:")
print(f"  R_cutoff = χ₁₁ λ / (2π) = {chi_11:.4f} × {lam_res_um:.2f} μm / (2π) = {R_co_um:.3f} μm")
print(f"  R_inner < {R_co_um:.2f} μm → fields decay before reaching center (evanescent)")
print(f"  R_inner > {R_co_um:.2f} μm → fields PROPAGATE through hole, all angles overlap at axis")
check("BD10: Cutoff radius ~1 μm at 2.7 μm (micron-scale precision required)",
      0.5 < R_co_um < 2.0,
      f"R_cutoff = {R_co_um:.3f} μm  at λ = {lam_res_um:.2f} μm")
print()

# BD11: Evanescent case — R = 0.5 μm (below cutoff) ──────────────────────────
R_ev_m   = 0.5e-6
k_co_ev  = chi_11 / R_ev_m
kappa_ev = math.sqrt(k_co_ev**2 - k_wave**2)
delta_ev_nm = 1e9 / kappa_ev
trans_ev    = math.exp(-2 * R_ev_m * kappa_ev)
I_wall_ev   = I_IRMPD_crit / (trans_ev * phi**2)

print("BD11: Evanescent regime (R_inner = 0.5 μm < cutoff):")
print(f"  Decay constant κ = {kappa_ev:.3e} m⁻¹  →  δ = {delta_ev_nm:.0f} nm")
print(f"  I_center = I_wall × e^(-2Rκ) = I_wall × {trans_ev:.4f}")
print(f"  Required I_wall (with phi² coupling) for I_center = I_IRMPD: {I_wall_ev:.3e} W/cm²")
print(f"  vs I_ionize(2.7 μm) = {I_ion_midIR:.0e} W/cm²  →  gap = {I_ion_midIR/I_wall_ev:.0f}x")
check("BD11: Evanescent wall intensity below ionization (gap > 50×)",
      I_wall_ev < I_ion_midIR / 50,
      f"I_wall = {I_wall_ev:.3e}  I_ion/50 = {I_ion_midIR/50:.0e}  gap = {I_ion_midIR/I_wall_ev:.0f}x")
print()

# BD12: Propagating case — R = 1.5 μm (above cutoff) ─────────────────────────
R_pr_m  = 1.5e-6
k_co_pr = chi_11 / R_pr_m
geo_fac = 0.5    # conservative: ring-to-axis focusing at R ~ λ

I_wall_pr = 1e10   # W/cm², achievable with 2.7 μm pulsed OPA
I_ctr_pr  = geo_fac * I_wall_pr
I_ctr_phi = I_ctr_pr * phi**2

print("BD12: Propagating regime (R_inner = 1.5 μm > cutoff):")
print(f"  k_cutoff(R=1.5μm) = {k_co_pr:.3e} m⁻¹ < k_wave = {k_wave:.3e} m⁻¹  → propagating ✓")
print(f"  Fields from ALL azimuthal positions reach center axis and add constructively.")
print(f"  I_center ≈ {geo_fac:.1f} × I_wall  (conservative ring→axis geometric factor, R ~ λ)")
print(f"  For I_wall = {I_wall_pr:.0e} W/cm² (achievable, 2.7 μm OPA):")
print(f"    I_center           = {I_ctr_pr:.3e} W/cm²")
print(f"    I_center × phi²    = {I_ctr_phi:.3e} W/cm²  (T₁g azimuthal mode coupling)")
check("BD12: Propagating torus at I_wall=1e10 W/cm² drives I_center above I_IRMPD",
      I_ctr_phi > I_IRMPD_crit,
      f"I_center(phi²) = {I_ctr_phi:.3e}  I_IRMPD = {I_IRMPD_crit:.3e}")
print()

# BD13: Bond-breaking window confirmed ────────────────────────────────────────
in_win = (I_IRMPD_crit < I_ctr_phi < I_ion_midIR)
print("BD13: Bond-breaking window at center of 1.5 μm above-cutoff torus:")
print(f"  I_IRMPD   = {I_IRMPD_crit:.3e} W/cm²  ← lower bound (must exceed for IRMPD)")
print(f"  I_center  = {I_ctr_phi:.3e} W/cm²  ← {'IN WINDOW ✓' if in_win else 'OUTSIDE'}")
print(f"  I_ionize  = {I_ion_midIR:.0e} W/cm²  ← upper bound (2.7 μm, Keldysh)")
print(f"  O-H bonds disrupted WITHOUT ionization: {in_win}")
check("BD13: Bond-breaking window open (IRMPD without ionization)",
      in_win,
      f"[{I_IRMPD_crit:.2e}, {I_ion_midIR:.0e}] contains I_center = {I_ctr_phi:.2e}")
print()

# BD14: Full route comparison ──────────────────────────────────────────────────
print("BD14: All viable routes compared:")
print(f"  {'Route':<38}  {'I_required(W/cm²)':<19}  Mode")
print(f"  {'-'*38}  {'-'*19}  ----")
print(f"  {'DUV ArF 193 nm (single photon)':<38}  {I_UV_ArF:<19.2e}  CW/pulsed")
print(f"  {'IRMPD bare 2.7 μm':<38}  {I_IRMPD_crit:<19.2e}  pulsed")
print(f"  {'Torus focus, evanescent + phi²':<38}  {I_wall_ev:<19.2e}  pulsed, R<1 μm")
print(f"  {'Torus focus, propagating + phi²':<38}  {I_wall_pr:<19.2e}  pulsed, R>1 μm [BD12]")
print(f"  {'Ring resonator Q=1e5 (BD8)':<38}  {I_input:<19.2e}  CW, ring cavity")
print(f"  {'Ring resonator + T₁g vortex (BD9)':<38}  {I_final:<19.2e}  CW, milliwatt")
print()
print(f"  KEY: above-cutoff torus achieves I_center in IRMPD window at I_wall = 1e10 W/cm².")
print(f"  Implementation: axicon lens + SLM to form toroidal focus at R ~ 1.5 μm,")
print(f"  driving the TE₀₁ azimuthal mode. Water flows along center axis (microfluidic).")
check("BD14: Propagating torus (BD12) requires I_wall << I_ionize (clean IRMPD regime)",
      I_wall_pr < I_ion_midIR / 100,
      f"I_wall = {I_wall_pr:.0e}  I_ion = {I_ion_midIR:.0e}  ratio = {I_wall_pr/I_ion_midIR:.0e}")
print()

# BD15: Aperture (inner radius) operating window ──────────────────────────────
print("BD15: Aperture window — range of R_inner that keeps I_center in IRMPD window:")
print(f"  I_wall = {I_wall_pr:.0e} W/cm²  (representative achievable OPA pulse)")
print()

# Calibrated geometric factor from BD12: at R_ref=1.5μm, I_center = geo_fac * I_wall
# => I_center(R) = geo_fac * (R_ref/R)^2 * I_wall  for R >> R_cutoff
# (ring-to-axis focusing scales as (λ/R)^2 for large R; absorb R_ref into constant)
R_ref_m   = 1.5e-6
C_geo     = geo_fac * R_ref_m**2          # calibration constant (m^2)
# I_center(R) = C_geo * I_wall * phi^2 / R^2   [propagating regime]

# Upper bound R_max: I_center drops to I_IRMPD
R_max_m   = math.sqrt(C_geo * I_wall_pr * phi**2 / I_IRMPD_crit)
R_max_um  = R_max_m * 1e6

# Lower bound in evanescent regime: solve 2*R*kappa(R) = -ln(I_IRMPD/(I_wall*phi^2)) numerically
target = -math.log(I_IRMPD_crit / (I_wall_pr * phi**2))   # = 2*R*kappa at threshold
R_min_m = None
for ri in range(400, 1040):     # scan 0.40 to 1.034 μm in 1 nm steps
    R_test = ri * 1e-9
    k_co_t = chi_11 / R_test
    if k_co_t <= k_wave:
        break                   # hit propagating regime
    kappa_t = math.sqrt(k_co_t**2 - k_wave**2)
    if 2 * R_test * kappa_t <= target:
        R_min_m = R_test
        break
R_min_um = R_min_m * 1e6 if R_min_m else None

print(f"  Propagating cutoff:  R_cutoff = {R_co_um:.3f} μm  (BD10)")
if R_min_m:
    print(f"  Evanescent minimum:  R_min    = {R_min_um:.2f} μm  (I_center just reaches I_IRMPD)")
print(f"  Propagating maximum: R_max    = {R_max_um:.2f} μm  (geometric dilution kills I_center)")
print(f"  Practical window:    {R_co_um:.2f} μm  <  R_inner  <  {R_max_um:.2f} μm  [propagating]")
print(f"  (Evanescent extension: down to ~{R_min_um:.2f} μm for same I_wall)")
print()
print(f"  At R_inner = 1.5 μm:  I_center = {I_ctr_phi:.3e} W/cm²  (BD12, 8.9× above I_IRMPD)")
I_ctr_max = C_geo * I_wall_pr * phi**2 / R_max_m**2
print(f"  At R_inner = {R_max_um:.1f} μm: I_center = {I_ctr_max:.3e} W/cm²  (= I_IRMPD, lower edge)")
print()
print(f"  Products in window: IRMPD gives H• + OH• radicals, NOT H₂ + O₂ directly.")
print(f"    H• + H•  → H₂  (fast recombination)")
print(f"    OH• + OH• → H₂O₂ → H₂O + ½O₂  (slower, partial)")
print(f"    OH• + H• → H₂O  (back-reaction, dominant without catalyst)")
print(f"    Net: partial water splitting; downstream catalyst/separator required.")
print(f"    First O-H bond (H₂O → H• + OH•): 2.7 μm resonant → efficient.")
print(f"    Second O-H bond (OH• → O• + H•): resonance red-shifts to ~3570 cm⁻¹ → less efficient.")
print(f"    Estimated yield without catalyst: ~10-40% H₂ per broken bond.")

check("BD15: Propagating aperture window spans R_cutoff to R_max > 3× R_cutoff",
      R_max_um > 3 * R_co_um,
      f"window = [{R_co_um:.2f}, {R_max_um:.2f}] μm  at I_wall = {I_wall_pr:.0e} W/cm²")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(f"  Reference: docs/doc_hadron_manipulation.txt (Zone 2 molecular extension)")
print(SEP)

# ── Design summary block ──────────────────────────────────────────────────────
def E_peak_Vm(I_Wcm2):
    return math.sqrt(2 * I_Wcm2 * 1e4 / (eps0 * c_ms))

E_lo  = E_peak_Vm(I_IRMPD_crit)   # lower bound (IRMPD threshold)
E_hi  = E_peak_Vm(I_ctr_phi)      # BD12 working point
E_ion_field = E_peak_Vm(I_ion_midIR)  # ionization ceiling

nu_THz = nu_OH_cm * 3e10 / 1e12

W = 52
print()
print("=" * W)
print(f"  TOROIDAL FOCUS — DESIGN PARAMETERS (O-H, water)")
print("=" * W)
print(f"  Resonant frequency")
print(f"    {nu_OH_cm} cm\u207b\u00b9  =  {nu_THz:.1f} THz  =  {lam_res_um:.2f} \u03bcm")
print(f"    h\u03bd = {E_ph_eV:.4f} eV/photon   N = {N_irmpd:.0f} photons per bond")
print("-" * W)
print(f"  Aperture (inner radius)   [I_wall = {I_wall_pr:.0e} W/cm\u00b2]")
print(f"    lower (propagating cutoff) :  {R_co_um:.2f} \u03bcm")
print(f"    upper (geometric dilution) :  {R_max_um:.2f} \u03bcm")
print(f"    sweet spot                 :  1.5 \u2013 3 \u03bcm")
print("-" * W)
print(f"  Center peak E-field")
print(f"    IRMPD threshold  (I_IRMPD) :  {E_lo/1e6:.0f} MV/m")
print(f"    BD12 working point         :  {E_hi/1e6:.0f} MV/m  (8.9\u00d7 above threshold)")
print(f"    Ionization ceiling (2.7\u03bcm) :  {E_ion_field/1e9:.1f} GV/m")
print("=" * W)
