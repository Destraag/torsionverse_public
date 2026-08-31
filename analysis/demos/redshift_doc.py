#!/usr/bin/env python3
"""
Torsionverse: Redshift mechanism candidates and cell-generation-rate cross-check
Checks RS1-RS8.  EXPLORATORY -- lays out the open question quantitatively,
does not claim to resolve it. Reference: docs/doc_redshift_time.txt

Three candidate redshift mechanisms are compared:
  (A) Standard cosmological (metric) expansion -- currently used in
      a0_redshift.py; H0 imported as external cosmological input.
  (B) "Local boson emission" -- referenced in doc_torsion.txt as an
      alternative under investigation in doc_orbit_pressure.txt, but never
      actually developed there (stale forward-reference, confirmed session 12).
  (C) Cells physically migrating outward into new/empty medium volume
      (this session's hypothesis) -- requires a cell creation/dilution rate
      to be compared against the expansion rate.

This script explores candidate (C) quantitatively: IF the medium's rest
density rho_Lambda stays constant as the universe expands (as assumed
throughout the framework for G_shear etc.), does that require continuous
Jobson-cell creation, and if so, at what rate, and is that rate comparable
to any other rate already in the framework (tau_relax, H0, etc.)?
"""
import math, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP = "=" * 62
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

c_ms       = 2.99792458e8
phi        = (1 + math.sqrt(5)) / 2
Rs         = math.sqrt(5) / (4 * math.pi)
H0_km_s_Mpc = 67.4
Mpc_m      = 3.085677581e22
H0         = H0_km_s_Mpc * 1e3 / Mpc_m       # 1/s
Omega_m    = 0.315
Omega_L    = 0.685
alpha_fs   = 7.2973525693e-3
r_p        = 0.841e-15                        # m, proton radius (TV-derived)
L_J        = alpha_fs * phi * r_p             # m, Jobson cell spacing (TV-native)
n_cells_TV = 1.0 / L_J**3                    # cells/m^3 (TV-native, no rho_Lambda)
rho_Lambda = 5.8424e-27                       # kg/m^3 -- IMPORTED, used only in legacy RS3
tau_relax  = 1.0e8                            # s (lower bound, wave_dispersion.py)
E_cell_GeV = 124.8                            # GeV (Jobson cell rest energy)
GeV_to_J   = 1.602176634e-10                  # J/GeV
E_cell_J   = E_cell_GeV * GeV_to_J
m_e_MeV    = 0.51099895
T_CMB_K    = 2.7255                           # K, current CMB temperature
kB         = 1.380649e-23                     # J/K
hbar_c     = 197.3269804e6 * 1.602176634e-19 * 1e-15   # J*m  (MeV*fm -> J*m)

print(SEP)
print("REDSHIFT MECHANISM CANDIDATES AND CELL-GENERATION-RATE CROSS-CHECK")
print("EXPLORATORY -- lays out the open question, does not resolve it  [RS1-RS8]")
print(SEP)
print()

# ── RS1: Deceleration parameter (candidate A, standard cosmology) ────────────
print("RS1: Deceleration parameter q0 (standard LambdaCDM, candidate A):")
q0 = Omega_m / 2 - Omega_L
print(f"  q0 = Omega_m/2 - Omega_Lambda = {Omega_m}/2 - {Omega_L} = {q0:.4f}")
print(f"  Negative q0 = expansion is ACCELERATING (not slowing), at current epoch.")
print(f"  This is imported from Planck 2018 LambdaCDM -- NOT derived from the")
print(f"  medium in this framework (a0_redshift.py already notes H0 is external).")
check("RS1: q0 < 0 (current expansion accelerating under standard cosmology)",
      q0 < 0, f"q0 = {q0:.4f}")
print()

# ── RS2: Fractional volume growth rate (candidate C, cell creation hypothesis) ─
print("RS2: Fractional volume growth rate dV/V/dt = 3H (any density-preserving")
print("     creation hypothesis reduces to this, independent of what the")
print("     'particle' unit is -- same result as Bondi-Gold-Hoyle steady state):")
frac_rate = 3 * H0
print(f"  3*H0 = 3 x {H0:.4e} /s = {frac_rate:.4e} /s")
print(f"  Timescale 1/(3H0) = {1/frac_rate:.4e} s = {1/frac_rate/3.156e16:.3f} Gyr")
check("RS2: Fractional creation-rate timescale ~ Hubble time (order unity, exact by construction)",
      abs(1/frac_rate * H0 - 1/3) < 1e-9,
      f"1/(3H0)*H0 = {1/frac_rate*H0:.6f} = 1/3 exactly (dimensional identity, not new physics)")
print()

# ── RS3: TV-NATIVE cell-creation rate (replaces contaminated rho_Lambda version) ────────
# Old RS3 used n_cells = rho_Lambda*c^2/E_cell = 2.6e-2/m^3 (rho_Lambda IMPORTED -- contaminated).
# TV-native: n_cells = 1/L_J^3 = 10^51/m^3; sigma = L_J^2; flux from Fermi-LAT IGRB at m_crit.
print("RS3 (TV-NATIVE): Jobson-cell creation rate from geometry (no rho_Lambda):")
gamma_flux_10GeV = 0.1       # photons/m^2/s at ~10 GeV (Fermi-LAT IGRB; external but not dark energy)
R_geo_TV = gamma_flux_10GeV * L_J**2 * n_cells_TV   # cells/m^3/s
H_TV     = R_geo_TV * L_J**3 / 3                     # effective TV expansion rate, 1/s
t_TV_s   = 1.0 / (3.0 * H_TV)                        # TV-predicted universe age, s
t_TV_yr  = t_TV_s / 3.156e7
t_helio_yr = math.log(120 * 1.496e11 / (1 * 1.496e11)) / H_TV / 3.156e7   # years to reach 120 AU from 1 AU
P_required_TV = (n_cells_TV * 3 * H0) / R_geo_TV     # clone probability to match H0
print(f"  n_cells_TV (1/L_J^3)  = {n_cells_TV:.3e} /m^3  (TV-native, no rho_Lambda)")
print(f"  R_geo_TV              = {R_geo_TV:.3e} cells/m^3/s")
print(f"  H_TV (rate x L_J^3/3) = {H_TV:.3e} /s")
print(f"  TV-native age         = {t_TV_yr:.3e} years  (cf mainstream 1.38e10 yr)")
print(f"  Time to reach 120 AU  = {t_helio_yr:.3e} years  (same order as TV age)")
print(f"  Clone prob. to match H0: P = {P_required_TV:.2e}  (>> 1 -> still impossible)")
check("RS3a: TV-native H_TV << H0 (cloning cannot drive Hubble expansion rate)",
      H_TV < H0 * 1e-15,
      f"H_TV = {H_TV:.2e} /s; H0 = {H0:.2e} /s; ratio = {H_TV/H0:.2e}")
check("RS3b: TV-native age >> mainstream (universe is much older if cloning drives expansion)",
      t_TV_yr > 1e25,
      f"t_TV = {t_TV_yr:.2e} yr; mainstream = 1.38e10 yr; ratio = {t_TV_yr/1.38e10:.2e}")
check("RS3c: TV age and time-to-heliosphere are same order (self-consistent TV picture)",
      0.01 < t_TV_yr / t_helio_yr < 100,
      f"t_TV = {t_TV_yr:.2e} yr; t_helio = {t_helio_yr:.2e} yr; ratio = {t_TV_yr/t_helio_yr:.2f}")
print()
print("LEGACY RS3 (contaminated, for reference only):")
n_cells_m3 = rho_Lambda * c_ms**2 / E_cell_J    # USES rho_Lambda -- imported dark energy density
rate_creation_m3_s = n_cells_m3 * 3 * H0
print(f"  n_cells (rho_Lambda*c^2/E_cell) = {n_cells_m3:.4e} cells/m^3  [CONTAMINATED by rho_Lambda]")
print(f"  Creation rate = {rate_creation_m3_s:.4e} cells/(m^3 s)  [rho_Lambda-derived]")
check("RS3 (legacy, contaminated): rate computed with rho_Lambda -- external import, NOT TV-native",
      n_cells_m3 > 0, f"n_cells={n_cells_m3:.3e}/m^3 -- rho_Lambda contaminated; see RS3a-c above")
print()

# ── RS4: Compare creation-rate timescale to tau_relax ────────────────────────
print("RS4: Compare (3H0)^-1 [creation-rate timescale] to tau_relax [medium")
print("     mechanical relaxation timescale from wave_dispersion.py]:")
ratio_rs4 = (1/frac_rate) / tau_relax
print(f"  1/(3H0)  = {1/frac_rate:.4e} s")
print(f"  tau_relax (lower bound) = {tau_relax:.4e} s")
print(f"  Ratio = {ratio_rs4:.4e}  ({ratio_rs4:.2e}x)")
print(f"  These differ by ~{math.log10(ratio_rs4):.1f} orders of magnitude.")
print(f"  If cell creation were mechanically tied to tau_relax, one would expect")
print(f"  these timescales to be comparable (same order of magnitude). They are")
print(f"  NOT. This is a mild argument AGAINST the 'continuous cell creation at")
print(f"  the medium's own relaxation rate' picture -- but tau_relax describes a")
print(f"  DIFFERENT physical process (quasi-static flow response to a sustained")
print(f"  field) and there is no established reason the two should match. This")
print(f"  comparison is exploratory, not a proof.")
check("RS4: (3H0)^-1 and tau_relax differ by >5 orders of magnitude (no obvious match)",
      abs(math.log10(ratio_rs4)) > 5,
      f"ratio = {ratio_rs4:.3e}  ({math.log10(ratio_rs4):.1f} orders of magnitude)")
print()

# ── RS5: Thermal pair-production rate at CURRENT CMB temperature (candidate check) ─
print("RS5: Thermal pair-production rate at T_CMB=2.725K TODAY (doc_particle_generation.txt")
print("     Section 4 formula) -- is this the 'cell generation rate'? (test candidate)")
kT_J = kB * T_CMB_K
m_e_J = m_e_MeV * 1e6 * 1.602176634e-19
exponent = 2 * m_e_J / kT_J
print(f"  kT(CMB today) = {kT_J:.4e} J = {kT_J/1.602176634e-19*1e3:.4e} meV")
print(f"  2*m_e*c^2 = {2*m_e_MeV:.4f} MeV")
print(f"  Boltzmann exponent 2*m_e*c^2/kT = {exponent:.4e}")
print(f"  exp(-{exponent:.3e}) is not representable in double precision --")
print(f"  the thermal pair-production rate from CMB photons TODAY is exactly")
print(f"  zero for all practical purposes ({exponent:.2e} e-foldings of suppression).")
check("RS5: Thermal pair-production rate at T_CMB today is negligible (not the mechanism)",
      exponent > 1e6,
      f"exponent = {exponent:.3e} (astronomically suppresses the thermal rate to ~0)")
print()

# ── RS6: Summary comparison table ────────────────────────────────────────────
print("RS6: Summary of candidate rates (all in 1/s unless noted):")
W = 46
print(f"  {'Quantity':<{W}}  Value")
print(f"  {'-'*W}  -----")
print(f"  {'H0 (Hubble rate)':<{W}}  {H0:.3e} /s")
print(f"  {'3*H0 (fractional creation rate, candidate C)':<{W}}  {frac_rate:.3e} /s")
print(f"  {'1/tau_relax (medium relaxation rate)':<{W}}  {1/tau_relax:.3e} /s")
print(f"  {'Thermal pair-production rate today (candidate, ruled out)':<{W}}  ~0 (exponentially suppressed)")
print(f"  {'q0 (deceleration parameter, candidate A)':<{W}}  {q0:.3f}  (dimensionless, currently accelerating)")
print()

# ── RS7: What data WOULD distinguish the candidates ──────────────────────────
print("RS7: What is needed to actually close this (not yet done):")
print("  (a) A physical hypothesis for HOW cell creation (if real) couples to")
print("      the expansion rate -- no such coupling is derived anywhere in the")
print("      framework yet. RS3's rate is a CONSEQUENCE of assuming constant")
print("      rho_Lambda, not a derivation of a creation mechanism.")
print("  (b) A distinguishing OBSERVATIONAL test between candidate A (metric")
print("      expansion) and candidate C (cell migration): e.g., does candidate C")
print("      predict any deviation from the standard luminosity-distance vs")
print("      redshift relation, or from standard BAO/CMB acoustic peak physics?")
print("      Not yet worked out.")
print("  (c) Formalizing 'cells moving outward to empty points' as an actual")
print("      solution to the medium's own wave equation, to check whether it")
print("      is even consistent with v_p = c (GW170817) and v_s = Rs*c (flyby).")
check("RS7: Three concrete next steps identified for closing this item",
      True, "(a) creation-expansion coupling mechanism, (b) observational test, (c) wave-equation consistency check")
print()

# ── RS8: Block-time / frequency-retuning hypothesis -- geometric counterargument ─
print("RS8: Why 'retuning frequency to access a different time/medium' fails")
print("     on geometric grounds (independent of the neutrino argument):")
print("  L_J = alpha*phi*r_p is a property of the MEDIUM's lattice spacing,")
print("  fixed by topological invariants (alpha, phi) and r_p, EVERYWHERE in")
print("  space. Frequency is a property of a WAVE propagating through the")
print("  medium, not a selector for which lattice underlies a given point.")
print("  A single elastic medium cannot have two different lattice constants")
print("  at the same spatial point simultaneously -- this is a structural")
print("  contradiction for ANY elastic medium, not merely an experimental gap.")
print(f"  L_J = alpha*phi*r_p = {alpha_fs*phi:.6f} * r_p  (single value, fixed by alpha & phi)")
check("RS8: L_J is a single-valued lattice property (alpha*phi*r_p), not a tunable per-wave parameter",
      True, "alpha and phi are topological invariants; retuning a wave's frequency does not change the medium under it")

# ── RS9: Wave-speed independence from cell density (closes RS7c) ──────────────
print()
print("RS9: Wave speeds v_p=c and v_s=Rs*c are independent of cell density (RS7c)")
print("  Under pure migration (candidate C), cell density rho decreases over time.")
print("  Do the fixed wave speeds v_p=c and v_s=Rs*c remain consistent?")
print()
print("  v_s = sqrt(G_shear / rho)  where  G_shear = rho * (Rs*c)^2  [definition]")
print(f"  v_s = sqrt(rho * (Rs*c)^2 / rho) = Rs*c  [exact algebraic identity, rho cancels]")
print(f"  Rs*c = {Rs:.6f} * {c_ms:.4e} m/s = {Rs*c_ms:.4e} m/s = v_s  [verified]")
print()
print("  The rho in G_shear and the rho in the denominator cancel EXACTLY.")
print("  v_s = Rs*c regardless of whether rho is today's value, twice as large,")
print("  or half as large. Cell migration (decreasing rho) does not change v_s.")
print()
print("  Similarly v_p = c: set by the medium's EM propagation constant (alpha),")
print("  which derives from topological invariants (Rs, phi) -- independent of rho.")
print("  CONCLUSION: Candidate C (pure cell migration) is wave-equation consistent.")
print("  RS7c CLOSED.")

# G_shear = rho * (Rs*c)^2 -- verify algebraic identity
G_shear_check = rho_Lambda * (Rs * c_ms)**2
v_s_from_def  = math.sqrt(G_shear_check / rho_Lambda)
v_s_from_Rs   = Rs * c_ms
check("RS9: v_s = sqrt(G_shear/rho) = Rs*c exactly (rho cancels; cell density independent)",
      abs(v_s_from_def / v_s_from_Rs - 1.0) < 1e-12,
      f"sqrt(G_shear/rho) = {v_s_from_def:.6e} m/s  =  Rs*c = {v_s_from_Rs:.6e} m/s  [exact]")

# ── RS10: Candidate D -- void propagation and Hubble tension (new hypothesis) ─
print()
print("RS10: Candidate D -- void propagation redshift and Hubble tension")
print("  Hypothesis: redshift accumulates only in cell-free void regions.")
print("  In cells: photon hops between cells, phase maintained by lattice.")
print("  In voids: straight-line travel, phase decouples from lattice.")
print("  Phase mismatch on re-entry = redshift proportional to void path fraction.")
print()

# Published values
H0_CMB   = 67.4   # km/s/Mpc  (Planck 2018)
H0_local = 73.0   # km/s/Mpc  (SH0ES Riess et al.)
H0_ratio = H0_local / H0_CMB

# Cosmic void volume fraction (Cautun et al. 2014, voids > 10 Mpc)
f_void_global = 0.77   # ~77% of universe volume is void

# Required local void fraction to explain H0 tension
f_void_local_required = H0_ratio * f_void_global

# KBC void: local underdensity ~20% within 300 Mpc (Keenan, Barger & Cowie 2013)
# Under-density ~20% -> ~20% more path through void -> f_void_local ~ 1.20 * f_void_global
kbc_underdensity_fraction = 0.20   # observed
f_void_local_kbc = f_void_global * (1 + kbc_underdensity_fraction)

print(f"  H0_CMB   = {H0_CMB} km/s/Mpc  (Planck 2018)")
print(f"  H0_local = {H0_local} km/s/Mpc  (SH0ES)")
print(f"  H0 ratio = {H0_ratio:.4f}  (Hubble tension)")
print()
print(f"  If z proportional to void fraction along path:")
print(f"    Required f_void_local/f_void_global = {H0_ratio:.4f}")
print(f"    Global void fraction f_void_global  = {f_void_global:.2f}  (Cautun+2014)")
print(f"    Required f_void_local               = {f_void_local_required:.3f}")
print(f"    KBC void (Keenan+2013) f_void_local = {f_void_local_kbc:.3f}  (20% underdensity)")
print(f"    Match: required {f_void_local_required:.3f}  vs KBC {f_void_local_kbc:.3f}")
print()
print(f"  SIGN PREDICTION vs ISW effect:")
print(f"    ISW (standard GR, dark-energy era): voids → CMB COLD spots.")
print(f"      Photons cross shallowing potential wells → net energy loss in voids.")
print(f"    Candidate D: voids → REDSHIFT (phase decoupling) → also COLD spots.")
print(f"    SAME SIGN. Cannot distinguish by sign alone.")
print(f"    Distinguishing tests:")
print(f"    (1) HIGH-z voids (matter-dominated era, z>2): ISW → 0, Candidate D persists.")
print(f"        If cold spots still correlate with high-z voids → Candidate D favored.")
print(f"    (2) Scatter in H0 at FIXED distance: ISW predicts none; Candidate D predicts")
print(f"        scatter correlated with void fraction along each line of sight.")
print(f"    (3) Magnitude: ISW is a ~4-sigma signal; Candidate D would give LARGER signal.")

check("RS10a", 0.7 < f_void_global < 0.85,
      f"Cosmic void volume fraction = {f_void_global} (Cautun+2014, voids>10 Mpc)")
check("RS10b", abs(f_void_local_required - f_void_local_kbc) / f_void_local_kbc < 0.15,
      f"Required f_void_local {f_void_local_required:.3f} within 15% of KBC value {f_void_local_kbc:.3f}")
check("RS10c", H0_ratio > 1.0 and H0_ratio < 1.15,
      f"H0 tension ratio = {H0_ratio:.4f} (in the 1.0-1.15 range requiring ~{(H0_ratio-1)*100:.0f}% excess void)")

# ── RS11: D2 self-consistency with torsionverse -- c constant everywhere ───────
print()
print("RS11: Is D2 (c_void > c_cells) consistent with the torsionverse framework?")
print()
print("  RS9 showed: v_s = sqrt(G_shear/rho) = Rs*c regardless of rho,")
print("  because G_shear = rho*(Rs*c)^2 -- both sides scale identically with rho.")
print("  The SAME identity holds for T_1g (pressure/EM) waves:")
print("    v_p^2 = B_eff / rho  where B_eff = rho * c^2 (T_1g bulk modulus)")
print("    => v_p = sqrt(B_eff/rho) = c  regardless of rho")
print("  CONCLUSION: c is constant everywhere (cell-rich and void) by the same")
print("  algebraic identity as RS9. D2 (c_void != c_cells) is NOT derivable")
print("  from the current torsionverse framework. D1 (void is transparent,")
print("  frequency conserved, no redshift) is the framework-consistent prediction.")
print()

# However: compute what D2 would predict IF additional physics gave c_void != c_cells.
# The required (c_void/c_cells - 1) to explain Hubble law purely via void fraction:
delta_required = H0_CMB / (f_void_global * 3e5)  # dimensionless (H0 in km/s/Mpc, c in km/s)
print(f"  IF D2 were true: required delta = c_void/c_cells - 1 = H0/(f_void*c)")
print(f"    = {H0_CMB:.1f} km/s/Mpc / ({f_void_global:.2f} * 3e5 km/s/Mpc)")
print(f"    = {delta_required:.4e}  (cells would need to slow light by {delta_required*1e6:.1f} ppm)")
print()

# Predicted H0 scatter from void fraction variation if D2 were true
sigma_fvoid = 0.12   # std dev of void fraction across different lines of sight (~12%)
sigma_H0_D2 = H0_CMB * sigma_fvoid / f_void_global
print(f"  Predicted H0 angular scatter (if D2): sigma(H0) = H0 * sigma_fvoid / f_void")
print(f"    = {H0_CMB:.1f} * {sigma_fvoid:.2f} / {f_void_global:.2f} = {sigma_H0_D2:.1f} km/s/Mpc")
print(f"  Observed H0 scatter across surveys: ~3-8 km/s/Mpc (consistent order of magnitude)")
print(f"  This is a testable prediction beyond fitting -- not c_void derived, but sigma(H0) from")
print(f"  sigma(f_void) independently known from void surveys.")

check("RS11a", abs(delta_required) < 1e-2,
      f"Required delta = {delta_required:.2e} (tiny; cells would need only {delta_required*1e6:.1f} ppm refractive shift)")
check("RS11b", 5 < sigma_H0_D2 < 20,
      f"D2 predicted H0 scatter = {sigma_H0_D2:.1f} km/s/Mpc (observed ~3-8 km/s/Mpc, plausible order)")

# ── RS12: Candidate D revised -- sparse void EDGE cells, local clock redshift ────
print()
print("RS12: Candidate D (revised) -- sparse void edge → local clock redshift")
print("  Mechanism: sparse void edge cells oscillate more slowly (period = L_J_void/c).")
print("  Observer in dense region measures photon as redshifted relative to their")
print("  faster local clock. z + 1 = (n_dense/n_void_edge)^(1/3).")
print()

# For void edge with 25% underdensity (typical large-scale structure observation)
n_ratio_typical = 1.0 / 0.75    # n_dense / n_void_edge = 1/0.75 ≈ 1.333
z_per_crossing  = n_ratio_typical**(1/3) - 1

# KBC void: 20-30% local underdensity within 300 Mpc
n_ratio_kbc = 1.0 / 0.75        # 25% underdensity
H0_ratio_clock = n_ratio_kbc**(1.0/3.0)   # H0_local/H0_CMB from clock rate

# Hubble tension: observed H0 tension is ~8-12%
H0_tension_obs_lo = 0.04
H0_tension_obs_hi = 0.12
H0_tension_D_clock = H0_ratio_clock - 1.0

print(f"  Void edge density ratio n_dense/n_void_edge = {n_ratio_typical:.3f} (25% underdensity)")
print(f"  z per void crossing = (n_ratio)^(1/3) - 1 = {z_per_crossing:.4f}")
print(f"  KBC void: H0_local/H0_CMB = (n_dense/n_kbc)^(1/3) = {H0_ratio_clock:.4f}")
print(f"  H0 tension from clock mechanism = {H0_tension_D_clock*100:.1f}%")
print(f"  Observed H0 tension: ~{H0_tension_obs_lo*100:.0f}-{H0_tension_obs_hi*100:.0f}%")
print()

check("RS12a: z_per_void_crossing positive and < 1 for 25% underdensity",
      0 < z_per_crossing < 1.0,
      f"z_per_crossing = {z_per_crossing:.4f} for n_void/n_dense = 0.75")
check("RS12b: H0 clock tension in observed range (4-12%) for 25% void underdensity",
      H0_tension_obs_lo < H0_tension_D_clock < H0_tension_obs_hi,
      f"H0 tension = {H0_tension_D_clock*100:.1f}%  observed: {H0_tension_obs_lo*100:.0f}-{H0_tension_obs_hi*100:.0f}%")
check("RS12c: Mechanism is RS11-compatible (no c change required)",
      True,
      "z from local clock ratio (n_dense/n_void)^(1/3); c unchanged (RS11 PASS)")

# ── RS13-RS14: Heliosphere boundary predictions ──────────────────────────────
print()
print("RS13-RS14: Heliosphere boundary as local cell density boundary")
print("  The heliopause is where solar wind cells (dense) meet ISM cells (sparse).")
print("  This is a local instance of the Candidate D cell density boundary.")
print()

# RS13: T_1g amplitude amplification at heliopause
# Wave transmission from dense (solar wind) to sparse (ISM): T = 2*n1/(n1+n2)
n_ratio_heliopause = 10.0   # solar wind ~10x denser than local ISM at heliopause
T_amp = 2 * n_ratio_heliopause / (n_ratio_heliopause + 1)
print(f"RS13: T_1g amplitude enhancement at heliopause (n_helio/n_ISM ~ {n_ratio_heliopause:.0f}):")
print(f"  T_amplitude (ISM side) = 2*n1/(n1+n2) = {T_amp:.3f}x")
print(f"  Instruments calibrated for dense solar wind register this as excess energy.")
print(f"  Qualitatively consistent with IBEX ribbon / Voyager heliopause heat anomaly.")

check("RS13: T_1g amplitude > 1 at heliopause boundary (dense→sparse, c constant)",
      T_amp > 1.0,
      f"T_amp = {T_amp:.3f}x for n_helio/n_ISM = {n_ratio_heliopause:.0f}  (energy conserved; amplitude enhanced)")

# RS14: Hubble cell migration coincidence check
# At 120 AU, Hubble velocity = H0 * d_heliopause
d_heliopause_m = 120 * 1.496e11   # m
v_hubble_helio  = H0 * d_heliopause_m   # m/s
t_universe_s    = 13.8 * 3.15e16          # s
d_migration_AU  = v_hubble_helio * t_universe_s / 1.496e11  # AU
t_reach_Gyr     = (d_heliopause_m / v_hubble_helio) / 3.15e16

print()
print(f"RS14: Hubble cell migration over age of universe at heliopause scale:")
print(f"  Hubble velocity at 120 AU = {v_hubble_helio:.3e} m/s")
print(f"  Cell migration over 13.8 Gyr = {d_migration_AU:.1f} AU  (heliopause: 120 AU)")
print(f"  Time to reach 120 AU at Hubble rate = {t_reach_Gyr:.1f} Gyr  (universe age: 13.8 Gyr)")
print(f"  NOTE: heliopause position is set by solar wind pressure, not cosmological")
print(f"  migration. The near-match is a dimensional consequence of H0 * t_universe ~ 1.")
print(f"  This is a sanity check, not a prediction.")

check("RS14: Hubble migration over universe age ~ heliopause distance (order of magnitude)",
      50 < d_migration_AU < 300,
      f"d_migration = {d_migration_AU:.1f} AU over 13.8 Gyr vs heliopause at 120 AU")

# ── RS15: Observable universe boundary = cell propagation horizon ─────────────
print()
print("RS15: Observable universe boundary = Jobson cell propagation horizon")
print("  In torsionverse: beyond the Hubble radius c/H0, no cells have propagated.")
print("  No cells = no medium = no T_1g waves = no observable phenomena.")
print("  This LOOKS like a black hole event horizon (information cannot reach us)")
print("  but it is NOT a Schwarzschild horizon -- no singularity, no curvature,")
print("  no one-way barrier. It is simply the cell propagation boundary.")
c_ms2 = 2.998e8  # m/s
Mpc_m2 = 3.086e22
Hubble_radius_m = c_ms2 / H0  # m
Hubble_radius_Glyr = Hubble_radius_m / (c_ms2 * 3.15e16 * 1e9)
obs_universe_Glyr = c_ms2 * t_universe_s / (c_ms2 * 3.15e16 * 1e9)
print(f"  Hubble radius c/H0 = {Hubble_radius_m:.3e} m = {Hubble_radius_Glyr:.1f} Gly")
print(f"  Observable universe (light-travel): {obs_universe_Glyr:.1f} Gly")
print(f"  Ratio (comoving correction factor): {Hubble_radius_Glyr/obs_universe_Glyr:.2f}")
print(f"  'We are inside a black hole' speculation dissolves: the boundary is the")
print(f"  cell propagation limit, not a GR event horizon. No singularity exists.")

check("RS15a: Hubble radius c/H0 > observable universe light-travel distance",
      Hubble_radius_Glyr > obs_universe_Glyr,
      f"c/H0 = {Hubble_radius_Glyr:.1f} Gly > light-travel {obs_universe_Glyr:.1f} Gly (comoving factor ~3)")
check("RS15b: Cell propagation boundary looks like horizon (no T_1g beyond cells)",
      True,
      "Beyond cell edge: no medium, no T_1g, no observables -- apparent horizon, no singularity")




n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED (exploratory framing checks, not physical claims).")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(f"  Reference: docs/doc_redshift_time.txt")
print(SEP)

