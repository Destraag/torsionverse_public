Analysis Scripts — Torsion Universe Framework
========================================================

This directory contains computational analysis scripts supporting the
theoretical work in whitepaper.txt and medium_chains.txt.

All scripts require only Python's standard library (no pip installs).
Run from the PROJECT ROOT (torsionverse/):

  ALPHA DERIVATION CHAIN:
    python analysis/alpha/wave_path_test.py
    python analysis/alpha/epsilon_search.py
    python analysis/alpha/hopf_stability.py
    python analysis/alpha/hopf_c4b_correction.py

  MEDIUM / SCALE TESTS:
    python analysis/medium/scale_check.py
    python analysis/medium/kappa_sensitivity.py
    python analysis/medium/medium_properties.py

  GRAVITY / ASTROPHYSICAL:
    python analysis/gravity/flyby_anomaly.py
    python analysis/gravity/pioneer_anomaly.py
    python analysis/gravity/mond_rotation.py

  COSMOLOGICAL:
    python analysis/cosmo/gw170817_check.py
    python analysis/cosmo/electroweak_check.py

  MISC:
    python analysis/misc/whitepaper_audit.py

FOLDER STRUCTURE:
  analysis/
    constants.py           Shared module — imported by all scripts in subdirs
    README.txt             This file
    alpha/                 C4a → C4b → wave path → stability derivation chain
    medium/                Rs, medium properties, scale tests
    gravity/               MOND, Pioneer, flyby anomaly, rotation curves
    cosmo/                 GW, electroweak, aurora, pion-nucleon sigma
    misc/                  Audit tools, exploratory scripts

IMPORT NOTE:
  Scripts in subdirectories import constants.py via:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
    from constants import *
  constants.py stays in analysis/ (parent of all subdirs).

========================================================
SHARED MODULE
========================================================

constants.py
  Shared physical constants for all analysis scripts.
  CODATA 2018 / PDG 2022 values, plus framework-derived quantities.
  Key exports: Rs_exact (= sqrt(5)/(4*pi)), alpha, c_m_s, H0_s, cH0,
  G, M_sun_kg, M_E_kg, R_E_m, omega_E, AU_m, m_b_GeV, sigma_piN_MeV,
  kappa_GeV_per_fm, a0_m_s2, R_nuclear, R_hadronic, R_galactic,
  R_s_mean, SEP.
  Import with: from constants import *

========================================================
SCRIPTS
========================================================

scale_check.py  [Parts 1-6 — core R_s derivation]
  Cross-scale dimensionless ratio analysis for the medium constitutive
  law hypothesis. Identifies R_s = sqrt(5)/(4*pi) = 0.17794 appearing
  at nuclear, hadronic, and galactic scales.

  Key findings:
    - R2/R1 = 2*pi EXACTLY (particle crossing-ring to full EM loop)
    - m_b = 4.18 GeV identified as hadronic saturation energy (0.8% off)
    - sigma_piN = 59.1 MeV identified as nuclear saturation energy (2.7% off)
    - R_s = sqrt(5)/(4*pi): three-scale cluster B mean 0.17846 +/- 0.00199
    - a0 = R_s * c * H0  (MOND formula derived, not fitted; 0.8% accuracy)
    - B meson radius prediction: r_B ~ 0.24 fm (testable with lattice QCD)
    - Sensitivity analysis: R_s stable to 0.5% under input variation
    - Cosmological evolution: R_s ratios constant if H0 alpha c unchanged

  Source: medium_chains.txt Sections 9-14; whitepaper.txt Appendix E

flyby_anomaly.py  [Parts A-B — spin-loading and rotation saturation]
  Torsion spin-loading prediction for the Earth flyby anomaly +
  wrapper closure / rotation saturation threshold analysis.

  Key findings:
    - Anderson K formula verified against six flyby events (RMS < 1 mm/s)
    - c_torsion implied from flyby data: c_torsion = R_s * c = 0.1779c
      (R_s appears as torsion propagation speed ratio -- 4th appearance)
    - Saturation condition: v_rot / v_esc >= R_s for full wrapper closure
    - Earth: 23% saturated; Jupiter/Saturn: fully saturated
    - K_max (full saturation) ~ 4.3x Anderson K
    - PREDICTIONS: Jupiter/Saturn flyby K ~ 4x Earth; Venus K ~ 0
    - CAVEAT: Polar geometry may shift c_torsion estimate; see Sec 13

  Source: medium_chains.txt Sections 13-14

pioneer_anomaly.py  [Pioneer check — open problem]
  Tests whether torsion mass-loading can explain the Pioneer anomaly
  (a_P = 8.74e-10 m/s^2, constant 20-70 AU).

  Key findings:
    - F1 (R_s * g_Newton):          FAIL -- thousands of times too large
    - F2 (R_s * sqrt(g * cH0)):     FAIL -- r-dependent, not constant
    - F3 (standard MOND):            FAIL -- negligible at Pioneer distances
    - F4 (kinematic wake):           OPEN -- needs medium density derivation
    - Descriptive: a_P = 7.22 * R_s * c * H0  (accurate, not derived)
    - PATH 1 (thermal recoil): Turyshev 2012 explains ~75%; may be complete
    - PATH 2 (torsion wake): Requires medium density from QCD + R_s
    - PATH 3 (modified MOND mu): Unlikely -- Pioneer in deep Newtonian regime
    - RECOMMENDED: Determine if thermal fully accounts for Pioneer first.

  Source: medium_chains.txt Section 14, Open Question Q1

mond_rotation.py  [MOND rotation curve test — C3 verification]
  Tests MOND formula with a0 = R_s * c * H0 (no free parameter) against
  galaxy rotation data. Four test parts:
    A. a0 framework vs SPARC/RAR measured a0 (153-galaxy fit)
    B. Baryonic Tully-Fisher for 12 galaxies from the literature
    C. BTF slope regression and gas-subsample statistics
    D. H0 tension: what H0 does the a0 measurement imply?
    E. Milky Way spot check (MOND borderline regime)

  Key findings:
    - a0 = R_s * c * H0 = 1.210e-10 m/s^2  (H0 = 70.0)
    - SPARC/RAR measured a0 = 1.20 +/- 0.02e-10 m/s^2
    - Difference: +0.85%  (0.51 sigma) -- STRONG SUPPORT for C3
    - Gas-dominated BTF geometric mean: ~14% above prediction (M/L uncertainty)
    - H0 implied by a0 measurement: 69.4 km/s/Mpc (between Planck and SH0ES)
    - Milky Way: v_mond = 200 km/s vs 232 km/s observed (-14%; M(<R0) uncertain)
    - NEXT: Download SPARC.dat (Lelli+2016) for definitive 153-galaxy test

  Source: whitepaper.txt Appendix E (C3), medium_chains.txt Section 12

========================================================
ALPHA DERIVATION SERIES  (C4a → C4b → wave path → C4c)
========================================================

These scripts form a progressive derivation chain for the fine-structure
constant alpha from first principles (Hopf geometry + medium constants).
Run in order for the full derivation narrative.

hopf_c4.py  [early C4a exploration — superseded]
  Initial Hopf geometry exploration leading to C4a.
  Superseded by hopf_c4_correction.py and c4a_candidates.py.

hopf_c4_phi_hit.py  [phi/pi combinations scan — superseded]
  Brute-force search for alpha from phi/pi combinations.
  Superseded by c4a_candidates.py.

c4a_numbers.py  [C4a reference numbers]
  Computes C4a = sqrt(5)*phi/(16*pi^3) and related quantities.
  Reference values: alpha_C4a, 1/alpha_C4a, gap from CODATA.

c4a_candidates.py  [C4a → C4b discovery]
  Systematic scan discovers C4b quadratic: 2*alpha^2-(4*pi^2/phi)*alpha+Rs=0.
  Shows C4b is 100x closer to CODATA than C4a (0.00056% vs 0.060%).
  Key: n=2 coefficient identified as winding number, not fitted.
  Source: whitepaper.txt Appendix D

alpha_precision_check.py  [C4b precision audit]
  Confirms C4b is 37,127 sigma from CODATA — cannot replace CODATA for
  precision work. G-2 prediction fails by 23,157 sigma with alpha_C4b.
  Purpose: keep the claim honest (C4b is a conjecture, not a replacement).

hopf_c4_correction.py  [C4a gap analysis — 4 parts]
  Four-part analysis of the C4a gap (-4.406e-6).
  Parts A-D: gap not Schwinger-scale; no clean first-order correction;
  prior art positioning; four rigorous-derivation criteria (D1-D4).
  Conclusion: Step D2 (Hopf linking integral) is the mathematical target.
  Source: whitepaper.txt Appendix D, c4a_theory.txt Parts 0-0b

hopf_c4b_correction.py  [C4b double-spin analysis — 4 parts]
  Parallel four-part analysis for C4b.
  Part C: coefficient scan over 13 geometric values confirms n=2 is unique
  best-fit integer: closest to n_exact=2.0187 by factor 54.
  Part D: double-spin hypothesis — figure-8 rotates around horizontal axis,
  making crossing ring a (1,2) torus knot with linking number 2.
  n_exact = 2.01869; residual = 0.01869.
  Source: c4a_theory.txt Part 0b, whitepaper.txt Appendix D Step 2b

c4b_residual_scale.py  [residual 0.01869 geometric search]
  Tests geometric candidates for n_exact - 2 = 0.01869.
  Near-hits: Rs/pi^2 (-3.5%), alpha*phi^2 (+2.2%). None definitive.
  Source: c4a_theory.txt Part 0b

c4b_residual_medium.py  [residual 0.01869 medium search]
  Tests solar system bodies' saturation fraction sat_frac = v_rot/v_esc/Rs.
  Sun: sat_frac_Sun = 0.01812 (-3.1% from residual).
  Only solar system body within 4%; all others 78-8000% off.
  Kickback hypothesis: local medium coupling efficiency sets fractional n.
  Source: c4a_theory.txt Part 0b, whitepaper.txt Appendix D

pulsar_double_spin.py  [cross-scale n=2 verification]
  PSR B1828-11 as macroscopic (1,2) torus knot.
  Modulation ratio P1/P2 = 1.996 (0.07 sigma from n=2).
  Beam traces (1,2) torus knot in angle-space per precession cycle.
  Two-scale n=2: quantum (C4b) and pulsar both select n=2 topology.
  Needs ~3x better timing precision to resolve n_exact vs n=2.
  Source: c4a_theory.txt Part 0b, whitepaper.txt Appendix D Step 2b

wave_path_test.py  [wave path mechanism — C4b residual closure]
  Tests whether resonant wave on (1,2) torus knot explains residual.
  Model: phi = 2*theta + epsilon*sin(k*theta) on Hopf torus.
  k=2 resonant wave at epsilon=0.11938 gives n_EM = n_exact = 2.01869 exactly.
  Direction verified: wave increases arc length → n_EM > 2 (correct direction).
  Architecture complete: topology (int 2) + geometry (frac 0.01869) = n_exact.
  Source: c4a_theory.txt Part 0c, whitepaper.txt Appendix D Step 2c

epsilon_search.py  [closed-form search for wave amplitude]
  High-precision (N=50000) determination of epsilon = 0.11937954.
  Best closed-form candidate: 3/(8*pi) = (p+q)/(4*R2) = 0.11936621 (0.011% off).
  Geometric reading: (1+2)/(4*2*pi) = sum of winding numbers / (4*major radius).
  Alpha error with 3/(8*pi): -6.7e-8% (vs C4b n=2: -0.000560%).
  Cross-scale: epsilon = sin(6.86 deg); Mercury inclination 7.00 deg (+2.1%).
  Source: c4a_theory.txt Part 0c, whitepaper.txt Appendix D Step 2c

hopf_stability.py  [stability analysis: is smooth (1,2) path stable?]
  Tests whether the smooth (1,2) torus knot is stable or unstable under
  k=2 resonant wave perturbations using the curvature energy (LIA proxy
  for Biot-Savart self-interaction).
  Method: computes d²E_curv/deps² at eps=0 via Frenet-Serret curvature.
  Key result: d²E_curv/deps² = -15.1 < 0  →  UNSTABLE.
  The smooth path is a local maximum of curvature energy.
  The k=2 wave arises SPONTANEOUSLY — not externally driven.
  Equilibrium eps* found by minimising E_eff = E_curv + lambda*(n_EM-n_exact)^2.
  Finding (2026-08-16): curvature minimum (eps~0.153) and winding minimum
  (eps~0.119) do not coincide → LIA is insufficient; full Biot-Savart
  self-energy must be minimised. This is a non-local integral problem,
  not a local PDE — sharpens the remaining step to a well-defined target.
  The remaining open calculation: derive lambda from the full Biot-Savart
  self-energy of a (1,2) torus knot current loop on the Hopf torus.
  This is a well-defined eigenvalue problem, not a free parameter search.
  Source: c4a_theory.txt Part 0c, whitepaper.txt Appendix D Step 2c

biot_savart_min.py  [Biot-Savart self-energy minimisation — route RULED OUT]
  Numerically evaluates Biot-Savart self-energy E_BS(eps) for the
  (1,2) torus knot with resonant wave at N=80 quadrature points.
  Finding: dL/deps|0 = -3.27 (arc length DECREASES), minimum at eps~0.3.
  E'_BS(0) = -31 (negative, drives eps upward), but N=80 too coarse for E''.
  Conclusion: Classical Biot-Savart route does NOT select eps=0.119.
  The Biot-Savart equilibrium requires unphysical core radius a~29*R1.
  Source: writhe_min.py (confirms), agenda.txt Item C

writhe_min.py  [Gauss writhe integral; Gap 1 ruled out; Gap 2 confirmed]
  Numerically evaluates the Gauss writhe integral (no regularisation needed)
  for the (1,2) torus knot with resonant wave amplitude eps in [0, 0.25].
  Finding: Wr(0) = -0.075 (negative, near-zero). Literature formula wrong
  for R1/R2=1/(2*pi). dWr/deps = 0 at eps~0.17, NOT at eps=0.119.
  GAP 1 CONCLUSION: Biot-Savart/writhe route ruled out as mechanism for
  selecting eps=0.119. Classical self-energy minimisation does not work.
  GAP 2 CONFIRMED: Rs = sqrt(p^2+q^2)/(4*pi) = geometric identity.
  Rs = norm of winding vector (1,2) divided by 4*pi. NOT empirical.
  Source: agenda.txt Item C, whitepaper.txt Appendix E status

hopf_linking_integral.py  [Gap 3: Hopf fiber linking integral — identities confirmed]
  Verifies that ALL THREE constants in the C4b quadratic derive from
  (p,q)=(1,2) and R2=2*pi via exact algebraic identities:
    phi = (1+sqrt(p^2+q^2))/2  [EXACT: golden ratio IS the (1,2) winding invariant]
    Rs  = sqrt(p^2+q^2)/(4*pi) [EXACT: matches writhe_min.py Gap 2 result]
    Q   = 2*R2^2/(1+sqrt(p^2+q^2)) = 4*pi^2/phi  [EXACT: coupling from Hopf scale]
  The C4b quadratic reduces to pure geometry:
    q*alpha^2 - [2*R2^2/(1+||(p,q)||)]*alpha + ||(p,q)||/(4*pi) = 0
  Inputs: (p,q)=(1,2), R2=2*pi. No free parameters whatsoever.
  phi, Rs, n=2 are ALL derived outputs, not assumed inputs.
  Alpha error from fully geometric formula: -0.000560% (same as C4b).
  The golden ratio appears in nature uniquely for (1,2) systems because
  ||w|| = sqrt(5) is the only winding vector whose norm gives phi = (1+sqrt5)/2.
  Cross-scale: other (p,q) pairs give phi-like constants, but none = phi.
  Gap 3 reduced to single calculation: show Q = CS[A] on S^3 for the
  (1,2) Hopf connection A (one differential-forms computation).
  Source: whitepaper.txt Appendix D Step 3, Appendix E status

========================================================
MEDIUM / ASTROPHYSICAL SERIES
========================================================

medium_properties.py
  Derives medium constitutive relations from Rs.
  Phase speeds, impedance, dispersion relations.

kappa_sensitivity.py
  Sensitivity analysis for the sigma_piN coupling constant kappa.
  Tests stability of nuclear-scale R_s identification.

a0_redshift.py
  Tests cosmological evolution of a0 = Rs*c*H0 with redshift.

aurora_test.py
  Aurora/magnetosphere test of the medium wave speed prediction.

gw170817_check.py
  Verifies v_gravitational = c from GW170817 (tensor wave speed).
  Checks v_shear = Rs*c prediction (scalar/longitudinal mode).

electroweak_check.py
  Tests medium model against electroweak scale (W/Z masses, Fermi constant).

sigma_piN.py
  Detailed analysis of sigma_piN = 59.1 MeV as nuclear saturation energy.
  Verifies R_s at nuclear scale.

sparc_btf.py
  SPARC dataset (Lelli+2016) baryonic Tully-Fisher analysis.
  Requires SPARC_Lelli2016c.mrt and SPARC_table1.dat (included).

pioneer_medium.py
  Medium density estimate from Pioneer anomaly constraint.
  Feeds into PATH 2 (kinematic wake) analysis.

new_ground.py
  Exploratory: new observational tests beyond the established conjectures.

wrapper_closure.py
  Rotation saturation / wrapper closure threshold analysis.
  Companion to flyby_anomaly.py Part B.

verify_constants.py
  Verifies constants.py values against CODATA/PDG sources.

whitepaper_audit.py
  Audits whitepaper.txt for internal consistency (numbers cross-check).

========================================================
FOLDER ORGANISATION (PROPOSED — not yet implemented)
========================================================

Current: all scripts flat in analysis/
Proposed split into subdirectories:

  analysis/
    alpha/           C4a/C4b/wave-path derivation chain
      c4a_numbers.py
      c4a_candidates.py
      alpha_precision_check.py
      hopf_c4.py
      hopf_c4_phi_hit.py
      hopf_c4_correction.py
      hopf_c4b_correction.py
      c4b_residual_scale.py
      c4b_residual_medium.py
      pulsar_double_spin.py
      wave_path_test.py
      epsilon_search.py
      hopf_stability.py

    medium/          Rs, medium properties, astrophysical tests
      constants.py
      medium_properties.py
      kappa_sensitivity.py
      scale_check.py
      verify_constants.py

    gravity/         MOND, Pioneer, flyby, rotation curves
      flyby_anomaly.py
      pioneer_anomaly.py
      pioneer_medium.py
      mond_rotation.py
      sparc_btf.py
      wrapper_closure.py
      a0_redshift.py

    cosmo/           GW, electroweak, cosmological
      gw170817_check.py
      electroweak_check.py
      aurora_test.py
      sigma_piN.py

    misc/
      new_ground.py
      whitepaper_audit.py
      scale_ratios.py
      scale_check.py

  NOTE: Reorganisation is a separate task. Scripts currently import
  from constants.py using relative paths; any move requires updating
  those imports. Discuss before implementing.
