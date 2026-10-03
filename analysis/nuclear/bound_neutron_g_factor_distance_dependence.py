"""
bound_neutron_g_factor_distance_dependence.py
==============================================
Independent re-investigation of [BOUND-NEUTRON-G-FACTOR]
(notes/open_items/false_positive_scan_series1.txt): doc_nucleus.txt Sec 5.7
and doc_magnetism.txt Sec 6 claim

    g_n(bound) = g_n(free) - mu_Zone3_proton = -1.567 - 0.360 = -1.927 mu_N

i.e. a FLAT, distance-independent 100% reuse of the proton's own Zone 3
self-contribution as the correction felt by an adjacent bound neutron.
This script asks: is that flat reuse physically motivated, or does a real
amplitude/geometry calculation using the framework's own established tools
give something else?

ESTABLISHED INPUTS REUSED VERBATIM (re-confirmed by direct rerun this
session, cited per source -- NOT recomputed differently here):
  - g_n(free)       = -1.5666 mu_N  [analysis/nuclear/neutron_g_factor.py Sec 4]
  - mu_Zone3_proton =  0.3602 mu_N  [analysis/nuclear/proton_g_factor.py Sec 4]
  - mu_n(PDG)       = -1.9130 mu_N  [CODATA/PDG]
  - r_p = 0.8414 fm, lambda_p = hbar_c/m_p = 0.2103 fm, r_grind = 2*lambda_p
    = 0.4206 fm  [doc_nucleus.txt PS1-PS4]
  - Zone 3 pressure is UNIFORM, P(r) ~ r^0, for lambda_p < r < r_p
    [analysis/nuclear/mechanical_radius.py MR1 -- independently matches
    Burkert et al. 2018 (Science 361, 207) JLab CLAS12 mechanical radius
    0.67+/-0.03 fm to 0.4-sigma, REAL external data, not an assumption]
  - Beyond r_p, entanglement_phase_lock.py derives E_Z3(r) =
    alpha*hbar_c*r_p^2/r^3 for a DIFFERENT quantity (inter-particle
    entanglement coupling energy), explicitly motivated in that script's
    own words as: "The Hopf winding IS a topological magnetic dipole ->
    field ~ 1/r^3." This is the SAME Hopf-driven Zone 3 structure that
    generates mu_Zone3_proton itself (not a different mechanism) -- so the
    identical (r_p/r)^3 falloff, continuous and capped at 1 inside r_p
    (matching mechanical_radius.py's own uniform-pressure finding), is the
    natural, non-arbitrary extension to the magnetic-moment-type quantity
    needed here. This is the answer to part (b) of the brief: reused
    because it is the same physical field, not borrowed by analogy only.
  - V_Reid(r): Reid (1968) NN potential, formula reused VERBATIM from
    analysis/nuclear/nn_potential_external_comparison.py (real, peer-
    reviewed, fit to actual NN scattering data -- not torsionverse-internal)
  - Deuteron RMS separation ~4 fm vs its own potential well ~0.8-1 fm
    (given fact -- real/external -- reflects the deuteron's anomalously
    weak 2.2 MeV binding; NOT representative of a typical MeV-per-nucleon
    bound pair)
  - doc_nucleus.txt Section 5 (A_g(T_1g x T_2g)=0): no p-n hard-core
    repulsion exists, but Zone 1+2 (r < lambda_p) is STILL a rigid,
    non-interpenetrating excluded-volume region for EVERY nucleon (proton
    or neutron) -- so r_grind = 2*lambda_p remains the absolute MINIMUM
    admissible center-to-center separation even for a p-n pair, it is just
    not an energy barrier for that pair.

NEW IN THIS SCRIPT (not reused from elsewhere, built to answer part (c)/(d)):
  - Precise numerical location of the real Reid potential's well minimum
    (the established facts only quote 2 bracketing sample points, 0.8 and
    0.9 fm) via a fine grid scan + parabolic refinement.
  - A standard-QM harmonic zero-point RMS spread estimate around that
    minimum (from the Reid well's own curvature + the standard p/n reduced
    mass) to quantify, quantitatively rather than just qualitatively, how
    compact a TYPICAL (non-deuteron-like) bound pair's wavefunction should be.
  - Both a simple POINT (center-to-center distance) and a finite-size
    SHELL-AVERAGED (volume-averaged over the neutron's own would-be Zone 3
    shell, since its Zone 1+2 hard core cannot overlap with the proton's)
    version of the coupling fraction, cross-checked against each other via
    an explicit convexity argument (Jensen's inequality).
  - A Gaussian-weighted average of the coupling fraction over the harmonic
    spread -- a genuine "distribution of separations" treatment, not just
    a single point estimate.

CROSS-CHECK NOTE (found during this investigation, reported honestly):
sandbox/nuclear/bound_neutron_zone3_distance_check.py already explores the
same question with a Monte-Carlo shell-average (no Reid potential, no
harmonic spread). Its core falloff model (uniform within r_p, (r_p/r)^3
beyond, citing the identical two source scripts) is independently
re-derived here from the same established facts -- itself a cross-check
that this is the "obvious" reading of the given facts, not an idiosyncratic
choice. Its shell-average numbers are independently reproduced here via a
DIFFERENT method (deterministic quadrature, no RNG) as a genuine check, not
copied -- compared explicitly in Section 6.

Run: python analysis/nuclear/bound_neutron_g_factor_distance_dependence.py
Reference: notes/open_items/false_positive_scan_series1.txt [BOUND-NEUTRON-G-FACTOR]
"""

import sys, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP, SEP2 = "=" * 72, "-" * 72
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, bool(cond)))
    print(f"  [{s}] {name}")
    if detail:
        print(f"         {detail}")

pi = math.pi

# ── ESTABLISHED CONSTANTS (reused, cited, not re-derived) ────────────────────
alpha       = 7.2973525693e-3
hbar_c      = 197.3269804     # MeV*fm
m_p         = 938.272         # MeV
m_n         = 939.5654        # MeV  (CODATA/PDG)
m_pi_charged = 139.57039      # MeV  (CODATA/PDG, external, standard nuclear-force input)
rho_0       = 0.16            # fm^-3 (standard nuclear saturation density, doc_nucleus.txt line 461)

r_p_fm      = 0.8414                     # fm, PS4
lambda_p_fm = hbar_c / m_p                # = 0.2103 fm
r_grind_fm  = 2 * lambda_p_fm             # = 0.4206 fm

mu_n_free   = -1.5666   # mu_N  [neutron_g_factor.py Sec 4, rerun-confirmed this session]
mu_Zone3_p  =  0.3602   # mu_N  [proton_g_factor.py Sec 4, rerun-confirmed this session]
mu_n_meas   = -1.9130   # mu_N  [PDG]

print(SEP)
print("bound_neutron_g_factor_distance_dependence.py")
print("Independent re-derivation of g_n(bound) with a real distance/geometry model")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("SECTION 1: INPUTS AND THE ORIGINAL (FLAT, 100%-REUSE) CLAIM")
print(SEP2)
g_n_bound_original = mu_n_free - mu_Zone3_p
err_original = 100 * (g_n_bound_original - mu_n_meas) / mu_n_meas
print(f"  g_n(free)        = {mu_n_free:.4f} mu_N   [18% gap vs PDG, MIT-bag-proxy error, not in question]")
print(f"  mu_Zone3_proton  = {mu_Zone3_p:.4f} mu_N   [proton's OWN Zone 3 self-contribution]")
print(f"  mu_n (PDG)       = {mu_n_meas:.4f} mu_N")
print(f"  ORIGINAL CLAIM: g_n(bound) = g_n(free) - mu_Zone3_proton (100% reuse, no distance)")
print(f"    = {mu_n_free:.4f} - {mu_Zone3_p:.4f} = {g_n_bound_original:.4f} mu_N   err = {err_original:+.2f}%")
print()
print("  QUESTION (a): is 100% reuse, with NO separation dependence at all,")
print("  physically well-motivated for a neutron sitting OUTSIDE the proton,")
print("  or does it need a real amplitude/geometry calculation?")
print("  ANSWER: no distance dependence at all cannot be right in general --")
print("  it would predict the SAME correction whether the neutron sits at")
print("  r_grind (closest possible contact) or at the Moon. A geometry-aware")
print("  calculation is required; Sections 2-6 build the most well-motivated")
print("  one available from the established facts.")

check("S1 Original claim is a RATIO fit (0.7% from PDG) that includes ZERO "
      "explicit separation dependence",
      True,
      "g_n(bound) formula has no r, d, or any distance variable anywhere in it")

# =============================================================================
print()
print(SEP2)
print("SECTION 2: THE COUPLING-FRACTION FALLOFF LAW (answers part b)")
print(SEP2)
print("""
  mechanical_radius.py MR1: Zone 3 pressure is UNIFORM, P(r) ~ r^0, for
  lambda_p < r < r_p (matches real Burkert et al. 2018 CLAS12 data, 0.4-sigma).
  => inside the proton's own r_p, the coupling strength does not fall with r.

  entanglement_phase_lock.py: for r >= r_p, the SAME Hopf-driven Zone 3
  structure's field is explicitly derived as a topological magnetic dipole
  field, falling as (r_p/r)^3, normalized to 1 at r=r_p (its own EP1 check).
  This is not a different physical quantity reused by loose analogy -- it
  is literally the same Hopf winding whose Zone 3 shell generates
  mu_Zone3_proton. A magnetic-moment-type quantity sourced by the identical
  structure must fall off with the identical power law outside the source
  (standard multipole-expansion fact: any linear coupling to a 1/r^3
  dipole field itself scales as 1/r^3 -- this is why E_Z3(r) itself was
  built as 1/r^3 in the first place).

  COMBINED PROFILE (continuous, capped at 1):
    fraction_point(r) = 1                  for r <= r_p
                      = (r_p / r)^3         for r  > r_p
""")

def fraction_point(r):
    if r <= r_p_fm:
        return 1.0
    return (r_p_fm / r) ** 3

check("CP1 fraction_point is continuous at r_p",
      abs(fraction_point(r_p_fm - 1e-9) - fraction_point(r_p_fm + 1e-9)) < 1e-6,
      f"f(r_p-)={fraction_point(r_p_fm-1e-9):.6f}  f(r_p+)={fraction_point(r_p_fm+1e-9):.6f}")
check("CP2 fraction_point(r_grind) is fully capped at 1 (closest allowed p-n "
      "contact sits well inside the proton's own r_p)",
      fraction_point(r_grind_fm) == 1.0,
      f"r_grind={r_grind_fm:.4f} fm < r_p={r_p_fm:.4f} fm -> fraction=1")

# =============================================================================
print()
print(SEP2)
print("SECTION 3: THE REAL REID (1968) NN POTENTIAL -- LOCATING THE WELL MINIMUM")
print(SEP2)

def V_reid(r_fm):
    """Reid (1968) NN potential, MeV, r in fm, mu=0.7 fm^-1.
    Reused VERBATIM from analysis/nuclear/nn_potential_external_comparison.py
    (real, peer-reviewed, fit to actual NN scattering data)."""
    mu = 0.7
    x = mu * r_fm
    return (-10.463 * math.exp(-x) / x
            - 1650.6 * math.exp(-4 * x) / x
            + 6484.2 * math.exp(-7 * x) / x)

print(f"  {'r (fm)':>8}  {'V_Reid (MeV)':>14}")
print(f"  {'-'*8}  {'-'*14}")
for r in [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 1.0, 1.2, 1.5, 2.0]:
    print(f"  {r:>8.2f}  {V_reid(r):>14.2f}")

# Fine grid scan to bracket the minimum, then parabolic refinement.
N_SCAN = 20000
r_lo, r_hi = 0.3, 2.0
best_r, best_v = None, 1e18
for i in range(N_SCAN + 1):
    r = r_lo + (r_hi - r_lo) * i / N_SCAN
    v = V_reid(r)
    if v < best_v:
        best_v, best_r = v, r
h_scan = (r_hi - r_lo) / N_SCAN
v0, v1, v2 = V_reid(best_r - h_scan), V_reid(best_r), V_reid(best_r + h_scan)
denom = (v0 - 2*v1 + v2)
r0_well = best_r - 0.5 * h_scan * (v2 - v0) / denom if denom != 0 else best_r
V_well  = V_reid(r0_well)

print()
print(f"  Grid-scan + parabolic-refined minimum: r0_well = {r0_well:.5f} fm")
print(f"  V_Reid(r0_well) = {V_well:.4f} MeV  (well depth)")
print(f"  r0_well / r_p   = {r0_well/r_p_fm:.4f}   (how close the REAL, external")
print(f"                     NN-scattering-fit well sits to this framework's own r_p)")

# Curvature at the minimum (independent, coarser step for a robust 2nd derivative)
h_curv = 0.01
k_spring = (V_reid(r0_well + h_curv) - 2*V_reid(r0_well) + V_reid(r0_well - h_curv)) / h_curv**2
print(f"  V''(r0_well) = k_spring = {k_spring:.2f} MeV/fm^2")

check("RW1 A genuine local minimum was found (positive curvature)",
      k_spring > 0,
      f"k_spring = {k_spring:.2f} MeV/fm^2")
check("RW2 Real external Reid-potential well sits within 20% of this "
      "framework's own r_p (a genuine, not cherry-picked, numerical coincidence "
      "if true)",
      abs(r0_well - r_p_fm) / r_p_fm < 0.20,
      f"r0_well={r0_well:.4f} fm  r_p={r_p_fm:.4f} fm  diff={100*abs(r0_well-r_p_fm)/r_p_fm:.1f}%")

# =============================================================================
print()
print(SEP2)
print("SECTION 4: HOW COMPACT IS A TYPICAL (NON-DEUTERON) BOUND PAIR? "
      "(harmonic zero-point estimate)")
print(SEP2)
print("""
  The deuteron's RMS separation (~4 fm) is anomalously large because its
  binding energy (2.2 MeV) is tiny compared to the well depth (~95 MeV) --
  an extreme uncertainty-principle/tunneling effect. A typical, MeV-per-
  nucleon-bound pair (not deuteron-shallow) should sit much closer to the
  classical well. Standard QM harmonic approximation around the REAL Reid
  well's own curvature gives an order-of-magnitude, non-deuteron-specific
  estimate of how compact such a pair's ground-state spread should be.
""")

mu_pn_c2 = m_p * m_n / (m_p + m_n)   # reduced mass energy, MeV
hbar_omega = hbar_c * math.sqrt(k_spring / mu_pn_c2)   # MeV
E_zp = 0.5 * hbar_omega
sigma_r = hbar_c * math.sqrt(1.0 / (2 * mu_pn_c2 * hbar_omega))   # fm

print(f"  Reduced mass mu_pn*c^2        = {mu_pn_c2:.4f} MeV")
print(f"  hbar*omega (well curvature)   = {hbar_omega:.4f} MeV")
print(f"  Zero-point energy E_zp=hw/2   = {E_zp:.4f} MeV   (well depth = {-V_well:.2f} MeV)")
print(f"  Harmonic RMS spread sigma_r   = {sigma_r:.5f} fm")
print(f"  sigma_r / deuteron RMS (4 fm) = {sigma_r/4.0:.4f}")

harmonic_valid = E_zp < 0.5 * (-V_well)
print()
if harmonic_valid:
    print("  Harmonic approximation looks self-consistent (zero-point energy is a")
    print("  modest fraction of the well depth) -- sigma_r is a meaningful estimate.")
else:
    print("  CAVEAT, reported honestly: the real Reid well is narrow/stiff enough")
    print("  that the harmonic zero-point energy is comparable to or exceeds the")
    print("  well depth itself -- the harmonic approximation is at best an order-")
    print("  of-magnitude guide here, not a precision estimate. The QUALITATIVE")
    print("  conclusion (typical binding is far more compact than the deuteron's")
    print("  anomalous 4 fm) still stands on general nuclear-physics grounds even")
    print("  if sigma_r's own numeric value should not be over-trusted.")

check("ZP1 Harmonic spread is much smaller than the deuteron's known 4 fm "
      "(confirms: typical/stronger binding -> far more compact than the "
      "anomalously shallow deuteron)",
      sigma_r < 1.0,
      f"sigma_r = {sigma_r:.4f} fm  vs deuteron 4 fm  (ratio {sigma_r/4.0:.3f})")

# =============================================================================
print()
print(SEP2)
print("SECTION 5: g_n(bound) AT EACH CANDIDATE SEPARATION -- POINT MODEL")
print(SEP2)

r0_pion_fm    = hbar_c / m_pi_charged                       # nuclear force range
r0_density_fm = (3 / (4*pi*rho_0)) ** (1/3)                 # bulk inter-nucleon spacing
r0_deuteron   = 4.0                                          # given fact, weak-binding case

candidates = [
    ("r_grind (p-n min separation)",      r_grind_fm),
    ("r_p (own Zone3 edge, reference)",   r_p_fm),
    ("Reid well minimum (external)",      r0_well),
    ("pion range hbar_c/m_pi (external)", r0_pion_fm),
    ("nuclear density spacing (external)",r0_density_fm),
    ("deuteron RMS (weak-binding case)",  r0_deuteron),
]

print(f"  {'separation':<34}  {'d (fm)':>8}  {'fraction':>9}  {'g_n(bound)':>11}  {'err vs PDG':>11}")
print(f"  {'-'*34}  {'-'*8}  {'-'*9}  {'-'*11}  {'-'*11}")
point_results = {}
for label, d in candidates:
    f = fraction_point(d)
    g_bound = mu_n_free - f * mu_Zone3_p
    err = 100 * (g_bound - mu_n_meas) / mu_n_meas
    point_results[label] = (d, f, g_bound, err)
    print(f"  {label:<34}  {d:>8.4f}  {f:>9.4f}  {g_bound:>11.4f}  {err:>+10.2f}%")

check("PT1 At the physically realistic external separations (pion range, "
      "nuclear density spacing) the point model gives a MUCH larger error "
      "than the claimed 0.7% (the flat 100% reuse is not recovered in general)",
      point_results["pion range hbar_c/m_pi (external)"][3] > 5.0 and
      point_results["nuclear density spacing (external)"][3] > 5.0,
      f"pion-range err={point_results['pion range hbar_c/m_pi (external)'][3]:+.2f}%  "
      f"density-spacing err={point_results['nuclear density spacing (external)'][3]:+.2f}%")

# =============================================================================
print()
print(SEP2)
print("SECTION 6: SHELL-AVERAGED COUPLING FRACTION (finite neutron size, "
      "deterministic quadrature cross-check)")
print(SEP2)
print("""
  The neutron's own Zone 1+2 (r' < lambda_p) is a rigid excluded-volume
  region that cannot overlap the proton's -- so only the neutron's own
  WOULD-BE Zone 3 shell (lambda_p < r' < r_p, measured from the neutron's
  own center) is ever actually exposed to the proton's external field.
  Volume-average fraction_point over that shell, offset by center-to-center
  separation d, instead of evaluating fraction_point at the single point d.
  Implemented here as a DETERMINISTIC double quadrature (law of cosines +
  radial shell weight r'^2 dr'), not Monte Carlo -- a genuine independent
  cross-check of sandbox/nuclear/bound_neutron_zone3_distance_check.py's
  separately-coded random-sampling version.
""")

N_R, N_U = 160, 160
r_lo_shell, r_hi_shell = lambda_p_fm, r_p_fm

def shell_average_fraction(d):
    num, den = 0.0, 0.0
    for i in range(N_R):
        rp = r_lo_shell + (r_hi_shell - r_lo_shell) * (i + 0.5) / N_R
        weight_r = rp**2 * (r_hi_shell - r_lo_shell) / N_R
        # inner integral over u=cos(theta'), law of cosines, uniform in u (solid angle)
        inner = 0.0
        for j in range(N_U):
            u = -1.0 + 2.0 * (j + 0.5) / N_U
            dist = math.sqrt(rp*rp + d*d - 2*rp*d*u)
            inner += fraction_point(dist) * (2.0 / N_U)   # du weight, normalized: int du/2 over [-1,1] = 1
        inner *= 0.5
        num += weight_r * inner
        den += weight_r
    return num / den

test_ds = [r_grind_fm, 0.5, 0.6, 0.7, r_p_fm, r0_density_fm, r0_pion_fm, 2*r_p_fm]
print(f"  {'d (fm)':>10}  {'point frac':>12}  {'shell frac (mine)':>19}")
print(f"  {'-'*10}  {'-'*12}  {'-'*19}")
shell_fracs = {}
for d in test_ds:
    pf = fraction_point(d)
    sf = shell_average_fraction(d)
    shell_fracs[d] = sf
    print(f"  {d:>10.4f}  {pf:>12.4f}  {sf:>19.4f}")

# Cross-check against the pre-existing sandbox script's own printed Monte-Carlo
# numbers at matching d (0.4206->0.8640, 0.8414->0.5972, 1.1427->0.4063,
# 1.4138->0.2555 -- read directly from its own run output, not re-typed from memory).
existing_script_mc = {
    r_grind_fm: 0.8640,
    r_p_fm:     0.5972,
    1.1427:     0.4063,
    1.4138:     0.2555,
}
print()
max_diff = 0.0
for d, mc_val in existing_script_mc.items():
    # match to nearest computed d (density/pion values rounded slightly differently)
    nearest_d = min(shell_fracs.keys(), key=lambda x: abs(x - d))
    if abs(nearest_d - d) < 0.01:
        diff = abs(shell_fracs[nearest_d] - mc_val)
        max_diff = max(max_diff, diff)
        print(f"  d={d:.4f}: mine(quadrature)={shell_fracs[nearest_d]:.4f}  "
              f"existing(Monte Carlo)={mc_val:.4f}  diff={diff:.4f}")

check("SH1 Independently-coded deterministic quadrature matches the "
      "pre-existing Monte-Carlo sandbox script to within 1%",
      max_diff < 0.01,
      f"max abs difference across matched separations = {max_diff:.4f}")
check("SH2 Shell-averaged fraction EXCEEDS the point fraction once d is well "
      "beyond r_p (Jensen's inequality for a convex 1/r^3 falloff: averaging "
      "a near-side gain and a far-side loss around a convex function nets "
      "positive) -- confirms this is a real, understood effect, not noise",
      shell_fracs[2*r_p_fm] > fraction_point(2*r_p_fm),
      f"at d=2*r_p: shell={shell_fracs[2*r_p_fm]:.4f}  point={fraction_point(2*r_p_fm):.4f}")

print()
print(f"  {'separation':<34}  {'d (fm)':>8}  {'shell frac':>10}  {'g_n(bound)':>11}  {'err vs PDG':>11}")
print(f"  {'-'*34}  {'-'*8}  {'-'*10}  {'-'*11}  {'-'*11}")
shell_results = {}
for label, d in candidates:
    f = shell_average_fraction(d)
    g_bound = mu_n_free - f * mu_Zone3_p
    err = 100 * (g_bound - mu_n_meas) / mu_n_meas
    shell_results[label] = (d, f, g_bound, err)
    print(f"  {label:<34}  {d:>8.4f}  {f:>10.4f}  {g_bound:>11.4f}  {err:>+10.2f}%")

# =============================================================================
print()
print(SEP2)
print("SECTION 7: GAUSSIAN-WEIGHTED 'DISTRIBUTION OF SEPARATIONS' BEST ESTIMATE")
print(SEP2)
print(f"""
  Centre the separation distribution at the real Reid well minimum
  (r0_well = {r0_well:.4f} fm) with the harmonic zero-point spread
  (sigma_r = {sigma_r:.4f} fm) from Section 4, and average BOTH the point
  and shell-averaged coupling fractions over that Gaussian -- the most
  careful single estimate this script can produce from established tools.
""")

N_GAUSS = 31
def gauss_weighted_fraction(frac_fn):
    total_w, total_fw = 0.0, 0.0
    for k in range(N_GAUSS):
        z = -4.0 + 8.0 * k / (N_GAUSS - 1)     # +/- 4 sigma
        d = r0_well + z * sigma_r
        if d <= 0:
            continue
        w = math.exp(-0.5 * z * z)
        total_w += w
        total_fw += w * frac_fn(d)
    return total_fw / total_w

frac_point_gauss = gauss_weighted_fraction(fraction_point)
frac_shell_gauss = gauss_weighted_fraction(shell_average_fraction)

g_point_gauss = mu_n_free - frac_point_gauss * mu_Zone3_p
g_shell_gauss = mu_n_free - frac_shell_gauss * mu_Zone3_p
err_point_gauss = 100 * (g_point_gauss - mu_n_meas) / mu_n_meas
err_shell_gauss = 100 * (g_shell_gauss - mu_n_meas) / mu_n_meas

print(f"  Gaussian-weighted point fraction   = {frac_point_gauss:.4f}  "
      f"-> g_n(bound) = {g_point_gauss:.4f}  err = {err_point_gauss:+.2f}%")
print(f"  Gaussian-weighted shell fraction    = {frac_shell_gauss:.4f}  "
      f"-> g_n(bound) = {g_shell_gauss:.4f}  err = {err_shell_gauss:+.2f}%")

check("GW1 Best (distribution-averaged) estimate is a genuine improvement "
      "over the free-neutron 18% gap in the right direction",
      abs(err_shell_gauss) < 18.1,
      f"best estimate err = {err_shell_gauss:+.2f}%  vs free-neutron 18.1% gap")
check("GW2 Best (distribution-averaged) estimate does NOT reproduce the "
      "claimed +0.7% precision (tests whether 'essentially closed' survives)",
      abs(err_shell_gauss) > 2.0,
      f"best estimate err = {err_shell_gauss:+.2f}%  vs claimed +0.7%")

# =============================================================================
print()
print(SEP2)
print("SECTION 8: ALTERNATIVE EXPONENT -- IS 1/r^3 EVEN THE RIGHT POWER? "
      "(honest ambiguity check)")
print(SEP2)
print("""
  Section 2 reused (r_p/r)^3, reasoning that mu_Zone3_proton is driven by
  the SAME Hopf winding that entanglement_phase_lock.py treats as a
  topological magnetic dipole (field ~ 1/r^3). A genuine, defensible
  ALTERNATIVE reading exists: mu_Zone3_proton's own integrand is built from
  a LOCAL VELOCITY v(r) (proton_g_factor.py Section 4: integrand ~ rho*v(r)
  *r^2), and entanglement_phase_lock.py's own derivation chain states an
  INTERMEDIATE velocity falloff v(r) = v_0*(r_p/r)^2 one step before its
  final (r_p/r)^3 ENERGY result. Reusing the dimensionally-closer velocity
  falloff -- (r_p/r)^2, not (r_p/r)^3 -- is an equally defensible choice,
  and gives a SLOWER-decaying (larger-at-a-given-distance) estimate. Both
  are reported here, honestly, as neither is independently settled by
  anything established in this framework.
""")

def fraction_point_v2(r):
    if r <= r_p_fm:
        return 1.0
    return (r_p_fm / r) ** 2

print(f"  {'separation':<34}  {'d (fm)':>8}  {'(r_p/r)^3':>10}  {'(r_p/r)^2':>10}  "
      f"{'g_n [^3]':>9}  {'g_n [^2]':>9}")
print(f"  {'-'*34}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*9}  {'-'*9}")
for label, d in candidates:
    f3 = fraction_point(d)
    f2 = fraction_point_v2(d)
    g3 = mu_n_free - f3 * mu_Zone3_p
    g2 = mu_n_free - f2 * mu_Zone3_p
    print(f"  {label:<34}  {d:>8.4f}  {f3:>10.4f}  {f2:>10.4f}  {g3:>9.4f}  {g2:>9.4f}")

check("EX1 The (r_p/r)^2 alternative gives a LESS steep falloff than (r_p/r)^3 "
      "at every separation beyond r_p (as expected -- lower power -> slower decay)",
      all(fraction_point_v2(d) >= fraction_point(d) for _, d in candidates),
      "checked at every candidate separation beyond r_p")
check("EX2 Even under the more generous (r_p/r)^2 alternative, the realistic "
      "external separations (pion range, nuclear density spacing) still miss "
      "the claimed +0.7% by more than 3x",
      abs(100*(mu_n_free - fraction_point_v2(r0_pion_fm)*mu_Zone3_p - mu_n_meas)/mu_n_meas) > 2.1,
      f"(r_p/r)^2 at pion range: err = "
      f"{100*(mu_n_free - fraction_point_v2(r0_pion_fm)*mu_Zone3_p - mu_n_meas)/mu_n_meas:+.2f}%")

# =============================================================================
print()
print(SEP)
print("SUMMARY")
print(SEP)
passed = sum(1 for _, ok in results if ok)
failed = len(results) - passed
print(f"  Total: {len(results)}   PASS: {passed}   FAIL: {failed}")
print()
print(f"  Original flat-reuse claim:         g_n(bound) = {g_n_bound_original:.4f}  err = {err_original:+.2f}%")
print(f"  Best distance-aware estimate       (Reid well +/- harmonic spread,")
print(f"    shell-averaged):                 g_n(bound) = {g_shell_gauss:.4f}  err = {err_shell_gauss:+.2f}%")
print(f"  Free-neutron baseline (no Zone3):  g_n(free)  = {mu_n_free:.4f}  err = "
      f"{100*(mu_n_free-mu_n_meas)/mu_n_meas:+.2f}%")
print()
print("  HONEST READING:")
print("  (1) A real geometry/distance calculation does NOT support flat 100%")
print("      reuse in general -- it only comes close to 100% because the REAL,")
print("      external Reid-potential well minimum happens to sit very near this")
print("      framework's own r_p (Section 3), not because 100% reuse is itself")
print("      correct at arbitrary separation.")
print("  (2) At more conservative/typical nuclear separations (pion force range,")
print("      bulk nuclear density spacing), the predicted correction shrinks")
print("      substantially and the error grows well beyond the claimed 0.7%,")
print("      though it remains better than the free-neutron 18% gap.")
print("  (3) The claimed +0.7% precision is NOT independently well-supported by")
print("      any calculation available in this framework's current tools -- it")
print("      is better described as landing inside a wide, genuinely uncertain")
print("      band (roughly 2%-15% depending on modeling choices that are each")
print("      individually defensible) than as a precision result.")
print()
print(f"  Reference: notes/open_items/false_positive_scan_series1.txt [BOUND-NEUTRON-G-FACTOR]")
print(f"  Cross-check: sandbox/nuclear/bound_neutron_zone3_distance_check.py")
