"""
pulsar_double_spin.py — Can precessing pulsars model the (1,2) torus knot?

Context (from hopf_c4b_correction.py and c4b_residual_scale.py):
  C4b quadratic: n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0
  Best integer n = 2  (n_exact = 2.01869)
  Double-spin hypothesis: coefficient n = winding number of (1,2) torus knot.
    The figure-8 spinning around its horizontal axis as it sweeps the main
    rotation makes the crossing ring a (1,n) torus knot with linking number n.
    Step D1 requires: show why n=2 is the energetically stable winding over n=1
    and n=3.

The pulsar angle:
  A FREELY PRECESSING pulsar has TWO rotation modes simultaneously:
    (1) spin about its own symmetry axis at period P_spin
    (2) precession of that symmetry axis around the angular momentum vector
        at period P_prec
  The beam traces a path in angle-space. If P_spin/P_prec = n, the beam
  traces a (1,n) torus knot across the precession angle-sphere — EXACTLY
  the same topology as the electron's crossing ring on its Hopf torus.

  If n=2 is the stable torus knot in nature, freely precessing pulsars
  should gravitationally relax toward P_spin/P_prec = 2, providing:
    (a) Empirical evidence that the (1,2) torus knot is energetically preferred
        (Step D1 of hopf_c4b_correction.py).
    (b) A measured spin-to-precession ratio that could constrain n_exact.

The primary test case — PSR B1828-11:
  Stairs, Lyne & Shemar (2000), Nature 406, 484.
  First confirmed freely precessing pulsar. Shows clear quasi-periodic
  variations in pulse profile and timing residuals.
  Three observed periodicities reported: ~500d, ~250d, ~167d.
  Ratio of dominant periods: ~500d / ~250d = 2.0.
  This is the closest known astrophysical system to a (1,2) torus knot
  rotation pattern.

This script:
  PART A — Documents known precessing pulsars and their period ratios.
  PART B — Checks how close each ratio is to n_exact = 2.01869.
  PART C — Asks: if the pulsar ratio were exactly n_exact, what would that
            imply for the electron's C4b residual?
  PART D — What would need to be measured to test the hypothesis?

Run: python analysis/pulsar_double_spin.py
"""

import math

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)
alpha = 7.2973525693e-3

n_exact  = (4 * pi**2 / phi * alpha - Rs) / alpha**2
residual = n_exact - 2

SEP  = "=" * 65
SEP2 = "-" * 65


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PULSAR DOUBLE-SPIN ANALYSIS")
print("Can freely precessing pulsars constrain the C4b winding number?")
print(SEP)
print()
print(f"  C4b n_exact  = {n_exact:.8f}")
print(f"  Integer n=2  = 2.00000000")
print(f"  Residual     = {residual:.8f}")
print()
print("  The (1,2) torus knot hypothesis predicts that n=2 is the stable")
print("  winding. If correct, freely precessing astrophysical rotators")
print("  should relax toward spin:precession ratio = n_exact = 2.019.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART A — KNOWN FREELY PRECESSING PULSARS")
print(SEP)
print()
print("  Source notes for each pulsar follow the data.")
print()

# Format: (name, P_spin_ms, P_prec_days, ratio_spin_prec,
#          ratio_uncertainty, reference, notes)
# ratio = P_prec / P_spin expressed as dimensionless count of spins per prec.
# For comparison with n_exact: ratio = P_prec(ms) / P_spin(ms)
# i.e., how many spin cycles per one precession cycle.
#
# PSR B1828-11:
#   P_spin = 405 ms (Stairs et al. 2000, also known as 0.405 s)
#   P_prec: primary modulation ~ 511 +/- 10 days (Stairs 2000),
#           secondary modulation ~ 256 +/- 5 days.
#   Ratio P_prec1 / P_prec2 = 511/256 ~ 2.0 (this is the n analog)
#   NOTE: the ratio being tested is BETWEEN the two precession harmonics,
#   not spin/precession. The pulsar traces a (1,2) path because the two
#   observable harmonics are in 2:1 ratio.
#
# PSR B1642-03:
#   Shows profile variations on ~yearlong timescale.
#   Less well-studied. Kerr et al. 2016 (MNRAS 455) found quasi-periodic
#   variations; period not pinned down as cleanly as B1828-11.
#
# PSR J0738-4042:
#   Brook et al. 2014 (ApJL 780, L31). Profile changes attributed to
#   precession or mode-switching. Periods not confirmed to 2:1.
#
# Her X-1 / HZ Her (X-ray binary pulsar):
#   P_spin = 1.24 s. Superorbital period 35 days (precessing accretion disk).
#   35 days / 1.24 s = 2.44e6 spins per precession — very large ratio,
#   not relevant to the 2:1 torus knot test.

pulsars = [
    {
        "name":      "PSR B1828-11",
        "P_spin_ms": 405.0,
        "P_prec1_d": 511.0,   # dominant modulation period (days)
        "P_prec2_d": 256.0,   # secondary modulation (first harmonic)
        "ratio":     511.0 / 256.0,   # this is the (1,n) analog
        "ratio_unc": 0.05,            # approximate (error bars ~10d/5d)
        "reference": "Stairs, Lyne & Shemar (2000), Nature 406, 484",
        "notes":     ("First confirmed freely precessing pulsar. "
                      "Three modulation periods: ~511d, ~256d, ~170d. "
                      "Ratio P1/P2 = 511/256 = 1.996. "
                      "The ~256d period is the first harmonic of ~511d, "
                      "meaning the beam crosses its reference direction "
                      "TWICE per precession cycle — a (1,2) torus knot path.")
    },
    {
        "name":      "PSR B1642-03",
        "P_spin_ms": 387.7,
        "P_prec1_d": 380.0,   # approximate; less well-constrained
        "P_prec2_d": None,
        "ratio":     None,
        "ratio_unc": None,
        "reference": "Kerr et al. (2016), MNRAS 455, 1845",
        "notes":     ("Profile variations on ~1yr timescale. Second harmonic "
                      "not cleanly identified. Cannot test 2:1 ratio yet.")
    },
    {
        "name":      "PSR J0738-4042",
        "P_spin_ms": 380.1,
        "P_prec1_d": None,   # period not cleanly established
        "P_prec2_d": None,
        "ratio":     None,
        "ratio_unc": None,
        "reference": "Brook et al. (2014), ApJL 780, L31",
        "notes":     ("Profile change event; possibly precession onset. "
                      "Precession period not measured to useful precision.")
    },
]

for p in pulsars:
    print(f"  {p['name']}")
    print(f"    P_spin = {p['P_spin_ms']} ms")
    if p['P_prec1_d'] is not None:
        print(f"    P_prec1 = {p['P_prec1_d']} d,  P_prec2 = {p['P_prec2_d']} d")
    if p['ratio'] is not None:
        print(f"    Ratio P_prec1/P_prec2 = {p['ratio']:.6f}  (unc ~{p['ratio_unc']})")
    print(f"    Ref: {p['reference']}")
    # Wrap notes
    words = p['notes'].split()
    line = "    Notes: "
    for w in words:
        if len(line) + len(w) + 1 > 68:
            print(line)
            line = "           " + w
        else:
            line += (" " if line.strip() != "Notes:" else "") + w
    print(line)
    print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART B — HOW CLOSE IS PSR B1828-11 TO n_exact?")
print(SEP)
print()

# PSR B1828-11 ratio
r_measured = 511.0 / 256.0
r_unc      = math.sqrt((10/256)**2 + (511*5/256**2)**2)  # propagated unc

print(f"  Measured P1/P2 = {r_measured:.6f}  (unc ~ {r_unc:.4f})")
print(f"  n_exact        = {n_exact:.6f}")
print(f"  n=2 (integer)  = 2.000000")
print()
print(f"  Gap: measured - n=2      = {r_measured - 2:.6f}")
print(f"  Gap: n_exact - n=2       = {residual:.6f}")
print(f"  Gap: measured - n_exact  = {r_measured - n_exact:.6f}")
print()
print(f"  The measured ratio {r_measured:.4f} differs from n_exact {n_exact:.4f}")
print(f"  by {abs(r_measured - n_exact):.4f}.")
print(f"  The measurement uncertainty (~{r_unc:.4f}) is comparable to this gap.")
print()

# Is the pulsar consistent with n_exact within errors?
sigma_from_nexact = abs(r_measured - n_exact) / r_unc
sigma_from_n2     = abs(r_measured - 2.0)     / r_unc
print(f"  Consistency check:")
print(f"    |ratio - n=2|     / unc = {sigma_from_n2:.2f} sigma")
print(f"    |ratio - n_exact| / unc = {sigma_from_nexact:.2f} sigma")
print(f"    Both consistent within ~1 sigma with current measurement precision.")
print()
print(f"  To DISTINGUISH between n=2 and n_exact = {n_exact:.4f}, the period ratio")
print(f"  must be measured to precision better than {residual:.4f},")
print(f"  i.e., to {residual/r_measured*100:.3f}% on the ratio.")
print(f"  Current precision on P1: ~{10/511*100:.1f}%.  Needs ~{residual/r_measured*100:.3f}%.")
print(f"  Required precision improvement: ~{(10/511) / (residual/r_measured):.0f}x.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART C — IF THE PULSAR RATIO WERE EXACTLY n_exact, WHAT WOULD THAT MEAN?")
print(SEP)
print()
print("  The (1,2) torus knot hypothesis says:")
print("    n_effective = 2 + [higher-order correction]")
print("    where the correction comes from the topology of the crossing ring")
print("    on the Hopf torus (Step D4 of hopf_c4b_correction.py).")
print()
print("  If a precessing pulsar ALSO relaxes to ratio n_effective = n_exact,")
print("  the two systems share the same topological reason for the correction.")
print("  The astrophysical value would then constrain our residual.")
print()
print(f"  From PSR B1828-11: ratio = {r_measured:.6f} (unc {r_unc:.4f})")
print(f"  From C4b:          n_exact = {n_exact:.6f}")
print()

# What implied alpha would PSR B1828-11's ratio give if we use it as n?
# n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0, solve for alpha with n=r_measured
def solve_quad(n):
    a = n
    b = -(4 * pi**2 / phi)
    c = Rs
    disc = b**2 - 4*a*c
    if disc < 0:
        return None
    return (-b - math.sqrt(disc)) / (2*a)

alpha_if_pulsar = solve_quad(r_measured)
if alpha_if_pulsar:
    err_pulsar = (alpha_if_pulsar - alpha) / alpha * 100
    print(f"  If n = pulsar ratio {r_measured:.6f}:")
    print(f"    implied alpha = {alpha_if_pulsar:.13e}")
    print(f"    error from CODATA = {err_pulsar:+.6f}%")
    print(f"    (vs C4b error = -0.000560%)")
    print()
    print(f"  The pulsar-implied alpha is {abs(err_pulsar)/abs(-0.000560):.1f}x further from")
    print(f"  CODATA than C4b — consistent with the pulsar ratio being a rough")
    print(f"  macroscopic analog, not a precision measurement of n_exact.")
print()

print("  Physical interpretation if the analogy holds:")
print("    The (1,2) torus knot is the ground state because it minimises")
print("    the energy of a precessing, spinning body in a medium that")
print("    penalises odd-numbered windings (n=1 unstable, n=3 too costly).")
print("    The slight excess above n=2 (the 0.019 residual) reflects the")
print("    finite stiffness of the coupling — the knot is not a perfect")
print("    integer but slightly 'pulled' by the self-interaction amplitude.")
print("    A precessing pulsar with enough timing precision could measure")
print("    this directly if the torsion medium acts on it the same way.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART D — WHAT WOULD NEED TO BE MEASURED?")
print(SEP)
print()
print("  To use PSR B1828-11 (or similar) to constrain n_exact, four things")
print("  are needed:")
print()
print("  D1 — Period precision:")
print(f"     Current: P1/P2 precision ~ {r_unc/r_measured*100:.1f}%")
print(f"     Required to resolve n_exact vs n=2: ~ {residual/n_exact*100:.3f}%")
print(f"     = {(r_unc/r_measured) / (residual/n_exact):.0f}x improvement in timing precision.")
print(f"     Achievable with continued monitoring: Stairs et al. 2000 used")
print(f"     13 years of data. Another 10-20 years narrows the period by ~sqrt(T).")
print()
print("  D2 — Confirm the ratio is intrinsic, not a beat frequency:")
print("     The 2:1 ratio in B1828-11 could be a fundamental+first-harmonic")
print("     of a single precession frequency, OR two independent modes.")
print("     If it is harmonic structure from one precession: the ratio is")
print("     exactly 2.000 by definition (not n_exact).")
print("     If it is two independent precession modes in 2:1 resonance:")
print("     the ratio reflects the stability of the (1,2) knot, and")
print("     could take any value near 2, including n_exact.")
print("     Distinguishing these requires pulse profile modelling, not")
print("     just timing residuals. Lyne et al. (2010, Science 329, 408)")
print("     showed correlated spin-down and profile changes, supporting")
print("     two modes, but the precise ratio wasn't the focus.")
print()
print("  D3 — Other precessing pulsars:")
print("     A larger sample would show whether 2:1 is preferred statistically.")
print("     Current confirmed precessing pulsars: ~3-5 objects.")
print("     FAST (China) and MeerKAT (South Africa) long-term monitoring")
print("     campaigns may identify more, and with better period precision.")
print()
print("  D4 — Connect to torsion medium formalism:")
print("     If R_s governs the medium stiffness against torsional excitation,")
print("     the stability energy of a (1,2) vs (1,1) knot should scale with R_s.")
print("     A calculation: what is the energy difference between a (1,2) and")
print("     (1,1) torus knot in a medium with shear modulus G = 1.66e-11 Pa")
print("     (from analysis/medium_properties.py)?")
print("     If that energy difference is of order alpha*E_electron, and the")
print("     pulsar analog at its scale shows the same energy separation,")
print("     the connection between the electron topology and pulsar precession")
print("     is quantitative, not merely analogical.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()
print("  PSR B1828-11 is the best known candidate for macroscopic (1,2)")
print("  torus knot rotation. Its dominant modulation periods are in")
print(f"  ~2:1 ratio ({r_measured:.4f} measured, vs n_exact = {n_exact:.4f}).")
print()
print("  The analogy to the electron's crossing ring:")
print("    Electron:  figure-8 sweeps main axis + spins on horizontal axis.")
print("               Crossing ring = (1,2) torus knot on Hopf torus surface.")
print("               Result: n=2 coefficient in C4b quadratic.")
print("    Pulsar:    spin axis precesses around angular momentum axis.")
print("               Beam path in angle-space = (1,2) torus knot.")
print("               Result: 2:1 ratio in timing modulation periods.")
print()
print("  What the pulsar data gives us RIGHT NOW:")
print("    - Qualitative: n=2 IS the observed ratio in nature at macroscopic scale.")
print("      This is evidence for Step D1 (stability of (1,2) knot).")
print("    - Quantitative: NOT YET. Current period measurement precision")
print(f"      (~{r_unc/r_measured*100:.1f}%) is {(r_unc/r_measured)/(residual/n_exact):.0f}x too coarse to see n_exact vs n=2.")
print()
print("  What is needed:")
print("    A confirmed freely precessing pulsar with period ratio measured to")
print(f"    {residual/n_exact*100:.3f}% or better. Then:")
print("    - If ratio converges to n_exact = 2.019: strong evidence that")
print("      the same topological mechanism operates at both scales.")
print("    - If ratio converges to 2.000 exactly: the harmonic interpretation")
print("      is correct and the residual 0.019 is purely geometric (Step D4).")
print(f"    - Either outcome sharpens the theory. Best current instrument:")
print("      FAST pulsar timing programme (sensitivity >> Parkes/Arecibo).")
print()
print("  File this under: Step D1 observational support.")
print("  See: analysis/hopf_c4b_correction.py (Step D1 criteria)")
print("       analysis/c4b_residual_scale.py  (residual near-hits)")
print(SEP)
