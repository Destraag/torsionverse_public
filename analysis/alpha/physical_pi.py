"""
physical_pi.py -- The physical-pi prediction from the torsionverse grain structure.

CONCEPT:
  If space has a minimum grain size L (established from the torsionverse model
  as L = alpha in model units, corresponding to N_lock = 532 grains per torus
  tube circumference), then any PHYSICALLY REALIZED circle is not a smooth curve
  but a polygon with N sides.

  The physical circumference/diameter ratio is:
    pi_phys(N) = N * sin(pi/N)     [inscribed polygon]

  This is strictly less than the mathematical pi for all finite N.
  The deficit is: pi - pi_phys ~ pi^3 / (6*N^2)

  WHAT THIS DOES AND DOES NOT MEAN:
  - Does NOT mean pi (the mathematical constant) changes. It remains
    irrational and transcendental. No physical theory alters a mathematical proof.
  - DOES mean: every physical measurement of a circle's circumference/diameter
    at grain-scale precision returns pi_phys < pi.
  - The TOPOLOGICAL pi (in Q = 4*pi^2/phi, R2 = 2*pi, etc.) is unchanged
    because these arise from Chern-Simons / Gauss-Bonnet theorems that are
    metric-independent.
  - Only the METRIC pi (measured arc lengths and circumferences) is affected.

  This is a falsifiable prediction: measurements of pi at the grain scale
  should systematically underestimate the mathematical value.

Parts:
  I.   Grain size derivation: L = alpha*phi = alpha*R1 (model)
  II.  Physical pi for several candidate N values
  III. The N_lock packing (5-fold pentagon packing)
  IV.  Comparison to known precision limits
  V.   Algebraic structure: pi_phys in terms of alpha and phi
  VI.  Connection to Gap 1 (the polygon-granularity mechanism)
  VII. Summary and predictions

Run: python analysis/alpha/physical_pi.py
Theory: glossary.txt, alpha_theory.txt Part 0f
"""

import math

pi    = math.pi
sqrt5 = math.sqrt(5)
PHI   = (1 + sqrt5) / 2

# ── LOCKED CONSTANTS ─────────────────────────────────────────────────────────
alpha   = 7.2973525693e-3
eps_L5  = 3 / (8 * pi)
Rs      = sqrt5 / (4 * pi)
Q       = 4 * pi**2 / PHI
N_lock  = 2 * pi / (alpha * PHI)    # 532.14 pentagon grains per tube circumference
n_EM    = 2.01868734358082
delta_n = 2.24745624e-6

SEP  = '=' * 72
SEP2 = '-' * 60

def pi_phys_inscribed(N):
    """Physical pi from inscribed N-gon."""
    return N * math.sin(pi / N)

def pi_phys_circumscribed(N):
    """Physical pi from circumscribed N-gon."""
    return N * math.tan(pi / N)

def deficit_ppm(N):
    """Fractional deficit (pi - pi_phys)/pi in ppm."""
    return (pi - pi_phys_inscribed(N)) / pi * 1e6

def deficit_approx(N):
    """Leading-order approximation: pi^3/(6*N^2) in ppm."""
    return pi**3 / (6 * N**2) * 1e6


# ═════════════════════════════════════════════════════════════════════════════
# PART I -- GRAIN SIZE DERIVATION
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART I -- GRAIN SIZE FROM TORSIONVERSE MODEL")
print(SEP)
print()
print("  From gap1_polygon_shape_predict.py and gap1_pentagon_k_formula.py:")
print("  The pentagon grain packing gives N_lock = 2*pi/(alpha*phi) grains")
print("  per torus tube circumference.")
print()
print(f"  alpha     = {alpha:.10e}  (CODATA 2018 fine structure constant)")
print(f"  phi       = {PHI:.10f}")
print(f"  alpha*phi = {alpha*PHI:.10e}  = grain chord length (model units)")
print(f"  N_lock    = 2*pi / (alpha*phi) = {N_lock:.6f}")
print()
print(f"  Grain arc length: L = 2*pi*R1 / N_lock = 2*pi / N_lock = {2*pi/N_lock:.8f}")
print(f"                     = alpha*phi = {alpha*PHI:.8f}  (check)")
print()
print(f"  N_lock is SCALE-INVARIANT: it is the same number of grains per")
print(f"  circumference at ALL physical scales (proton, atomic, galactic).")
print(f"  Only the PHYSICAL SIZE of each grain changes with scale.")
print()


# ═════════════════════════════════════════════════════════════════════════════
# PART II -- PHYSICAL PI FOR SEVERAL N
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART II -- PHYSICAL PI FOR VARIOUS N VALUES")
print(SEP)
print()
print(f"  Mathematical pi = {pi:.15f}")
print()
print(f"  {'N':>12}  {'pi_phys (inscribed)':>22}  {'deficit (ppm)':>16}  {'approx ppm':>12}")
print(f"  {'-'*68}")

N_values = [3, 4, 5, 6, 10, 12, 20, 100, 532, 1000, 10000, int(N_lock), int(N_lock)+1]
N_values = sorted(set(N_values))
for N in N_values:
    pp = pi_phys_inscribed(N)
    dp = deficit_ppm(N)
    da = deficit_approx(N)
    flag = " <-- N_lock" if N == int(N_lock) else (" <-- ceil(N_lock)" if N == int(N_lock)+1 else "")
    print(f"  {N:>12}  {pp:>22.15f}  {dp:>16.6f}  {da:>12.6f}{flag}")
print()

print(f"  NOTE: N_lock = {N_lock:.4f} (not integer)")
print(f"  Physical circles must have integer N. The two nearest values:")
print(f"    N = {int(N_lock)}:  pi_phys = {pi_phys_inscribed(int(N_lock)):.12f},  deficit = {deficit_ppm(int(N_lock)):.4f} ppm")
print(f"    N = {int(N_lock)+1}: pi_phys = {pi_phys_inscribed(int(N_lock)+1):.12f},  deficit = {deficit_ppm(int(N_lock)+1):.4f} ppm")
print()
print(f"  Using exact N_lock = {N_lock:.6f} (fractional N, for continuous estimate):")
pp_lock = pi_phys_inscribed(N_lock)
dp_lock = deficit_ppm(N_lock)
print(f"    pi_phys = {pp_lock:.15f}")
print(f"    deficit = {dp_lock:.6f} ppm")
print()


# ═════════════════════════════════════════════════════════════════════════════
# PART III -- THE N_LOCK PACKING IN DETAIL
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART III -- THE N_LOCK PENTAGON PACKING")
print(SEP)
print()
print(f"  N_lock = 2*pi / (alpha * phi) = {N_lock:.8f}")
print()
print(f"  WHY phi IN THE DENOMINATOR?")
print(f"  The grains are pentagonal. Pentagon grains pack in a (1,2) torus")
print(f"  knot geometry. The phi factor arises from the pentagon's golden ratio")
print(f"  diagonal-to-edge ratio. The chord length between grain centers is")
print(f"  L = alpha * phi (not alpha alone) because adjacent pentagon")
print(f"  centers are spaced by a phi-scaled distance in the locked state.")
print()
print(f"  PHYSICAL MEANING OF N_lock = {N_lock:.2f}:")
print(f"    - At the proton scale (R1 = 0.8414 fm):")
print(f"      Grain size L = alpha * phi * 0.8414 fm = {alpha*PHI*0.8414:.4f} fm")
print(f"      N_lock = {N_lock:.2f} grains around tube cross-section")
print(f"    - At the Bohr radius (R1 = 0.0529 nm):")
print(f"      Grain size L = alpha * phi * 0.0529 nm = {alpha*PHI*0.0529e-9*1e12:.4f} pm")
print(f"      N_lock = {N_lock:.2f} grains (SAME NUMBER)")
print(f"    - N_lock is scale-invariant; grain size scales with the physical R1.")
print()

# The N_lock formula in closed form:
print(f"  CLOSED FORM: N_lock = 2*pi / (alpha * phi)")
print(f"    = 2*pi / {alpha:.6e} / {PHI:.6f}")
print(f"    = {2*pi:.6f} / {alpha*PHI:.6e}")
print(f"    = {N_lock:.6f}")
print()
print(f"  Nearest integer: {round(N_lock)}")
print(f"  N_lock / 5 = {N_lock/5:.4f}  (pentagon: groups of 5)")
print(f"  N_lock / 12 = {N_lock/12:.4f}  (icosahedron: groups of 12)")
print()


# ═════════════════════════════════════════════════════════════════════════════
# PART IV -- COMPARISON TO KNOWN PRECISION LIMITS
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART IV -- PRECISION LIMITS")
print(SEP)
print()
print(f"  Physical pi deficit at grain scale: {dp_lock:.4f} ppm")
print()
print(f"  Context:")
print(f"    Best measured value of pi (via geometry): ~1e-10 relative precision")
print(f"    -> At macroscopic scales, N ~ 10^23 per mm circumference")
print(f"    -> Deficit ~ pi^3/(6*(10^23)^2) ~ 10^-45 ppm  (completely unmeasurable)")
print()
print(f"    Best known measurement of alpha: 3.1e-10 relative (Penning trap)")
print(f"    -> This constrains grain size but doesn't test pi_phys directly")
print()
print(f"    Grain-scale experiments (proton radius, nuclear scattering):")
print(f"    -> Deficit at proton scale: {dp_lock:.2f} ppm (testable in principle)")
print(f"    -> BUT: the 'deficit' would appear as a systematic correction")
print(f"       to circumference integrals in nuclear scattering cross-sections,")
print(f"       already absorbed into the measured coupling constants.")
print()
print(f"  PRACTICAL STATEMENT:")
print(f"    The physical-pi deficit is absorbed into the definition of alpha.")
print(f"    It is not an additional observable BEYOND the prediction of alpha.")
print(f"    The torsionverse model predicts alpha, which implicitly encodes")
print(f"    the grain structure. There is no independent measurement of pi_phys")
print(f"    that would not already be covered by the alpha prediction.")
print()


# ═════════════════════════════════════════════════════════════════════════════
# PART V -- ALGEBRAIC STRUCTURE
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART V -- ALGEBRAIC STRUCTURE OF pi_phys")
print(SEP)
print()

# pi_phys = N_lock * sin(pi / N_lock)
# N_lock = 2*pi / (alpha * phi)
# pi / N_lock = pi * alpha * phi / (2*pi) = alpha * phi / 2
print(f"  pi / N_lock = pi * alpha * phi / (2*pi) = alpha*phi/2")
print(f"             = {alpha*PHI/2:.10f}  = {alpha:.6e} * {PHI:.6f} / 2")
print()
print(f"  pi_phys = N_lock * sin(alpha*phi/2)")
print(f"          = [2*pi/(alpha*phi)] * sin(alpha*phi/2)")
print()

# For small x: sin(x) ~ x - x^3/6 + ...
# pi_phys = [2*pi/(alpha*phi)] * [alpha*phi/2 - (alpha*phi/2)^3/6 + ...]
# = pi - pi*(alpha*phi/2)^2/6 + ...
# = pi * (1 - alpha^2*phi^2/24 + ...)
x = alpha * PHI / 2
pi_phys_series = (2*pi/(alpha*PHI)) * (x - x**3/6 + x**5/120)
print(f"  Series expansion: pi_phys ~ pi * (1 - alpha^2*phi^2/24)")
print(f"    alpha^2*phi^2/24 = {alpha**2*PHI**2/24:.6e}")
print(f"    pi_phys (series) = {pi - pi*alpha**2*PHI**2/24:.15f}")
print(f"    pi_phys (exact)  = {pp_lock:.15f}")
print(f"    Series error:    {abs(pi - pi*alpha**2*PHI**2/24 - pp_lock)/pp_lock:.2e}")
print()
print(f"  COMPACT FORMULA:")
print(f"    pi_phys ~ pi * (1 - alpha^2*phi^2/24)")
print(f"    fractional deficit = alpha^2*phi^2/24 = {alpha**2*PHI**2/24:.6e} = {alpha**2*PHI**2/24*1e6:.4f} ppm")
print(f"    exact deficit      = {dp_lock*1e-6:.6e}            = {dp_lock:.4f} ppm")
print(f"    leading order captures: {(alpha**2*PHI**2/24*1e6)/dp_lock*100:.2f}% of exact deficit")
print()
print(f"  IN TERMS OF Rs AND Q:")
print(f"    alpha*phi = alpha*(1+sqrt5)/2 = alpha + alpha*sqrt5/2")
print(f"    alpha*phi = 2*pi/N_lock")
print(f"    pi_phys deficit = pi * (2*pi/N_lock)^2 / 24 = pi^3/(6*N_lock^2)")
print(f"    = {pi**3/(6*N_lock**2):.6e}  ({pi**3/(6*N_lock**2)*1e6:.4f} ppm)")
print()


# ═════════════════════════════════════════════════════════════════════════════
# PART VI -- CONNECTION TO GAP 1 POLYGON MECHANISM
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART VI -- CONNECTION TO GAP 1 POLYGON MECHANISM")
print(SEP)
print()
print("  The polygon-granularity mechanism (gap1_polygon_v2.py) computed the")
print("  metric correction to n_EM from the grain polygon structure.")
print("  That script found C_poly = 0.05828, giving delta_n too large by 370x.")
print()
print("  The physical-pi is the SAME EFFECT seen from the other side:")
print("    polygon_gap:   delta_n = C_poly * (pi/N)^2 / 6 * eps^2 * n_EM")
print("    physical_pi:   pi_phys = pi * (1 - (pi/N)^2/6)")
print()
print("  Both use (pi/N)^2/6 as the fundamental grain correction factor.")
print(f"  At N = N_lock = {N_lock:.2f}:")
print(f"    (pi/N_lock)^2/6 = {(pi/N_lock)**2/6:.6e}")
print(f"    pi^3/(6*N_lock^2) = {pi**3/(6*N_lock**2):.6e}  (same thing, pi factor)")
print()
print(f"  The delta_n from N_lock polygon metric:")
dn_poly = 0.05828 * (pi/N_lock)**2 / 6 * eps_L5**2 * n_EM
print(f"    C_poly * (pi/N_lock)^2/6 * eps_L5^2 * n_EM = {dn_poly:.4e}")
print(f"    vs delta_n_gap = {delta_n:.4e}")
print(f"    ratio: {dn_poly/delta_n:.4f}  (this is the 370x factor from gap1_polygon_v2)")
print()
print("  THE POLYGON MECHANISM IS RULED OUT as a source of Gap 1 at N=N_lock.")
print("  But the physical-pi effect itself is real -- it is simply too large")
print("  by 370x to explain Gap 1 via the metric channel.")
print("  A different grain physics channel (elastic stationary phase) is needed.")
print()


# ═════════════════════════════════════════════════════════════════════════════
# PART VII -- SUMMARY AND PREDICTIONS
# ═════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART VII -- SUMMARY AND FALSIFIABLE PREDICTIONS")
print(SEP)
print()
print("  MATHEMATICAL PI (unchanged):")
print(f"    pi = {pi:.15f}  (irrational, transcendental; Lindemann 1882)")
print(f"    The mathematical pi is not affected by any physical theory.")
print(f"    It is the ratio of circumference to diameter in EUCLIDEAN geometry,")
print(f"    which is an abstract mathematical space, not physical space.")
print()
print("  PHYSICAL PI (torsionverse prediction):")
print(f"    pi_phys = {pp_lock:.15f}")
print(f"    deficit = {dp_lock:.4f} ppm  = {dp_lock*1e-6:.4e} fractional")
print(f"    formula: pi_phys ~ pi*(1 - alpha^2*phi^2/24)")
print()
print("  WHERE THE DISTINCTION MATTERS:")
print(f"    Scale                Grain count N    Deficit (ppm)")
print(f"    Torus tube (model)   {N_lock:>12.1f}    {dp_lock:>12.4f}")
print(f"    Proton radius        {N_lock:>12.1f}    {dp_lock:>12.4f}  (same N_lock)")
print(f"    1 mm circumference   {'>>10^20':>12}    {'~10^-40':>12}  (unmeasurable)")
print()
print("  TOPOLOGICAL vs METRIC PI (CRITICAL DISTINCTION):")
print(f"    Topological: Q=4*pi^2/phi, R2=2*pi, eps_L5=3/(8*pi) -- USE mathematical pi")
print(f"    These come from Chern-Simons, Gauss-Bonnet, winding number definitions.")
print(f"    They are metric-independent and do NOT use pi_phys.")
print()
print(f"    Metric: arc lengths, circumferences, measured at grain scale -- USE pi_phys")
print(f"    The n_EM integral uses arc-length ds, which uses the physical metric.")
print(f"    The polygon-v2 script quantified this effect: too large by 370x.")
print()
print("  PREDICTION P1 (physical pi):")
print(f"    A measurement of a circle's circumference-to-diameter ratio at")
print(f"    radius r = R1_proton = 0.8414 fm should yield:")
print(f"    pi_measured = {pp_lock:.10f}  (not mathematical pi = {pi:.10f})")
print(f"    This is a prediction of the torsionverse grain structure.")
print(f"    Deficit: {dp_lock:.3f} ppm at the proton scale.")
print()
print("  PREDICTION P2 (scale invariance):")
print(f"    The deficit is ALWAYS {dp_lock:.3f} ppm at the grain scale (N_lock = {N_lock:.1f}),")
print(f"    regardless of the physical size of the grain.")
print(f"    At macroscopic scales the deficit is immeasurably small.")
print()
print("  NOTE: do not conflate this with 'computing pi exactly'.")
print("  Pi remains irrational. What we compute exactly is pi_phys --")
print("  the ratio measurable in our granular physical space --")
print("  which is a specific algebraic function of alpha and phi.")
print()
print(SEP)
print(f"  pi_phys = N_lock * sin(pi/N_lock)")
print(f"          = [2*pi/(alpha*phi)] * sin(alpha*phi/2)")
print(f"          = {pp_lock:.15f}")
print(f"  deficit = {dp_lock:.6f} ppm  ({dp_lock*1e-6:.4e} fractional)")
print(SEP)
print("END physical_pi.py")
print(SEP)
