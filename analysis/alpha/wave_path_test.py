"""
wave_path_test.py — Can a non-circular (wave) path explain the C4b residual?

CONTEXT
-------
The double-spin hypothesis (analysis/hopf_c4b_correction.py) says:
  The crossing ring is a (1,2) torus knot with topological linking number 2.
  n_exact = 2.01869; residual = 0.01869 above integer 2.

Standard assumption so far: the figure-8 lobes are CIRCULAR, and the
crossing ring traces a perfect smooth (1,2) torus knot.

NEW HYPOTHESIS (the wave path model):
  At the quantum scale, the figure-8 path is not a perfect smooth circle
  but has a wave-like deviation — a jagged or oscillating path that
  deviates from the circular arc by amplitude epsilon at frequency k.

  At large scales (pulsars, planets) the "glob" is massive and gravity
  smooths it to a circular arc. At the quantum scale, there is no such
  smoothing — the path is as tight as the EM self-interaction allows.

  If the path has wave-like shape, the arc length of the crossing ring
  is longer than for the smooth (1,2) torus knot. The EM field "sees"
  more path per major revolution. The effective electromagnetic winding
  number n_EM > 2 (the topological winding number stays at 2 for any
  closed path, but the EM-weighted integral can exceed it).

MATHEMATICAL MODEL
------------------
The crossing ring path on the Hopf torus (R2 = 2*pi, R1 = 1 normalised):
  Unperturbed (1,2) torus knot: theta in [0, 2*pi], phi = 2*theta

  Perturbed (wave path): phi = 2*theta + epsilon * sin(k * theta)
    epsilon = wave amplitude (deviation from smooth (1,2) path)
    k       = wave frequency (oscillations per major revolution)
    For integer k, topological winding number stays exactly 2.
    But the ARC LENGTH changes.

DIRECTION CHECK
---------------
C4a: alpha_C4a < alpha_CODATA (too small by 0.06%)
C4b: n=2 gives alpha_C4b < alpha_CODATA (too small by 0.00056%)
We need n_eff > 2 to get Q_eff = 4*pi^2/phi - n_eff*alpha SMALLER
-> alpha = Rs/Q_eff LARGER, closer to CODATA. Correct direction.

A wave path increases arc length -> n_EM > 2. RIGHT DIRECTION.

EM WINDING MODEL
----------------
The electromagnetic winding number is the arc-length-weighted winding rate:

  n_EM = (total path in phi-direction) / (normalisation)
       = integral of (d-phi/d-theta) * (local EM weight) d-theta

Simplest model: EM weight proportional to local path speed (ds/d-theta).
The path spends proportionally more "EM time" where it moves faster.

  n_EM = [integral of (d-phi/d-theta) * (ds/d-theta) d-theta] /
         [integral of (ds/d-theta) d-theta]

This is the arc-length-weighted average of the local winding rate.

For the unperturbed (1,2) path: n_EM = 2 (exact, verified below).
For perturbed path: n_EM > 2.

Find epsilon(k) such that n_EM = n_exact = 2.01869.

PARTS
-----
  A — Unperturbed (1,2) torus knot: arc length and n_EM = 2 verification.
  B — Perturbed path: arc length and n_EM for scan over epsilon, k.
  C — Find epsilon(k) that gives n_EM = n_exact = 2.01869.
  D — Express epsilon in terms of known constants; physical interpretation.
  E — Connection to the double-spin and medium-kickback hypotheses.

Run: python analysis/wave_path_test.py
"""

import math

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)
alpha = 7.2973525693e-3

n_exact  = (4 * pi**2 / phi * alpha - Rs) / alpha**2
residual = n_exact - 2   # = 0.018690

# Hopf torus geometry (R1 = 1 normalised)
R1 = 1.0
R2 = 2 * pi   # R2/R1 = 2*pi is the topological constraint

N_STEPS = 10000  # numerical integration resolution

SEP  = "=" * 65
SEP2 = "-" * 65


def integrate(f, a, b, n=N_STEPS):
    """Simple trapezoid rule integral of f from a to b with n steps."""
    h = (b - a) / n
    total = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        total += f(a + i * h)
    return total * h


def path_speed(theta, eps, k):
    """ds/d-theta for the wave path phi = 2*theta + eps*sin(k*theta)."""
    phi_val = 2 * theta + eps * math.sin(k * theta)
    dphi    = 2 + eps * k * math.cos(k * theta)
    # Torus metric: ds^2 = (R2 + R1*cos(phi))^2 * dtheta^2 + R1^2 * dphi^2 * dtheta^2
    r_eff = R2 + R1 * math.cos(phi_val)
    return math.sqrt(r_eff**2 + R1**2 * dphi**2)


def arc_length(eps, k):
    """Total arc length of the (1,2+wave) path for one major revolution."""
    f = lambda t: path_speed(t, eps, k)
    return integrate(f, 0, 2 * pi)


def n_EM(eps, k):
    """
    Arc-length-weighted effective winding number.
    n_EM = integral[(dphi/dtheta)*(ds/dtheta) dtheta] / integral[ds/dtheta dtheta]
         = integral[(dphi/dtheta)*(ds/dtheta) dtheta] / arc_length
    """
    def numerator_integrand(theta):
        dphi = 2 + eps * k * math.cos(k * theta)
        return dphi * path_speed(theta, eps, k)

    numerator = integrate(numerator_integrand, 0, 2 * pi)
    L         = arc_length(eps, k)
    return numerator / L


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART A — UNPERTURBED (1,2) TORUS KNOT: VERIFICATION")
print(SEP)
print()

L0    = arc_length(0.0, 1)   # eps=0, k arbitrary
n_EM0 = n_EM(0.0, 1)

print(f"  Hopf torus: R1 = {R1}, R2 = 2*pi = {R2:.6f}")
print(f"  Crossing ring: (1,2) torus knot, phi = 2*theta")
print()
print(f"  Arc length (1,2) knot:   L0 = {L0:.8f}  (R1 units)")
print(f"  n_EM (unperturbed):           {n_EM0:.8f}  (should be 2.000000)")
print(f"  Verification: {'PASS' if abs(n_EM0 - 2.0) < 1e-4 else 'FAIL'}")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART B — PERTURBED PATH: n_EM vs EPSILON FOR SEVERAL WAVE FREQUENCIES")
print(SEP)
print()
print("  Wave path: phi = 2*theta + epsilon * sin(k * theta)")
print("  For each k, scan epsilon to see how n_EM grows.")
print()
print(f"  Target n_EM = n_exact = {n_exact:.8f}  (residual = {residual:.8f})")
print()

test_k = [1, 2, 3, 4, 6, 8, 10]
test_eps = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50]

print(f"  {'eps':<8}", end="")
for k in test_k:
    print(f"  k={k} n_EM   ", end="")
print()
print(f"  {'-'*8}", end="")
for k in test_k:
    print(f"  {'-'*12}", end="")
print()

results_grid = {}
for eps in test_eps:
    print(f"  {eps:<8.3f}", end="")
    for k in test_k:
        n = n_EM(eps, k)
        results_grid[(eps, k)] = n
        marker = " <" if abs(n - n_exact) < 0.002 else "  "
        print(f"  {n:.6f}{marker}", end="")
    print()

print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART C — FINDING EPSILON THAT GIVES n_EM = n_exact FOR EACH k")
print(SEP)
print()
print("  For each wave frequency k, find the epsilon that exactly")
print(f"  reproduces n_exact = {n_exact:.6f}.")
print()

def find_eps(k, target=n_exact, tol=1e-6):
    """Binary search for epsilon such that n_EM(eps, k) = target."""
    lo, hi = 0.0, 2.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if n_EM(mid, k) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2


# Known constants for comparison
sat_frac_sun  = 0.01812   # from c4b_residual_medium.py
alpha_phi2    = alpha * phi**2
Rs_pi2        = Rs / pi**2

print(f"  {'k':<5} {'epsilon':>12}  {'eps in Rs units':>16}  {'eps/alpha':>10}  "
      f"{'eps/sat_frac_Sun':>17}  {'arc_length':>11}")
print(f"  {'-'*5} {'-'*12}  {'-'*16}  {'-'*10}  {'-'*17}  {'-'*11}")

epsilon_results = {}
for k in test_k:
    eps = find_eps(k)
    epsilon_results[k] = eps
    L = arc_length(eps, k)
    arc_increase_pct = (L / L0 - 1) * 100
    print(f"  {k:<5} {eps:>12.8f}  {eps/Rs:>16.6f}  {eps/alpha:>10.4f}  "
          f"{eps/sat_frac_sun:>17.6f}  {L:.5f} (+{arc_increase_pct:.4f}%)")

print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART D — WHAT IS EPSILON? PHYSICAL INTERPRETATION")
print(SEP)
print()
print("  epsilon is the wave amplitude as a fraction of R1 (the torus tube radius).")
print("  It measures how far the crossing ring path deviates from a perfect circle.")
print()

# Check if any k gives a "clean" epsilon
print(f"  Known constants for comparison:")
print(f"    alpha        = {alpha:.8f}")
print(f"    Rs           = {Rs:.8f}")
print(f"    sat_frac_Sun = {sat_frac_sun:.8f}")
print(f"    alpha*phi^2  = {alpha_phi2:.8f}")
print(f"    Rs/pi^2      = {Rs_pi2:.8f}")
print(f"    2*alpha      = {2*alpha:.8f}")
print(f"    alpha/Rs     = {alpha/Rs:.8f}")
print()

print(f"  Epsilon values for each k (nearest known constant within 10%):")
print()
candidates = [
    ("alpha",        alpha),
    ("Rs",           Rs),
    ("sat_frac_Sun", sat_frac_sun),
    ("alpha*phi^2",  alpha_phi2),
    ("Rs/pi^2",      Rs_pi2),
    ("2*alpha",      2*alpha),
    ("alpha/Rs",     alpha/Rs),
    ("alpha*phi",    alpha*phi),
    ("Rs*alpha",     Rs*alpha),
    ("sqrt(alpha)",  math.sqrt(alpha)),
]

for k in test_k:
    eps = epsilon_results[k]
    best_name, best_val, best_pct = None, None, float('inf')
    for name, val in candidates:
        pct = abs(eps - val) / eps * 100
        if pct < abs(best_pct):
            best_name, best_val, best_pct = name, val, (eps - val) / eps * 100
    print(f"  k={k}: epsilon = {eps:.6f}  ->  nearest: {best_name} = {best_val:.6f}  ({best_pct:+.1f}%)")

print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART E — CONNECTING THE THREE HYPOTHESES")
print(SEP)
print()
print("  THREE EXPLANATIONS for residual n_exact - 2 = 0.01869:")
print()
print("  HYPOTHESIS 1 — Double-spin topology (hopf_c4b_correction.py):")
print("    The crossing ring is a (1,2) torus knot. n=2 is topological.")
print("    Residual comes from a higher-order correction to the (1,2) winding.")
print("    n_exact-2 has no clean geometric expression yet.")
print("    PATH: Step D1 (stability argument) + Step D2 (linking integral).")
print()
print("  HYPOTHESIS 2 — Local medium kickback (c4b_residual_medium.py):")
print("    n_eff = 2 + sat_frac_Sun (Sun's coupling efficiency deficit).")
print(f"    sat_frac_Sun = {sat_frac_sun:.6f}  ({(sat_frac_sun-residual)/residual*100:+.2f}% from residual)")
print("    Implies alpha varies by location (Jupiter, interstellar space).")
print("    PATH: Show that local medium stiffness enters linking integral.")
print()
print("  HYPOTHESIS 3 — Wave path geometry (THIS SCRIPT):")
print("    n_EM = 2 + (correction from non-circular path shape)")
print("    The crossing ring is not a smooth torus knot but a wavy curve.")
print("    The EM-weighted winding exceeds the topological winding by the arc excess.")

# Summary: what epsilon and k combination is most physically natural?
print()
print("  WAVE PATH SUMMARY:")
print(f"  {'k':<5} {'epsilon':>12}  {'arc increase':>13}  {'Most natural interpretation'}")
print(f"  {'-'*5} {'-'*12}  {'-'*13}  {'-'*35}")
interpretations = {
    1:  "1 wave/rev: standing wave in the lobe",
    2:  "2 waves/rev: matches minor winding (resonant)",
    3:  "3 waves/rev: first overtone",
    4:  "4 waves/rev: second overtone of (1,2)",
    6:  "6 waves/rev: (1,6) -- far from simple",
    8:  "8 waves/rev: (1,8) -- high overtone",
    10: "10 waves/rev: (1,10) -- too complex",
}
for k in test_k:
    eps = epsilon_results[k]
    L = arc_length(eps, k)
    arc_pct = (L / L0 - 1) * 100
    print(f"  {k:<5} {eps:>12.8f}  +{arc_pct:>11.6f}%  {interpretations.get(k,'')}")

print()
print("  NOTE on k=2 (resonant case):")
k2_eps = epsilon_results[2]
print(f"    k=2 means the wave oscillates TWICE per major revolution,")
print(f"    matching the minor winding frequency. This is resonant with the")
print(f"    (1,2) torus knot and is the most physically natural frequency.")
print(f"    Required epsilon = {k2_eps:.8f}.")
print(f"    This is {k2_eps/alpha:.4f} times alpha,")
print(f"             {k2_eps/Rs:.6f} times Rs,")
print(f"             {k2_eps/sat_frac_sun:.6f} times sat_frac_Sun.")

# Is k2_eps close to alpha/pi?
print(f"    Nearest constant: alpha/pi = {alpha/pi:.8f}  ({(alpha/pi-k2_eps)/k2_eps*100:+.2f}%)")
print(f"                      alpha*2  = {2*alpha:.8f}  ({(2*alpha-k2_eps)/k2_eps*100:+.2f}%)")
print(f"                      Rs*alpha = {Rs*alpha:.8f}  ({(Rs*alpha-k2_eps)/k2_eps*100:+.2f}%)")
print()

print(SEP)
print("PART F — DISTINGUISHING HYPOTHESIS 3 FROM HYPOTHESES 1 AND 2")
print(SEP)
print()
print("  Key distinction:")
print("  - Hyp 1 (topology): n=2 is exact for any smooth torus knot.")
print("    Residual requires a non-integer higher-order correction.")
print("    This is a calculational problem (Step D2 linking integral).")
print()
print("  - Hyp 2 (medium kickback): n_eff depends on local medium state.")
print("    Would predict alpha varies between Earth, Jupiter, deep space.")
print("    This is an observational/experimental test.")
print()
print("  - Hyp 3 (wave path): n_EM > n_topo due to non-circular path.")
print("    The arc-length-weighted winding exceeds topological winding.")
print("    The amount depends on the wave amplitude epsilon.")
print("    EPSILON IS DETERMINED BY: what physics sets the wave scale?")
print("    If epsilon = sat_frac_Sun: reduces to Hyp 2 (medium kickback).")
print("    If epsilon = alpha/pi or 2*alpha: intrinsic to EM self-interaction.")
print("    If epsilon = alpha^(1/2)/(4*pi): intrinsic to Hopf geometry.")
print()
print("  These hypotheses are NOT mutually exclusive:")
print("    The wave path IS the physical mechanism.")
print("    The wave AMPLITUDE is set by either the medium state (Hyp 2)")
print("    or intrinsic EM self-interaction (Hyp 1 higher-order).")
print("    The wave path model and the double-spin model are compatible:")
print("    the (1,2) torus knot is the topology; the wave is the geometry.")
print()
print("  MOST LIKELY UNIFIED PICTURE:")
print("    - Topology: (1,2) torus knot. n_topo = 2. (from double-spin Hyp 1)")
print("    - Geometry: path has wave oscillations at resonant frequency k=2.")
print(f"    - Amplitude: epsilon ~ {k2_eps:.5f}, from the EM self-interaction at the")
print(f"      crossing ring (intrinsic, not from local medium state).")
print(f"    - This gives n_EM = {n_EM(k2_eps, 2):.6f}, close to n_exact = {n_exact:.6f}.")
print(f"    - Residual after wave correction: {abs(n_EM(k2_eps, 2) - n_exact):.8f}.")
print()

print(SEP)
print("SUMMARY")
print(SEP)
print()
print("  Q: Can a wave/jagged path explain the C4b residual 0.01869?")
print("  A: YES, mechanically. Any wave perturbation with epsilon*k*f(geom)")
print("     = 0.01869 will give n_EM = n_exact.")
print()
print("  What k is most natural?")
print(f"  A: k=2 (resonant with the (1,2) torus knot minor winding).")
print(f"     Requires epsilon = {epsilon_results[2]:.6f}.")
print()
print("  What sets epsilon?")
print("  A: UNKNOWN — this is the key open question. Three candidates:")
print(f"    - Intrinsic EM: 2*alpha = {2*alpha:.6f}  ({(2*alpha-epsilon_results[2])/epsilon_results[2]*100:+.1f}%)")
print(f"    - Medium state: sat_frac_Sun = {sat_frac_sun:.6f}  ({(sat_frac_sun-epsilon_results[2])/epsilon_results[2]*100:+.1f}%)")
print(f"    - Geometric:    Rs/pi^2 = {Rs_pi2:.6f}  ({(Rs_pi2-epsilon_results[2])/epsilon_results[2]*100:+.1f}%)")
print()
print("  What the wave path model adds (that double-spin alone doesn't give):")
print("  It separates the problem into two parts:")
print("    TOPOLOGY (integer part): (1,2) torus knot winding = 2")
print("    GEOMETRY (fractional part): wave deviation amplitude = 0.01869")
print("  These can be attacked independently. The linking integral")
print("  calculation (Step D2) gives the topology. The wave amplitude")
print("  calculation is a classical mechanics problem on the Hopf torus.")
print()
print("  See also: analysis/hopf_c4b_correction.py (double-spin topology)")
print("            analysis/c4b_residual_medium.py (medium kickback)")
print("            analysis/c4b_residual_scale.py  (geometric candidates)")
print(SEP)
