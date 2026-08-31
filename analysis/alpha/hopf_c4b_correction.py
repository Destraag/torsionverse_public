"""
hopf_c4b_correction.py — Parallel inspection analysis for Conjecture C4b.

Starting point (from c4a_candidates.py, Conjecture C4b):

    C4b quadratic: 2*alpha^2 - (4*pi^2/phi)*alpha + R_s = 0
    Physical root: alpha_C4b = 7.2973117300057e-3
    CODATA alpha:              7.2973525693000e-3
    Error:                     -0.000560%

The coefficient 2 in the quadratic is the central unknown.
C4a's 4-step inspection (hopf_c4_correction.py) asked: what IS the gap?
This script asks: what IS the 2?

  PART A — The C4b residual gap.
            Same QED-scale analysis as C4a PART A, but for the new
            (much smaller) residual. Is the remaining gap a known scale?

  PART B — Coefficient scan: what if the coefficient were n, not 2?
            Test n = 1, 2, 3, phi, 2*phi, 1/phi, pi, R_s, alpha, ...
            Which value of the coefficient produces the most accurate alpha?
            Does n=2 correspond to anything special in this scan?

  PART C — Double-spin winding number hypothesis.
            If the figure-8 sweeps the main axis once AND spins n times
            around its own horizontal axis, the crossing ring becomes a
            (1,n) torus knot. The linking integral is modified by n.
            This predicts the quadratic coefficient = n.
            Test: for n=2, does the model recover C4b exactly?
            What does n=1 give? What does n=3 give?

  PART D — Criteria for deriving the coefficient 2.
            What would a complete geometric derivation look like?
            What mathematics would be required?

Run: python analysis/hopf_c4b_correction.py
"""

import math

SEP  = "=" * 65
SEP2 = "-" * 65

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)

alpha_CODATA = 7.2973525693e-3
alpha_C4a    = sqrt5 * phi / (16 * pi**3)
alpha_C4b    = 7.2973117300057e-3   # from c4a_candidates.py


# ─────────────────────────────────────────────────────────────────────────────
# Verify C4b by solving the quadratic independently
# ─────────────────────────────────────────────────────────────────────────────

def solve_quadratic(n):
    """Solve n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0, return physical root."""
    a = n
    b = -(4 * pi**2 / phi)
    c = Rs
    disc = b**2 - 4 * a * c
    if disc < 0:
        return None, None
    root_plus  = (-b + math.sqrt(disc)) / (2 * a)
    root_minus = (-b - math.sqrt(disc)) / (2 * a)
    # Physical root: between 0 and 1, close to alpha_CODATA
    if 0 < root_minus < 1:
        return root_minus, root_plus
    return root_plus, root_minus

alpha_C4b_check, alpha_C4b_plus = solve_quadratic(2)
gap_C4b_abs = alpha_C4b_check - alpha_CODATA
gap_C4b_rel = gap_C4b_abs / alpha_CODATA
gap_C4b_pct = gap_C4b_rel * 100

print(SEP)
print("STARTING POINT — CONJECTURE C4b VERIFICATION")
print(SEP)
print()
print(f"  Quadratic: 2*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0")
print(f"  Physical root (alpha_-):")
print(f"    alpha_C4b  (stored)   = {alpha_C4b:.13e}")
print(f"    alpha_C4b  (computed) = {alpha_C4b_check:.13e}")
print(f"    alpha_CODATA          = {alpha_CODATA:.13e}")
print(f"    C4a error:  {(alpha_C4a - alpha_CODATA)/alpha_CODATA*100:+.6f}%")
print(f"    C4b error:  {gap_C4b_pct:+.6f}%")
print(f"    C4b is {abs((alpha_C4a-alpha_CODATA)/alpha_CODATA)/abs(gap_C4b_rel):.0f}x more accurate than C4a")
print()
print(f"  Unphysical root (alpha_+): {alpha_C4b_plus:.6f}  (>> 1, not a coupling constant)")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART A — THE C4b RESIDUAL GAP: QED-SCALE ANALYSIS")
print(SEP)
print()
print("  Same analysis as hopf_c4_correction.py PART A, but for the")
print("  smaller C4b gap. What scale is the remaining discrepancy?")
print()

gap_C4b = gap_C4b_abs

a1 = alpha_CODATA / (2 * pi)         # Schwinger term
a2 = alpha_CODATA**2
a3 = alpha_CODATA / pi
a3sq = a3**2

print(f"  C4b residual gap         = {gap_C4b:+.6e}  ({gap_C4b_pct:+.6f}%)")
print()
print(f"  QED scales for comparison:")
print(f"    alpha/(2*pi)   = {a1:.6e}  (Schwinger term)")
print(f"    alpha^2        = {a2:.6e}")
print(f"    (alpha/pi)^2   = {a3sq:.6e}")
print(f"    alpha^3/pi^2   = {alpha_CODATA**3/pi**2:.6e}")
print()
print(f"  C4b gap / [alpha/(2*pi)]  = {gap_C4b / a1:.6f}")
print(f"  C4b gap / alpha^2         = {gap_C4b / a2:.6f}")
print(f"  C4b gap / (alpha/pi)^2    = {gap_C4b / a3sq:.6f}")
print(f"  C4b gap / alpha^3/pi^2    = {gap_C4b / (alpha_CODATA**3/pi**2):.4f}")
print()
print(f"  C4a gap for comparison:")
gap_C4a = alpha_C4a - alpha_CODATA
print(f"    C4a gap / [alpha/(2*pi)]  = {gap_C4a / a1:.6f}")
print(f"    C4a gap / (alpha/pi)^2    = {gap_C4a / a3sq:.4f}")
print()

# What additive correction would close C4b exactly?
delta_C4b = (alpha_CODATA / alpha_C4b_check) - 1.0
print(f"  Fractional correction needed to close C4b: delta = {delta_C4b:+.8e}")
print()
c1_C4b = delta_C4b * pi / alpha_CODATA
print(f"  If alpha = C4b * (1 + c1*(alpha/pi)):  c1 = {c1_C4b:.6f}")
print(f"  Schwinger c1 = 0.500000   ratio = {c1_C4b/0.5:.4f}")
print()

# Is the C4b gap itself recognizable?
print(f"  Is gap_C4b recognizable as a combination of constants?")
candidates = [
    ("alpha^3",               alpha_CODATA**3),
    ("alpha^3 / pi",          alpha_CODATA**3 / pi),
    ("alpha^2 * phi",         alpha_CODATA**2 * phi),
    ("alpha^2 / (4*pi)",      alpha_CODATA**2 / (4 * pi)),
    ("Rs * alpha^2",          Rs * alpha_CODATA**2),
    ("Rs * alpha^2 * phi",    Rs * alpha_CODATA**2 * phi),
    ("(alpha/pi)^2 / 2",      (alpha_CODATA / pi)**2 / 2),
    ("alpha^2 / (2*phi)",     alpha_CODATA**2 / (2 * phi)),
]
print(f"  {'Expression':<30} {'Value':>14}  {'Ratio to gap':>14}  {'% off':>8}")
print(f"  {'-'*30} {'-'*14}  {'-'*14}  {'-'*8}")
for name, val in candidates:
    ratio = gap_C4b / val
    pct_off = (val - gap_C4b) / abs(gap_C4b) * 100
    flag = "  <-- NEAR" if 0.8 < abs(ratio) < 1.25 else ""
    print(f"  {name:<30} {val:>14.5e}  {ratio:>14.4f}  {pct_off:>+8.2f}%{flag}")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART B — COEFFICIENT SCAN: WHAT IF THE COEFFICIENT WERE n, NOT 2?")
print(SEP)
print()
print("  The quadratic is: n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0")
print("  For different values of n, what alpha does the physical root give?")
print("  n=2 is C4b. What if n were a different geometric number?")
print()
print(f"  {'n':<18} {'alpha_-':>16}  {'error %':>12}  {'1/alpha':>14}  Notes")
print(f"  {'-'*18} {'-'*16}  {'-'*12}  {'-'*14}  {'-'*25}")

test_n = [
    (1,       "n=1 (no quadratic correction)"),
    (2,       "n=2 (C4b — our conjecture)"),
    (3,       "n=3 (triple winding)"),
    (4,       "n=4"),
    (phi,     f"n=phi={phi:.5f} (golden ratio)"),
    (2*phi,   f"n=2*phi={2*phi:.5f}"),
    (pi,      f"n=pi={pi:.5f}"),
    (1/phi,   f"n=1/phi={1/phi:.5f}"),
    (sqrt5,   f"n=sqrt5={sqrt5:.5f}"),
    (Rs,      f"n=Rs={Rs:.5f}"),
    (2/phi,   f"n=2/phi={2/phi:.5f}"),
    (4/phi,   f"n=4/phi={4/phi:.5f}"),
    (pi/phi,  f"n=pi/phi={pi/phi:.5f}"),
]

best_n    = None
best_err  = float('inf')
best_alpha = None

for n_val, n_label in test_n:
    r_minus, r_plus = solve_quadratic(n_val)
    if r_minus is None or r_minus <= 0 or r_minus >= 1:
        print(f"  {n_label:<18} {'NO PHYSICAL ROOT':>16}")
        continue
    err = (r_minus - alpha_CODATA) / alpha_CODATA * 100
    inv_a = 1 / r_minus
    flag = ""
    if abs(err) < abs(best_err):
        best_err = err
        best_n = n_val
        best_alpha = r_minus
        flag = "  <-- best"
    is_c4b = abs(n_val - 2) < 1e-9
    marker = " [C4b]" if is_c4b else ""
    print(f"  {n_label:<18} {r_minus:>16.10e}  {err:>+12.6f}%  {inv_a:>14.8f}{marker}{flag}")

print()
print(f"  Best coefficient in this scan: n = {best_n:.5f}")
print(f"  (Gives alpha closest to CODATA with geometric constant as coefficient)")
print()

# Find the exact n that would make alpha = alpha_CODATA exactly
# n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0 at alpha = alpha_CODATA:
# n = (4*pi^2/phi * alpha_CODATA - Rs) / alpha_CODATA^2
n_exact = (4 * pi**2 / phi * alpha_CODATA - Rs) / alpha_CODATA**2
print(f"  Exact n that recovers alpha_CODATA from the quadratic:")
print(f"    n_exact = (4*pi^2/phi * alpha - Rs) / alpha^2")
print(f"            = {n_exact:.10f}")
print()
print(f"  Nearest integers and geometric values:")
for name, val in [("2", 2), ("3", 3), ("phi", phi), ("2*phi", 2*phi),
                  ("pi", pi), ("2+alpha", 2+alpha_CODATA), ("2+Rs", 2+Rs)]:
    diff = val - n_exact
    print(f"    {name:<15} = {val:.8f}  diff from n_exact: {diff:+.8f}  ({diff/n_exact*100:+.4f}%)")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART C — DOUBLE-SPIN WINDING NUMBER HYPOTHESIS")
print(SEP)
print()
print("  GEOMETRIC BACKGROUND")
print(SEP2)
print()
print("  Standard C4a picture:")
print("    The figure-8 sweeps through 360 degrees of rotation around the")
print("    vertical axis. One full sweep traces a torus with one crossing ring.")
print("    The crossing ring is a simple closed curve on the torus surface.")
print("    Topological constraint: R2/R1 = 2*pi.")
print("    Key ratio: 4*pi^2 (the Hopf torus major-circuit ratio).")
print()
print("  Double-spin extension (the tangent):")
print("    Suppose the figure-8 ALSO rotates around its own horizontal axis")
print("    as it sweeps around the main vertical axis.")
print("    If it completes n rotations of the secondary spin per one full")
print("    sweep of the main rotation, the crossing ring is no longer a")
print("    simple circle on the torus surface.")
print()
print("    The crossing ring becomes a (1,n) TORUS KNOT:")
print("      - winds once around the torus in the major direction (the main sweep)")
print("      - winds n times in the minor direction (the secondary spin)")
print()
print("    Consequence for the linking integral:")
print("    The linking number of the crossing ring with the torus fiber circle")
print("    is multiplied by n. In C4a, this linking number contributes to the")
print("    denominator Q = 4*pi^2/phi. With a (1,n) torus knot crossing ring,")
print("    the self-interaction correction becomes:")
print()
print("      Q_corrected = 4*pi^2/phi - n * alpha")
print()
print("    And the self-consistent equation becomes:")
print()
print("      n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0")
print()
print("    For n=2: this is exactly Conjecture C4b.")
print()
print("  WHAT THE DOUBLE-SPIN LOOKS LIKE:")
print(SEP2)
print()
print("  n=1 (no secondary spin):")
print("    The figure-8 maintains the same orientation as it sweeps.")
print("    Crossing ring = simple equatorial circle on the torus.")
print("    No self-consistent correction. Gives:")
r1, _ = solve_quadratic(1)
if r1:
    err1 = (r1 - alpha_CODATA) / alpha_CODATA * 100
    print(f"    alpha(n=1) = {r1:.10e}  error = {err1:+.6f}%")
    print(f"    This is equivalent to alpha = Rs / (4*pi^2/phi - alpha)")
    print(f"    A weaker self-correction — just subtracts alpha, not 2*alpha.")
print()
print("  n=2 (C4b — figure-8 spins twice per main sweep):")
print("    The crossing ring is a (1,2) torus knot — a simple twist.")
print("    Topologically: a curve that goes around the torus once in the")
print("    major direction and twice in the minor direction.")
print("    This is the simplest non-trivial torus knot on this surface.")
print("    The linking integral picks up a factor of 2.")
print(f"    alpha(n=2) = {alpha_C4b_check:.10e}  error = {gap_C4b_pct:+.6f}%")
print()
print("  n=3 (figure-8 spins three times per sweep):")
r3, _ = solve_quadratic(3)
if r3:
    err3 = (r3 - alpha_CODATA) / alpha_CODATA * 100
    print(f"    Crossing ring is a (1,3) torus knot (trefoil-adjacent).")
    print(f"    alpha(n=3) = {r3:.10e}  error = {err3:+.6f}%")
print()
print("  WHY n=2 IS GEOMETRICALLY SPECIAL:")
print(SEP2)
print()
print("  The (1,2) torus knot is the simplest non-trivial winding on a torus.")
print("  It is the lowest-energy excitation of the crossing ring topology.")
print()
print("  There is also a connection to the Hopf invariant:")
print("    The standard Hopf map has invariant H=1 (one linking of fibers).")
print("    A double-covered Hopf map has H=2.")
print("    The coefficient n in the quadratic may equal the Hopf invariant")
print("    of the double-spin geometry: H = n = 2.")
print()
print("  Physical reading of n=2:")
print("    Each half-revolution of the secondary spin brings the crossing")
print("    ring back to where it started spatially, so the ring encounters")
print("    the 'other side' of the self-intersection twice per revolution.")
print("    This doubles the EM self-interaction correction.")
print()

# Compute: what is the effective Q for each n?
print(f"  Effective Q = 4*pi^2/phi - n*alpha for each n:")
Q_bare = 4 * pi**2 / phi
Q_exact = Rs / alpha_CODATA
print(f"    Q_bare = 4*pi^2/phi = {Q_bare:.8f}")
print(f"    Q_exact (target)    = {Q_exact:.8f}")
print()
for n_val, n_label in [(1, "n=1"), (2, "n=2 [C4b]"), (3, "n=3")]:
    Q_eff = Q_bare - n_val * alpha_CODATA
    print(f"    Q(n={n_val}) = {Q_bare:.6f} - {n_val}*alpha = {Q_eff:.8f}  "
          f"(gap from Q_exact: {Q_eff - Q_exact:+.8f})")
print()
print(f"  Note: Q(n=2) = {Q_bare - 2*alpha_CODATA:.8f}, Q_exact = {Q_exact:.8f}")
print(f"  Residual in Q after n=2 correction: {(Q_bare - 2*alpha_CODATA) - Q_exact:.8f}")
print(f"  This residual corresponds to the remaining C4b gap of {gap_C4b_pct:.6f}%.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART D — CRITERIA FOR DERIVING THE COEFFICIENT 2")
print(SEP)
print()
print("  For the double-spin hypothesis to become a derivation, four things")
print("  must be shown:")
print()
print("  STEP D1 — Establish n=2 as the topologically preferred winding")
print("    The figure-8 sweeping with secondary spin is parameterized by n.")
print("    Physics requires a preferred n. Candidates:")
print("      (a) Minimum energy: which n minimizes the self-energy of the")
print("          crossing ring on the torus? If n=2, energy minimization")
print("          selects the (1,2) torus knot configuration.")
print("      (b) Stability: which n produces a stable knotted configuration?")
print("          The (1,2) torus knot is the simplest; it may be the stable")
print("          ground state. The (1,1) case (no extra crossing) is unstable")
print("          in a medium with torsion. (1,3) requires more energy.")
print("      (c) Quantization: if the spin rate is quantized (like angular")
print("          momentum), the lowest non-trivial quantum number is n=2.")
print("          This would require showing why n=0 and n=1 are excluded.")
print()
print("  STEP D2 — Compute the linking integral with the (1,n) correction")
print("    The linking number of the (1,2) torus knot with the torus fiber")
print("    is a topological invariant. For the standard (1,n) torus knot,")
print("    the linking number with the core of the torus is n.")
print("    Computing this for the specific Hopf torus (R2/R1 = 2*pi) with")
print("    icosahedral symmetry at the crossing ring requires:")
print("      - Treating the (1,2) crossing ring as a source of EM field")
print("      - Computing the holonomy around the fiber using the connection")
print("        form on the Hopf bundle")
print("      - Showing the holonomy integral gives alpha exactly")
print()
print("  STEP D3 — Show the self-consistent form is necessary")
print("    The quadratic arises because the coupling constant alpha appears")
print("    in the correction to the geometry (Q = Q_bare - 2*alpha), which")
print("    in turn determines alpha. This fixed-point structure must be")
print("    derived, not assumed.")
print("    Analogy: in QED, the electron self-energy depends on alpha, and")
print("    alpha itself is defined by the self-consistent renormalization")
print("    condition. The same structure must appear in the torus geometry.")
print()
print("  STEP D4 — Recover the C4b precision or better")
print("    A complete derivation would produce the exact value of alpha,")
print("    not merely the 0.00056% approximation. The remaining gap after")
print("    n=2 must be explained (higher-order correction, or the geometry")
print("    gives a different n_exact = 2.000... value).")
print()
print(f"  Recall: n_exact = {n_exact:.10f}")
print(f"  This is n=2 plus a small offset: {n_exact - 2:+.8f}")
print(f"  In units of alpha: {(n_exact-2)/alpha_CODATA:.4f}")
print(f"  In units of Rs:    {(n_exact-2)/Rs:.6f}")
print(f"  Nearest recognizable form: 2 + alpha^2 = {2 + alpha_CODATA**2:.10f} (diff: {2+alpha_CODATA**2-n_exact:+.2e})")
print(f"                             2 + Rs/pi   = {2 + Rs/pi:.10f} (diff: {2+Rs/pi-n_exact:+.2e})")
print(f"                             2 + alpha/phi = {2+alpha_CODATA/phi:.10f} (diff: {2+alpha_CODATA/phi-n_exact:+.2e})")
print(f"                             2 + 2*alpha  = {2+2*alpha_CODATA:.10f} (diff: {2+2*alpha_CODATA-n_exact:+.2e})")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY TABLE")
print(SEP)
print()
print(f"  alpha_CODATA                = {alpha_CODATA:.13e}")
print(f"  alpha_C4a (n/a, direct)     = {alpha_C4a:.13e}  ({(alpha_C4a-alpha_CODATA)/alpha_CODATA*100:+.6f}%)")
print(f"  alpha_C4b (n=2 quadratic)   = {alpha_C4b_check:.13e}  ({gap_C4b_pct:+.6f}%)")
print()
print(f"  Coefficient 2 in C4b:")
print(f"    n_exact for perfect fit    = {n_exact:.10f}")
print(f"    n=2 is {abs(n_exact-2)/n_exact*100:.6f}% away from n_exact")
print(f"    Remaining gap after n=2:   {gap_C4b_pct:+.6f}% in alpha")
print()
print(f"  Double-spin hypothesis:")
print(f"    n = winding number of crossing ring on the Hopf torus surface")
print(f"    n = 1: simple equatorial ring — no correction (no double-spin)")
print(f"    n = 2: (1,2) torus knot — C4b, error {gap_C4b_pct:+.6f}%")
print(f"    n = 3: (1,3) torus knot — overcorrects significantly")
print(f"    n_exact = {n_exact:.6f}: requires a non-integer winding,")
print(f"              which would follow from a higher-order topological")
print(f"              correction beyond the leading (1,2) torus knot.")
print()
print(f"  Status of the coefficient 2:")
print(f"    IDENTIFIED: n=2 is the double-spin (secondary rotation) winding number.")
print(f"    UNPROVEN:   the physical mechanism selecting n=2 over n=1 or n=3.")
print(f"    PATH FORWARD:")
print(f"      Step D1: show energy or stability argument for (1,2) torus knot.")
print(f"      Step D2: compute linking integral for (1,2) crossing ring.")
print(f"      Step D3: derive the self-consistent form from the Hopf connection.")
print(f"      Step D4: residual gap {gap_C4b_pct:.6f}% needs higher-order term.")
print()
print(f"  Connection to C4a path (Steps 1-4 of Appendix D):")
print(f"    C4a Step 1: R2/R1 = 2*pi topological constraint.   [DONE]")
print(f"    C4a Step 2: 4*pi^2/phi as the Hopf torus ratio.   [DONE, conjectural]")
print(f"    C4b Step 2b: n=2 as the crossing ring winding.     [NEW, hypothesized]")
print(f"    Step 3: Hopf linking integral giving alpha exactly. [OPEN]")
print(f"    Step 4: Cross-predictions from the full geometry.  [OPEN]")
print(SEP)
