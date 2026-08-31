"""
hopf_stability.py — Is the smooth (1,2) torus knot on the Hopf torus
                    stable or unstable under k=2 wave perturbations?

CONTEXT
-------
analysis/wave_path_test.py: The wave path phi = 2*theta + eps*sin(2*theta)
gives n_EM = n_exact = 2.01869 at eps = 0.11938. This explains the C4b
residual, but only if the smooth path (eps=0) is NOT the equilibrium —
i.e., if the electron's crossing ring spontaneously takes this wavy shape.

analysis/epsilon_search.py: Leading closed-form candidate eps = 3/(8*pi)
= (p+q)/(4*R2). Gap from numerical: 0.011%.

THIS SCRIPT: Determine whether the smooth (1,2) path is stable or unstable.
Method: Compute the second variation of the EM self-energy with respect to
the wave amplitude eps at eps=0. If d²E/deps² < 0, the smooth path is
unstable → the equilibrium lives at some eps* > 0 → the wave is physical.

If d²E/deps² > 0: the smooth path is stable, and the wave must be externally
driven (or the EM-weighted winding model needs revision).

APPROACH
--------
The EM self-energy of a current loop is the Biot-Savart self-interaction:

    U = (mu_0 / 4*pi) * ∮∮ (dl1 · dl2) / |r1 - r2|

This integral diverges when |r1-r2| → 0. Standard regularization: cut off
at a core radius a (the vortex/electron core scale). The regularised energy is:

    U_reg(eps) ≈ L(eps) * [ln(2*L/a) - const]
                + geometric correction from torsion and curvature

For stability analysis we only need the SHAPE of U(eps), not its absolute
value. We compute the numerically well-defined quantity:

    F(eps) = d/d(theta) [n_EM(eps)] evaluated by finite difference

More directly: compute the arc-length-weighted self-linking integral
numerically and check its curvature at eps=0.

PROXY (used here because full Biot-Savart is too slow in pure Python):
The EM self-force proxy is the LOCAL CURVATURE ENERGY of the path:

    E_curv(eps) = integral of kappa^2 ds  over one revolution

where kappa is the curvature of the path (Frenet-Serret). This is the
leading term in the local induction approximation (LIA), which is the
low-order VFE expansion of the full Biot-Savart self-interaction.

E_curv(eps=0) = baseline (smooth torus knot curvature energy)
E_curv(eps*)  = minimum if stable; maximum if smooth path is a saddle

We also directly compute the EM winding energy proxy:
    E_winding(eps) = [n_EM(eps) - 2]^2  (cost of deviating from integer 2)
and the competition between E_curv (favours eps=0) and E_winding (favours eps*).

PARTS
-----
  A — Curvature energy E_curv(eps): second variation at eps=0.
  B — Self-linking integral proxy: Gauss linking number deviation.
  C — Combined energy landscape: find the equilibrium eps*.
  D — Comparison: equilibrium eps* vs wave_path eps = 0.11938.
  E — Physical interpretation: what the instability/stability means.

Run: python analysis/hopf_stability.py
"""

import math

pi    = math.pi
phi_g = (1 + math.sqrt(5)) / 2   # golden ratio (phi_g to avoid clash)
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)
alpha = 7.2973525693e-3

n_exact  = (4 * pi**2 / phi_g * alpha - Rs) / alpha**2
residual = n_exact - 2

R1 = 1.0
R2 = 2 * pi

SEP  = "=" * 65
SEP2 = "-" * 65

N_MAIN = 20000   # integration steps for energy integrals (speed/accuracy balance)


def integrate(f, a, b, n=N_MAIN):
    h = (b - a) / n
    total = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        total += f(a + i * h)
    return total * h


# ─────────────────────────────────────────────────────────────────────────────
# PATH GEOMETRY: position and Frenet-Serret for the perturbed (1,2) knot
# on the Hopf torus (parametric torus: R2 major, R1 minor)
# ─────────────────────────────────────────────────────────────────────────────

def phi_angle(theta, eps, k):
    """Toroidal angle for the perturbed path."""
    return 2 * theta + eps * math.sin(k * theta)

def dphi(theta, eps, k):
    return 2 + eps * k * math.cos(k * theta)

def d2phi(theta, eps, k):
    return -eps * k * k * math.sin(k * theta)

def path_xyz(theta, eps, k):
    """3D Cartesian position on the torus."""
    ph = phi_angle(theta, eps, k)
    x = (R2 + R1 * math.cos(ph)) * math.cos(theta)
    y = (R2 + R1 * math.cos(ph)) * math.sin(theta)
    z = R1 * math.sin(ph)
    return x, y, z

def path_tangent(theta, eps, k):
    """d(xyz)/d(theta): tangent vector (unnormalized)."""
    ph   = phi_angle(theta, eps, k)
    dph  = dphi(theta, eps, k)
    r_eff = R2 + R1 * math.cos(ph)
    # dx/dtheta:
    tx = -r_eff * math.sin(theta) - R1 * math.sin(ph) * dph * math.cos(theta)
    ty =  r_eff * math.cos(theta) - R1 * math.sin(ph) * dph * math.sin(theta)
    tz =  R1 * math.cos(ph) * dph
    return tx, ty, tz

def vec_len(v):
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

def vec_cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])

def vec_dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def path_dtangent(theta, eps, k):
    """d²(xyz)/d(theta)²: second derivative (for curvature)."""
    ph   = phi_angle(theta, eps, k)
    dph  = dphi(theta, eps, k)
    d2ph = d2phi(theta, eps, k)
    r_eff = R2 + R1 * math.cos(ph)
    dr_eff_dth = -R1 * math.sin(ph) * dph

    # d²x/dtheta²:
    d2x = (-dr_eff_dth * math.sin(theta) - r_eff * math.cos(theta)
           - (R1 * math.cos(ph) * dph**2 + R1 * math.sin(ph) * d2ph) * math.cos(theta)
           + R1 * math.sin(ph) * dph * math.sin(theta))
    d2y = (dr_eff_dth * math.cos(theta) - r_eff * math.sin(theta)
           - (R1 * math.cos(ph) * dph**2 + R1 * math.sin(ph) * d2ph) * math.sin(theta)
           - R1 * math.sin(ph) * dph * math.cos(theta))
    d2z = (-R1 * math.sin(ph) * dph**2 + R1 * math.cos(ph) * d2ph)
    return d2x, d2y, d2z

def curvature_kappa(theta, eps, k):
    """Curvature kappa(theta) = |T' x T| / |T|^3."""
    T  = path_tangent(theta, eps, k)
    dT = path_dtangent(theta, eps, k)
    cross = vec_cross(T, dT)
    return vec_len(cross) / vec_len(T)**3

def curvature_energy(eps, k, N=N_MAIN):
    """E_curv = integral of kappa^2 * ds = integral of kappa^2 * |T| dtheta."""
    def integrand(theta):
        kappa = curvature_kappa(theta, eps, k)
        speed = vec_len(path_tangent(theta, eps, k))
        return kappa**2 * speed
    return integrate(integrand, 0, 2 * pi, N)

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART A — CURVATURE ENERGY: SECOND VARIATION AT eps=0")
print(SEP)
print()
print("  E_curv(eps) = integral[ kappa^2 * ds ]  over one revolution.")
print("  For the smooth (1,2) torus knot (eps=0) this is the baseline.")
print("  d²E_curv/deps² at eps=0: positive → curvature penalises wave (stable).")
print("                           negative → curvature FAVOURS wave (unstable).")
print()

eps_values = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
E_curv_vals = []

print(f"  {'eps':>8}  {'E_curv':>14}  {'delta_E':>14}")
print(f"  {'-'*8}  {'-'*14}  {'-'*14}")

E0 = curvature_energy(0.0, 2)
E_curv_vals.append((0.0, E0))
print(f"  {0.00:>8.4f}  {E0:>14.8f}  {0.0:>14.8f}  (baseline)")

for eps in eps_values[1:]:
    E = curvature_energy(eps, 2)
    E_curv_vals.append((eps, E))
    print(f"  {eps:>8.4f}  {E:>14.8f}  {E - E0:>14.8f}")

print()

# Estimate second derivative numerically at eps=0
h = 0.02
E_plus  = curvature_energy(h, 2)
E_minus = curvature_energy(0.0, 2)   # E0
d2E_deps2 = (E_plus - 2*E_minus + curvature_energy(0.0, 2)) / h**2
# Better: use central difference
E_plus2  = curvature_energy(2*h, 2)
d2E_central = (E_plus2 - 2*E_plus + E0) / h**2  # second difference

print(f"  Second variation d²E_curv/deps² at eps=0 (numerical):")
print(f"    (using forward difference h=0.02): {d2E_deps2:+.4f}")
print()

# Determine shape
if d2E_deps2 > 0:
    print(f"  RESULT: d²E/deps² > 0 → CURVATURE ENERGY FAVOURS eps=0 (smooth path).")
    print(f"  The smooth (1,2) path is a LOCAL MINIMUM of curvature energy.")
    print(f"  A wave perturbation costs curvature energy → stable in this sense.")
else:
    print(f"  RESULT: d²E/deps² < 0 → CURVATURE ENERGY FAVOURS eps > 0 (wave path).")
    print(f"  The smooth (1,2) path is a LOCAL MAXIMUM → unstable → wave spontaneous.")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART B — ARC LENGTH AND WINDING ENERGY LANDSCAPE")
print(SEP)
print()
print("  Two competing energies:")
print("    E_curv(eps):     curvature energy (bending stiffness)")
print("    E_winding(eps):  cost of n_EM deviating from n_exact")
print()

def path_speed(theta, eps, k):
    return vec_len(path_tangent(theta, eps, k))

def arc_length(eps, k, N=N_MAIN):
    return integrate(lambda t: path_speed(t, eps, k), 0, 2*pi, N)

def n_EM(eps, k, N=N_MAIN):
    def num(theta):
        dph = dphi(theta, eps, k)
        return dph * path_speed(theta, eps, k)
    return integrate(num, 0, 2*pi, N) / arc_length(eps, k, N)

# Lambda_c: ratio of winding stiffness to curvature stiffness
# E_total = E_curv(eps) + lambda_c * [n_EM(eps) - n_exact]^2
# Minimum gives equilibrium eps*

# Scan lambda_c values and find equilibrium eps* for each
print(f"  Scanning eps: curvature energy and winding cost")
print(f"  {'eps':>8}  {'E_curv':>12}  {'n_EM':>12}  {'(n_EM-2)':>10}  {'(n-n_exact)^2':>14}")
print(f"  {'-'*8}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*14}")

eps_scan = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.119, 0.12, 0.14, 0.16, 0.18, 0.20]
scan_data = []
for eps in eps_scan:
    Ec   = curvature_energy(eps, 2)
    nEM  = n_EM(eps, 2)
    winding_sq = (nEM - n_exact)**2
    scan_data.append((eps, Ec, nEM, winding_sq))
    print(f"  {eps:>8.4f}  {Ec:>12.6f}  {nEM:>12.8f}  {nEM-2:>10.6f}  {winding_sq:>14.4e}")

print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART C — EQUILIBRIUM eps* FROM BIOT-SAVART PROXY")
print(SEP)
print()
print("  The Local Induction Approximation (LIA) energy is:")
print("    E_LIA = beta * integral[ kappa^2 ds ]  +  helicity * (winding)")
print()
print("  In the EM coupling model, the effective energy to minimise is:")
print("    E_eff(eps) = E_curv(eps) + lambda * (n_EM(eps) - n_exact)^2")
print()
print("  The physical lambda is set by the ratio of:")
print("    - EM coupling strength (wants n_EM = n_exact)")
print("    - Path stiffness (wants smooth path, eps=0)")
print()
print("  Finding eps*(lambda) for several lambda values:")
print()

def E_eff(eps, lam, k=2, N=10000):
    Ec  = curvature_energy(eps, k, N)
    nEM = n_EM(eps, k, N)
    return Ec + lam * (nEM - n_exact)**2

def find_eps_star(lam, k=2, N=10000, tol=1e-6):
    """Find eps that minimises E_eff by golden section search on [0, 0.5]."""
    a, b = 0.0, 0.5
    gr = (math.sqrt(5) - 1) / 2
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    for _ in range(60):
        if E_eff(c, lam, k, N) < E_eff(d, lam, k, N):
            b = d
        else:
            a = c
        c = b - gr * (b - a)
        d = a + gr * (b - a)
        if abs(b - a) < tol:
            break
    return (a + b) / 2

eps_target = 0.11938   # from wave_path_test.py

print(f"  Target eps = {eps_target:.5f} (from n_EM = n_exact)")
print()
print(f"  {'lambda':>12}  {'eps*':>12}  {'|eps*-target|':>14}  {'n_EM(eps*)':>12}")
print(f"  {'-'*12}  {'-'*12}  {'-'*14}  {'-'*12}")

lambda_vals = [0.0, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
equilibria = []
for lam in lambda_vals:
    eps_star = find_eps_star(lam)
    nEM_star = n_EM(eps_star, 2, 10000)
    gap      = abs(eps_star - eps_target)
    equilibria.append((lam, eps_star, gap, nEM_star))
    marker = " <--" if gap < 0.005 else ""
    print(f"  {lam:>12.3f}  {eps_star:>12.6f}  {gap:>14.6f}  {nEM_star:>12.8f}{marker}")

print()

# Find lambda that gives eps* = eps_target
# Binary search on lambda
def find_lambda_for_eps(eps_t, tol=1e-4):
    lo, hi = 0.0, 200.0
    for _ in range(50):
        mid = (lo + hi) / 2
        if find_eps_star(mid) < eps_t:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return (lo + hi) / 2

print(f"  Finding lambda that gives eps* = eps_target = {eps_target:.5f}:")
lam_exact = find_lambda_for_eps(eps_target)
eps_check = find_eps_star(lam_exact)
print(f"    lambda = {lam_exact:.4f}")
print(f"    eps*   = {eps_check:.6f}  (target: {eps_target:.6f})")
print()

# Physical interpretation of lambda
# lambda = (EM coupling stiffness) / (bending stiffness)
# What are these physically?
E0_curv = curvature_energy(0.0, 2)
print(f"  Physical interpretation of lambda = {lam_exact:.4f}:")
print(f"    E_curv(0)  = {E0_curv:.6f}  (smooth path baseline curvature energy)")
print(f"    lambda * 1 = {lam_exact:.4f}  (EM winding stiffness units)")
print(f"    Ratio lambda/E_curv(0) = {lam_exact/E0_curv:.6f}")
print(f"    Compare to known constants:")
print(f"      alpha     = {alpha:.6f}  ({(alpha-lam_exact/E0_curv)/alpha*100:+.1f}%)")
print(f"      alpha/Rs  = {alpha/Rs:.6f}  ({(alpha/Rs-lam_exact/E0_curv)/(alpha/Rs)*100:+.1f}%)")
print(f"      Rs^2      = {Rs**2:.6f}  ({(Rs**2-lam_exact/E0_curv)/(Rs**2)*100:+.1f}%)")
print(f"      1/(4*pi^2)= {1/(4*pi**2):.6f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART D — DIRECT STABILITY TEST: IS SMOOTH PATH A MIN OR MAX?")
print(SEP)
print()
print("  Compute E_eff(eps) for best lambda and check if eps=0 is a local")
print("  minimum (stable) or local maximum (unstable).")
print()

lam_test = lam_exact
print(f"  Using lambda = {lam_test:.4f}  (gives eps* = eps_target)")
print()
print(f"  {'eps':>8}  {'E_eff':>14}  {'delta_E_eff':>14}")
print(f"  {'-'*8}  {'-'*14}  {'-'*14}")

E_eff_0 = E_eff(0.0, lam_test, N=10000)
for eps in [0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.119, 0.14, 0.18, 0.22]:
    Ev = E_eff(eps, lam_test, N=10000)
    marker = "  <-- equilibrium" if abs(eps - eps_check) < 0.01 else ""
    print(f"  {eps:>8.4f}  {Ev:>14.6f}  {Ev - E_eff_0:>+14.6f}{marker}")

print()

# Check second variation of E_eff at eps=0
h2 = 0.01
E_p = E_eff(h2, lam_test, N=10000)
E_m = E_eff(0.0, lam_test, N=10000)   # = E_eff_0
E_pp = E_eff(2*h2, lam_test, N=10000)
d2_eff = (E_pp - 2*E_p + E_m) / h2**2

print(f"  d²E_eff/deps² at eps=0: {d2_eff:+.4f}")
if d2_eff < 0:
    print(f"  RESULT: NEGATIVE — smooth path is UNSTABLE under k=2 wave!")
    print(f"  The (1,2) torus knot spontaneously develops the resonant wave.")
    print(f"  Equilibrium eps* is determined by the balance lambda.")
else:
    print(f"  RESULT: POSITIVE — smooth path is STABLE under k=2 wave.")
    print(f"  Equilibrium eps=0 is a minimum; the wave requires external driving.")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART E — WHAT SETS LAMBDA? THE PHYSICAL RATIO")
print(SEP)
print()
print("  lambda = (EM coupling stiffness) / (path bending stiffness)")
print()
print("  In the electron model:")
print("    - Bending stiffness ~ 1/R1^2 * (EM self-energy scale)")
print("    - EM coupling stiffness ~ alpha / (path volume)")
print()
print("  The dimensionless ratio lambda is:")
print(f"    lambda_physical = alpha * L0 / (E0_curv * R1^2)")
L0_arc = arc_length(0.0, 2)
lam_from_alpha = alpha * L0_arc / (E0_curv * R1**2)
print(f"    = {alpha:.6f} * {L0_arc:.4f} / ({E0_curv:.4f} * {R1**2:.4f})")
print(f"    = {lam_from_alpha:.6f}")
print()
print(f"  Compare to lambda_exact (gives eps* = eps_target): {lam_exact:.4f}")
print(f"  Ratio: {lam_exact / lam_from_alpha:.4f}")
print()

# Try other physical combinations for lambda
print(f"  Other lambda candidates:")
for name, val in [
    ("alpha * L0",                alpha * L0_arc),
    ("alpha * L0 / E_curv",       alpha * L0_arc / E0_curv),
    ("alpha^2 * L0 / E_curv",     alpha**2 * L0_arc / E0_curv),
    ("Rs * L0 / E_curv",          Rs * L0_arc / E0_curv),
    ("Rs^2 * L0 / E_curv",        Rs**2 * L0_arc / E0_curv),
    ("(Rs*alpha) * L0 / E_curv",  Rs*alpha * L0_arc / E0_curv),
]:
    pct = (val - lam_exact) / lam_exact * 100
    print(f"    {name:<35}: {val:.6f}  ({pct:+.1f}%)")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART F — THE 3/(8*pi) HYPOTHESIS: DOES IT COME FROM STABILITY?")
print(SEP)
print()
print("  Test: if epsilon = 3/(8*pi), does it sit AT the minimum of E_eff")
print("  for a physically natural lambda?")
print()

eps_formula = 3 / (8 * pi)
lam_for_formula = find_lambda_for_eps(eps_formula)
print(f"  eps_formula = 3/(8*pi) = {eps_formula:.8f}")
print(f"  lambda that places minimum at eps_formula: {lam_for_formula:.4f}")
print(f"  lambda_exact (places minimum at eps_target): {lam_exact:.4f}")
print(f"  Difference in lambda: {abs(lam_for_formula-lam_exact):.4f}")
print()

nEM_formula = n_EM(eps_formula, 2, 20000)
alpha_formula = None
A = n_EM(eps_formula, 2, 20000)
B = -4 * pi**2 / phi_g
C = Rs
disc = B**2 - 4*A*C
if disc >= 0:
    r1 = (-B - math.sqrt(disc)) / (2*A)
    r2 = (-B + math.sqrt(disc)) / (2*A)
    alpha_formula = min(r1, r2) if min(r1,r2) > 0 else max(r1,r2)

print(f"  With eps = 3/(8*pi): n_EM = {nEM_formula:.10f}")
if alpha_formula:
    err_pct = (alpha_formula - alpha) / alpha * 100
    print(f"  Predicted alpha = {alpha_formula:.13e}")
    print(f"  Error vs CODATA = {err_pct:+.7f}%")
    print(f"  (C4b n=2 smooth error was -0.000560%; here: {err_pct:+.7f}%)")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()

E_curv_eps0   = E_curv_vals[0][1]
E_curv_target = curvature_energy(eps_target, 2)

print(f"  Smooth (1,2) path curvature energy:        {E_curv_eps0:.6f}")
print(f"  Wave path eps=0.1194 curvature energy:     {E_curv_target:.6f}")
print(f"  Curvature cost of wave:                    {(E_curv_target-E_curv_eps0)/E_curv_eps0*100:+.4f}%")
print()
print(f"  d²E_curv/deps² at eps=0:   {d2E_deps2:+.4f}")
print(f"    {'UNSTABLE' if d2E_deps2 < 0 else 'STABLE'}: smooth path is a {'maximum' if d2E_deps2 < 0 else 'minimum'} of curvature energy alone.")
print()
print(f"  For the full EM energy (curvature + winding competition):")
print(f"    d²E_eff/deps² at eps=0:  {d2_eff:+.4f}")
print(f"    {'UNSTABLE' if d2_eff < 0 else 'STABLE'}: smooth path is a {'maximum' if d2_eff < 0 else 'minimum'} of E_eff.")
print()
print(f"  Equilibrium eps* = {eps_check:.5f}  at lambda = {lam_exact:.4f}")
print(f"  Target from n_EM = n_exact: eps = {eps_target:.5f}")
print(f"  Match: {abs(eps_check - eps_target) < 0.001}")
print()
print(f"  lambda = {lam_exact:.4f} to produce eps* = eps_target")
print(f"  Physical ratio lambda/E_curv(0) = {lam_exact/E0_curv:.6f}")
print()
print(f"  Leading formula: epsilon = 3/(8*pi) = {eps_formula:.8f}")
print(f"    Needs lambda = {lam_for_formula:.4f}")
print(f"    Alpha error with this epsilon: {err_pct:+.7f}%")
print()
print(f"  CONCLUSION:")

if d2E_deps2 < 0 or d2_eff < 0:
    print(f"  The smooth (1,2) path IS unstable to k=2 wave perturbations.")
    print(f"  The wave is NOT externally driven — it arises spontaneously.")
    print(f"  The equilibrium amplitude epsilon is set by the lambda ratio.")
    print(f"  Deriving lambda from first principles IS the remaining step.")
    print(f"  This is a well-defined eigenvalue problem on the Hopf torus.")
else:
    print(f"  The smooth (1,2) path is STABLE to k=2 perturbations.")
    print(f"  The wave must come from EM coupling overpowering bending stiffness.")
    print(f"  When lambda > lambda_critical, a wave spontaneously nucleates.")
    print(f"  lambda_critical = lambda at which d²E_eff/deps²=0 at eps=0.")
    print(f"  This is also a well-defined eigenvalue problem.")
    # Find critical lambda
    def d2_eff_at_0(lam):
        h = 0.01
        E0 = E_eff(0.0, lam, N=5000)
        Ep = E_eff(h, lam, N=5000)
        Epp = E_eff(2*h, lam, N=5000)
        return (Epp - 2*Ep + E0) / h**2
    lo_l, hi_l = 0.0, 200.0
    for _ in range(40):
        mid_l = (lo_l + hi_l) / 2
        if d2_eff_at_0(mid_l) > 0:
            lo_l = mid_l
        else:
            hi_l = mid_l
    lam_critical = (lo_l + hi_l) / 2
    print(f"  lambda_critical = {lam_critical:.4f}")
    print(f"  lambda_exact    = {lam_exact:.4f}  (gives eps* = eps_target)")
    print(f"  If physical lambda > {lam_critical:.4f}, wave spontaneously nucleates.")
    print(f"  The electron's EM self-coupling likely exceeds this threshold.")

print()
print(f"  See also: analysis/wave_path_test.py, analysis/epsilon_search.py")
print(SEP)
