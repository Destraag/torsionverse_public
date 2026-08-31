"""
biot_savart_min.py — Gap 1 and Gap 2: Biot-Savart self-energy minimisation

BACKGROUND
----------
From analysis/alpha/hopf_stability.py:
  - The smooth (1,2) torus knot is UNSTABLE (d2E/deps2 = -15.1 < 0).
  - The LIA (local induction approximation) curvature energy alone places
    the minimum at eps ~ 0.153, not the target eps = 0.1194.
  - The FULL non-local Biot-Savart self-energy is the correct functional.

TWO ROUTES
----------
Route A — ARC LENGTH (thin-wire limit, pure geometry, fast):
  In the Biot-Savart thin-wire limit:
    E_self(eps) ≈ L(eps) * [log(2L/a) - C]
  Since log(2L/a) >> 1 and varies slowly with eps, the minimum of E_self
  tracks the minimum of the arc length L(eps).
  This reduces Gap 1 to: where does L(eps) have its minimum?
  No double integral required. Pure 1D geometry.

Route B — FULL BIOT-SAVART (double integral, non-local correction):
  E_BS(eps) = integral integral T1.T2 / sqrt(|r1-r2|^2 + a^2) dth1 dth2
  Slower (N^2 evaluations) but gives the non-local correction and the
  Gap 2 energy scale (potentially Rs).

GAP 1 TEST:
  Is the arc length minimum at eps = 3/(8*pi) = (p+q)/(4*R2) = 0.11937?
  If yes: Gap 1 closes from pure geometry of the Hopf torus + (1,2) knot.

GAP 2 TEST:
  Does the non-local Biot-Savart energy scale equal Rs = sqrt(5)/(4*pi)?
  Is Rs = sqrt(p^2+q^2)/(2*R2) the natural self-energy scale of the (p,q)
  torus knot on the Hopf torus?

Run: python analysis/alpha/biot_savart_min.py
"""

import math

pi    = math.pi
phi_g = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)
alpha = 7.2973525693e-3

n_exact  = (4 * pi**2 / phi_g * alpha - Rs) / alpha**2
residual = n_exact - 2

R1 = 1.0
R2 = 2 * pi   # Hopf torus constraint: R2/R1 = 2*pi
p, q = 1, 2   # (p,q) torus knot winding numbers

eps_target  = 0.1193795395   # from epsilon_search.py (N=50000)
eps_formula = 3 / (8 * pi)   # = (p+q)/(4*R2) = 3/(8*pi)
eps_Rs      = Rs * 2 / 3     # Formula B from epsilon_search.py

N_ARC = 20000   # integration steps for arc length (1D, fast)
N_BS  = 80      # grid size for Biot-Savart double integral (N^2 evaluations)

SEP  = "=" * 65
SEP2 = "-" * 65

def integrate(f, a, b, n):
    h = (b - a) / n
    total = 0.5 * (f(a) + f(b))
    for i in range(1, n):
        total += f(a + i * h)
    return total * h

# ─────────────────────────────────────────────────────────────────────────────
# PATH GEOMETRY
# ─────────────────────────────────────────────────────────────────────────────

def phi_val(theta, eps):
    return q * theta + eps * math.sin(q * theta)

def dphi_val(theta, eps):
    return q + eps * q * math.cos(q * theta)

def path_xyz(theta, eps):
    ph = phi_val(theta, eps)
    r_eff = R2 + R1 * math.cos(ph)
    return (r_eff * math.cos(theta),
            r_eff * math.sin(theta),
            R1 * math.sin(ph))

def tangent_xyz(theta, eps):
    ph   = phi_val(theta, eps)
    dph  = dphi_val(theta, eps)
    r_eff = R2 + R1 * math.cos(ph)
    dr_dth = -R1 * math.sin(ph) * dph
    tx = dr_dth * math.cos(theta) - r_eff * math.sin(theta)
    ty = dr_dth * math.sin(theta) + r_eff * math.cos(theta)
    tz = R1 * math.cos(ph) * dph
    return (tx, ty, tz)

def speed(theta, eps):
    t = tangent_xyz(theta, eps)
    return math.sqrt(t[0]**2 + t[1]**2 + t[2]**2)

def arc_length(eps, N=N_ARC):
    return integrate(lambda t: speed(t, eps), 0, 2*pi, N)

def dL_deps(eps, N=N_ARC, h=1e-5):
    """Numerical derivative of arc length with respect to eps."""
    return (arc_length(eps+h, N) - arc_length(eps-h, N)) / (2*h)


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART A — ARC LENGTH L(eps): SCAN AND MINIMUM")
print(SEP)
print()
print("  Thin-wire Biot-Savart: E_self ~ L(eps) * log(2L/a)")
print("  Minimum of E_self ~ minimum of arc length L(eps).")
print()
print("  Scanning L(eps) for the k=2 perturbed (1,2) torus knot:")
print()

L0 = arc_length(0.0)
print(f"  L(0) = {L0:.8f}  (smooth (1,2) torus knot baseline)")
print()

eps_scan = [0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.11, 0.12,
            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.20, 0.22, 0.25]
L_vals = []

print(f"  {'eps':>8}  {'L(eps)':>14}  {'L-L0':>12}  {'(L-L0)/L0 %':>12}")
print(f"  {'-'*8}  {'-'*14}  {'-'*12}  {'-'*12}")

for eps in eps_scan:
    L = arc_length(eps)
    L_vals.append((eps, L))
    print(f"  {eps:>8.4f}  {L:>14.8f}  {L-L0:>+12.8f}  {(L-L0)/L0*100:>+12.6f}%")

print()

# Find minimum by golden section search
gr = (math.sqrt(5) - 1) / 2
a_gs, b_gs = 0.0, 0.3
c_gs = b_gs - gr * (b_gs - a_gs)
d_gs = a_gs + gr * (b_gs - a_gs)
for _ in range(60):
    if arc_length(c_gs) < arc_length(d_gs):
        b_gs = d_gs
    else:
        a_gs = c_gs
    c_gs = b_gs - gr * (b_gs - a_gs)
    d_gs = a_gs + gr * (b_gs - a_gs)
    if abs(b_gs - a_gs) < 1e-9:
        break
eps_min_L = (a_gs + b_gs) / 2
L_min     = arc_length(eps_min_L)

print(f"  Arc length MINIMUM (golden section search):")
print(f"    eps_min(L)   = {eps_min_L:.10f}")
print(f"    L_min        = {L_min:.10f}")
print(f"    L0 - L_min   = {L0-L_min:.10f}  ({(L0-L_min)/L0*100:.6f}%)")
print()
print(f"  Comparison to candidates:")
print(f"    eps_min(L)   = {eps_min_L:.10f}")
print(f"    eps_formula  = 3/(8*pi) = {eps_formula:.10f}  diff = {(eps_formula-eps_min_L)/eps_min_L*100:+.5f}%")
print(f"    eps_target   = {eps_target:.10f}  diff = {(eps_target-eps_min_L)/eps_min_L*100:+.5f}%")
print(f"    eps_Rs       = Rs*2/3   = {eps_Rs:.10f}  diff = {(eps_Rs-eps_min_L)/eps_min_L*100:+.5f}%")
print()

gap1_closed = abs(eps_min_L - eps_formula) / eps_formula * 100
print(f"  GAP 1 STATUS (arc length route):")
if gap1_closed < 0.1:
    print(f"  STRONG MATCH: Arc length minimum at 3/(8*pi) to {gap1_closed:.4f}%.")
    print(f"  If this holds analytically: epsilon = (p+q)/(4*R2) is proven from")
    print(f"  pure Hopf torus geometry. No free parameters.")
elif gap1_closed < 1.0:
    print(f"  NEAR MATCH ({gap1_closed:.3f}%): Arc length minimum close to 3/(8*pi).")
    print(f"  Non-local Biot-Savart correction may close the remaining gap.")
else:
    print(f"  NO MATCH ({gap1_closed:.2f}%): Arc length minimum is NOT at 3/(8*pi).")
    print(f"  Full non-local Biot-Savart required to locate true minimum.")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART B — ARC LENGTH DERIVATIVES AT eps=0")
print(SEP)
print()
print("  Analytic Taylor expansion: L(eps) = L0 + L1*eps + L2*eps^2 + ...")
print("  By the eps -> -eps symmetry analysis:")
print("    If L(eps) is NOT even: dL/deps|_0 != 0 -> L has non-zero linear term")
print("    If L(eps) IS even: dL/deps|_0 = 0 -> minimum condition is 2nd order")
print()

h_deriv = 1e-4
L_plus  = arc_length(+h_deriv)
L_zero  = arc_length(0.0)
L_minus = arc_length(-h_deriv)

dL_0  = (L_plus - L_minus) / (2 * h_deriv)
d2L_0 = (L_plus - 2*L_zero + L_minus) / h_deriv**2

print(f"  dL/deps|_0  = {dL_0:.6f}")
print(f"  d2L/deps2|_0 = {d2L_0:.4f}")
print()

if abs(dL_0) < 1e-3:
    print(f"  dL/deps ~ 0: L(eps) is approximately EVEN in eps.")
    print(f"  Minimum condition: d2L/deps2 = 0 (at the inflection point).")
    print(f"  This occurs at eps beyond the initial curvature.")
else:
    print(f"  dL/deps != 0: L(eps) has a non-zero linear term.")
    print(f"  L DECREASES as eps increases from 0 (dL/deps < 0)." if dL_0 < 0 else
          f"  L INCREASES as eps increases from 0.")
    print(f"  The minimum of L (where dL/deps = 0) is at eps > 0.")

# Find where dL/deps = 0 numerically (the minimum of L)
def dL_num(eps, h=1e-4):
    return (arc_length(eps+h) - arc_length(eps-h)) / (2*h)

# Binary search for zero of dL/deps
lo_d, hi_d = 0.0, 0.3
if dL_num(0.0) * dL_num(0.3) < 0:  # sign change exists
    for _ in range(50):
        mid_d = (lo_d + hi_d) / 2
        if dL_num(mid_d) > 0:
            hi_d = mid_d
        else:
            lo_d = mid_d
        if hi_d - lo_d < 1e-8:
            break
    eps_deriv_zero = (lo_d + hi_d) / 2
    print()
    print(f"  dL/deps = 0 at eps = {eps_deriv_zero:.10f}  (numerical)")
    print(f"  This is the true arc length minimum (cross-check with Part A).")
    print(f"  Matches golden section: {abs(eps_deriv_zero - eps_min_L) < 1e-5}")
else:
    print(f"  No sign change in dL/deps on [0, 0.3] — minimum may be outside range.")

print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART C — FORMULA SEARCH: WHAT IS eps_min(L) IN CLOSED FORM?")
print(SEP)
print()
print(f"  eps_min(L) = {eps_min_L:.10f}  (from golden section, N=20000)")
print()

# Test candidate formulas involving p=1, q=2, R2=2*pi, R1=1, sqrt5, phi
print(f"  Test formulas (p=1, q=2, R2=2*pi, R1=1):")
print()

candidates = [
    ("(p+q)/(4*R2)",         (p+q)/(4*R2)),
    ("(p+q)/(4*R2*R1)",      (p+q)/(4*R2*R1)),
    ("3/(8*pi)",             3/(8*pi)),
    ("sqrt(p^2+q^2)/(4*pi)", sqrt5/(4*pi)),
    ("sqrt(p^2+q^2)/(2*R2)", sqrt5/(2*R2)),
    ("Rs",                   Rs),
    ("Rs*2/3",               Rs*2/3),
    ("q/(4*pi)",             q/(4*pi)),
    ("q*R1/(2*R2)",          q*R1/(2*R2)),
    ("p*q/(4*pi)",           p*q/(4*pi)),
    ("(q/R2)^(1/2)",         math.sqrt(q/R2)),
    ("q/(2*R2-1)",           q/(2*R2-1)),
    ("1/(pi+1/q)",           1/(pi+1/q)),
    ("q/(R2*pi)",            q/(R2*pi)),
    ("sqrt(q/(R2^2+q^2))",   math.sqrt(q/(R2**2+q**2))),
    ("q*R1/R2",              q*R1/R2),
    ("q/(2*(R2+q))",         q/(2*(R2+q))),
    ("q^2/(4*pi*R2)",        q**2/(4*pi*R2)),
    ("phi_g/(4*pi+phi_g)",   phi_g/(4*pi+phi_g)),
    ("(phi_g-1)/R2",         (phi_g-1)/R2),
    ("1/(2*pi+pi/phi_g)",    1/(2*pi+pi/phi_g)),
    ("Rs*phi_g/(2*pi)",      Rs*phi_g/(2*pi)),
    ("Rs/phi_g",             Rs/phi_g),
    ("p/(4*pi-q)",           p/(4*pi-q)),
    ("sqrt5/(2*R2+1)",       sqrt5/(2*R2+1)),
    ("sqrt5/(4*pi+2)",       sqrt5/(4*pi+2)),
]

print(f"  {'Formula':<32} {'Value':>14}  {'% diff from eps_min':>20}")
print(f"  {'-'*32} {'-'*14}  {'-'*20}")

hits = []
for name, val in candidates:
    pct = (val - eps_min_L) / eps_min_L * 100
    hits.append((abs(pct), pct, name, val))

hits.sort()
for _, pct, name, val in hits[:15]:
    marker = "  ***" if abs(pct) < 0.1 else ("  **" if abs(pct) < 0.5 else ("  *" if abs(pct) < 2.0 else ""))
    print(f"  {name:<32} {val:>14.10f}  {pct:>+19.5f}%{marker}")

print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART D — FULL BIOT-SAVART DOUBLE INTEGRAL (Route B)")
print(SEP)
print()
print(f"  E_BS(eps) = integral integral [T1.T2] / sqrt(|r1-r2|^2 + a^2) dth1 dth2")
print(f"  Using N={N_BS}x{N_BS} grid ({N_BS**2} evaluations per eps), a=R1*0.2")
print()

a_core = R1 * 0.2   # regularisation core radius

def dot3(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def dist3_sq(r1, r2):
    return (r1[0]-r2[0])**2 + (r1[1]-r2[1])**2 + (r1[2]-r2[2])**2

def E_BS(eps, N=N_BS, a=a_core):
    """Regularised Biot-Savart self-energy (double trapezoid sum, periodic)."""
    h = 2 * pi / N
    total = 0.0
    for i in range(N):
        th1 = i * h
        r1  = path_xyz(th1, eps)
        t1  = tangent_xyz(th1, eps)
        for j in range(N):
            th2 = j * h
            r2  = path_xyz(th2, eps)
            t2  = tangent_xyz(th2, eps)
            d2  = dist3_sq(r1, r2) + a**2
            kern = dot3(t1, t2) / math.sqrt(d2)
            total += kern
    return total * h * h

print(f"  Computing E_BS at eps=0, h, 2h (h=0.05) for Taylor expansion...")
print(f"  (This may take a minute for N={N_BS})")
print()

h_bs = 0.05
E0_bs  = E_BS(0.0)
Eh_bs  = E_BS(h_bs)
E2h_bs = E_BS(2*h_bs)

print(f"  E_BS(0)    = {E0_bs:.6f}")
print(f"  E_BS(h)    = {Eh_bs:.6f}   (h=0.05)")
print(f"  E_BS(2h)   = {E2h_bs:.6f}   (2h=0.10)")
print()

# Taylor: E(eps) = E0 + a2*eps^2 + a4*eps^4 + ...
# E''(0) = 2*a2 = (Eh - 2*E0 + E_{-h}) / h^2 = (Eh - E0) / h^2 (using symmetry E(-h)=E0?)
# Actually E(-eps) ≠ E(eps) in general (shown in Part B).
# So use: E(h) ≈ E(0) + E'(0)*h + E''(0)/2*h^2
# and: E(2h) ≈ E(0) + E'(0)*2h + E''(0)/2*(2h)^2
# From these two: E'(0) = (4*E(h) - E(2h) - 3*E(0)) / (2*h)
#                 E''(0) = (E(2h) - 2*E(h) + E(0)) / h^2

dE_BS_0  = (4*Eh_bs - E2h_bs - 3*E0_bs) / (2*h_bs)
d2E_BS_0 = (E2h_bs - 2*Eh_bs + E0_bs) / h_bs**2

print(f"  E'_BS(0)  = {dE_BS_0:.6f}   (linear term)")
print(f"  E''_BS(0) = {d2E_BS_0:.4f}  (stability: {'UNSTABLE' if d2E_BS_0 < 0 else 'STABLE'})")
print()

# Find minimum of E_BS via scan
print(f"  Scanning E_BS(eps) for full minimum:")
eps_bs_scan = [0.00, 0.05, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20]
E_bs_vals = [(0.0, E0_bs), (h_bs, Eh_bs), (2*h_bs, E2h_bs)]
print(f"  {'eps':>8}  {'E_BS':>14}  {'E_BS - E_BS(0)':>16}")
print(f"  {'-'*8}  {'-'*14}  {'-'*16}")
print(f"  {0.00:>8.4f}  {E0_bs:>14.4f}  {0.0:>+16.4f}")
print(f"  {h_bs:>8.4f}  {Eh_bs:>14.4f}  {Eh_bs-E0_bs:>+16.4f}")
print(f"  {2*h_bs:>8.4f}  {E2h_bs:>14.4f}  {E2h_bs-E0_bs:>+16.4f}")

prev_E = E2h_bs
for eps in [0.12, 0.14, 0.16, 0.18, 0.20]:
    Ev = E_BS(eps)
    E_bs_vals.append((eps, Ev))
    print(f"  {eps:>8.4f}  {Ev:>14.4f}  {Ev-E0_bs:>+16.4f}")

# Find minimum
min_eps_bs = min(E_bs_vals, key=lambda x: x[1])[0]
print()
print(f"  E_BS minimum in scan: eps ≈ {min_eps_bs:.4f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART E — Rs FROM SELF-ENERGY SCALE (GAP 2)")
print(SEP)
print()
print("  The Biot-Savart self-energy for a thin wire:")
print("    E_BS(eps) = L(eps) * log(2*L(eps) / (a*e)) + E_nonlocal(eps)")
print()
print("  The non-local part E_nonlocal is a-independent.")
print("  Extract it by subtraction: E_nonlocal = E_BS - L*log(2L/a)")
print()

import math as _math
e_euler = _math.exp(1)

def E_local_approx(eps, a=a_core):
    """Thin-wire approximation for the local (log-divergent) part."""
    L = arc_length(eps)
    return L * _math.log(2 * L / (a * e_euler))

E_nonlocal_0 = E0_bs - E_local_approx(0.0)
print(f"  E_BS(0)         = {E0_bs:.6f}")
print(f"  E_local(0)      = {E_local_approx(0.0):.6f}  [L0*log(2L0/(a*e))]")
print(f"  E_nonlocal(0)   = {E_nonlocal_0:.6f}")
print()

# Normalise by L0 and R2 combinations
print(f"  Normalised non-local energy scale E_nonlocal / (geometric quantities):")
print()
norm_candidates = [
    ("L0",                L0),
    ("R2",                R2),
    ("R2^2",              R2**2),
    ("L0*R2",             L0*R2),
    ("L0/(2*pi)",         L0/(2*pi)),
    ("4*pi^2",            4*pi**2),
    ("L0^2/(4*pi^2)",     L0**2/(4*pi**2)),
]
for name, norm in norm_candidates:
    ratio = E_nonlocal_0 / norm if norm != 0 else float('nan')
    # Check if ratio ≈ Rs, Rs^2, sqrt5, phi, etc.
    rs_match  = abs(ratio - Rs)/Rs*100 if ratio > 0 else 999
    rs2_match = abs(ratio - Rs**2)/Rs**2*100 if ratio > 0 else 999
    s5_match  = abs(ratio - sqrt5/4)/sqrt5*100/4 if ratio > 0 else 999
    best = min([(rs_match,'Rs'), (rs2_match,'Rs^2'), (abs(ratio-1/pi)/(1/pi)*100,'1/pi'), (abs(ratio-1/(2*pi))/(1/(2*pi))*100,'1/(2pi)')])
    print(f"  E_nl / {name:<22} = {ratio:>10.6f}   (nearest: {best[1]}, {best[0]:.1f}% off)")

print()

# Check if E_nonlocal relates to Rs directly
print(f"  Key check: E_nonlocal(0) / (L0 * Rs):")
ratio_L_Rs = E_nonlocal_0 / (L0 * Rs)
print(f"    = {ratio_L_Rs:.6f}")
print(f"    Compare to pi = {pi:.6f}  ({(ratio_L_Rs-pi)/pi*100:+.2f}%)")
print(f"    Compare to 2  = 2.000000  ({(ratio_L_Rs-2)/2*100:+.2f}%)")
print(f"    Compare to e  = {e_euler:.6f}  ({(ratio_L_Rs-e_euler)/e_euler*100:+.2f}%)")
print()

# Check E_nonlocal at eps=eps_target vs eps=0
E_nl_target = E_BS(eps_target) - E_local_approx(eps_target)
print(f"  E_nonlocal(eps_target=0.1194) = {E_nl_target:.6f}")
print(f"  E_nonlocal(0)                 = {E_nonlocal_0:.6f}")
print(f"  Delta E_nonlocal              = {E_nl_target - E_nonlocal_0:+.6f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()
print(f"  ARC LENGTH MINIMUM (thin-wire Biot-Savart limit):")
print(f"    eps_min(L)   = {eps_min_L:.10f}  [from golden section, N=20000]")
print(f"    3/(8*pi)     = {eps_formula:.10f}  [= (p+q)/(4*R2)]")
print(f"    eps_target   = {eps_target:.10f}  [from wave_path_test.py]")
print(f"    Diff(eps_min vs 3/(8*pi)): {(eps_min_L-eps_formula)/eps_formula*100:+.5f}%")
print()
print(f"  FULL BIOT-SAVART (N={N_BS}):")
print(f"    E''_BS(0) = {d2E_BS_0:.4f}  ({'UNSTABLE' if d2E_BS_0 < 0 else 'STABLE'})")
print(f"    Minimum scan: eps ≈ {min_eps_bs:.4f}")
print()
print(f"  INTERPRETATION:")
if abs(eps_min_L - eps_formula)/eps_formula < 0.002:
    print(f"  ROUTE A CONFIRMS GAP 1:")
    print(f"  The arc length minimum of the k=2 perturbed (1,2) torus knot on the")
    print(f"  Hopf torus is at eps = 3/(8*pi) = (p+q)/(4*R2) to high precision.")
    print(f"  This is a GEOMETRIC RESULT: no free parameters, no Biot-Savart needed.")
    print(f"  The formula epsilon = (p+q)/(4*R2) is determined by Hopf geometry + knot.")
    print()
    print(f"  ANALYTIC PROOF NEEDED:")
    print(f"  Show analytically that dL/deps = 0 at eps = (p+q)/(4*R2) for any (p,q)")
    print(f"  torus knot on the Hopf torus with R2 = 2*pi.")
    print(f"  This is a 1D integral identity — tractable via residue calculus or")
    print(f"  by reducing to a standard elliptic integral.")
else:
    print(f"  Route A: arc length minimum at {eps_min_L:.5f}, NOT 3/(8*pi) = {eps_formula:.5f}.")
    print(f"  Non-local Biot-Savart correction is needed to explain the difference.")
    print(f"  Route B minimum is at eps ≈ {min_eps_bs:.4f}.")
print()
print(f"  GAP 2 (Rs from self-energy scale): see Part E above.")
print(f"  E_nonlocal(0) / (L0 * Rs) = {ratio_L_Rs:.4f}")
print(f"  If this ratio has a clean geometric value, Rs is proven from (p,q) + Hopf.")
print()
print(f"  See also: analysis/alpha/hopf_stability.py  (instability proof)")
print(f"            analysis/alpha/epsilon_search.py   (3/(8*pi) candidate)")
print(f"            analysis/alpha/wave_path_test.py   (n_EM model)")
print(SEP)
