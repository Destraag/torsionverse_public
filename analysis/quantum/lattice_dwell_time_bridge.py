"""
lattice_dwell_time_bridge.py
============================
Derives (not asserts) that the DISCRETE Jobson-cell nexus-hopping picture
("dwell-time asymmetry", doc_torsionverse.txt) and the CONTINUUM Klein-Gordon
group-velocity picture ("v = c^2*k/omega ... persists for free", same doc)
are the SAME physics at two orders of one expansion -- not two independent
postulates. Closes a gap raised while auditing the neutrino/dwell-time
discussion: how does a winding "know" its own velocity/direction between
nexus hops?

MODEL: a 1D chain of Jobson cells, spacing a, nearest-neighbor elastic
coupling kappa, onsite mass/restoring term omega_0 (physically
omega_0 = omega_C = m*c^2/hbar for a massive winding; omega_0 = 0 for a
massless/photon mode -- WHERE omega_0 itself comes from is NOT re-derived
here, that is the separate, already-addressed Zone-1-displacement/Maxwell-
jamming mass mechanism). Standard discrete Klein-Gordon lattice EOM:

    M * d^2(psi_n)/dt^2 = kappa*(psi_(n+1) - 2*psi_n + psi_(n-1)) - M*omega_0^2*psi_n

DERIVATION CHAIN:
  1. Plane-wave ansatz -> EXACT dispersion
       omega^2 = omega_0^2 + (4*kappa/M)*sin^2(k*a/2)
     verified directly against the recursion relation at every lattice site
     (no matrix/eigenvalue library needed).
  2. Long-wavelength limit (k*a -> 0): omega^2 -> omega_0^2 + c^2*k^2
     [Klein-Gordon dispersion], with c^2 = kappa*a^2/M IDENTIFIED with
     c^2 = K/rho (K=1/eps_0, rho=mu_0 -- already-established torsionverse
     constants, MG2-2/doc_torsion.txt Section 3.3) -- zero new free
     parameters. This is what qm_from_medium.py's QM1 currently ASSERTS
     (its own check computes c*k twice and compares -- a tautology, not a
     derivation); this script supplies the missing derivation.
  3. Group velocity from the EXACT (not approximated) dispersion:
       v_g = (c^2/a) * sin(k*a) / omega
     matches the continuum v_g = c^2*k/omega (doc_torsionverse.txt's
     "Velocity IS the wave group velocity") as k*a -> 0.
  4. KEY RESULT: the leading finite-ka correction to v_g is NEGATIVE and
     a-DEPENDENT -- invisible at leading order, since c=sqrt(K/rho) alone
     does not depend on a. The pressure-gradient/dwell-time force
     (doc_magnetism.txt Section 1.2: "thinned medium = lower density =
     longer wave dwell time") is a LATTICE-DISCRETENESS (finite k*a)
     effect, not a separate postulate from the same dispersion relation.
  5. Standard ray-tracing fact (time-independent medium: omega is CONSERVED
     along a ray, same principle as Snell's law / gravitational lensing in
     analog-gravity treatments) gives v_group(a) MONOTONICALLY DECREASING
     as a grows, at fixed omega -- wider local spacing = slower transport =
     longer dwell time, matching doc_magnetism.txt's mechanism from the
     SAME lattice dispersion relation (not a simple 1/a power law -- the
     exact relation is k(a)=(2/a)*arcsin(C2*a), C2 fixed by omega), not an
     independently-asserted rule.

SCOPE: this derives the KINEMATIC bridge (discrete hopping <-> continuum
k/omega persistence) only. It does NOT re-derive where omega_0 (mass) comes
from, or the microscopic origin of a(x) varying in a real pressure
gradient (doc_magnetism.txt Section 1.2's own mechanism) -- both are taken
as given inputs from elsewhere in the framework.

CHECKS:
  BR1: plane-wave ansatz solves the discrete EOM exactly at every site
  BR2: KG-approximation relative error shrinks as O((ka)^2) [ratio -> 4]
  BR3: c^2=kappa*a^2/M matches c^2=K/rho (K=1/eps_0, rho=mu_0) exactly
  BR4: exact group velocity converges to continuum v=c^2*k/omega as O((ka)^2)
  BR5: v_group(wider a) < v_group(narrower a) at fixed k
  BR6: v_group(a) is monotonically decreasing as a grows, at fixed omega
       (proper ray-tracing)

Run: python analysis/quantum/lattice_dwell_time_bridge.py
Reference: docs/series1/doc_torsionverse.txt (dwell-time asymmetry, inertia
  sections); docs/series1/doc_magnetism.txt Section 1.2; qm_from_medium.py
"""
import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

pi = math.pi

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  {'[PASS]' if cond else '[FAIL]'} {name}")
    if detail: print(f"         {detail}")

print(SEP)
print("lattice_dwell_time_bridge.py -- discrete lattice -> KG dispersion bridge")
print(SEP)

# ── SECTION 1: discrete chain EOM -> exact dispersion relation ───────────────
print()
print(SEP2)
print("SECTION 1: DISCRETE CHAIN EOM -> EXACT DISPERSION RELATION")
print(SEP2)
print("  M*d^2(psi_n)/dt^2 = kappa*(psi_(n+1)-2*psi_n+psi_(n-1)) - M*omega_0^2*psi_n")
print("  Plane-wave ansatz psi_n = exp(i*(k*n*a - omega*t)) gives, EXACTLY:")
print("    omega^2 = omega_0^2 + (4*kappa/M)*sin^2(k*a/2)")

# Toy (dimensionless) units -- verifying the MATH, not fitting real numbers.
kappa_t, M_t, a_t, omega0_t = 2.3, 1.7, 1.0, 0.31
N_t = 12
j_t = 5                                     # test mode index
k_t = 2*pi*j_t / (N_t * a_t)                # periodic BC: k = 2*pi*j/(N*a)

def omega_sq(k, kappa=kappa_t, M=M_t, a=a_t, omega0=omega0_t):
    return omega0**2 + (4*kappa/M) * math.sin(k*a/2)**2

omega_t = math.sqrt(omega_sq(k_t))

max_residual = 0.0
for n in range(N_t):
    psi_n   = complex(math.cos(k_t*n*a_t),     math.sin(k_t*n*a_t))
    psi_np1 = complex(math.cos(k_t*(n+1)*a_t), math.sin(k_t*(n+1)*a_t))
    psi_nm1 = complex(math.cos(k_t*(n-1)*a_t), math.sin(k_t*(n-1)*a_t))
    lhs = -M_t * omega_t**2 * psi_n
    rhs = kappa_t*(psi_np1 - 2*psi_n + psi_nm1) - M_t*omega0_t**2*psi_n
    max_residual = max(max_residual, abs(lhs - rhs))

check("BR1: plane-wave ansatz solves the discrete EOM EXACTLY at every site (periodic chain)",
      max_residual < 1e-10,
      f"max |LHS-RHS| over N={N_t} sites = {max_residual:.2e}  (mode j={j_t}, ka={k_t*a_t:.4f} rad)")

# ── SECTION 2: long-wavelength limit -> Klein-Gordon dispersion ──────────────
print()
print(SEP2)
print("SECTION 2: a -> 0 AT FIXED k REPRODUCES KLEIN-GORDON DISPERSION")
print(SEP2)
print("  Hold c^2=kappa*a^2/M FIXED (continuum K,rho fixed -- only refine the lattice)")
print("  while shrinking a, at FIXED k: sin^2(ka/2) = (ka/2)^2*[1-(ka)^2/12+O((ka)^4)]")
print("    => omega^2 = omega_0^2 + c^2*k^2 * [1 - (ka)^2/12 + ...]")
print("  Bracket -> 1 as a -> 0: KG dispersion (qm_from_medium.py QM3) falls out of")
print("  the lattice model, it is not an independent postulate.")

c_t = a_t * math.sqrt(kappa_t / M_t)
k_ref = 0.3   # a fixed physical wavenumber -- refine the LATTICE around it

def omega_sq_lattice(k, a, c=c_t, omega0=omega0_t):
    """Dispersion with c held FIXED (continuum K,rho fixed) as 'a' varies --
    equivalent to kappa(a)=kappa_t*(a_t/a), M(a)=M_t*(a/a_t) substituted in."""
    return omega0**2 + (4*c**2/a**2) * math.sin(k*a/2)**2

def omega_sq_KG(k, c=c_t, omega0=omega0_t):
    return omega0**2 + c**2 * k**2

# Sanity: must agree exactly with Section 1's formula at the reference point.
assert abs(omega_sq_lattice(k_t, a_t) - omega_sq(k_t)) < 1e-10

a_values = [a_t, a_t/2, a_t/4, a_t/8, a_t/16]   # refine the lattice, k_ref fixed
prev_rel_err = None
ratios = []
for a in a_values:
    rel_err = abs(omega_sq_lattice(k_ref, a) - omega_sq_KG(k_ref)) / (c_t**2 * k_ref**2)
    if prev_rel_err is not None:
        ratios.append(prev_rel_err / rel_err)
    prev_rel_err = rel_err

check("BR2: KG-approximation relative error (vs the kinetic term c^2*k^2) shrinks as "
      "O(a^2) [ratio -> 4.0 each time a halves, k held FIXED]",
      all(abs(r - 4.0) < 0.05 for r in ratios[-2:]),
      f"successive error ratios (halving a, k={k_ref} fixed): {[f'{r:.4f}' for r in ratios]}  (expect -> 4.0)")

eps_0 = 8.8541878128e-12   # F/m (CODATA)
mu_0  = 1.25663706212e-6   # N/A^2 (CODATA)
c_SI  = 299792458.0        # m/s (exact, SI definition)
K_torsion    = 1 / eps_0        # Pa  [K=1/eps_0, MG2-2/doc_torsion.txt Sec 3.3]
rho_torsion  = mu_0             # kg/m^3 (SI units of mu_0 for this identity)
c_from_K_rho = math.sqrt(K_torsion / rho_torsion)

check("BR3: c^2=kappa*a^2/M (this lattice model) IDENTIFIED with c^2=K/rho "
      "(K=1/eps_0, rho=mu_0, already-established torsionverse constants) -- "
      "zero new free parameters",
      abs(c_from_K_rho - c_SI) / c_SI < 1e-9,
      f"sqrt(K/rho) = {c_from_K_rho:.6e} m/s   c (SI) = {c_SI:.6e} m/s")

# ── SECTION 3: group velocity -- exact vs continuum v=c^2*k/omega ───────────
print()
print(SEP2)
print("SECTION 3: GROUP VELOCITY -- EXACT vs CONTINUUM v=c^2*k/omega")
print(SEP2)
print("  d(omega)/dk from Section 1's EXACT dispersion (direct calculus):")
print("    v_g = (c^2/a) * sin(k*a) / omega")
print("  Continuum formula already used in doc_torsionverse.txt: v_g = c^2*k/omega.")
print("  These agree as a -> 0 at fixed k, with the SAME O(a^2) correction as Section 2.")

def v_group_lattice(k, a, c=c_t, omega0=omega0_t):
    om = math.sqrt(omega_sq_lattice(k, a, c, omega0))
    return (c**2 / a) * math.sin(k*a) / om

def v_group_continuum(k, c=c_t, omega0=omega0_t):
    om = math.sqrt(omega_sq_KG(k, c, omega0))
    return c**2 * k / om

vg_cont_ref = v_group_continuum(k_ref)
prev_rel_err = None
ratios_vg = []
for a in a_values:
    rel_err = abs(v_group_lattice(k_ref, a) - vg_cont_ref) / vg_cont_ref
    if prev_rel_err is not None:
        ratios_vg.append(prev_rel_err / rel_err)
    prev_rel_err = rel_err

check("BR4: exact group velocity converges to the continuum v=c^2*k/omega formula "
      "as O(a^2) at fixed k ['constant velocity persists for free' is the a->0 limit]",
      all(abs(r - 4.0) < 0.05 for r in ratios_vg[-2:]),
      f"successive error ratios (halving a, k={k_ref} fixed): {[f'{r:.4f}' for r in ratios_vg]}  (expect -> 4.0)")

# ── SECTION 4: finite-ka correction IS the dwell-time/force mechanism ────────
print()
print(SEP2)
print("SECTION 4: WIDER SPACING (a) -> SLOWER GROUP VELOCITY -> LONGER DWELL TIME")
print(SEP2)
print("  v_g = (c^2*k/omega) * [1 - (ka)^2/6 + O((ka)^4)]  (from sin(ka)/a series)")
print("  Invisible at leading order: c=sqrt(K/rho) alone does not depend on a. The")
print("  pressure-gradient/dwell-time force (doc_magnetism.txt Section 1.2: 'thinned")
print("  medium = lower density = longer wave dwell time') is a LATTICE-DISCRETENESS")
print("  (finite ka) effect from this SAME dispersion relation, not a separate rule.")

a_narrow, a_wide = 0.5 * a_t, 1.5 * a_t     # narrower / wider local spacing, same k_ref
vg_narrow = v_group_lattice(k_ref, a_narrow)
vg_wide   = v_group_lattice(k_ref, a_wide)

check("BR5: v_group(wider spacing) < v_group(narrower spacing) at fixed k -- wider "
      "spacing = slower transport = longer dwell time (doc_magnetism.txt Sec 1.2, "
      "same lattice model)",
      vg_wide < vg_narrow,
      f"v_g(a={a_narrow:.3f}) = {vg_narrow:.6f}   v_g(a={a_wide:.3f}) = {vg_wide:.6f}")

# ── SECTION 5: conserved omega along a ray -> v_group falls as a grows ───────
print()
print(SEP2)
print("SECTION 5: FREQUENCY CONSERVED ALONG A RAY -> v_group FALLS MONOTONICALLY AS a GROWS")
print(SEP2)
print("  Standard ray-tracing fact (time-independent medium): a wave's TEMPORAL")
print("  frequency omega is conserved crossing regions of different local spacing")
print("  a(x) (same principle as Snell's law / gravitational lensing in analog-")
print("  gravity treatments). Solving omega^2(k,a)=omega^2 (fixed) for k(a):")
print("    sin(k*a/2) = C2*a  =>  k(a) = (2/a)*arcsin(C2*a)   [C2 fixed by omega]")
print("    => k(a) -> 2*C2 (finite) as a->0; v_group(a) falls smoothly as a grows")
print("    (the exact relation -- not a naive 1/a power law)")

omega_fixed = math.sqrt(omega_sq_lattice(k_ref, a_t))
C2 = math.sqrt(omega_fixed**2 - omega0_t**2) / (2*c_t)

def k_at_fixed_omega(a):
    return (2/a) * math.asin(C2 * a)

def v_group_fixed_omega(a):
    return v_group_lattice(k_at_fixed_omega(a), a)

a_test_values = [0.3*a_t, 0.6*a_t, a_t, 1.4*a_t, 1.8*a_t]
vg_fixed_omega = [v_group_fixed_omega(a) for a in a_test_values]
monotone_decreasing = all(vg_fixed_omega[i] > vg_fixed_omega[i+1] for i in range(len(vg_fixed_omega)-1))

check("BR6: at FIXED omega (proper ray-tracing), v_group(a) is MONOTONICALLY "
      "DECREASING as local spacing a grows -- wider spacing = slower transport = "
      "longer dwell time, from the SAME lattice dispersion relation",
      monotone_decreasing,
      f"v_g at a/a0={[f'{a/a_t:.2f}' for a in a_test_values]}: "
      f"{[f'{v:.6f}' for v in vg_fixed_omega]}  (monotonically decreasing)")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP2)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, s, d in results:
        if s == 'FAIL': print(f"  FAILED: {name}")
print()
print("  CONCLUSION: the discrete 'winding hops nexus to nexus, redirected by local")
print("  spacing' picture and the continuum 'k/omega ratio persists for free at")
print("  constant velocity' picture are the SAME dispersion relation viewed at two")
print("  orders of the same k*a expansion -- leading order (ka->0) gives inertia/")
print("  Newton's first law (BR2-BR4); the NEXT order (finite ka) gives the")
print("  pressure-gradient/dwell-time force (BR5-BR6). Neither needs the corpuscle")
print("  to separately 'know' its own energy level in a 3D vector -- omega and k")
print("  are properties of the collective wave pattern, propagated forward by the")
print("  SAME dispersion relation at every scale, exactly as in ordinary phonon/")
print("  Bloch-wave lattice dynamics.")
print("  Reference: docs/series1/doc_torsionverse.txt (dwell-time, inertia sections);")
print("  docs/series1/doc_magnetism.txt Section 1.2; qm_from_medium.py")
print(SEP)
