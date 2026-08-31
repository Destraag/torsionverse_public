"""
constants.py — Single source of truth for all locked constants in the
torsionverse / alpha derivation project.

USAGE IN ANY SCRIPT:
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
    from analysis.alpha.constants import *

    OR from the project root:
    from analysis.alpha.constants import *

PURPOSE:
  All scripts previously defined constants inline, making it easy to use
  stale or approximated values by accident. This module contains every
  locked value with its source script and the date it was confirmed.

  RULE: If a constant appears here it is LOCKED — do not redefine it in a
  script unless deliberately testing a variation. If you update a constant,
  update the version date and source below.

  RULE: When computing a derived quantity in a script, import it from here
  rather than retyping a number. If the derived quantity is not here, add it.

Last updated: 2026-08-18 (session 5)
Confirmed by: gap1_cgeo_analytic.py, gap1_triangle_contact.py,
              gap1_lagrangian_elastic.py, gap1_anisotropic_compliance.py,
              gap1_vertex_bending_modulus.py, gap1_boussinesq_vertex.py
"""

import math

# ── MATHEMATICAL CONSTANTS ────────────────────────────────────────────────────
pi    = math.pi
sqrt2 = math.sqrt(2)
sqrt3 = math.sqrt(3)
sqrt5 = math.sqrt(5)
PHI   = (1 + sqrt5) / 2          # golden ratio = 1.6180339887...


# ── CODATA 2018 FUNDAMENTAL CONSTANTS ────────────────────────────────────────
# Source: NIST CODATA 2018
alpha = 7.2973525693e-3           # fine structure constant (CODATA 2018)


# ── TORUS GEOMETRY ────────────────────────────────────────────────────────────
# Source: alpha_theory.txt, hopf_linking_integral.py
R2    = 2 * pi                    # Hopf torus major radius (model units)
R1    = 1.0                       # Hopf torus minor radius (model units)

# (p,q) = (1,2) torus knot — these emerge from the winding; do NOT assume them
p_wind = 1
q_wind = 2


# ── WAVE GEOMETRY ─────────────────────────────────────────────────────────────
# Source: gap1_cgeo_analytic.py (confirmed), alpha_theory.txt
eps_L5   = 3 / (8 * pi)          # = 0.11936620731892  topological lock point
gamma_w  = 2 * eps_L5             # = 3/(4*pi)  wave shear strain amplitude


# ── GRAIN GEOMETRY ────────────────────────────────────────────────────────────
# Source: gap1_gj5_phi_identity.py, gap1_cgeo_analytic.py
gj5       = 1 / (2 * PHI**2)     # = 0.19098300562505  pentagon jamming threshold
                                   # = 1 - cos(pi/5) = a_W/2  (ALL equivalent, proven)
gamma_c   = gj5 / gamma_w         # = 0.79998774324288  normalised jamming threshold
                                   # NOTE: NOT 4/5 exactly; irrational

# Grain edge length = alpha (model units — grain size scales with EM coupling)
L_grain   = alpha                  # = 7.2973525693e-3

# Icosahedral vertex geometry
Omega_ico = pi / 3                 # angular deficit at icosahedral vertex (= 2pi - 5*pi/3)
alpha_c   = math.asin(5 / 6)      # cone half-opening angle at vertex = 56.4427 deg


# ── WAVE COUPLING CONSTANTS ───────────────────────────────────────────────────
# Source: gap1_cgeo_analytic.py (Richardson-confirmed, locked to 14 sig figs)
dn   = 0.16857744391041            # d(n_EM)/d(eps) at eps_L5
d2n  = 0.19763679211711            # d^2(n_EM)/d(eps)^2 at eps_L5
I_el = 0.15614610339308            # elastic integral (exact closed form)
C_geo = 10.33418281379304          # geometric coupling (alpha-independent)

# Cross-check identities (should all be ~equal, verified):
# C_geo * I_el * tan(pi/5) == d2n/dn == 1.17237981...
# C_geo = d2n / (dn * I_el * tan(pi/5))


# ── HOPF/C4b STRUCTURE CONSTANTS ─────────────────────────────────────────────
# Source: hopf_linking_integral.py, writhe_min.py
Rs  = sqrt5 / (4 * pi)            # = 0.17794064290479  Hopf coupling constant
Q   = 4 * pi**2 / PHI             # = 24.39900390158560  Hopf suppression factor
n_EM   = 2.01868734358082         # EM-weighted winding number at eps_L5 (Richardson-confirmed)
n_exact = 2.01868959103706        # n required to give alpha_CODATA from C4b


# ── GAP 1 ─────────────────────────────────────────────────────────────────────
# Source: gap1_richardson.py (confirmed NOT numerical artifact)
delta_n   = 2.24745624e-6         # n_exact - n_EM  (Gap 1 winding number gap)
delta_eps = 1.33318918e-5         # eps_num - eps_L5 (Gap 1 in epsilon units)
f_frac    = 1.11688995e-4         # delta_eps / eps_L5 (fractional gap)

# Canonical Gap 1 identity (0.0084% accurate):
# delta_n = K * alpha^2 * d2n   where K = 4*eps_L5/sqrt5
K_gap = 4 * eps_L5 / sqrt5       # = 3/(2*pi*sqrt5)  canonical prefactor


# ── CONTACT STIFFNESS k_n ─────────────────────────────────────────────────────
# IMPORTANT: TWO VALUES EXIST — use the locked identity, not the HCP fit.
#
# k_n_LOCKED: from the exact identity k_n = 2*C_geo*delta_n*alpha^3
#   Source: gap1_exact_prefactor.py, confirmed gap1_cgeo_analytic.py
#   This is the CORRECT value to use in any first-principles comparison.
#
# k_n_HCP: from the phenomenological HCP breathing-room fit (sqrt(3)-alpha)/2 * alpha^5
#   Source: gap1_triangle_contact.py (0.006% fit)
#   This is a PHENOMENOLOGICAL approximation — DO NOT use in exact derivations.
#
# WHY THEY DIFFER: k_n_LOCKED = 1.80507e-11; k_n_HCP = 1.784e-11 (1.18% off)
#   The 1.18% difference propagates into any formula using k_n_HCP.
#   Always use k_n_LOCKED for new calculations.

k_n_LOCKED = 2 * C_geo * delta_n * alpha**3   # = 1.80506718e-11  (AUTHORITATIVE)
k_n_HCP    = 1.784e-11                         # STALE — phenomenological only

# Alias for convenience (always the locked value):
k_n = k_n_LOCKED


# ── ANISOTROPIC COMPLIANCE (session 5, 2026-08-18) ───────────────────────────
# Source: gap1_anisotropic_compliance.py (commit 673dd5d)
# The mechanism that closes Gap 1: jammed-arc vertex compliance s_jam ~ 0.283*alpha
frac_d2n_jam  = -0.0406020036     # jammed-arc fraction of d2n (N=2M, gap1_frac_d2n_precision.py, ±1.6e-5)
frac_dn_jam   = -0.1547287739     # jammed-arc fraction of dn  (N=2M, gap1_frac_d2n_precision.py, ±9.5e-7)
s_jam         = 2.0647518466e-3   # required anisotropic compliance (= 0.283*alpha)
k_vertex_needed = k_n / s_jam     # = 8.640e-9  required vertex spring constant


# ── MEDIUM ELASTIC CONSTANTS ──────────────────────────────────────────────────
# Source: alpha_theory.txt (torsion medium definition), gap1_boussinesq_vertex.py
# The torsion medium has v_s = Rs*c, which defines G = Rs^2 in model units (c=1).
# HOWEVER: for grain contact mechanics, model units set G_medium = 1 (shear modulus).
# All grain elastic calculations use these values.
G_medium = 1.0                     # shear modulus (model units)
nu_medium = (1 - 2*Rs**2) / (2*(1 - Rs**2))  # exact: (1-2Rs²)/(2(1-Rs²)) = 0.483651...
                                   # Source: whitepaper.txt; Rs = sqrt5/(4*pi)
                                   # Simplification: 1-nu = 1/(2*(1-Rs²))
                                   # Previously stored as 0.4837 (4 sig figs; diff = 0.010%)
E_medium = 2 * G_medium * (1 + nu_medium)   # Young's modulus
E_r_medium = E_medium / (1 - nu_medium**2)  # = 3.7148  reduced modulus

# Plane-strain correction factor for cone-bending formula:
# k_cone_3D = k_cone_2D * 1/(1-nu)
plane_strain_factor = 1 / (1 - nu_medium)   # = 1.93686


# ── LOBKOVSKY-WITTEN CONE VERTEX (session 5, 2026-08-18) ─────────────────────
# Source: gap1_vertex_bending_modulus.py (commit a394853),
#         gap1_boussinesq_vertex.py (commit 62cd984)
# The icosahedral grain vertex bends like a thin elastic cone.
# Full formula (with correct E and plane-strain):
#   k_vertex = E_medium * L^3 * Omega^2 / (12*(1-nu^2)) / (24*pi*log(a/r_v)) * 1/(1-nu)
# Numerical result:
k_vertex_LW  = (E_medium * L_grain**3 * Omega_ico**2
               / (12 * (1 - nu_medium**2))
               / (24 * pi * math.log((L_grain*sqrt3/2) / (L_grain/sqrt3)))
               * plane_strain_factor)    # = 8.7137e-9

# ── COOPERATIVE VERTEX STIFFENING (session 5, 2026-08-19) ────────────────────
# Source: gap1_phi_coupling.py (committed this session)
#
# The isolated-grain LW formula misses two cooperative back-reaction paths
# that activate when the icosahedral shell goes taut at eps = eps_L5:
#
#   PATH 1 — 5 direct edges (A→B_i):
#     k₁ = (5-√5)/2 · k_n  = (√5/PHI) · k_n
#     Geometry: cos²(edge-axis angle) = 1/(√5·PHI)  EXACT
#     Derivation: δL_i = u·cos(α_c); F_z = k_n·u·cos²(α_c); sum 5 edges
#
#   PATH 2 — 5 face midpoints (A→M_i, M_i = midpt of far edge B_iB_{i+1}):
#     k₂ = (√5-2) · k_n  = (1/PHI³) · k_n
#     Geometry: cos²(midpt-axis angle) = (√5-2)/5  EXACT
#     Derivation: face acts rigid (k_LW >> k_n); midpoint pulled toward A
#
#   COMBINED (exact algebraic identity PHI = √5/PHI + 1/PHI³):
#     k_vertex_eff = k_vertex_LW + PHI · k_n
#
#   Residual vs k_vertex_needed: -0.26%  (0.36 sigma from frac_d2n_jam precision)
#
# NOTE: sin(alpha_c) = 5/6 stored in constants IS THE CONE METRIC FACTOR
#   (sector angle 5π/3 out of 2π), NOT the 3D edge-axis angle.
#   The true 3D edge-axis cos²= (1-1/√5)/2 = 1/(√5·PHI).

k_vertex_eff = k_vertex_LW + PHI * k_n_LOCKED   # = 8.7428e-9  closes Gap 1 to 0.26%


# ── TOPOLOGICAL CONSTANTS ─────────────────────────────────────────────────────
N_lock  = 2 * pi / (alpha * PHI)  # = 532.14  grains per tube circumference
a_W     = 1 / PHI**2              # = 0.38197  icosahedral acceptance window edge
                                   # gj5 = a_W/2  (proven)


# ── CONVENIENCE VERIFICATION ─────────────────────────────────────────────────
def verify_constants():
    """Run internal consistency checks. Call this from any script to confirm
    the constants module loaded correctly."""
    tol = 1e-8
    checks = [
        ("gj5 = 1/(2*PHI^2)",      abs(gj5 - 1/(2*PHI**2)) < tol),
        ("gj5 = a_W/2",            abs(gj5 - a_W/2) < tol),
        ("gj5 = 1-cos(pi/5)",      abs(gj5 - (1-math.cos(pi/5))) < 1e-12),
        ("eps_L5 = 3/(8*pi)",      abs(eps_L5 - 3/(8*pi)) < tol),
        ("gamma_c = gj5/gamma_w",  abs(gamma_c - gj5/gamma_w) < tol),
        ("C_geo = d2n/(dn*I_el*tan(pi/5))",
                                   abs(C_geo - d2n/(dn*I_el*math.tan(pi/5))) < 1e-6),
        ("k_n = 2*C_geo*delta_n*a3",
                                   abs(k_n - 2*C_geo*delta_n*alpha**3) < 1e-20),
        ("delta_n = K*a2*d2n ±0.01%",
                                   abs(delta_n/(K_gap*alpha**2*d2n) - 1) < 1e-4),
        ("Rs = sqrt5/(4*pi)",      abs(Rs - sqrt5/(4*pi)) < tol),
        ("Q = 4*pi^2/PHI",         abs(Q - 4*pi**2/PHI) < tol),
    ]
    all_ok = True
    for name, result in checks:
        status = "OK" if result else "FAIL"
        if not result:
            all_ok = False
        print(f"  [{status}] {name}")
    return all_ok


if __name__ == "__main__":
    print("constants.py — locked values for torsionverse alpha derivation")
    print("=" * 60)
    print()
    print(f"  alpha          = {alpha:.13e}")
    print(f"  eps_L5         = {eps_L5:.14f}")
    print(f"  gj5            = {gj5:.14f}")
    print(f"  gamma_c        = {gamma_c:.14f}")
    print(f"  dn             = {dn:.14f}")
    print(f"  d2n            = {d2n:.14f}")
    print(f"  I_el           = {I_el:.14f}")
    print(f"  C_geo          = {C_geo:.14f}")
    print(f"  delta_n        = {delta_n:.14e}")
    print(f"  k_n (LOCKED)   = {k_n:.10e}")
    print(f"  k_n (HCP stale)= {k_n_HCP:.10e}   *** DO NOT USE FOR EXACT WORK ***")
    print(f"  s_jam          = {s_jam:.10e}")
    print(f"  k_vertex_needed= {k_vertex_needed:.10e}")
    print(f"  k_vertex_LW    = {k_vertex_LW:.10e}")
    print(f"  E_medium       = {E_medium:.8f}")
    print(f"  nu_medium      = {nu_medium:.4f}")
    print(f"  plane_strain   = {plane_strain_factor:.8f}")
    print()
    print("Internal consistency checks:")
    ok = verify_constants()
    print()
    if ok:
        print("  All checks passed.")
    else:
        print("  WARNING: one or more checks FAILED.")
