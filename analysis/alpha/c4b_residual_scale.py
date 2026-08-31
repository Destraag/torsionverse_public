"""
c4b_residual_scale.py — Can any cross-scale constant explain the n_exact residual?

Context:
  C4b quadratic: n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0
  The physical root with n=2 gives alpha_C4b, error -0.000560% from CODATA.

  The exact n that recovers alpha_CODATA is:
    n_exact = (4*pi^2/phi * alpha_CODATA - Rs) / alpha_CODATA^2
            = 2.018689591...

  So n=2 is off by a residual:
    residual = n_exact - 2 = 0.018689591...

  Question: does this residual correspond to any constant appearing at
  the other scales in the cross-scale ratio table, or to any simple
  combination of {alpha, phi, pi, Rs, R3, R5}?

  If a scale constant matches, the double-spin hypothesis gains support:
  the (1,2) torus knot is the leading term, and the residual winding
  correction (the 0.0187 beyond integer 2) comes from a higher-scale
  coupling. If no scale constant matches, n=2 is simply the best integer
  and the residual is a higher-order topological correction to be derived.

Run: python analysis/c4b_residual_scale.py
"""

import math

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)

alpha       = 7.2973525693e-3
me_MeV      = 0.51099895
mp_MeV      = 938.27208816
r_proton_fm = 0.8414
kappa_GeV_per_fm = 0.9
nuclear_binding_MeV = 8.0
a0_m_s2     = 1.2e-10
c_m_s       = 2.99792458e8
H0_km_s_Mpc = 70.0
Mpc_in_m    = 3.085677581e22
H0_s        = H0_km_s_Mpc * 1e3 / Mpc_in_m
cH0         = c_m_s * H0_s

SEP  = "=" * 65
SEP2 = "-" * 65

# ─────────────────────────────────────────────────────────────────────────────
# Compute n_exact and residual
# ─────────────────────────────────────────────────────────────────────────────

n_exact  = (4 * pi**2 / phi * alpha - Rs) / alpha**2
residual = n_exact - 2       # = 0.018689591...

print(SEP)
print("C4B RESIDUAL SCALE ANALYSIS")
print("What fills the gap between n=2 and n_exact=2.0187?")
print(SEP)
print()
print(f"  n_exact  = {n_exact:.10f}")
print(f"  n=2      = 2.0000000000")
print(f"  residual = {residual:.10f}")
print(f"  (This is what a higher-order correction to C4b must equal.)")
print()

# ─────────────────────────────────────────────────────────────────────────────
# The five cross-scale ratios from scale_check.py
# ─────────────────────────────────────────────────────────────────────────────

R1 = alpha / (2 * pi)                                        # particle/EM
R2 = alpha                                                    # full EM
R3 = nuclear_binding_MeV / mp_MeV                            # nuclear binding
R4 = (kappa_GeV_per_fm * r_proton_fm) / (mp_MeV / 1000.0)   # hadronic
R5 = a0_m_s2 / cH0                                           # galactic

print(SEP)
print("PART 1 — CROSS-SCALE CONSTANTS AS CANDIDATES")
print(SEP)
print()
print(f"  {'Expression':<35} {'Value':>12}  {'Ratio to residual':>18}  {'% off':>8}")
print(f"  {'-'*35} {'-'*12}  {'-'*18}  {'-'*8}")

cross_scale = [
    ("R1 = alpha/(2*pi)",         R1),
    ("R2 = alpha",                R2),
    ("R3 = nuclear binding/mp",   R3),
    ("R4 = hadronic (kappa*rp/mp)",R4),
    ("R5 = galactic (a0/cH0)",    R5),
    ("R1 * 2*pi  [= alpha]",      R1 * 2 * pi),
    ("R3 / alpha",                R3 / alpha),
    ("R5 * 2*pi",                 R5 * 2 * pi),
    ("R5 / R1",                   R5 / R1),
    ("R4 / 2*pi",                 R4 / (2 * pi)),
]

for name, val in cross_scale:
    ratio = val / residual
    pct   = (val - residual) / residual * 100
    flag  = "  <-- NEAR" if abs(pct) < 10 else ""
    print(f"  {name:<35} {val:>12.6f}  {ratio:>18.6f}  {pct:>+8.2f}%{flag}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# Simple combinations of {alpha, phi, pi, Rs}
# ─────────────────────────────────────────────────────────────────────────────

print(SEP)
print("PART 2 — SIMPLE GEOMETRIC COMBINATIONS")
print(SEP)
print()
print(f"  {'Expression':<35} {'Value':>12}  {'% off from residual':>20}")
print(f"  {'-'*35} {'-'*12}  {'-'*20}")

geometric = [
    ("2*pi*alpha",                2 * pi * alpha),
    ("pi*alpha",                  pi * alpha),
    ("phi*alpha",                 phi * alpha),
    ("2*phi*alpha",               2 * phi * alpha),
    ("4*alpha",                   4 * alpha),
    ("alpha/Rs",                  alpha / Rs),
    ("Rs*alpha",                  Rs * alpha),
    ("Rs*alpha*phi",              Rs * alpha * phi),
    ("2*Rs",                      2 * Rs),
    ("Rs/pi^2",                   Rs / pi**2),
    ("Rs/pi",                     Rs / pi),
    ("sqrt5/pi^2",                sqrt5 / pi**2),
    ("alpha^(1/2)",               alpha**0.5),
    ("alpha^(1/2)/pi",            alpha**0.5 / pi),
    ("alpha^(1/2)/(2*pi)",        alpha**0.5 / (2 * pi)),
    ("1/(4*pi^2/phi)",            phi / (4 * pi**2)),    # = 1/Q_bare
    ("alpha * 1/phi^2",           alpha / phi**2),
    ("alpha * phi^2",             alpha * phi**2),
    ("(alpha/pi)^(1/2)",          (alpha / pi)**0.5),
    ("Rs^2/alpha",                Rs**2 / alpha),
    ("2*alpha^2/Rs",              2 * alpha**2 / Rs),
]

near_hits = []
for name, val in geometric:
    if val <= 0:
        continue
    pct = (val - residual) / residual * 100
    flag = "  <-- NEAR" if abs(pct) < 10 else ""
    if abs(pct) < 5:
        near_hits.append((name, val, pct))
    print(f"  {name:<35} {val:>12.6f}  {pct:>+20.3f}%{flag}")
print()

# ─────────────────────────────────────────────────────────────────────────────
# Near-hit summary
# ─────────────────────────────────────────────────────────────────────────────

print(SEP)
print("PART 3 — NEAR HIT SUMMARY")
print(SEP)
print()
print(f"  residual = n_exact - 2 = {residual:.8f}")
print()

if near_hits:
    print(f"  The following expressions are within 5% of the residual:")
    print()
    for name, val, pct in near_hits:
        print(f"    {name:<35} = {val:.8f}  ({pct:+.3f}%)")
    print()
else:
    print(f"  No expression within 5% of the residual found.")
    print(f"  Best candidates:")
    all_vals = cross_scale + geometric
    sorted_vals = sorted(
        [(n, v) for n, v in all_vals if v > 0],
        key=lambda x: abs(x[1] - residual) / residual
    )
    for name, val in sorted_vals[:6]:
        pct = (val - residual) / residual * 100
        print(f"    {name:<35} = {val:.8f}  ({pct:+.3f}%)")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Interpretation
# ─────────────────────────────────────────────────────────────────────────────

print(SEP)
print("PART 4 — INTERPRETATION")
print(SEP)
print()
print("  The residual 0.01869 is what separates n=2 from n_exact.")
print("  If a scale constant matches it, the interpretation is:")
print("    n_effective = 2 + [scale constant]")
print("    meaning the (1,2) torus knot leading term is corrected by a")
print("    cross-scale coupling -- the secondary winding is slightly more")
print("    than integer 2 due to a physical mixing with the matched scale.")
print()
print("  If no scale constant matches (no near hit found), the interpretation is:")
print("    n=2 is the leading integer winding number.")
print("    The residual 0.01869 is a higher-order geometric correction")
print("    (not from another physical scale) that a full linking integral")
print("    calculation (Step D2/Step 3) would produce analytically.")
print("    This is the more conservative and likely scenario.")
print()

# Is the residual itself expressible as n*alpha for small integer n?
print("  Residual / alpha:")
print(f"    residual / alpha = {residual / alpha:.6f}")
print(f"    Nearest integer: {round(residual/alpha)}")
print(f"    (An integer here would mean residual = k*alpha for integer k)")
print()
print("  Residual / Rs:")
print(f"    residual / Rs = {residual / Rs:.6f}")
print(f"    Nearest simple fraction: {round(residual/Rs*10)/10:.1f}")
print()
print("  Residual in 1/alpha terms (i.e., what it means for 1/alpha prediction):")
inv_alpha_C4b = 1 / (7.2973117300057e-3)
inv_alpha_CODATA = 1 / alpha
print(f"    1/alpha_C4b     = {inv_alpha_C4b:.8f}")
print(f"    1/alpha_CODATA  = {inv_alpha_CODATA:.8f}")
print(f"    gap in 1/alpha  = {inv_alpha_C4b - inv_alpha_CODATA:+.8f}")
print(f"    This gap = {(inv_alpha_C4b - inv_alpha_CODATA):.6f} in 1/alpha units.")
print(f"    Note: 1/(2*pi) = {1/(2*pi):.6f}  (Hopf suppression scale)")
print(f"    Ratio: gap / (1/2pi) = {(inv_alpha_C4b - inv_alpha_CODATA) * 2 * pi:.6f}")
print()
print(SEP)
print("CONCLUSION")
print(SEP)
print()
if near_hits:
    print("  At least one scale constant is within 5% of the residual.")
    print("  See near-hit summary above for the candidate(s).")
    print("  This DOES NOT prove a cross-scale coupling -- it identifies")
    print("  a candidate worth testing by a forward derivation.")
else:
    print("  No scale constant within 5% found.")
    print("  The residual 0.01869 does not obviously come from any known")
    print("  physical scale. It is consistent with a higher-order topological")
    print("  correction intrinsic to the (1,2) torus knot geometry.")
    print()
    print("  Bottom line for the double-spin hypothesis:")
    print("    n=2 (integer) is correct as the LEADING winding number.")
    print("    The residual correction should come from Step D2 (linking")
    print("    integral calculation), not from another physical scale.")
    print("    Deep-diving the double-spin geometry is worth doing on its")
    print("    own terms, independent of the cross-scale comparison.")
print(SEP)
