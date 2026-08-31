"""
mechanical_radius.py
====================
Computes the proton mechanical radius from the Zone 3 spinning cell pressure
profile and compares to the measured value from gravitational form factors.

MEASUREMENT:
  Burkert et al. 2018 (Science 361, 207) -- Jefferson Lab CLAS12
  r_mech = 0.67 +/- 0.03 fm  (mechanical radius from GPDs)
  This is SMALLER than the charge radius r_p = 0.8414 fm.

THEORY (torsion medium Zone 3):
  Zone 3: lambda_p < r < r_p  (spinning co-rotation shell)
  The pressure reduction from spinning cells varies as P(r) ~ r^n.
  The mechanical RMS radius:
    <r^2> = int_{lambda_p}^{r_p} r^2 * P(r) * r^2 dr / int P(r) * r^2 dr
  Four physically motivated models:
    n = -2  (angular momentum conservation: v ~ 1/r, P ~ rho*v^2 ~ 1/r^2)
    n =  0  (uniform pressure in Zone 3: Bernoulli from shear wave G = Rs^2*K)
    n = +1  (Couette-like intermediate flow)
    n = +2  (rigid body rotation: v ~ r, P ~ rho*v^2 ~ r^2)

Run: python analysis/nuclear/mechanical_radius.py
Reference: docs/doc_nucleus.txt, ESSENTIALLY CLOSED section
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, phi, r_p, hbar_c, L_J as L_J_fm

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

pi      = math.pi
Rs      = math.sqrt(5) / (4*pi)
r_p_fm  = r_p * 1e15            # metres -> fm
m_p     = 938.272               # MeV
lambda_p_fm = hbar_c / m_p      # proton Compton wavelength (fm) = 0.2103 fm

# Measured mechanical radius (Burkert et al. 2018, Science 361, 207)
r_mech_measured    = 0.67   # fm  (best value)
r_mech_uncertainty = 0.03   # fm  (1-sigma)

# ── SECTION 1: Zone 3 geometry ────────────────────────────────────────────────
print(SEP)
print("SECTION 1: ZONE 3 GEOMETRY")
print(SEP2)
print(f"  lambda_p = {lambda_p_fm:.4f} fm  (Zone 1/2 boundary, jamming scale)")
print(f"  r_p      = {r_p_fm:.4f} fm  (Zone 3/4 boundary, charge radius)")
print(f"  Zone 3 width = r_p - lambda_p = {r_p_fm - lambda_p_fm:.4f} fm")
print(f"  r_p/lambda_p = {r_p_fm/lambda_p_fm:.4f}  (should be ~4)")
print()
print(f"  Measured mechanical radius: {r_mech_measured} +/- {r_mech_uncertainty} fm")
print(f"  (Burkert et al. 2018, Jefferson Lab CLAS12, gravitational form factor)")

# ── SECTION 2: Pressure profile models ───────────────────────────────────────
print()
print(SEP)
print("SECTION 2: MECHANICAL RADIUS FOR FOUR PRESSURE PROFILES")
print(SEP2)

def mechanical_radius(n, r_lo, r_hi):
    """
    Compute pressure-weighted RMS radius for P(r) ~ r^n in [r_lo, r_hi].
    <r^2> = int r^2 * r^n * r^2 dr / int r^n * r^2 dr
           = int r^(n+4) dr / int r^(n+2) dr
    """
    if abs(n + 5) < 1e-10:  # singularity at n = -5
        return None
    if abs(n + 3) < 1e-10:  # singularity at n = -3
        return None

    def intpow(r, p):
        if abs(p + 1) < 1e-10:
            return math.log(r)
        return r**(p+1) / (p+1)

    numerator   = intpow(r_hi, n+4) - intpow(r_lo, n+4)
    denominator = intpow(r_hi, n+2) - intpow(r_lo, n+2)

    if denominator == 0 or numerator/denominator < 0:
        return None
    return math.sqrt(numerator / denominator)

models = [
    (-2, "Angular momentum conservation  v~1/r, P~1/r^2"),
    ( 0, "Uniform pressure (shear wave G=Rs^2*K)        "),
    (+1, "Couette-like intermediate flow  P~r            "),
    (+2, "Rigid body rotation             v~r, P~r^2    "),
]

print(f"  {'n':>4}  {'Physical model':<46}  {'r_mech (fm)':>12}  {'vs measured':>12}")
print(f"  {'-'*4}  {'-'*46}  {'-'*12}  {'-'*12}")

best_n, best_r, best_model = None, None, None
min_dist = 1e10
for n, label in models:
    r_mech = mechanical_radius(n, lambda_p_fm, r_p_fm)
    if r_mech is None:
        print(f"  {n:>4}  {label}  {'N/A':>12}")
        continue
    dist = abs(r_mech - r_mech_measured)
    sigma = dist / r_mech_uncertainty
    match = f"{sigma:.1f}σ"
    print(f"  {n:>4}  {label}  {r_mech:>12.4f}  {match:>12}")
    if dist < min_dist:
        min_dist, best_n, best_r, best_model = dist, n, r_mech, label

print()
print(f"  Best match: n = {best_n}  ({best_model.strip()})")
print(f"              r_mech = {best_r:.4f} fm  ({min_dist/r_mech_uncertainty:.1f}σ from measured)")

# ── SECTION 3: Find n that gives exactly r_mech = 0.67 fm ────────────────────
print()
print(SEP)
print("SECTION 3: PRESSURE PROFILE THAT MATCHES r_mech = 0.67 fm EXACTLY")
print(SEP2)

# Binary search for n
n_lo, n_hi = -2.0, 2.0
for _ in range(100):
    n_mid = (n_lo + n_hi) / 2
    r = mechanical_radius(n_mid, lambda_p_fm, r_p_fm)
    if r is None: break
    if r < r_mech_measured:
        n_lo = n_mid
    else:
        n_hi = n_mid

n_exact = (n_lo + n_hi) / 2
r_exact = mechanical_radius(n_exact, lambda_p_fm, r_p_fm)

print(f"  For r_mech = {r_mech_measured} fm exactly:")
print(f"    Pressure profile P(r) ~ r^n  with  n = {n_exact:.4f}")
print(f"    Verification: r_mech({n_exact:.4f}) = {r_exact:.5f} fm")
print()

# Physical interpretation of n_exact
print(f"  Physical interpretation:")
print(f"    n = 0.0:  uniform pressure  [P = const]")
print(f"    n = 0.1:  very mild radial increase  [slightly rigid-body-like]")
print(f"    n_exact = {n_exact:.3f}  [between uniform and mild angular momentum increase]")
print()

# Shear modulus connection
G_over_K = Rs**2
P_shear_norm = G_over_K  # G/K = Rs^2 is the shear to bulk modulus ratio
print(f"  Torsion medium shear wave connection:")
print(f"    G/K = Rs^2 = {G_over_K:.6f}")
print(f"    The uniform pressure model (n=0) corresponds to constant G/K = Rs^2")
print(f"    in Zone 3 -- consistent with the shear wave (v_s = Rs*c) being the")
print(f"    dominant coupling mechanism in the spinning Zone 3 shell.")

# ── SECTION 4: Key checks ─────────────────────────────────────────────────────
print()
print(SEP)
print("SECTION 4: CHECKS")
print(SEP2)

r_mech_uniform = mechanical_radius(0, lambda_p_fm, r_p_fm)
r_mech_angmom  = mechanical_radius(-2, lambda_p_fm, r_p_fm)

check("MR1 Uniform pressure (n=0) gives r_mech within 1-sigma of measurement",
      abs(r_mech_uniform - r_mech_measured) <= r_mech_uncertainty,
      f"r_mech(n=0) = {r_mech_uniform:.4f} fm  measured = {r_mech_measured} +/- {r_mech_uncertainty} fm")

check("MR2 r_mech < r_p  (mechanical radius < charge radius)",
      r_mech_uniform < r_p_fm,
      f"r_mech = {r_mech_uniform:.4f} fm < r_p = {r_p_fm:.4f} fm")

check("MR3 r_mech > lambda_p  (mechanical radius > jamming scale)",
      r_mech_uniform > lambda_p_fm,
      f"r_mech = {r_mech_uniform:.4f} fm > lambda_p = {lambda_p_fm:.4f} fm")

check("MR4 n_exact is close to 0 (uniform pressure is near-exact model)",
      abs(n_exact) < 0.5,
      f"n_exact = {n_exact:.4f}  (|n| < 0.5 means near-uniform)")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP)
passed = sum(1 for _,s,_ in results if s=="PASS")
failed = sum(1 for _,s,_ in results if s=="FAIL")
print(f"  Total checks: {len(results)}   PASS: {passed}   FAIL: {failed}")
print()
print(f"  Zone 3 boundary: lambda_p = {lambda_p_fm:.4f} fm  to  r_p = {r_p_fm:.4f} fm")
print(f"  Uniform pressure model (P = const):  r_mech = {r_mech_uniform:.4f} fm")
print(f"  Measured (Burkert 2018):             r_mech = {r_mech_measured} +/- {r_mech_uncertainty} fm")
print(f"  Agreement: {abs(r_mech_uniform-r_mech_measured)/r_mech_uncertainty:.1f}σ")
print()
print(f"  PHYSICAL CONCLUSION:")
print(f"  The Zone 3 spinning cells have approximately UNIFORM pressure,")
print(f"  consistent with the shear modulus G = Rs^2*K being the constant")
print(f"  coupling in the zone. This uniform G is the torsion medium's")
print(f"  shear wave (v_s = Rs*c) maintaining constant pressure throughout")
print(f"  the co-rotating shell. The mechanical radius ({r_mech_uniform:.3f} fm)")
print(f"  is smaller than the charge radius ({r_p_fm:.4f} fm) because the")
print(f"  pressure is concentrated in Zone 3, not at r_p itself.")
print()
print(f"  Reference: docs/doc_nucleus.txt")
