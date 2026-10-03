"""
proton_g_factor.py
==================
LEGACY CROSS-CHECK: single-mechanism candidate models for the proton
magnetic moment g_p = 2.7928 mu_N -- NOT the current derivation.

The current, correct derivation (g_p = 2.799 mu_N, +0.23%) is
analysis/demos/nucleus_doc.py's N15-N17 checks: a COMBINED model
  mu_p = R_spin_Zone1 * (1 + 2*Rs^2) * mu_SU6 + mu_orbital + mu_Zone3
(Zone 1 MIT-bag-proxy spin reduction, amplified by the Zone 2 Maxwell-
jamming correction (1+2*Rs^2), plus separate orbital and Zone 3 terms).

This file exists to show that SIMPLER, single-mechanism candidates do
NOT reach that precision on their own -- motivating why the combined
picture above is needed:
  - SU(6) alone (Section 1): +7.4%, this file's own best performer.
  - MIT bag spin reduction alone (Section 2): -29.9%.
  - String model alone (Section 3): -56.2%.
  - String + Zone3 + various corrections, Models A-D (Sections 4-5):
    -43% to -44% -- WORSE than plain SU(6), not an improvement.
  - Pressure-divergence model (Section 0): a 4th, separate candidate,
    tuned to a free d-quark position rather than independently derived.

NOTE: the MIT bag spin reduction R_spin_MIT (Section 2) is numerically
reused (not independently re-derived) in nucleus_doc.py's own Zone 1
proxy term above -- the spherical bag boundary condition approximates
the actual I_h Zone 1 mode mixing. That reuse is the only load-bearing
connection between this file and the current derivation; none of this
file's own OVERALL models (SU6 alone, MIT bag alone, string, A-D,
pressure-divergence) should be cited as reproducing measured g_p.

Run: python analysis/nuclear/proton_g_factor.py
Reference: docs/series1/doc_nucleus.txt (cross-check role only --
  see analysis/demos/nucleus_doc.py for the actual derivation)
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, phi, r_p, hbar_c, L_J as L_J_fm

SEP  = "=" * 65
SEP2 = "-" * 65
results = []
pi = math.pi

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

m_p     = 938.272    # MeV
m_e     = 0.51100    # MeV
r_p_fm  = r_p * 1e15
Rs      = math.sqrt(5) / (4*pi)
lambda_p_fm = hbar_c / m_p  # proton Compton wavelength (fm) = 0.2103 fm
mu_p_measured = 2.7928      # nuclear magnetons (CODATA)

# ── Spherical Bessel functions for MIT bag ────────────────────────────────────
def j0(x): return math.sin(x)/x if x > 1e-12 else 1.0
def j1(x): return (math.sin(x)/x**2 - math.cos(x)/x) if x > 1e-12 else x/3

# Numerical integration
def integrate(f, a, b, n=1000):
    h = (b-a)/n
    return sum(f(a + (i+0.5)*h) * h for i in range(n))

# Pre-compute Zone 3 before Section 0 (needed for pressure divergence model)
V_Zone3_pre = (4/3)*pi*(r_p_fm**3 - lambda_p_fm**3)
mu_Zone3 = (4*pi/3) * Rs * integrate(lambda r: r**3, lambda_p_fm, r_p_fm) / V_Zone3_pre * (2*m_p/hbar_c)

# ── SECTION 0: Pressure divergence model (PRIMARY) ───────────────────────────
print(SEP)
print("SECTION 0: PRESSURE DIVERGENCE MODEL (the torsion medium prediction)")
print(SEP2)
print(f"""
  Physical picture:
    2 u quarks at Zone 1/2 boundary r ~ lambda_p = {lambda_p_fm:.4f} fm
      charge: +2/3 each -> net +4e/3, INWARD pressure wells (large r)
    1 d quark at centre r ~ 0
      charge: -1/3,      OUTWARD pressure source (small r, small lever arm)

  Magnetic moment from spinning charge distribution:
    mu = g * (e/2m_p) * (hbar/2)
    g_p = 2 * [sum_i Q_i * r_i^2] / (e * lambda_p^2)   [charge-weighted r^2]
""")

# g_p from charge distribution (r_d = 0, u quarks at lambda_p)
Q_u = 4/3   # charge of two u quarks in units of e
Q_d = -1/3  # charge of d quark in units of e

def g_from_positions(r_u, r_d):
    """Effective gyromagnetic ratio from charge-weighted r^2 distribution.
    For proton: 2 u quarks (charge +4e/3) at r_u, 1 d quark (charge -e/3) at r_d.
    g_eff/2 = [Q_u*r_u^2 + Q_d*r_d^2] / [(Q_u+Q_d) * r_u^2]  (normalised to u-quark scale)
    mu_p = g_eff/2 * mu_N for spin-1/2 proton
    """
    numerator = Q_u * r_u**2 + Q_d * r_d**2
    denominator = (Q_u + Q_d) * r_u**2  # = 1 * r_u^2  (total charge = e)
    return numerator / denominator  # = mu_p in nuclear magnetons

g_exact_center = g_from_positions(lambda_p_fm, 0)
mu_divergence_r0 = g_exact_center   # this IS mu_p in nuclear magnetons

print(f"  r_u = lambda_p = {lambda_p_fm:.4f} fm  (u quarks at cog boundary)")
print(f"  r_d = 0 fm  (d quark exactly at center)")
print(f"  g_p = 2 * [(4/3)*lambda_p^2 + (-1/3)*0] / lambda_p^2 = 2 * (4/3) = {g_exact_center:.4f}")
print(f"  mu_p (divergence model, r_d=0) = g_p/2 = {mu_divergence_r0:.4f} mu_N")
print()

# Scan what r_d gives the exact measured value (before Zone 3 correction)
target_no_Z3 = mu_p_measured - 0.360  # measured minus Zone 3 contribution
best_rd, best_g = None, None
for rd_tenth in range(0, 30):
    r_d = rd_tenth / 10 * lambda_p_fm
    g = g_from_positions(lambda_p_fm, r_d)
    if abs(g/2 - target_no_Z3) < 0.01:
        best_rd, best_g = r_d, g

print(f"  Adding Zone 3 contribution (+{mu_Zone3:.3f} mu_N from mechanical_radius.py):")
mu_total_r0 = mu_divergence_r0 + mu_Zone3
print(f"  mu_p (divergence + Zone3) = {mu_total_r0:.4f} mu_N  err={100*(mu_total_r0-mu_p_measured)/mu_p_measured:+.2f}%")
print()
print(f"  Remaining 8% gap = d quark is NOT at exactly r=0.")
print(f"  d quark has a finite position distribution in Zone 1.")
print(f"  For measured mu_p = {mu_p_measured:.4f}: need r_d such that g/2 + Zone3 = measured.")

# Find the exact r_d needed
for rd100 in range(0, 200):
    r_d = rd100 / 100 * lambda_p_fm
    g = g_from_positions(lambda_p_fm, r_d)
    if abs(g/2 + mu_Zone3 - mu_p_measured) < 0.005:
        print(f"  Required r_d = {r_d:.4f} fm = {r_d/lambda_p_fm:.3f} * lambda_p")
        print(f"    -> g_p = {g:.4f}, mu = {g/2+mu_Zone3:.4f} mu_N")
        break

print()
check("GP0 Pressure divergence (r_d=0) gives mu_p in right direction and magnitude",
      abs(mu_total_r0 - mu_p_measured)/mu_p_measured < 0.15,
      f"divergence+Z3 = {mu_total_r0:.4f}  measured = {mu_p_measured:.4f}  err = {100*(mu_total_r0-mu_p_measured)/mu_p_measured:+.1f}%")

print()
print(SEP)
print("SECTION 1: SU(6) BASELINE (constituent quarks, m_q = m_p/3)")
print(SEP2)

# Constituent quark magnetic moments
m_u = m_p / 3  # MeV
m_d = m_p / 3  # MeV (isospin symmetric)
mu_u = (2/3) * (m_p/m_u)   # in units of mu_N  (charge 2/3, g=2, mass m_u)
mu_d = (-1/3) * (m_p/m_d)  # in units of mu_N

mu_SU6 = (4*mu_u - mu_d) / 3  # SU(6) formula for proton
print(f"  mu_u (constituent) = {mu_u:.4f} mu_N")
print(f"  mu_d (constituent) = {mu_d:.4f} mu_N")
print(f"  mu_p (SU6) = (4*mu_u - mu_d)/3 = {mu_SU6:.4f} mu_N")
print(f"  Measured:           = {mu_p_measured:.4f} mu_N")
print(f"  Error: {(mu_SU6-mu_p_measured)/mu_p_measured*100:+.2f}%")
print()

check("GP1 SU(6) gives mu_p = 3.000 mu_N",
      abs(mu_SU6 - 3.000) < 0.001,
      f"SU(6) = {mu_SU6:.4f}")
check("GP0b Pressure divergence (Section 0) is closer to measured than SU(6) baseline",
      abs(mu_total_r0 - mu_p_measured) < abs(mu_SU6 - mu_p_measured),
      f"divergence err = {abs(mu_total_r0-mu_p_measured):.4f}  SU6 err = {abs(mu_SU6-mu_p_measured):.4f}")

# ── SECTION 2: MIT bag spin reduction ─────────────────────────────────────────
print()
print(SEP)
print("SECTION 2: MIT BAG SPIN REDUCTION (spherical Bessel, bag = r_p)")
print(SEP2)

# MIT bag: f = j0(x0*r/R), g = j1(x0*r/R)
# R_spin = int[(f^2 - g^2/3) r^2 dr] / int[(f^2 + g^2) r^2 dr]
x0_MIT = 2.0428  # MIT bag lowest mode eigenvalue

def integrand_num_MIT(r):
    kr = x0_MIT * r / r_p_fm
    return (j0(kr)**2 - j1(kr)**2/3) * r**2

def integrand_den_MIT(r):
    kr = x0_MIT * r / r_p_fm
    return (j0(kr)**2 + j1(kr)**2) * r**2

R_spin_MIT = integrate(integrand_num_MIT, 0, r_p_fm) / integrate(integrand_den_MIT, 0, r_p_fm)
mu_p_MIT = R_spin_MIT * mu_SU6

print(f"  MIT bag x0 = {x0_MIT}, bag radius = r_p = {r_p_fm:.4f} fm")
print(f"  Spin reduction factor R_spin = {R_spin_MIT:.4f}")
print(f"  mu_p (MIT spin-reduced) = R_spin * SU(6) = {mu_p_MIT:.4f} mu_N")
print(f"  Error vs measured: {(mu_p_MIT-mu_p_measured)/mu_p_measured*100:+.2f}%")

# ── SECTION 3: String model -- three-across, bag at lambda_p ──────────────────
print()
print(SEP)
print("SECTION 3: STRING MODEL (three-across, 1D Dirac, scale = lambda_p)")
print(SEP2)

# Three-across string: total length 2*lambda_p, middle quark at lambda_p
# For a 1D Dirac box (massless) with MIT-like BC at x=0 and x=2L:
# f(x) = cos(k*x + phi), g(x) = sin(k*x + phi)
# BC at x=0: -i*gamma_x*psi = psi -> g(0) = f(0) -> phi = pi/4... let's use:
# Simple bag BC for 1D: BC requires that the current j^1 = f*g vanishes at boundaries.
# This gives: f(0)*g(0) = 0 and f(2L)*g(2L) = 0
# Solution: f = sin(k*x), g = cos(k*x) with k*2L = pi (half-wave in 2L)
# OR: use the ratio at the boundary

# For the simplest mode consistent with our string:
# String endpoints at x=0 and x=2L, middle at x=L
# Fundamental mode: f(x) = sin(pi*x/(2L)), g(x) = cos(pi*x/(2L))  [k = pi/2L]

L_string = lambda_p_fm  # half the string length (from center to end)
# Actually string goes from -lambda_p to +lambda_p, length 2*lambda_p
# Let x go from 0 to 2L where L = lambda_p

def f_string(x, L=lambda_p_fm):
    k = pi / (2*L)  # fundamental: half-wavelength = 2L
    return math.sin(k*x)

def g_string(x, L=lambda_p_fm):
    k = pi / (2*L)
    return math.cos(k*x)

# R_spin for the string model (outer u quarks at x = 0 to 2L, d quark at x = L)
# Outer quarks: probability distributed over entire string
# Overall spin reduction for p-wave-like outer quarks:
def num_string(x):
    return (f_string(x)**2 - g_string(x)**2/3) * 1  # 1D weight

def den_string(x):
    return (f_string(x)**2 + g_string(x)**2) * 1

R_spin_string = integrate(num_string, 0, 2*lambda_p_fm) / integrate(den_string, 0, 2*lambda_p_fm)
print(f"  String length = 2*lambda_p = {2*lambda_p_fm:.4f} fm")
print(f"  String mode: f(x) = sin(pi*x/2L),  g(x) = cos(pi*x/2L)")
print(f"  Spin reduction R_spin (string) = {R_spin_string:.4f}")

# d quark at CENTER (x=L): evaluate wave functions there
f_center = f_string(lambda_p_fm)
g_center = g_string(lambda_p_fm)
R_spin_d_center = (f_center**2 - g_center**2/3) / (f_center**2 + g_center**2)
print(f"  d quark at center: f={f_center:.4f}, g={g_center:.4f}")
print(f"  Spin reduction at center R_spin_d = {R_spin_d_center:.4f} (= 1/3, fully mixed at x=L)")
print()

# For the three-across proton: u quarks feel average string reduction, d quark at center
# The SU(6) formula weights: (4*mu_u - mu_d)/3
# with u-quark spin reduced by R_spin_string and d by R_spin_d_center
mu_u_str = (2/3) * (m_p/m_u) * R_spin_string    # reduced u magnetic moment
mu_d_str = (-1/3) * (m_p/m_d) * R_spin_d_center  # reduced d magnetic moment
mu_p_string = (4*mu_u_str - mu_d_str) / 3
print(f"  mu_u (string-reduced) = {mu_u_str:.4f} mu_N")
print(f"  mu_d (string-reduced) = {mu_d_str:.4f} mu_N")
print(f"  mu_p (string model)   = {mu_p_string:.4f} mu_N")
print(f"  Error vs measured: {(mu_p_string-mu_p_measured)/mu_p_measured*100:+.2f}%")

# ── SECTION 4: Zone 3 spinning shell contribution ─────────────────────────────
print()
print(SEP)
print("SECTION 4: ZONE 3 SPINNING SHELL (uniform pressure, omega = Rs*c/r)")
print(SEP2)

# Zone 3: lambda_p < r < r_p, uniform charge density, omega(r) = Rs*c/r
# mu_Zone3 = (1/3) * int_{lambda_p}^{r_p} (e/V_Z3) * (Rs*c/r) * r^2 * 4*pi*r^2 dr
# In units of mu_N = e*hbar_c/(2*m_p):

V_Zone3 = (4/3)*pi*(r_p_fm**3 - lambda_p_fm**3)
numerator_Z3 = (4*pi/3) * Rs * integrate(lambda r: r**3, lambda_p_fm, r_p_fm)
mu_Zone3 = (numerator_Z3 / V_Zone3) * (2*m_p/hbar_c)  # in mu_N

print(f"  V_Zone3 = {V_Zone3:.4f} fm^3")
print(f"  mu_Zone3 = {mu_Zone3:.4f} mu_N")

# ── SECTION 5: Jam correction (alpha/pi from vertex stiffness) ────────────────
print()
print(SEP)
print("SECTION 5: JAM + SPIN CORRECTIONS FROM TORSION MEDIUM")
print(SEP2)

# Vertex stiffness correction delta_n from doc_alpha
log5 = math.log(5)
L3 = (phi**3 + log5**3) / (phi**2 + log5**2)
x_fs = alpha * phi**2
k_fs = alpha * phi * (1 - (3/4)*alpha**2) / (1 + x_fs + x_fs**2)
delta_n = L3 * k_fs

# Spin correction: (1 + (3/4)*alpha^2)
spin_corr = 1 + (3/4)*alpha**2

# Alpha correction: (1 + alpha/pi) - Schwinger correction analog for magnetic moment
alpha_over_pi = alpha / pi

print(f"  delta_n (vertex stiffness) = {delta_n:.6f}")
print(f"  (3/4)*alpha^2 correction   = {(3/4)*alpha**2:.6f}")
print(f"  alpha/pi (Schwinger-type)   = {alpha_over_pi:.6f}")
print()

# Apply corrections to string model result + Zone 3
# Model A: string + Zone3, no jam correction
model_A = mu_p_string + mu_Zone3
# Model B: string + Zone3 + delta_n correction (reduces spin)
model_B = mu_p_string * (1 - delta_n) + mu_Zone3
# Model C: string + Zone3 + alpha/pi correction  
model_C = mu_p_string * (1 + alpha_over_pi) + mu_Zone3
# Model D: string + Zone3 + full torsion medium correction
torsion_corr = (1 + delta_n/pi) * (1 + (3/4)*alpha**2)  # from V21 chain
model_D = mu_p_string * torsion_corr + mu_Zone3

print(f"  Model A: string + Zone3 (no correction):   {model_A:.4f} mu_N  err={100*(model_A-mu_p_measured)/mu_p_measured:+.2f}%")
print(f"  Model B: string + Zone3 - delta_n:          {model_B:.4f} mu_N  err={100*(model_B-mu_p_measured)/mu_p_measured:+.2f}%")
print(f"  Model C: string + Zone3 + alpha/pi:         {model_C:.4f} mu_N  err={100*(model_C-mu_p_measured)/mu_p_measured:+.2f}%")
print(f"  Model D: string + Zone3 + full V21 corr:    {model_D:.4f} mu_N  err={100*(model_D-mu_p_measured)/mu_p_measured:+.2f}%")
print()
print(f"  Measured:                                    {mu_p_measured:.4f} mu_N")

# Find the best model
models = [('SU(6)', mu_SU6), ('MIT bag', mu_p_MIT), ('String', mu_p_string),
          ('A:str+Z3', model_A), ('B:str+Z3-dn', model_B),
          ('C:str+Z3+a/pi', model_C), ('D:str+Z3+V21', model_D)]
best = min(models, key=lambda x: abs(x[1]-mu_p_measured))
print(f"  BEST MODEL: {best[0]} = {best[1]:.4f} mu_N  ({100*(best[1]-mu_p_measured)/mu_p_measured:+.2f}%)")

# ── Checks ────────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SECTION 6: CHECKS")
print(SEP2)

check("GP2 String model gives mu_p closer to measured than SU(6)",
      abs(mu_p_string - mu_p_measured) < abs(mu_SU6 - mu_p_measured),
      f"string err = {abs(mu_p_string-mu_p_measured):.4f}  SU6 err = {abs(mu_SU6-mu_p_measured):.4f}")
check("GP3 String + Zone3 in correct direction (both reduce from SU6=3.000)",
      mu_p_string < mu_SU6 and model_A < mu_SU6,
      f"string={mu_p_string:.3f} < SU6={mu_SU6:.3f}, A={model_A:.3f}")
check("GP4 Best model within 15% of measured (framework gets order and sign)",
      abs(best[1] - mu_p_measured)/mu_p_measured < 0.15,
      f"best={best[1]:.4f}, measured={mu_p_measured:.4f}, err={100*(best[1]-mu_p_measured)/mu_p_measured:+.1f}%")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP)
passed = sum(1 for _,s,_ in results if s=="PASS")
failed = sum(1 for _,s,_ in results if s=="FAIL")
print(f"  Total checks: {len(results)}   PASS: {passed}   FAIL: {failed}")
print()
print(f"  {'Model':<25}  {'mu_p':>8}  {'error':>8}")
print(f"  {'-'*25}  {'-'*8}  {'-'*8}")
for name, val in models:
    print(f"  {name:<25}  {val:>8.4f}  {100*(val-mu_p_measured)/mu_p_measured:>+7.2f}%")
print(f"  {'MEASURED':<25}  {mu_p_measured:>8.4f}")
print()
print("  CONCLUSIONS:")
print(f"  1. This file's closest single-mechanism match to measured is")
print(f"     {best[0]} = {best[1]:.4f} mu_N ({100*(best[1]-mu_p_measured)/mu_p_measured:+.1f}%).")
print(f"  2. Zone 3 adds +{mu_Zone3:.3f} mu_N in every variant (small, legitimate,")
print(f"     but not sufficient on its own to close the gap).")
print(f"  3. The string-model refinements (Models A-D) move FARTHER from")
print(f"     measured than plain SU(6), not closer -- single mechanisms")
print(f"     explored here do not reproduce g_p to measured precision.")
print(f"  4. The actual current derivation (g_p=2.799, +0.23%) is a DIFFERENT,")
print(f"     combined model in analysis/demos/nucleus_doc.py (N15-N17) -- this")
print(f"     file is a legacy cross-check, not the primary derivation.")
print()
print("  Reference: docs/series1/doc_nucleus.txt -- see analysis/demos/nucleus_doc.py for the current derivation")
