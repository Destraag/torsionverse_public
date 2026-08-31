"""
alpha_cg_correction.py
=======================
Investigates whether the new CG/I_h properties (from doc_higgs Section 5a)
explain the missing accuracy in the alpha quadratic derivation.

THE OPEN STEP IN alpha DERIVATION:
  n_exact = 2 + delta_n  where delta_n = 0.01869 (empirical, from CODATA alpha)
  delta_n = L3(phi, log5) * k_n/k_eff  [from alpha derivation structure]
  k_n/k_eff = delta_n/L3 = 0.01869/1.6138 = 0.01158
  This ratio is NOT yet derived from first principles (the "open step").

QUESTION: Does k_n/k_eff = alpha*phi?
  alpha*phi = 7.2974e-3 * 1.6180 = 0.01181
  vs empirical: 0.01158
  Discrepancy: (0.01181 - 0.01158)/0.01158 = 2.0%

Further: does the T_1g x T_1g CG structure give a cleaner answer?
  T_1g x T_1g = A_g + T_1g + H_g  (1+3+5 = 9 modes)
  Various CG-derived ratios tested against k_n/k_eff = 0.01158.

Run: python analysis/alpha/alpha_cg_correction.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3   # CODATA-2018
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4*pi)
log5  = math.log(5)

SEP  = "=" * 70
SEP2 = "-" * 70

# ── Known values from alpha derivation ───────────────────────────────────────
n       = 2               # linking number (exact integer)
n_exact = 2.01869         # from CODATA alpha (empirical)
delta_n = n_exact - n     # = 0.01869

# L3(phi, log5): the Born-weighted mean
L3 = (phi**3 + log5**3) / (phi**2 + log5**2)

# Empirical ratio that needs to be derived
kn_keff_empirical = delta_n / L3

# Alpha quadratic: residual after correction
alpha_codata = 7.2973525693e-3
Q = 4*pi**2/phi
alpha_quad_exact = (Q - math.sqrt(Q**2 - 4*n*Rs)) / (2*n)
delta_alpha = alpha_quad_exact - alpha_codata

print(SEP)
print("ALPHA DERIVATION: MISSING ACCURACY FROM CG STRUCTURE")
print(SEP)
print()
print(f"  Known quantities:")
print(f"    phi   = {phi:.10f}")
print(f"    log5  = {log5:.10f}")
print(f"    L3(phi,log5) = (phi^3+log5^3)/(phi^2+log5^2) = {L3:.10f}")
print(f"    Rs    = sqrt(5)/(4*pi) = {Rs:.10f}")
print()
print(f"  Open step:")
print(f"    delta_n = n_exact - 2 = {delta_n:.5f}")
print(f"    delta_n = L3 * k_n/k_eff  =>  k_n/k_eff = {kn_keff_empirical:.8f}")
print()
print(f"  Alpha residual with n_exact:")
print(f"    alpha (quadratic, n=n_exact) = {alpha_quad_exact:.10e}")
print(f"    alpha (CODATA)               = {alpha_codata:.10e}")
print(f"    delta_alpha / alpha          = {delta_alpha/alpha_codata:.2e}  ({delta_alpha/alpha_codata*100:.6f}%)")
print()

# ── CANDIDATE 1: k_n/k_eff = alpha*phi ───────────────────────────────────────
print(SEP)
print("CANDIDATE 1  k_n/k_eff = alpha*phi")
print(SEP2)
kn_c1 = alpha * phi
delta_n_c1 = L3 * kn_c1
residual_c1 = (delta_n_c1 - delta_n) / delta_n
print(f"  alpha*phi = {kn_c1:.8f}")
print(f"  delta_n from this: L3 * alpha * phi = {delta_n_c1:.8f}")
print(f"  vs empirical delta_n:                 {delta_n:.8f}")
print(f"  Residual: {residual_c1*100:+.4f}%")
print()

# ── CANDIDATE 2: k_n/k_eff = alpha*phi*(1 - correction) ─────────────────────
print(SEP)
print("CANDIDATE 2  Fibonacci correction to candidate 1")
print(SEP2)
# alpha*phi^2 = alpha*(phi+1) = alpha*phi + alpha
# Perhaps k_n/k_eff = alpha*phi / (1 + alpha*phi)?
kn_c2a = alpha*phi / (1 + alpha*phi)
delta_n_c2a = L3 * kn_c2a
print(f"  alpha*phi/(1+alpha*phi) = {kn_c2a:.8f}")
print(f"  delta_n:                  {delta_n_c2a:.8f}  residual: {(delta_n_c2a-delta_n)/delta_n*100:+.4f}%")

# alpha/phi?
kn_c2b = alpha/phi
delta_n_c2b = L3 * kn_c2b
print(f"  alpha/phi              = {kn_c2b:.8f}")
print(f"  delta_n:                 {delta_n_c2b:.8f}  residual: {(delta_n_c2b-delta_n)/delta_n*100:+.4f}%")

# alpha*(phi-1) = alpha/phi (since phi-1 = 1/phi)
# same as above

# alpha*sqrt(5)?  sqrt(5)=2*phi-1
kn_c2c = alpha*sqrt5
delta_n_c2c = L3 * kn_c2c
print(f"  alpha*sqrt(5)          = {kn_c2c:.8f}")
print(f"  delta_n:                 {delta_n_c2c:.8f}  residual: {(delta_n_c2c-delta_n)/delta_n*100:+.4f}%")

# Rs/phi^2?
kn_c2d = Rs / phi**2
delta_n_c2d = L3 * kn_c2d
print(f"  Rs/phi^2               = {kn_c2d:.8f}")
print(f"  delta_n:                 {delta_n_c2d:.8f}  residual: {(delta_n_c2d-delta_n)/delta_n*100:+.4f}%")
print()

# ── CANDIDATE 3: CG-derived ratio ────────────────────────────────────────────
print(SEP)
print("CANDIDATE 3  CG-derived ratios from T_1g x T_1g = A_g + T_1g + H_g")
print(SEP2)
print()
# T_1g x T_1g = A_g(1) + T_1g(3) + H_g(5) = 9 modes
# Possible ratios:
ratios = {
    "dim(T_1g)/dim(T1xT1)": 3/9,
    "dim(A_g)/dim(T1xT1)":  1/9,
    "dim(H_g)/dim(T1xT1)":  5/9,
    "dim(T_1g)/dim(T_2g+H)": 3/(3+5),
    "chi(T_1g,C5)/dim":     phi/3,
    "chi(A_g,C5)/chi(T1xT1,C5)": 1/phi**2,
    "1/(1+phi+phi^2)":      1/(1+phi+phi**2),
    "phi^2/(1+phi+phi^2)":  phi**2/(1+phi+phi**2),
}
for name, ratio in ratios.items():
    dn_pred = L3 * ratio
    pct = (dn_pred - delta_n)/delta_n * 100
    marker = " <-- CLOSEST" if abs(pct) < 5 else ""
    print(f"  {name:<32} = {ratio:.6f}  delta_n={dn_pred:.6f}  {pct:+.2f}%{marker}")
print()

# ── CANDIDATE 4: The alpha*phi^2 connection ───────────────────────────────────
print(SEP)
print("CANDIDATE 4  Connection to alpha^2*phi^2 (vev correction term)")
print(SEP2)
print()
# The vev correction is alpha^2*phi^2. Could delta_n involve alpha^2*phi^2?
# delta_n/L3 = k_n/k_eff = 0.01158
# alpha^2*phi^2 = 1.394e-4 -- much smaller
# But: alpha^2*phi^2 / alpha = alpha*phi^2 = 0.01910
kn_c4a = alpha*phi**2
print(f"  alpha*phi^2 = alpha*(phi+1) = {kn_c4a:.8f}")
dn_c4a = L3*kn_c4a
print(f"  delta_n = L3*alpha*phi^2 = {dn_c4a:.8f}  residual: {(dn_c4a-delta_n)/delta_n*100:+.4f}%")
print()

# What about: alpha*phi = 0.01181, alpha*phi^2 = 0.01910
# Both bracket k_n/k_eff = 0.01158
# Linear interpolation: 0.01158 = (1-x)*alpha*phi + x*alpha/phi where x is some weight
x_interp = (kn_keff_empirical - alpha*phi) / (alpha/phi - alpha*phi)
print(f"  Linear interp: k_n/k_eff = (1-x)*alpha*phi + x*alpha/phi")
print(f"    x = {x_interp:.4f}  (x=0 gives alpha*phi, x=1 gives alpha/phi)")
print()

# ── CANDIDATE 5: Does L3 itself involve the CG structure? ────────────────────
print(SEP)
print("CANDIDATE 5  Does L3(phi, log5) come from CG character sums?")
print(SEP2)
print()
print(f"  L3(phi, log5) = (phi^3+log5^3)/(phi^2+log5^2) = {L3:.8f}")
print(f"  phi  = {phi:.8f}  (chi(T_1g,C_5) = chi(E_1/2,C_5))")
print(f"  log5 = {log5:.8f}")
print(f"  Difference phi - log5 = {phi-log5:.8f}")
print()
# phi and log5 are VERY close:
print(f"  The proximity phi ≈ log5 is remarkable:")
print(f"    phi - log5 = {phi-log5:.6f}  ({(phi-log5)/phi*100:.4f}% fractional)")
print()
# log5 from icosahedral counting?
# The icosahedron has 12 vertices, 30 edges, 20 faces
# C_5 rotation elements: 12+12=24 out of 60 total
# log(5): 5-fold symmetry -> number of edges per vertex = 5
print(f"  Geometric origin of log5:")
print(f"    Each vertex has 5 nearest neighbors (5-fold coordination)")
print(f"    The Born weight log-sum over all 5 nearest-neighbor interactions:")
print(f"    sum_k log(k) for k=1..5? = log(120) ≠ log5")
print(f"    Or: log of vertex coordination number = log(5) = {log5:.6f}  [EXACT]")
print()
print(f"  INSIGHT: log5 = log(vertex coordination number of icosahedron)")
print(f"  phi = golden ratio = chi(T_1g,C_5)")
print(f"  L3 is the weighted mean of chi(T_1g,C_5) and log(coordination).")
print()

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP2)
print(f"  Target: k_n/k_eff = {kn_keff_empirical:.8f}")
print(f"  Best candidate: alpha*phi = {alpha*phi:.8f}  (2.0% off)")
print()
print(f"  The 2.0% gap between alpha*phi and k_n/k_eff is equivalent to:")
print(f"    delta_n_predicted = L3*alpha*phi = {L3*alpha*phi:.6f}")
print(f"    delta_n_empirical =               {delta_n:.6f}")
print(f"    Mass correction: delta_n determines n_exact which determines alpha.")
print(f"    2% error in delta_n -> {2*0.02/137:.2e} fractional error in alpha")
print(f"    = {2*0.02/137*7.297e-3:.2e} absolute = {2*0.02:.0f}x smaller than current residual")
print()
print(f"  CONCLUSION: k_n/k_eff = alpha*phi is a strong CANDIDATE but not exact.")
print(f"  If correct, it would provide a first-principles derivation of delta_n.")
print(f"  The remaining 2% gap may come from:")
print(f"    - Higher-order Fibonacci correction: alpha*phi*(1-alpha*phi^2)?")
print(f"    - L3 correction: using phi exactly vs L3(phi,log5)")
print(f"    - A CG coefficient I haven't yet identified")
print()
# Check: does alpha*phi^2 / (1+alpha*phi) = k_n/k_eff?
kn_test = alpha*phi / (1 + alpha*phi**2)
print(f"  Test: alpha*phi/(1+alpha*phi^2) = {kn_test:.8f}  target={kn_keff_empirical:.8f}")
print(f"  Residual: {(kn_test-kn_keff_empirical)/kn_keff_empirical*100:+.4f}%  <-- better!")
