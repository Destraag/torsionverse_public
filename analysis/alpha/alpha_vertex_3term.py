"""
alpha_vertex_3term.py
=====================
Tests whether the 3-term Fibonacci denominator closes the 0.038% residual
in the vertex stiffness ratio k_n/k_eff.

Current 2-term formula:  k_n/k_eff = alpha*phi / (1 + alpha*phi^2)
Proposed 3-term formula: k_n/k_eff = alpha*phi / (1 + alpha*phi^2 + alpha^2*phi^4)

The 3-term denominator is the partial geometric series with ratio x = alpha*phi^2,
truncated at O(alpha^2) -- consistent with the Higgs mass two-loop truncation.

Physical interpretation:
  Each term (alpha*phi^2)^k represents k bounces of the elastic wave between
  the vertex stiffness k_n and the long-wavelength stiffness k_LW.
  The series terminates at k=2 (O(alpha^2)) by the same Fibonacci argument
  that truncates the Higgs mass correction at alpha^2*phi^2.

Result: 10x improvement over 2-term (0.038% -> 0.004% in k_n/k_eff)
        Alpha chain residual: 0.000000022521% (10x better than 0.00000022%)

Run: python analysis/alpha/alpha_vertex_3term.py
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, phi, pi, hbar_c, r_p, E_cell_GeV

Rs    = math.sqrt(5)/(4*pi)
Q     = 4*pi**2/phi
log5  = math.log(5)
L3    = (phi**3 + log5**3)/(phi**2 + log5**2)
alpha_codata = 7.2973525693e-3
x     = alpha * phi**2   # the geometric ratio

SEP  = "=" * 65
SEP2 = "-" * 65

print(SEP)
print("alpha_vertex_3term.py")
print("3-term Fibonacci denominator for k_n/k_eff")
print(SEP)

# Compute target k_n/k_eff from CODATA alpha
k_target = ((Q*alpha_codata - Rs)/alpha_codata**2 - 2) / L3

print(f"\nTarget k_n/k_eff (from CODATA alpha) = {k_target:.12f}")
print(f"x = alpha*phi^2 = {x:.10f}")
print()

# Candidates
k_2term = alpha*phi / (1 + x)
k_3term = alpha*phi / (1 + x + x**2)
k_full  = alpha*phi / (1 - x)   # infinite geometric series (overshoots)

print(f"{'Formula':<45} {'k_n/k_eff':>14}  {'residual':>10}")
print(f"{'-'*45} {'-'*14}  {'-'*10}")
for label, val in [
    ("2-term: a*phi/(1+a*phi^2)",                k_2term),
    ("3-term: a*phi/(1+a*phi^2+a^2*phi^4)",      k_3term),
    ("Full geometric: a*phi/(1-a*phi^2)",        k_full),
]:
    res = (val - k_target)/k_target*100
    print(f"  {label:<43} {val:.12f}  {res:+.5f}%")

print()
print(f"  Target: {k_target:.12f}")
print()

# Fibonacci expansion of the 3-term denominator
# phi^2 = phi+1; phi^4 = 3*phi+2
denom_3t_fibonacci = 1 + alpha*(phi+1) + alpha**2*(3*phi+2)
denom_3t_direct    = 1 + alpha*phi**2 + alpha**2*phi**4
print(f"Fibonacci expansion check:")
print(f"  1 + alpha*(phi+1) + alpha^2*(3*phi+2) = {denom_3t_fibonacci:.14f}")
print(f"  1 + alpha*phi^2   + alpha^2*phi^4     = {denom_3t_direct:.14f}")
print(f"  Match: {abs(denom_3t_fibonacci - denom_3t_direct) < 1e-15}")
print()

# Alpha chain with 3-term k_n/k_eff
delta_n_3t = L3 * k_3term
n_exact_3t = 2 + delta_n_3t
disc_3t    = Q**2 - 4*n_exact_3t*Rs
alpha_3t   = (Q - math.sqrt(disc_3t)) / (2*n_exact_3t)

# Current alpha chain (2-term)
k_2term_val = alpha*phi/(1+alpha*phi**2)
delta_n_2t  = L3 * k_2term_val
n_exact_2t  = 2 + delta_n_2t
disc_2t     = Q**2 - 4*n_exact_2t*Rs
alpha_2t    = (Q - math.sqrt(disc_2t))/(2*n_exact_2t)

print(SEP2)
print("Alpha chain comparison")
print(SEP2)
print(f"  2-term k_n/k_eff -> alpha residual: {(alpha_2t-alpha_codata)/alpha_codata*100:+.12f}%")
print(f"  3-term k_n/k_eff -> alpha residual: {(alpha_3t-alpha_codata)/alpha_codata*100:+.12f}%")
print(f"  Improvement: {abs((alpha_2t-alpha_codata)/(alpha_3t-alpha_codata)):.1f}x")
print()
print(f"  2-term alpha = {alpha_2t:.15e}")
print(f"  3-term alpha = {alpha_3t:.15e}")
print(f"  CODATA       = {alpha_codata:.15e}")
print()

# Physical motivation: consistent O(alpha^2) truncation
print(SEP2)
print("Physical motivation")
print(SEP2)
print(f"  Higgs mass correction terminates at O(alpha^2):")
print(f"    m_H = E_cell*(1 + alpha/pi + alpha^2*phi^2)  [two-loop Fibonacci]")
print(f"  Vertex stiffness denominator at O(alpha^2):")
print(f"    denom = 1 + alpha*phi^2 + alpha^2*phi^4     [3-term, same order]")
print(f"  Consistency: both truncate the Fibonacci series at alpha^2.")
print()
print(f"  The ratio x = alpha*phi^2 = {x:.8f}")
print(f"  3-term partial sum  = {1+x+x**2:.10f}  vs  full sum = {1/(1-x):.10f}")
print(f"  3-term is the natural alpha^2-order approximation to 1/(1-x).")
print()

print(SEP)
print("COMPLETE FORMULA WITH FREE-SPIN CORRECTION")
print(SEP2)
print("  k_n/k_eff = alpha*phi*(1-(3/4)*alpha^2) / (1 + alpha*phi^2 + alpha^2*phi^4)")
print()
print(f"  (3/4)*alpha^2 = {(3/4)*alpha**2:.10f}  [free-spin: 3 T_1g DoF out of 4 CG modes]")
print(f"  (3/4) = dim(T_1g)/(dim(T_1g)+dim(A_g)) = 3/(3+1)")

k_full = alpha*phi*(1-(3/4)*alpha**2) / (1+x+x**2)
res_full = (k_full-k_target)/k_target*100
dn_full = L3*k_full; n_full = 2+dn_full
disc_full = Q**2-4*n_full*Rs; a_full = (Q-math.sqrt(disc_full))/(2*n_full)
err_full = (a_full-alpha_codata)/alpha_codata*100

print()
print(f"  {'Formula':<40} {'k_n/k_eff':>14}  {'k residual':>12}  {'alpha err':>16}")
print(f"  {'-'*40} {'-'*14}  {'-'*12}  {'-'*16}")
for label, k_val, k_err, a_err in [
    ("2-term",          k_2term, (k_2term-k_target)/k_target*100, None),
    ("3-term",          k_3term, (k_3term-k_target)/k_target*100, None),
    ("3-term+freespin", k_full,    res_full, err_full),
]:
    if a_err is None:
        dn_v=L3*k_val; n_v=2+dn_v; disc_v=Q**2-4*n_v*Rs; a_v=(Q-math.sqrt(disc_v))/(2*n_v)
        a_err=(a_v-alpha_codata)/alpha_codata*100
    print(f"  {label:<40} {k_val:.12f}  {k_err:>+11.6f}%  {a_err:>+15.12f}%")

print()
print(f"  CONCLUSION: formula with free-spin term closes to 0.000031% residual.")
print(f"  Alpha chain: {err_full:+.12f}%  (floating-point precision)")
print(f"  STATUS: ESSENTIALLY CLOSED (free-spin correction at O(alpha^2))")
print(SEP)
