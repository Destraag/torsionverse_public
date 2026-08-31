"""
c4a_candidates.py — Evaluate candidate formulas identified in c4a_theory.txt.

Tests two structural candidates raised in PART 7, QUESTION 1:
  Candidate A: alpha = 2*sqrt(5) / (4*pi*(51-sqrt(5)))
               i.e., Q = (51-sqrt5)/2 = 24 + 1/phi^2
  Candidate B: quadratic equation (shown to be unphysical in the doc,
               but verified here for completeness)

Also checks whether (51-sqrt5)/2 has a plausible geometric origin.

Run: python analysis/c4a_candidates.py
"""

import math

SEP  = "=" * 65
pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)

alpha_CODATA = 7.2973525693e-3
alpha_C4a    = sqrt5 * phi / (16 * pi**3)

print(SEP)
print("CANDIDATE A: Q = (51 - sqrt5) / 2")
print(SEP)
print()

Q_A = (51 - sqrt5) / 2
alpha_A = Rs / Q_A

err_A = (alpha_A - alpha_CODATA) / alpha_CODATA * 100

print(f"  Q_A = (51 - sqrt5)/2 = {Q_A:.10f}")
print(f"  Q_exact (R_s/alpha)  = {Rs/alpha_CODATA:.10f}")
print(f"  Q_A - Q_exact        = {Q_A - Rs/alpha_CODATA:+.8f}")
print()
print(f"  alpha_A = R_s / Q_A  = {alpha_A:.13e}")
print(f"  alpha_CODATA         = {alpha_CODATA:.13e}")
print(f"  error                = {err_A:+.6f}%")
print(f"  1/alpha_A            = {1/alpha_A:.8f}")
print(f"  1/alpha_CODATA       = {1/alpha_CODATA:.8f}")
print()

# What is (51-sqrt5)/2 geometrically?
print(f"  Structural notes on Q_A = (51-sqrt5)/2:")
print(f"    51 = 3 * 17 = ??? (no obvious topological meaning)")
print(f"    51 - sqrt5 = {51 - sqrt5:.8f}")
print(f"    (51-sqrt5)/2 = 24 + (3-sqrt5)/2 = 24 + 1/phi^2")
print(f"      where 1/phi^2 = (3-sqrt5)/2 = {1/phi**2:.8f}  (check: {(3-sqrt5)/2:.8f})")
print(f"    So Q_A = 24 + 1/phi^2")
print(f"    24 = 4! = order of S_4 = vertices of 24-cell (4D polytope)")
print(f"    1/phi^2 = 2 - phi = phi^(-2) = {1/phi**2:.8f}")
print(f"    phi^2 = phi + 1 = {phi**2:.8f}")
print()

# compare to C4a
print(f"  Comparison:")
print(f"    C4a error:   {(alpha_C4a-alpha_CODATA)/alpha_CODATA*100:+.6f}%")
print(f"    Cand A error: {err_A:+.6f}%")
print(f"    Cand A is {'better' if abs(err_A) < abs((alpha_C4a-alpha_CODATA)/alpha_CODATA*100) else 'worse'} than C4a by factor {abs(err_A)/abs((alpha_C4a-alpha_CODATA)/alpha_CODATA*100):.2f}")
print()

print(SEP)
print("CANDIDATE B: Quadratic (verification of doc analysis)")
print(SEP)
print()
print("  Equation: 2*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0")
print()

a_coeff = 2.0
b_coeff = -(4*pi**2/phi)
c_coeff = Rs

discriminant = b_coeff**2 - 4*a_coeff*c_coeff
print(f"  a = {a_coeff}")
print(f"  b = -(4*pi^2/phi) = {b_coeff:.8f}")
print(f"  c = Rs = {c_coeff:.8f}")
print(f"  discriminant = b^2 - 4ac = {discriminant:.8f}")
print(f"  sqrt(discriminant) = {math.sqrt(discriminant):.8f}")
print()

alpha_quad_plus  = (-b_coeff + math.sqrt(discriminant)) / (2*a_coeff)
alpha_quad_minus = (-b_coeff - math.sqrt(discriminant)) / (2*a_coeff)

print(f"  alpha_+ = {alpha_quad_plus:.10e}  (unphysical: >>1)")
print(f"  alpha_- = {alpha_quad_minus:.13e}")
print(f"  CODATA   = {alpha_CODATA:.13e}")
err_quad = (alpha_quad_minus - alpha_CODATA) / alpha_CODATA * 100
print(f"  error    = {err_quad:+.6f}%")
print(f"  1/alpha_- = {1/alpha_quad_minus:.8f}  (CODATA: {1/alpha_CODATA:.8f})")
print()
if abs(err_quad) < abs((alpha_C4a-alpha_CODATA)/alpha_CODATA*100):
    print(f"  *** QUADRATIC IS MORE ACCURATE THAN C4a ({abs(err_quad):.4f}% vs {abs((alpha_C4a-alpha_CODATA)/alpha_CODATA*100):.4f}%) ***")
    print(f"  Physical interpretation: alpha satisfies a self-consistent equation")
    print(f"  alpha*(4*pi^2/phi - 2*alpha) = R_s, meaning the geometric denominator")
    print(f"  Q = 4*pi^2/phi - 2*alpha is itself modified by the EM coupling.")
    print(f"  Analogous to QED self-energy: the bare Hopf ratio 4*pi^2/phi is")
    print(f"  corrected by -2*alpha (the leading self-interaction term).")
else:
    print(f"  Conclusion: quadratic with Q=4*pi^2/phi - 2*alpha does not improve on C4a.")
print()

print(SEP)
print("FURTHER CANDIDATES — SYSTEMATIC Q SEARCH NEAR 24.384")
print(SEP)
print()
print("  Looking for geometric expressions Q with |Q - R_s/alpha| < 0.01")
print()

target_Q = Rs / alpha_CODATA
print(f"  Target Q = R_s/alpha = {target_Q:.10f}")
print()

# Generate candidates from combinations of small integers, pi, phi, sqrt5
candidates = []

# Form: a + b/phi^n or a + b*phi^n or similar
for a in range(20, 30):
    for b in [-2,-1,1,2]:
        for n in [1,2,3]:
            val = a + b/phi**n
            if abs(val - target_Q) < 0.02:
                name = f"{a} + {b}/phi^{n}"
                candidates.append((name, val))
            val2 = a + b*phi**n
            if abs(val2 - target_Q) < 0.02:
                name2 = f"{a} + {b}*phi^{n}"
                candidates.append((name2, val2))

# Form: n*pi^m or n/pi^m
for n in range(1, 80):
    for m in [1, 2]:
        val = n / pi**m
        if abs(val - target_Q) < 0.02:
            name = f"{n}/pi^{m}"
            candidates.append((name, val))
        val2 = n * pi / m
        if abs(val2 - target_Q) < 0.02:
            name2 = f"{n}*pi/{m}"
            candidates.append((name2, val2))

# Form: n*sqrt5 + m
for n in [-3,-2,-1,1,2,3]:
    for m in range(15, 30):
        val = n*sqrt5 + m
        if abs(val - target_Q) < 0.02:
            name = f"{n}*sqrt5+{m}"
            candidates.append((name, val))

# Form involving Rs
for n in [1,2,4,8]:
    val = target_Q  # placeholder
    # Rs-based
    for k in range(130,145):
        val = k * Rs
        if abs(val - target_Q) < 0.02:
            name = f"{k}*Rs"
            candidates.append((name, val))

# Remove duplicates and sort by error
seen = set()
unique = []
for name, val in candidates:
    err = abs(val - target_Q)
    key = round(val, 8)
    if key not in seen:
        seen.add(key)
        unique.append((err, name, val))
unique.sort()

print(f"  {'Expression':<30} {'Value':>14}  {'Error from target':>18}")
print(f"  {'-'*30} {'-'*14}  {'-'*18}")
for err, name, val in unique[:20]:
    alpha_implied = Rs / val if val != 0 else 0
    alpha_err = (alpha_implied - alpha_CODATA)/alpha_CODATA*100 if val != 0 else 999
    marker = " <-- BETTER THAN C4a" if abs(alpha_err) < abs((alpha_C4a-alpha_CODATA)/alpha_CODATA*100) else ""
    print(f"  {name:<30} {val:>14.8f}  {err:>+14.8f}  alpha_err={alpha_err:+.4f}%{marker}")
print()

print(SEP)
print("SUMMARY")
print(SEP)
print()
print(f"  C4a:        Q = 4*pi^2/phi     = {4*pi**2/phi:.8f}  alpha_err = {(alpha_C4a-alpha_CODATA)/alpha_CODATA*100:+.6f}%")
print(f"  Cand A:     Q = 24 + 1/phi^2   = {24+1/phi**2:.8f}  alpha_err = {err_A:+.6f}%")
print()
if abs(err_A) < abs((alpha_C4a-alpha_CODATA)/alpha_CODATA*100):
    print(f"  Candidate A is MORE ACCURATE than C4a.")
    print(f"  However: '24 + 1/phi^2' is less geometrically natural than '4*pi^2/phi'.")
    print(f"  The integer 24 = 4! requires explanation. The sum form suggests")
    print(f"  two separate geometric contributions rather than one unified expression.")
    print(f"  C4a's '4*pi^2/phi' is a single ratio from one geometric object.")
    print(f"  Candidate A may be a better numerical approximation but is")
    print(f"  not obviously a better physical derivation candidate.")
else:
    print(f"  C4a remains more accurate than Candidate A.")
print()
print(f"  Both require Step 3 (Hopf linking integral) to become derivations.")
