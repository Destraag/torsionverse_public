"""
higgs_pq_uniqueness_exhaustive.py
=================================
Turns doc_jobson_cell.txt Section 7's "(2) Corroborating scan (not exhaustive):
of the 8 small (p,q) pairs checked..." into a genuinely EXHAUSTIVE proof over
the full infinite (p,q) space, using jobson_cell_doc.py's own closed-form
E_cell(p,q) = 2*pi*hbar_c / (alpha*phi_pq*r_p) formula, where
phi_pq = (1+sqrt(p^2+q^2))/2.

KEY FACT: phi_pq is STRICTLY INCREASING in p^2+q^2, so E_cell(p,q) = K/phi_pq
(K a p,q-independent constant) is STRICTLY, MONOTONICALLY DECREASING in
p^2+q^2. This turns "check infinitely many (p,q)" into "find the finite norm
bound where E_cell first drops below the 125 GeV window, exhaustively
enumerate every coprime pair inside that bound (a small finite lattice-point
count), and invoke monotonicity for everything beyond it" -- a complete,
rigorous proof, not a spot-check.

CHECKS:
  PQ1: E_cell(p,q) is exactly K/phi_pq (algebraic form matches numeric eval)
  PQ2: phi_pq is strictly increasing in p^2+q^2 (monotonicity, sampled)
  PQ3: Norm bound N_max where E_cell(p,q) permanently drops below 120 GeV
  PQ4: Every coprime (p,q), p<=q, with p^2+q^2 <= N_max is enumerated
  PQ5: (1,2) is the ONLY scalar (n=p*q even) pair in the window [120,130] GeV
       among ALL pairs with p^2+q^2 <= N_max
  PQ6: No pair beyond N_max can re-enter the window (monotonicity + bound)

Run: python analysis/higgs/higgs_pq_uniqueness_exhaustive.py
Reference: docs/series1/doc_jobson_cell.txt Section 7 argument (2);
           analysis/demos/jobson_cell_doc.py E_cell_pq() (J20/J21)
"""
import sys, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 70
SEP2 = "-" * 70
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

pi     = math.pi
alpha  = 7.2973525693e-3        # CODATA 2018
r_p_fm = 0.8414                 # fm, CODATA 2018 (matches jobson_cell_doc.py chain)
hbar_c = 197.3269804            # MeV*fm

def phi_pq(p, q):
    return (1 + math.sqrt(p**2 + q**2)) / 2

def E_cell_pq(p, q):
    """Same closed form as jobson_cell_doc.py's E_cell_pq()."""
    L_J_pq = alpha * phi_pq(p, q) * r_p_fm   # fm (r_p_fm already in fm, no unit conversion needed)
    return 2 * pi * hbar_c / L_J_pq / 1000   # GeV

# K = the p,q-independent constant: E_cell(p,q) = K / phi_pq(p,q)
K = E_cell_pq(1, 1) * phi_pq(1, 1)

print(SEP)
print("SECTION 1: CLOSED-FORM ALGEBRAIC IDENTITY")
print(SEP2)
print(f"  E_cell(p,q) = K / phi_pq(p,q),  K = {K:.6f} GeV  (p,q-independent)")
for p, q in [(1,1), (1,2), (2,3), (3,5)]:
    lhs = E_cell_pq(p, q)
    rhs = K / phi_pq(p, q)
    print(f"    ({p},{q}): E_cell={lhs:.4f}  K/phi_pq={rhs:.4f}")

check("PQ1: E_cell(p,q) = K/phi_pq(p,q) exactly, for all tested pairs",
      all(abs(E_cell_pq(p, q) - K/phi_pq(p, q)) < 1e-9 for p, q in
          [(1,1),(1,2),(1,3),(2,3),(1,4),(2,5),(3,5),(1,5),(2,7),(5,8),(10,17)]),
      "checked 11 pairs spanning small to large norm")

print()
print(SEP)
print("SECTION 2: MONOTONICITY (phi_pq strictly increasing in p^2+q^2)")
print(SEP2)
norms = [1, 2, 5, 8, 10, 13, 17, 25, 50, 100, 1000, 100000]
phis = [(1 + math.sqrt(n)) / 2 for n in norms]
strictly_increasing = all(phis[i] < phis[i+1] for i in range(len(phis)-1))
print(f"  norms tested: {norms}")
print(f"  phi_pq values: {[f'{p:.4f}' for p in phis]}")
print(f"  Algebraically obvious: d(phi_pq)/d(norm) = 1/(4*sqrt(norm)) > 0 for all norm>0.")
check("PQ2: phi_pq strictly increasing in p^2+q^2 (hence E_cell strictly decreasing)",
      strictly_increasing,
      "phi_pq = (1+sqrt(N))/2 has positive derivative in N for all N>0 -- proven algebraically, sampled numerically")

print()
print(SEP)
print("SECTION 3: FINITE NORM BOUND FOR THE 125 GeV WINDOW [120, 130] GeV")
print(SEP2)
WINDOW_LO, WINDOW_HI = 120.0, 130.0
# Since E_cell is strictly decreasing in norm, find the norm where E_cell = WINDOW_LO.
# E_cell = K/phi_pq = WINDOW_LO  =>  phi_pq = K/WINDOW_LO  =>  norm = (2*phi_pq-1)^2
phi_at_lo = K / WINDOW_LO
N_max = (2 * phi_at_lo - 1) ** 2
print(f"  Solve E_cell(norm) = {WINDOW_LO} GeV for norm (p^2+q^2):")
print(f"    phi_pq required = K/{WINDOW_LO} = {phi_at_lo:.6f}")
print(f"    norm = (2*phi_pq-1)^2 = {N_max:.6f}")
print(f"  For ANY (p,q) with p^2+q^2 > {N_max:.4f}: E_cell(p,q) < {WINDOW_LO} GeV, PERMANENTLY")
print(f"  (by PQ2 monotonicity -- it can never increase back into the window).")
print(f"  => Only coprime (p,q) with p^2+q^2 <= {math.floor(N_max)} need to be checked. FINITE SET.")

check("PQ3: N_max is finite and small (bounds the search to a tiny lattice region)",
      0 < N_max < 20,
      f"N_max = {N_max:.4f} -- only integer p^2+q^2 up to {math.floor(N_max)} matter")

print()
print(SEP)
print("SECTION 4: EXHAUSTIVE ENUMERATION OF EVERY COPRIME (p,q), p<=q, p^2+q^2<=N_max")
print(SEP2)

def gcd(a, b):
    while b: a, b = b, a % b
    return a

N_MAX_INT = math.floor(N_max)
all_pairs = []
for total in range(1, N_MAX_INT + 1):
    for p in range(1, int(math.isqrt(total)) + 2):
        q2 = total - p*p
        if q2 < p*p:  # enforce p<=q (avoid double-counting the same knot)
            continue
        q = math.isqrt(q2)
        if q*q != q2:
            continue
        if gcd(p, q) != 1:
            continue
        all_pairs.append((p, q))

print(f"  Every coprime (p,q) with p<=q and p^2+q^2 <= {N_MAX_INT}:")
scalars_in_window = []
for p, q in sorted(set(all_pairs)):
    E = E_cell_pq(p, q)
    n = p * q
    typ = "scalar" if n % 2 == 0 else "vector"
    in_window = WINDOW_LO <= E <= WINDOW_HI
    marker = "  <-- IN WINDOW" if in_window else ""
    print(f"    ({p},{q}): p^2+q^2={p*p+q*q}, n={n}, E_cell={E:.3f} GeV [{typ}]{marker}")
    if typ == "scalar" and in_window:
        scalars_in_window.append((p, q, E))

check("PQ4: exhaustive enumeration found >=1 pair (sanity: search space non-empty)",
      len(all_pairs) >= 1, f"{len(set(all_pairs))} distinct coprime pairs with p^2+q^2<={N_MAX_INT}")

check("PQ5: (1,2) is the ONLY scalar pair in [120,130] GeV among ALL pairs with p^2+q^2<=N_max",
      len(scalars_in_window) == 1 and scalars_in_window[0][:2] == (1, 2),
      f"scalars in window: {scalars_in_window}")

print()
print(SEP)
print("SECTION 5: TAIL BEYOND N_max (monotonicity closes the infinite remainder)")
print(SEP2)
# Spot-check several large-norm pairs beyond N_max to confirm they stay below the window
large_pairs = [(1, 8), (3, 7), (5, 12), (10, 20), (100, 101)]
tail_below_window = []
for p, q in large_pairs:
    E = E_cell_pq(p, q)
    tail_below_window.append(E < WINDOW_LO)
    print(f"    ({p},{q}): p^2+q^2={p*p+q*q}, E_cell={E:.4f} GeV  (< {WINDOW_LO}? {E < WINDOW_LO})")

check("PQ6: sampled large-norm pairs beyond N_max all fall below the window (confirms monotonic tail)",
      all(tail_below_window),
      f"{sum(tail_below_window)}/{len(tail_below_window)} sampled pairs below {WINDOW_LO} GeV")

print()
print(SEP)
print("CONCLUSION")
print(SEP2)
print(f"  E_cell(p,q) = K/phi_pq(p,q) is a closed-form, strictly monotonically")
print(f"  decreasing function of p^2+q^2 (PQ1-PQ2). This means the infinite (p,q)")
print(f"  space splits into: a FINITE region (p^2+q^2 <= {N_MAX_INT}) that can be")
print(f"  exhaustively enumerated by brute force (PQ3-PQ5), and an infinite tail")
print(f"  (p^2+q^2 > {N_MAX_INT}) that is PROVABLY, PERMANENTLY below the 125 GeV")
print(f"  window by monotonicity alone (PQ6) -- no further checking is needed there.")
print(f"  (1,2) is the unique scalar pair in the window across the ENTIRE infinite")
print(f"  (p,q) space, not just the {N_MAX_INT if False else 8} pairs previously spot-checked.")
print(f"  This upgrades doc_jobson_cell.txt Section 7 argument (2) from a")
print(f"  'corroborating, not exhaustive' scan to a genuinely exhaustive proof.")

print()
print(SEP)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"SUMMARY: {passed}/{len(results)} PASS, {failed} FAIL")
if failed == 0:
    print("ALL CHECKS PASSED.")
print(SEP)
