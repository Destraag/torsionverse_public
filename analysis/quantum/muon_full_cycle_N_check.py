"""
muon_full_cycle_N_check.py

Tests whether using the muon's FULL periodic cycle length (not just one
6-vertex circuit) as N in the ring self-consistency correction closes the
gap that N<=30 could not (muon_ring_selfconsistency_check.py). The muon
does not repeat the same 6 vertices -- per muon_orbit_count.py, the
minimal I_h-symmetric pattern is Orbit A: 10 distinct 6-vertex circuits,
visiting all 12 vertices exactly 5 times each (10*6 = 60 = 12*5) before
the full symmetric pattern repeats.

APPROACH:
  1. Solve INVERSELY: what N does eff = phi/(1+N*alpha) need to exactly
     reproduce each existing candidate eff_mu value?
  2. Compare those required-N values to the actual, already-established
     framework numbers: 6 (one circuit), 60 (full Orbit A cycle), 70
     (total distinct circuits), 30 (total edges), 12 (total vertices).

Run: python analysis/quantum/muon_full_cycle_N_check.py
"""
import math

SEP = "=" * 70
phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3

print(SEP)
print("SOLVING FOR N: eff = phi/(1+N*alpha) = target")
print(SEP)
print()

targets = {
    "bipyramid (used in formula)": (9 - math.sqrt(5)) / 8,
    "PDG-inverted exact": 0.856308161568,
    "trivial real-path collapse": 1.0,
}

framework_numbers = {
    "6 (one muon circuit)": 6,
    "12 (all vertices)": 12,
    "20 (all faces)": 20,
    "30 (all edges)": 30,
    "60 (Orbit A full cycle: 10 circuits x 6)": 60,
    "70 (total distinct muon circuits)": 70,
}

for name, target in targets.items():
    # phi/(1+N*alpha) = target  =>  N = (phi/target - 1)/alpha
    N_needed = (phi / target - 1) / alpha
    print(f"  Target eff = {target:.6f}  ({name})")
    print(f"    N needed = {N_needed:.4f}")
    closest_name, closest_val = min(framework_numbers.items(), key=lambda kv: abs(kv[1] - N_needed))
    print(f"    Closest framework number: {closest_name}  (diff = {N_needed - closest_val:+.4f})")
    print()

print(SEP)
print("DIRECT TEST: N=60 (full Orbit A cycle) in each candidate formula form")
print(SEP)
N = 60
eff_a = phi / (1 + (N - 1) * alpha)
eff_c = phi / (1 + N * alpha)
print(f"  N=60:  eff_a (N-1) = {eff_a:.6f}   eff_c (N) = {eff_c:.6f}")
print(f"  Compare to bipyramid=0.845492, PDG-exact=0.856308")
print()

# Also test N=60 in the FULL mass formula directly, not just as isolated eff
def mass_formula(eff_mu):
    log5 = math.log(5)
    Rs2 = (math.sqrt(5) / (4 * math.pi)) ** 2
    m_p = 938.272046
    poly = 5 * math.tan(math.pi / 5)
    L3 = (eff_mu**3 + log5**3) / (eff_mu**2 + log5**2)
    x = alpha * eff_mu**2
    k = alpha * eff_mu * (1 - 0.75 * alpha**2) / (1 + x + x**2)
    dn = L3 * k
    base = 2 * math.pi * alpha * (2 / math.sqrt(5)) * phi**2 * m_p
    corr = 1 + Rs2 + 2 * alpha
    return base * (1 + dn / poly) * corr

m_mu_pdg = 105.6583755
for label, eff in [("eff_a (N=60, N-1 form)", eff_a), ("eff_c (N=60, N form)", eff_c)]:
    m_pred = mass_formula(eff)
    err = (m_pred - m_mu_pdg) / m_mu_pdg * 100
    print(f"  {label}: eff={eff:.6f}  ->  m_mu = {m_pred:.4f} MeV  (err = {err:+.4f}%)")
