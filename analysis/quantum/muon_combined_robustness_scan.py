"""
muon_combined_robustness_scan.py

Robustness check for eff_mu = (1/phi) + delta*(mass_ratio)^p, found earlier
at delta=alpha, p=2/3 giving +0.0046% error (comparable to bipyramid's
-0.0029%). Scans a grid around that point to see whether it's a genuine,
gently-varying minimum (robust) or a narrow, fragile coincidence (not
robust) -- requested explicitly before treating it as a real lead.

Run: python analysis/quantum/muon_combined_robustness_scan.py
"""
import math

phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
log5 = math.log(5)
Rs2 = (math.sqrt(5) / (4 * math.pi)) ** 2
m_p = 938.272046
m_e_pdg = 0.51099895
m_mu_pdg = 105.6583755
poly = 5 * math.tan(math.pi / 5)
mass_ratio = m_mu_pdg / m_e_pdg

def mass_formula(eff_mu):
    L3 = (eff_mu**3 + log5**3) / (eff_mu**2 + log5**2)
    x = alpha * eff_mu**2
    k = alpha * eff_mu * (1 - 0.75 * alpha**2) / (1 + x + x**2)
    dn = L3 * k
    base = 2 * math.pi * alpha * (2 / math.sqrt(5)) * phi**2 * m_p
    corr = 1 + Rs2 + 2 * alpha
    return base * (1 + dn / poly) * corr

def err_for(delta, p):
    eff = 1/phi + delta * mass_ratio**p
    if eff <= 0:
        return None
    return (mass_formula(eff) - m_mu_pdg) / m_mu_pdg * 100

print("=" * 74)
print("ROBUSTNESS SCAN: eff_mu = 1/phi + delta*(mass_ratio)^p")
print("=" * 74)
print()
print("PART A: vary p (exponent), delta fixed at alpha")
print(f"  {'p':>8s} {'mass_ratio^p':>14s} {'eff_mu':>10s} {'error %':>10s}")
for p in [0.55, 0.60, 0.6667, 0.70, 0.75, 0.80]:
    eff = 1/phi + alpha * mass_ratio**p
    err = err_for(alpha, p)
    print(f"  {p:8.4f} {mass_ratio**p:14.4f} {eff:10.6f} {err:+10.4f}")
print()

print("PART B: vary delta (as multiple of alpha), p fixed at 2/3")
print(f"  {'delta/alpha':>12s} {'eff_mu':>10s} {'error %':>10s}")
for mult in [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]:
    err = err_for(alpha * mult, 2/3)
    eff = 1/phi + alpha*mult * mass_ratio**(2/3)
    print(f"  {mult:12.2f} {eff:10.6f} {err:+10.4f}")
print()

print("PART C: 2D grid -- how large is the 'good' region (|error| < 0.01%)?")
print(f"  {'p':>8s}", end="")
mults = [0.90, 0.95, 1.00, 1.05, 1.10]
for m in mults:
    print(f" {'delta='+str(m)+'a':>12s}", end="")
print()
count_good = 0
count_total = 0
for p in [0.55, 0.60, 0.6667, 0.70, 0.75, 0.80]:
    print(f"  {p:8.4f}", end="")
    for mult in mults:
        err = err_for(alpha * mult, p)
        count_total += 1
        marker = "*" if err is not None and abs(err) < 0.01 else " "
        if marker == "*":
            count_good += 1
        print(f" {err:+11.4f}{marker}", end="")
    print()

print()
print(f"Grid cells with |error| < 0.01%: {count_good}/{count_total}")
print()
print("=" * 74)
print(f"Reference: bipyramid error = {(mass_formula((9-math.sqrt(5))/8)-m_mu_pdg)/m_mu_pdg*100:+.4f}%")
print(f"Reference: PDG-inverted exact eff needed = 0.856308161568")
