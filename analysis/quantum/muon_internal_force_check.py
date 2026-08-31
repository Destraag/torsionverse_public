#!/usr/bin/env python3
"""
muon_internal_force_check.py

Internal consistency check (session 12): the VERIFIED muon mass formula
(lepton_mass.py LM6-LM8, -0.003% accuracy) uses eff_mu = (9-sqrt5)/8, built
from a ratio of TWO DIFFERENT deflection magnitudes (bipyramid apex vs
equatorial). But the framework's OWN later check (LM4b) shows the REAL
icosahedral zigzag path has ALL 6 deflections uniformly = cos(72 deg) --
no apex/equator split at all. This script recomputes the force-ratio using
the REAL path's deflection values instead of the bipyramid's, and checks
whether the resulting mass prediction is still plausible.

Reference: analysis/quantum/lepton_mass.py (Section 2-3, LM3-LM8).
"""
import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 66
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
log5 = math.log(5)
Rs2 = 5 / (16*pi**2)
m_p = 938.272046
m_mu_pdg = 105.6583755

print(SEP)
print("INTERNAL CONSISTENCY CHECK: eff_mu FROM REAL PATH vs BIPYRAMID")
print(SEP)

# ── Real icosahedral zigzag path (reproduces LM4b exactly) ──────────────────
r_e_ico = 1/(2*math.sin(pi/5))
r_plane_ico = 2*r_e_ico*math.sin(pi/10)
z_u_ico = math.sqrt(1 - r_plane_ico**2)/2
h_top_ico = z_u_ico + math.sqrt(1 - r_e_ico**2)
top_ico = (0, 0, h_top_ico)
bot_ico = (0, 0, -h_top_ico)
upper_ico = [(r_e_ico*math.cos(2*pi*k/5), r_e_ico*math.sin(2*pi*k/5), z_u_ico) for k in range(5)]
lower_ico = [(r_e_ico*math.cos(2*pi*k/5+pi/5), r_e_ico*math.sin(2*pi*k/5+pi/5), -z_u_ico) for k in range(5)]

def sub(a, b): return tuple(a[i]-b[i] for i in range(3))
def norm(v): return math.sqrt(sum(c*c for c in v))
def unit(v):
    n = norm(v); return tuple(c/n for c in v)
def dot(a, b): return sum(a[i]*b[i] for i in range(3))

def deflect(vin, v, vout):
    d_in = unit(sub(v, vin)); d_out = unit(sub(vout, v))
    return dot(d_in, d_out)

path_ico = [top_ico, upper_ico[0], lower_ico[0], bot_ico, lower_ico[2], upper_ico[2], top_ico]
deflections = [deflect(path_ico[i-1], path_ico[i], path_ico[i+1]) for i in range(1, 6)]

print(f"\n  Real path deflections (all 5 interior vertices): {[round(d,6) for d in deflections]}")
check("F1: real icosahedral path has ALL deflections = cos(72 deg) = 1/(2*phi), reproducing LM4b",
      all(abs(d - 1/(2*phi)) < 1e-8 for d in deflections),
      f"all = {1/(2*phi):.6f}")

# ── eff_mu recomputed from the REAL path (no apex/equator split exists) ─────
print()
print("RECOMPUTING eff_mu FROM THE REAL PATH (not the bipyramid)")
print("-"*66)
cos_apex_real = deflections[0]   # any of them -- all identical
cos_eq_real = deflections[0]     # SAME value -- there is no second angle type
ratio_real = abs(cos_apex_real) / abs(cos_eq_real)
eff_mu_real = (1 + ratio_real) / 2

print(f"  Real path has only ONE deflection type: {cos_apex_real:.8f}")
print(f"  ratio = |cos_apex|/|cos_eq| = {ratio_real:.8f}  (trivially 1, since both are the same number)")
print(f"  eff_mu_real = (1+ratio)/2 = {eff_mu_real:.8f}")
print(f"  ORIGINAL (bipyramid) eff_mu = (9-sqrt5)/8 = {(9-math.sqrt(5))/8:.8f}")

check("F2: eff_mu collapses to a TRIVIAL value (1.0) when using the real uniform-deflection path",
      abs(eff_mu_real - 1.0) < 1e-9, f"eff_mu_real = {eff_mu_real:.8f}")

# ── Recompute the muon mass using eff_mu_real instead ───────────────────────
print()
print("RECOMPUTING MUON MASS WITH eff_mu_real INSTEAD OF THE BIPYRAMID VALUE")
print("-"*66)

def mass_from_eff(eff_mu, poly_norm):
    L3 = (eff_mu**3 + log5**3) / (eff_mu**2 + log5**2)
    x = alpha * eff_mu**2
    k = alpha * eff_mu * (1 - 0.75*alpha**2) / (1 + x + x**2)
    dn = L3 * k
    base = 2*pi*alpha*(2/math.sqrt(5))*phi**2*m_p
    corr = 1 + Rs2 + 2*alpha
    return base * (1 + dn/poly_norm) * corr

eff_mu_bip = (9 - math.sqrt(5)) / 8
poly_mu = 5*math.tan(pi/5)   # same polygon normalization used in lepton_mass.py

m_mu_bip  = mass_from_eff(eff_mu_bip, poly_mu)
m_mu_real = mass_from_eff(eff_mu_real, poly_mu)

err_bip  = (m_mu_bip  - m_mu_pdg) / m_mu_pdg * 100
err_real = (m_mu_real - m_mu_pdg) / m_mu_pdg * 100

print(f"  Bipyramid eff_mu = {eff_mu_bip:.8f}  ->  m_mu = {m_mu_bip:.4f} MeV  (err = {err_bip:+.4f}%)")
print(f"  Real-path eff_mu = {eff_mu_real:.8f}  ->  m_mu = {m_mu_real:.4f} MeV  (err = {err_real:+.4f}%)")
print(f"  PDG m_mu = {m_mu_pdg} MeV")

plausible = abs(err_real) < 1.0   # within 1% -- generous "still plausible" bar
check("F3: using the REAL path's eff_mu still gives a mass within 1% of PDG (plausible)",
      plausible, f"err_real = {err_real:+.4f}%  (bar: |err| < 1%)")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
for name, status, detail in results:
    print(f"  {status}: {name}")
print(SEP)
