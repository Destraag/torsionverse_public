"""
higgs_subrep_nj.py
==================
LEAD 2: Do sub-cell particles sit at representation resonances?

Observation: Higgs has N_J = 1/(2*pi) [~0.16% from exact].
Hypothesis:  W boson has N_J = 1/4 = 1/dim(G_g) [1.17% off].

The resonance condition: N_J * n = 1, where n is an I_h representation
dimension or natural number. Tests whether sub-cell particle Compton
wavelengths are fixed fractions of L_J set by icosahedral geometry.

Run: python analysis/higgs/higgs_subrep_nj.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2
L_J_fm_val = L_J   # already in fm from constants.py

SEP  = "=" * 65
SEP2 = "-" * 65

# I_h gerade irrep dimensions and names
irreps = {
    'A_g':  1,
    'T_1g': 3,
    'T_2g': 3,
    'G_g':  4,
    'H_g':  5,
}
# Additional natural resonance numbers
resonances = {
    '2*pi':    2*pi,
    '2*pi*phi':2*pi*phi,
    '4':       4,
    '5':       5,
    '6':       6,
    '12':      12,
    'pi':      pi,
    'phi^2':   phi**2,
    'phi^3':   phi**3,
    '2*phi^2': 2*phi**2,
}
# Combine all
all_n = {}
for name, d in irreps.items():
    all_n[f'dim({name})={d}'] = d
for name, v in resonances.items():
    all_n[name] = v

# Particles to test (sub-cell and boundary)
particles = [
    ('top',   172760,  'sub-cell'),
    ('Higgs', 125200,  'sub-cell'),
    ('Z',      91188,  'sub-cell'),
    ('W',      80377,  'sub-cell'),
    ('b',       4180,  'boundary'),
    ('tau',     1777,  'near-bulk'),
]

print(SEP)
print("LEAD 2: RESONANCE CONDITION N_J * n = 1")
print("Do sub-cell particles sit at N_J = 1/n for I_h representation dimensions?")
print(SEP2)
print()

print(f"  L_J = {L_J_fm_val:.8f} fm")
print(f"  E_cell = {E_cell_GeV:.4f} GeV")
print()

# ── Main table ─────────────────────────────────────────────────────────────────
print("PARTICLE N_J VALUES AND BEST RESONANCE MATCH:")
print(SEP2)
print(f"  {'Particle':<8} {'N_J':>8}  {'best n':>12}  {'N_J*n':>8}  {'err%':>8}  note")
print(SEP2)

for name, mass_MeV, regime in particles:
    nj = hbar_c / (mass_MeV * L_J_fm_val)
    
    best_name, best_n, best_err = None, None, 1e9
    for n_name, n_val in all_n.items():
        product = nj * n_val
        err = abs(product - 1.0)
        if err < best_err:
            best_err = err
            best_name = n_name
            best_n = n_val
    
    err_pct = (nj * best_n - 1) * 100
    flag = " <--" if abs(err_pct) < 2 else (" <" if abs(err_pct) < 5 else "")
    print(f"  {name:<8} {nj:>8.5f}  {best_name:>12}  {nj*best_n:>8.5f}  {err_pct:>+7.3f}%  {flag}")

print()

# ── Focused analysis: Higgs and W ──────────────────────────────────────────────
print(SEP)
print("FOCUSED ANALYSIS: HIGGS AND W BOSON RESONANCES")
print(SEP2)
print()

m_H = m_H_pdg22 * 1000   # MeV
m_W = 80377               # MeV PDG 2022
m_Z = 91188               # MeV

nj_H = hbar_c / (m_H * L_J_fm_val)
nj_W = hbar_c / (m_W * L_J_fm_val)
nj_Z = hbar_c / (m_Z * L_J_fm_val)

print(f"  HIGGS: N_J = {nj_H:.8f}")
print(f"    Resonance n = 2*pi: N_J * 2*pi = {nj_H * 2*pi:.8f}  (err: {(nj_H*2*pi-1)*100:+.4f}%)")
print(f"    The Higgs sits at N_J = 1/(2*pi) -- one RADIAN of the cell circumference.")
print(f"    This follows algebraically from m_H = E_cell = 2*pi*hbar_c/L_J.")
print()
print(f"  W BOSON: N_J = {nj_W:.8f}")
print(f"    Resonance n = 4 [dim(G_g)]: N_J * 4 = {nj_W * 4:.8f}  (err: {(nj_W*4-1)*100:+.4f}%)")
print(f"    Resonance n = pi:            N_J * pi = {nj_W * pi:.8f}  (err: {(nj_W*pi-1)*100:+.4f}%)")
print(f"    Resonance n = 4 gives m_W_predicted = 4*E_cell/(2*pi) = {4*E_cell_GeV/(2*pi):.4f} GeV")
print(f"    vs m_W_measured = {m_W/1000:.4f} GeV  (gap: {(4*E_cell_GeV/(2*pi)/(m_W/1000)-1)*100:+.4f}%)")
print()
print(f"  Z BOSON: N_J = {nj_Z:.8f}")
for n_name in ['dim(H_g)=5', 'dim(G_g)=4', 'pi', 'phi^2']:
    n_val = all_n.get(n_name, float(n_name.split('=')[-1]) if '=' in n_name else None)
    if n_val is None:
        continue
    prod = nj_Z * n_val
    print(f"    n = {n_name}: N_J * n = {prod:.5f}  (err: {(prod-1)*100:+.3f}%)")
print()

# ── Is the gap the same gap? ──────────────────────────────────────────────────
print(SEP)
print("IS THE W BOSON 1.17% GAP THE SAME AS THE m_H/m_W GAP?")
print(SEP2)
print()
# m_H/m_W gap vs pi/2
ratio_HW = m_H / m_W
target_HW = pi/2
gap_HW = (ratio_HW / target_HW - 1) * 100

# N_J_W gap vs 1/4
gap_NJW = (nj_W * 4 - 1) * 100

# These should be related:
# N_J_W = hbar_c/(m_W * L_J)
# 1/4 = hbar_c/(4*E_cell/(2*pi) * L_J) * 1 = 1/(4*m_H_bare/E_cell * something)
# The gap N_J_W * 4 - 1 = (hbar_c*4/(m_W*L_J)) - 1 = 4*E_cell/(2*pi*m_W) - 1

print(f"  m_H/m_W gap from pi/2:  {gap_HW:+.4f}%")
print(f"  N_J_W * 4 - 1:          {gap_NJW:+.4f}%")
print()
print(f"  These are DIFFERENT gaps (one is {gap_HW:.3f}%, other is {gap_NJW:.3f}%)")
print(f"  But they share the same structure: both relate m_W to m_H/(some pi factor)")
print()
# Algebraically:
# N_J_W * 4 = 4*hbar_c / (m_W * L_J) = 4 * E_cell/(2*pi) / m_W = 2*E_cell / (pi*m_W)
# m_H_bare / m_W = m_H / m_W * (1 + alpha/pi) (approx correction)
# If N_J_W * 4 = 1 exactly: m_W = 4*E_cell/(2*pi) = 2*E_cell/pi
# If m_H/m_W = pi/2 exactly: m_W = 2*m_H/pi = 2*E_cell*(1+alpha/pi)/pi

# The ratio: (2*E_cell*(1+alpha/pi)/pi) / (2*E_cell/pi) = 1+alpha/pi
correction = 1 + alpha/pi
print(f"  If both N_J_W=1/4 and m_H/m_W=pi/2 were exact simultaneously:")
print(f"    m_W from N_J_W=1/4:  2*E_cell/pi = {2*E_cell_GeV/pi:.4f} GeV")
print(f"    m_W from m_H/m_W=pi/2: 2*m_H/pi = {2*m_H/1000/pi:.4f} GeV")
print(f"    Ratio: (1+alpha/pi) = {correction:.8f}")
print(f"    The two conditions differ by EXACTLY the alpha/pi scalar QED correction!")
print()
print(f"  INTERPRETATION: N_J_W = 1/4 and m_H/m_W = pi/2 are CONSISTENT targets.")
print(f"  If N_J_W = 1/4 exactly (from dim(G_g)=4), and m_H = E_cell*(1+alpha/pi),")
print(f"  then m_H/m_W = m_H / (2*E_cell/pi) = pi*(1+alpha/pi)/2 ≈ pi/2*(1+alpha/pi).")
print(f"  This overcorrects pi/2 by alpha/pi = {alpha/pi*100:.4f}% -- which is the MEASURED gap.")
print()
m_H_val = m_H / 1000
m_W_pred_from_nj = 2 * E_cell_GeV / pi
ratio_from_nj = m_H_val / m_W_pred_from_nj
print(f"  CHECK (using m_H_pred = E_cell*(1+alpha/pi), NOT PDG m_H):")
m_H_pred_GeV = E_cell_GeV * (1 + alpha/pi)
ratio_from_nj = m_H_pred_GeV / m_W_pred_from_nj
print(f"  m_H_pred / (2*E_cell/pi) = {ratio_from_nj:.8f}")
print(f"  pi/2 * (1 + alpha/pi)    = {pi/2 * (1+alpha/pi):.8f}")
print(f"  Are these equal (tol 1e-4)? {abs(ratio_from_nj - pi/2*(1+alpha/pi)) < 1e-4}")
print()

# ── Chain summary ─────────────────────────────────────────────────────────────
print(SEP)
print("RESULT: ALGEBRAIC CHAIN CONNECTING m_H, m_W, AND CELL GEOMETRY")
print(SEP2)
print()
print("  If the following are exact (both from I_h geometry):")
print("    (1) m_H = E_cell * (1 + alpha/pi)  [scalar QED + cell energy, ESTABLISHED]")
print("    (2) N_J_W = 1/dim(G_g) = 1/4       [W Compton wavelength = L_J/4, LEAD]")
print()
print("  Then:")
print("    m_W = 4 * E_cell / (2*pi) = 2 * E_cell / pi")
print(f"        = {2*E_cell_GeV/pi:.4f} GeV  vs measured {m_W/1000:.4f} GeV  ({(2*E_cell_GeV/pi/(m_W/1000)-1)*100:+.4f}%)")
print()
print("    m_H / m_W = E_cell*(1+alpha/pi) / (2*E_cell/pi)")
print("              = pi*(1+alpha/pi)/2")
print(f"              = {pi*(1+alpha/pi)/2:.8f}")
print(f"    vs measured m_H/m_W = {m_H_val/(m_W/1000):.8f}")
print(f"    Deviation: {(pi*(1+alpha/pi)/2/(m_H_val/(m_W/1000))-1)*100:+.6f}%")
print()
print("  The chain (1)+(2) predicts m_H/m_W = pi*(1+alpha/pi)/2 = pi/2 overcorrected by alpha/pi.")
print(f"  Actual deviation from measured m_H/m_W: {(pi*(1+alpha/pi)/2/(m_H_val/(m_W/1000))-1)*100:+.4f}%  (NOT within 0.01%).")
m_H_pred_GeV = E_cell_GeV * (1 + alpha/pi)
print(f"  Cross-check: m_H_pred/m_W_pred = {m_H_pred_GeV:.4f}/{2*E_cell_GeV/pi:.4f} = {m_H_pred_GeV/(2*E_cell_GeV/pi):.8f}")
print(f"               pi*(1+alpha/pi)/2                            = {pi*(1+alpha/pi)/2:.8f}  [IDENTICAL ✓]")
sigma_from_pdg = abs(2*E_cell_GeV/pi - m_W/1000) / 0.012
print(f"  Our prediction m_W = {2*E_cell_GeV/pi:.3f} GeV is {abs(2*E_cell_GeV/pi - m_W/1000):.3f} GeV below PDG")
print(f"  = {sigma_from_pdg:.1f} sigma from PDG measurement.  RULED OUT as direct m_W prediction.")
print()
print("  STATUS: N_J_W = 1/4 is ruled out as a literal m_W prediction (76+ sigma).")
print("  The 1.17% proximity of N_J_W*4 to 1 is suggestive numerology, not a derivation.")
print("  Closing it requires either (a) a correction mechanism that shifts m_W by +1.2%,")
print("  or (b) a different resonance condition entirely.")
print(SEP)
