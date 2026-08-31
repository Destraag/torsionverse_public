"""
higgs_resonance_series.py
=========================
Checks whether SM particle masses correspond to specific N_J values
derivable from icosahedral geometry (lattice resonances).

If the torsion medium creates preferred coupling scales, SM particles
might sit at N_J values set by icosahedral numbers, not arbitrary.

Run: python analysis/higgs/higgs_resonance_series.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2
Rs   = math.sqrt(5) / (4*pi)
hbar_c_fm = hbar_c   # MeV*fm

SEP  = "=" * 65
SEP2 = "-" * 65

# ── Icosahedral structural numbers ────────────────────────────────────────────
# These are "natural" values that might set preferred N_J positions
ih_numbers = {
    '1/(2*pi)':         1/(2*pi),
    '1/pi':             1/pi,
    '1/(pi*phi)':       1/(pi*phi),
    '1/(2*pi*phi)':     1/(2*pi*phi),
    '1/phi':            1/phi,
    '1/phi^2':          1/phi**2,
    '1/phi^3':          1/phi**3,
    'Rs':               Rs,
    'Rs/phi':           Rs/phi,
    '2*Rs':             2*Rs,
    'Rs^2':             Rs**2,
    'alpha/pi':         alpha/pi,
    '1':                1.0,
    'phi':              phi,
    'pi':               pi,
    'phi^2':            phi**2,
    '2*pi':             2*pi,
    '4*pi':             4*pi,
    'N_lock/(2*pi^2)':  N_lock/(2*pi**2),
    '12 (vertices)':    12.0,
    '20 (faces)':       20.0,
    '30 (edges)':       30.0,
    'phi^(1/7)':        phi**(1/7),
}

print(SEP)
print("N_J RESONANCE SERIES -- ICOSAHEDRAL LATTICE RESONANCES")
print(SEP2)
print(f"  L_J = {L_J*1e15:.6f} fm")
print(f"  E_cell = {E_cell_GeV:.4f} GeV")
print()

# ── Known SM particle masses ──────────────────────────────────────────────────
particles = [
    ('top',    172760,   'm_top (PDG 2022)'),
    ('Higgs',  125200,   'm_H (PDG 2022)'),
    ('Z',       91188,   'm_Z (PDG 2022)'),
    ('W',       80377,   'm_W (PDG 2022)'),
    ('b',        4180,   'm_b (MSbar)'),
    ('tau',      1777,   'm_tau'),
    ('c',        1270,   'm_c (MSbar)'),
    ('proton',    938.3, 'm_p'),
    ('muon',      105.7, 'm_mu'),
    ('e',         0.511, 'm_e'),
]

# Compute N_J for each particle
print("SM PARTICLE N_J VALUES:")
print(f"  {'Particle':<10} {'mass (MeV)':>12} {'N_J':>12}   {'Nearest ih#':>16}   {'error%':>8}")
print(SEP2)

for name, mass_MeV, label in particles:
    nj = hbar_c_fm / (mass_MeV * L_J)   # L_J from constants.py is already in fm
    
    # Find nearest icosahedral number
    nearest_name = min(ih_numbers.keys(), key=lambda k: abs(nj - ih_numbers[k]))
    nearest_val  = ih_numbers[nearest_name]
    err_pct = (nj/nearest_val - 1)*100
    print(f"  {name:<10} {mass_MeV:>12.1f} {nj:>12.5f}   {nearest_name:>16}={nearest_val:.5f}   {err_pct:>+8.3f}%")

print()

# ── Series check: N_J = k/(2*pi) for integer k ───────────────────────────────
print(SEP)
print("SERIES A: N_J = k/(2*pi), mass = E_cell * k")
print(SEP2)
print(f"  {'k':>3}  {'N_J':>8}  {'mass (GeV)':>12}  closest particle         err%")
print(SEP2)
for k in range(1, 8):
    nj_target = k / (2*pi)
    m_GeV = k * E_cell_GeV
    # find closest known particle
    closest = min(particles, key=lambda p: abs(p[1]/1000 - m_GeV))
    err = (closest[1]/1000/m_GeV - 1)*100
    marker = " <--" if abs(err) < 5 else ""
    print(f"  {k:>3}  {nj_target:>8.5f}  {m_GeV:>12.4f}  {closest[0]:<25} {err:>+6.2f}%{marker}")
print()

# ── Series check: N_J = k/phi^n for various n ────────────────────────────────
print(SEP)
print("SERIES B: N_J = 1/(k*phi^n) -- does phi power give known masses?")
print(SEP2)
print(f"  {'n':>3}  {'k':>3}  {'N_J':>8}  {'mass (GeV)':>12}  closest particle         err%")
print(SEP2)
for n in range(-2, 5):
    for k in [1, 2, 3, 4, 5]:
        nj_target = 1.0 / (k * phi**n) if n >= 0 else k * phi**(-n)
        if nj_target > 100 or nj_target < 0.001:
            continue
        m_GeV = hbar_c_fm / (nj_target * L_J * 1000)
        closest = min(particles, key=lambda p: abs(p[1]/1000 - m_GeV))
        err = (closest[1]/1000/m_GeV - 1)*100
        if abs(err) < 3.0:   # only show close matches
            marker = " <-- CLOSE" if abs(err) < 1.0 else " <"
            print(f"  {n:>3}  {k:>3}  {nj_target:>8.5f}  {m_GeV:>12.4f}  {closest[0]:<25} {err:>+6.2f}%{marker}")
print()

# ── Series check: N_J from I_h representation dimensions ─────────────────────
print(SEP)
print("SERIES C: N_J from I_h representation dimensions")
print(SEP2)
ih_dims    = [1, 3, 3, 4, 5]       # gerade dim
ih_dim_names = ['A_g', 'T_1g', 'T_2g', 'G_g', 'H_g']
print(f"  {'Rep':>5}  {'dim':>4}  {'N_J=1/dim':>10}  {'mass (GeV)':>12}  closest particle     err%")
print(SEP2)
for rep_name, d in zip(ih_dim_names, ih_dims):
    nj_target = 1.0 / d
    m_GeV = hbar_c_fm / (nj_target * L_J * 1000)
    closest = min(particles, key=lambda p: abs(p[1]/1000 - m_GeV))
    err = (closest[1]/1000/m_GeV - 1)*100
    marker = " <--" if abs(err) < 3 else ""
    print(f"  {rep_name:>5}  {d:>4}  {nj_target:>10.5f}  {m_GeV:>12.4f}  {closest[0]:<25} {err:>+6.2f}%{marker}")
print()

# ── Series D: N_J from vertex geometry distances ──────────────────────────────
print(SEP)
print("SERIES D: N_J from icosahedral distances (in L_J units)")
print(SEP2)
# Key distances in the icosahedron (as multiples of L_J)
distances = {
    'edge = L_J':            1.0,
    'phi*L_J (edge-share)':  phi,
    'phi^2*L_J':             phi**2,
    'sqrt(phi+2)*L_J':       math.sqrt(phi+2),
    'phi^2/sqrt(3)*L_J':     phi**2/math.sqrt(3),
    '2*phi*L_J':             2*phi,
    'pi*L_J':                pi,
    '2*pi*L_J':              2*pi,
}
print(f"  {'Distance':>30}  {'N_J = R/L_J':>12}  {'mass (GeV)':>12}  closest   err%")
print(SEP2)
for desc, d_ratio in distances.items():
    # N_J = lambda_Compton / L_J = lambda_Compton_as_multiple_of_L_J / 1
    # If N_J = 1/d_ratio then m = d_ratio * E_cell / (2*pi)
    nj_target = 1.0 / d_ratio
    m_GeV = hbar_c_fm / (nj_target * L_J * 1000)
    closest = min(particles, key=lambda p: abs(p[1]/1000 - m_GeV))
    err = (closest[1]/1000/m_GeV - 1)*100
    marker = " <--" if abs(err) < 5 else ""
    print(f"  {desc:>30}  {nj_target:>12.5f}  {m_GeV:>12.4f}  {closest[0]:<10} {err:>+6.2f}%{marker}")
print()

# ── The b-quark focused check ────────────────────────────────────────────────
print(SEP)
print("B-QUARK FOCUSED CHECK (N_J_b = 4.75)")
print(SEP2)
m_b = 4180   # MeV
L_J_fm_val = L_J  # already in fm
N_J_b = hbar_c_fm / (m_b * L_J_fm_val)
N_lock_val = 2*pi / (alpha*phi)
print(f"  N_J_b = {N_J_b:.6f}")
print()
print(f"  Checking icosahedral combinations:")
candidates = {
    'N_lock / 120': N_lock_val / 120,
    'N_lock / |I_h|': N_lock_val / 120,
    'N_lock / (4*pi^2)': N_lock_val / (4*pi**2),
    '5*Rs/alpha': 5*Rs/alpha,
    'phi^3': phi**3,
    'phi^3 - 1/phi': phi**3 - 1/phi,
    '5/phi^(1/3)': 5/phi**(1/3),
    '3*pi/2': 3*pi/2,
    '2*pi - Rs': 2*pi - Rs,
    '(5+phi^2)/phi': (5+phi**2)/phi,
}
for name, val in candidates.items():
    err = (N_J_b / val - 1) * 100
    marker = " <--" if abs(err) < 2 else ""
    print(f"    {name:<30} = {val:.6f}  ({err:>+6.3f}%){marker}")
print()

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()
print("  Series A (k/(2*pi)): N_J = k/(2*pi) gives m = k*E_cell.")
print("  Only k=1 matches a known particle (Higgs). Others don't match SM.")
print()
print("  Series B (phi powers): certain combinations within 2-3% of particles,")
print("  but none clean enough to claim geometric origin.")
print()
print("  Series C (I_h dimensions): 1/dim gives masses not matching SM particles.")
print()
print("  Series D (icosahedral distances): phi-related N_J values don't map")
print("  cleanly to known SM particle masses.")
print()
print("  B-quark: N_J_b = 4.75 is NOT cleanly explained by any icosahedral")
print("  combination within 2%. The b-quark mass is not yet geometrically derived.")
print()
print("  CONCLUSION: No clean resonance series found that predicts multiple SM")
print("  particle masses from icosahedral geometry. The Higgs (k=1 in Series A)")
print("  remains the only particle with a geometric N_J derivation.")
print("  Fermion masses require a deeper theory beyond simple resonance series.")
print(SEP)
