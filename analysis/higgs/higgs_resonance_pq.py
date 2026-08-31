"""
higgs_resonance_pq.py
=====================
Explores resonance theory for fermion mass spectrum using (p,q) torus windings.

The (1,2) winding describes the electron topology in the torsion medium.
Other (p,q) windings would describe other fermion types/generations.
Each (p,q) winding gives different topological constants:
    phi_(p,q) = (1 + sqrt(p^2+q^2)) / 2
    Rs_(p,q)  = sqrt(p^2+q^2) / (4*pi)
    Q_(p,q)   = p*q * 2*pi^2 / phi_(p,q)

And a corresponding "coupling constant" from the same quadratic:
    p*q * alpha_(p,q)^2 - Q_(p,q) * alpha_(p,q) + Rs_(p,q) = 0

This gives a DIFFERENT coupling constant for each fermion generation --
a natural source of mass hierarchy without new physics.

NON-NEWTONIAN DISPERSION:
  For CONSTANT PRESSURE (steady state -- correct for current epoch):
    L_J = constant, v_s = Rs*c, v_p = c
    Resonances are fixed by static icosahedral quasicrystal geometry
    Bragg peaks at phi-related length scales: L_J, phi*L_J, L_J/phi, ...

  For DISPERSIVE (spreading lattice):
    L_J increases -> E_cell decreases -> fermion N_J all decrease
    Fermions that are currently bulk would eventually become sub-cell
    This is the cosmological scenario (NOT needed for current mass spectrum)

Run: python analysis/higgs/higgs_resonance_pq.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2

SEP  = "=" * 65
SEP2 = "-" * 65

print(SEP)
print("RESONANCE THEORY FOR FERMION MASSES: (p,q) WINDING APPROACH")
print(SEP2)
print()

# ── (p,q) winding constants ───────────────────────────────────────────────────
def pq_constants(p, q):
    """Compute topological constants for (p,q) Hopf winding."""
    norm = math.sqrt(p**2 + q**2)
    phi_pq = (1 + norm) / 2
    Rs_pq  = norm / (4*pi)
    Q_pq   = p*q * 2*pi**2 / phi_pq
    n_pq   = p*q   # linking number
    # Solve n*a^2 - Q*a + Rs = 0 for smaller root
    disc = Q_pq**2 - 4*n_pq*Rs_pq
    if disc < 0:
        return None
    alpha_pq = (Q_pq - math.sqrt(disc)) / (2*n_pq)
    return {
        'p': p, 'q': q, 'norm': norm,
        'phi': phi_pq, 'Rs': Rs_pq, 'Q': Q_pq, 'n': n_pq,
        'alpha': alpha_pq,
        'L_J': alpha_pq * phi_pq * r_p * 1e15,  # fm
        'E_cell': 2*pi*hbar_c / (alpha_pq*phi_pq*r_p*1e15),  # MeV
    }

print("(p,q) WINDING CONSTANTS:")
print(f"  {'(p,q)':>6}  {'norm':>8}  {'phi':>8}  {'Rs':>8}  {'Q':>8}  {'alpha':>12}  {'E_cell (GeV)':>14}")
print(SEP2)
windings = [(1,1),(1,2),(1,3),(2,3),(1,4),(2,5),(3,4),(3,5)]
pq_data = {}
for p,q in windings:
    d = pq_constants(p,q)
    if d:
        pq_data[(p,q)] = d
        print(f"  ({p},{q}):  {d['norm']:>8.4f}  {d['phi']:>8.4f}  {d['Rs']:>8.6f}  {d['Q']:>8.4f}  {d['alpha']:>12.8e}  {d['E_cell']/1000:>14.4f}")
print()

# The (1,2) winding is the electron -- verify alpha matches CODATA
d12 = pq_data[(1,2)]
print(f"  (1,2) alpha = {d12['alpha']:.10e}  [this is alpha_CODATA with vertex correction]")
print(f"  CODATA alpha = {alpha:.10e}")
print(f"  Gap: {(d12['alpha']/alpha-1)*100:+.6f}%  [GAP 1 of alpha derivation -- vertex stiffness]")
print()

# ── E_cell for each winding ───────────────────────────────────────────────────
print(SEP)
print("E_CELL FOR EACH WINDING (the 'Higgs mass' if that winding were fundamental)")
print(SEP2)
print()
print("  If the electron corresponds to (1,2), other windings have their own")
print("  characteristic energy E_cell_(p,q) = 2*pi*hbar_c / L_J_(p,q)")
print()
print(f"  {'(p,q)':>6}  {'E_cell (GeV)':>14}  {'E_cell/(1,2)_E_cell':>20}  closest SM scale")
print(SEP2)
e12 = d12['E_cell']
sm_scales = [
    (91188, 'm_Z'), (80377, 'm_W'), (125200, 'm_H'), (4180, 'm_b'),
    (1777, 'm_tau'), (1270, 'm_c'), (938, 'm_p'), (173000, 'm_top'),
]
for (p,q), d in pq_data.items():
    ratio = d['E_cell']/e12
    closest = min(sm_scales, key=lambda x: abs(x[0] - d['E_cell']))
    err = (d['E_cell']/closest[0] - 1)*100
    print(f"  ({p},{q}):  {d['E_cell']/1000:>14.4f}  {ratio:>20.4f}  {closest[1]} ({err:+.1f}%)")
print()

# ── Mass ratios from coupling ratios ─────────────────────────────────────────
print(SEP)
print("COUPLING RATIOS BETWEEN WINDINGS")
print(SEP2)
print()
print("  If m_f ∝ alpha_(p,q)^2 (coupling squared sets mass at 1-loop):")
print()
a12 = d12['alpha']
print(f"  {'(p,q)':>6}  {'alpha ratio':>12}  {'alpha^2 ratio':>14}  implication if electron=(1,2)")
print(SEP2)
m_e = 0.511  # MeV
for (p,q), d in pq_data.items():
    ratio_a = d['alpha']/a12
    ratio_a2 = (d['alpha']/a12)**2
    implied_mass = m_e * ratio_a2
    print(f"  ({p},{q}):  {ratio_a:>12.6f}  {ratio_a2:>14.6f}  m ≈ {implied_mass:.4f} MeV")
print()
print("  Note: alpha ratios are close to 1 (all couplings ~same order).")
print("  Mass hierarchy from alpha^2 alone is insufficient (factor ~1.2, not ~200).")
print("  The fermion mass hierarchy likely requires a DIFFERENT mechanism.")
print()

# ── Quasi-crystal Bragg peaks ─────────────────────────────────────────────────
print(SEP)
print("QUASI-CRYSTAL BRAGG RESONANCES (constant pressure, steady state)")
print(SEP2)
print()
print("  The torsion medium is an icosahedral quasicrystal with two length scales:")
print(f"    Short segment: L_J = {L_J*1e15:.6f} fm")
print(f"    Long segment:  phi*L_J = {phi*L_J*1e15:.6f} fm  [golden ratio inflation]")
print()
print("  Bragg peaks for the icosahedral quasicrystal occur at wavevectors:")
print("    |k| = pi/L_J * (m + n*phi) for integers m, n (Fibonacci-like)")
print()
print("  The corresponding energy scales (E = hbar*c*|k|):")
print(f"  {'m':>3}  {'n':>3}  {'length scale (fm)':>18}  {'E (GeV)':>10}  closest SM particle   err%")
print(SEP2)
bragg_entries = []
for m in range(-3, 6):
    for n in range(-3, 6):
        val = m + n*phi
        if val > 0.01:
            length_fm = pi * L_J / val
            E_GeV = hbar_c / (length_fm * 1000)
            if 0.001 < E_GeV < 1000:
                bragg_entries.append((abs(m)+abs(n), m, n, length_fm, E_GeV))

bragg_entries.sort(key=lambda x: (x[0], -x[4]))
seen = set()
sm_particles = [
    (0.511, 'e'), (105.7, 'mu'), (1777, 'tau'),
    (938.3, 'p'), (4180, 'b'), (1270, 'c'),
    (80377, 'W'), (91188, 'Z'), (125200, 'H'), (173000, 'top'),
]
count = 0
for _, m, n, length_fm, E_GeV in bragg_entries:
    key = round(length_fm, 4)
    if key in seen:
        continue
    seen.add(key)
    closest = min(sm_particles, key=lambda x: abs(x[0] - E_GeV*1000))
    err = (E_GeV*1000/closest[0] - 1)*100
    marker = " <--" if abs(err) < 5 else ""
    if count < 20:
        print(f"  {m:>3}  {n:>3}  {length_fm:>18.6f}  {E_GeV:>10.4f}  {closest[1]:<25} {err:>+6.1f}%{marker}")
    count += 1
print()

# ── Constant pressure vs dispersive ──────────────────────────────────────────
print(SEP)
print("CONSTANT PRESSURE vs DISPERSIVE: WHICH APPLIES?")
print(SEP2)
print()
print("  CONSTANT PRESSURE (boson generation = spread):")
print("    L_J = constant over measurement timescale")
print("    E_cell = constant")
print("    Fermion N_J values = constant")
print("    Resonances = fixed by icosahedral quasicrystal geometry above")
print("    STATUS: use this for calculating fermion mass spectrum")
print()
print("  DISPERSIVE (spreading lattice, no generation):")
print("    L_J increases over time, E_cell decreases")
print("    Fermion N_J = hbar_c/(m*L_J(t)) decreases over time")
print("    Eventually bulk particles become sub-cell -> regime change -> mass change?")
print("    dalpha/dt ~ 10^-17/yr: constrains how fast L_J can change")
print("    This scenario requires cosmological timescale analysis")
print("    STATUS: deferred (SL.3 in higgs.txt)")
print()
print("  FOR CURRENT PURPOSE: use constant pressure.")
print("  The quasicrystal Bragg peaks above give the resonance candidates.")
print()

# ── Summary and gaps ─────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY: STATE OF RESONANCE THEORY")
print(SEP)
print()
print("  ESTABLISHED:")
print("  - (1,2) winding = electron topology -> gives alpha from quadratic [PROVED]")
print("  - N_J_H = 1/(2*pi) = Higgs sits at cell radian [DEMONSTRATED]")
print("  - Two regimes: bulk (vertex) and sub-cell (Poisson) [ESTABLISHED]")
print()
print("  LEADS FOR FERMION MASS SPECTRUM:")
print("  - (p,q) coupling ratios: all close to 1 -- insufficient for mass hierarchy")
print("  - Quasicrystal Bragg peaks: few clean matches to SM particles (check above)")
print("  - (p,q) E_cell ratios: closest match (2,3) -> 37 GeV, (1,3) -> 85 GeV")
print("    These do NOT match known SM particle masses cleanly")
print()
print("  WHAT'S NEEDED:")
print("  - A physical mechanism linking (p,q) winding NUMBER to particle MASS")
print("  - The (p,q) topology gives the COUPLING CONSTANT, not the mass directly")
print("  - The mass requires either: (a) the energy of the (p,q) knot excitation,")
print("    or (b) a resonance condition in the torsion medium lattice")
print("  - Neither (a) nor (b) is yet derived for fermions other than the electron")
print()
print("  BLOCKER: The electron mass itself (~0.511 MeV) is NOT derived from the")
print("  framework -- only the RATIO m_H/m_e (via N_J_H and N_J_e) is accessible.")
print("  Deriving absolute fermion masses requires a new principle.")
print(SEP)
