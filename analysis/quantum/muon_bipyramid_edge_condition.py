#!/usr/bin/env python3
"""
muon_bipyramid_edge_condition.py

Follow-up to muon_internal_force_check.py (F1-F3) and the "free muon
reverts to bipyramid" conclusion (doc_leptons.txt Sections 2.3, 4.3;
notes/koide_investigation_history.txt Part 1): the bipyramid used for
eff_mu is stated
to be the Johnson solid J13, defined by ALL EDGES EQUAL, which forces
h_t/r_e = 1/phi. This script checks that claim directly instead of just
asserting it.

  BR1: generalizes the bipyramid's two deflection angles (apex,
       equatorial) to a free height ratio x = h_t/r_e (r_e = 1), and
       verifies the formulas reduce to the documented cos_apex=1/(2phi),
       cos_eq=-1/sqrt5, eff_mu=(9-sqrt5)/8 exactly at x=1/phi.
  BR2: scans eff_mu(x) and the resulting predicted muon mass error across
       x, and checks that the x minimizing |error| does NOT coincide with
       1/phi -- i.e. 1/phi is not secretly a hidden best-fit to the
       measured mass (which would be circular, no different in kind from
       candidate (2) in doc_leptons.txt Section 4.3, which IS an explicit
       inversion against the measured mass). This is a deliberate
       negative-expectation check, run to rule out a circular
       justification, not to find one.
  BR3: verifies that x=1/phi is EXACTLY the equal-edge-length condition
       (apex-to-equator edge length = equator-to-equator edge length) for
       this specific 2-apex/1-ring path topology -- the actual,
       non-mass-fit geometric reason the "regular" bipyramid (J13) is
       used, previously only asserted in doc_leptons.txt prose ("all
       equal edges") without a step-by-step script derivation.

INTERPRETATION (context for why BR3 matters, not itself a numeric check):
  The STRUCTURAL G32 mode (doc_jobson_cell.txt Section 7.5) is fully
  pinned to the resting cell's lattice: both vertex POSITION (exact
  icosahedron coordinates) and edge LENGTH (all edges = L_J) are fixed by
  the surrounding 20-face closure requirement. The FREE muon (1 circuit,
  2 corpuscles, same section) is no longer forced to sit at those exact
  lattice positions -- there is no face-closure left to satisfy -- but
  there is no independent reason for its individual hop length to change
  either, since that is a property of the corpuscle/medium interaction,
  not of the lattice arrangement. "Equal-edge-length, free vertex
  position" is exactly the bipyramid's defining condition for this path
  topology, and BR3 shows it forces x=1/phi with no additional input.

  This resolves WHY x=1/phi specifically, GIVEN a 2-apex/1-ring path. It
  does NOT resolve why that simpler (4-vertex) topology is preferred over
  the icosahedron's own 2-apex/2-ring (6-vertex, antiprism-twisted)
  topology used by the structural G32 mode -- that remains open (see
  notes/koide_investigation_history.txt Part 1).

Reference: analysis/quantum/muon_internal_force_check.py (F1-F3),
           analysis/quantum/lepton_mass.py (LM4b, mass formula).
"""
import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP = "=" * 66
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
cos36 = phi / 2   # = cos(36 deg) exactly, r_e=1

print(SEP)
print("BR1 -- GENERALIZED BIPYRAMID DEFLECTION ANGLES vs HEIGHT RATIO x")
print(SEP)

def cos_apex(x):
    """Apex deflection cosine, pentagonal bipyramid, r_e=1, height ratio x."""
    return (cos36 - x*x) / (1 + x*x)

def cos_eq(x):
    """Equatorial deflection cosine, same bipyramid."""
    return (x*x - 1) / (1 + x*x)

def eff_mu(x):
    ratio = abs(cos_apex(x)) / abs(cos_eq(x))
    return (1 + ratio) / 2

x0 = 1 / phi
print(f"\n  At x0 = 1/phi = {x0:.8f}:")
print(f"    cos_apex(x0) = {cos_apex(x0):.8f}   (doc: 1/(2phi)  = {1/(2*phi):.8f})")
print(f"    cos_eq(x0)   = {cos_eq(x0):.8f}   (doc: -1/sqrt5  = {-1/math.sqrt(5):.8f})")
print(f"    eff_mu(x0)   = {eff_mu(x0):.8f}   (doc: (9-sqrt5)/8 = {(9-math.sqrt(5))/8:.8f})")

check("BR1a: cos_apex(1/phi) matches documented 1/(2phi)",
      abs(cos_apex(x0) - 1/(2*phi)) < 1e-10)
check("BR1b: cos_eq(1/phi) matches documented -1/sqrt5",
      abs(cos_eq(x0) - (-1/math.sqrt(5))) < 1e-10)
check("BR1c: eff_mu(1/phi) matches documented (9-sqrt5)/8",
      abs(eff_mu(x0) - (9-math.sqrt(5))/8) < 1e-10)


print()
print(SEP)
print("BR2 -- DOES MINIMIZING MASS ERROR SELECT x=1/phi? (expect: NO)")
print(SEP)

def mass_from_eff(eff, poly_norm):
    L3 = (eff**3 + log5**3) / (eff**2 + log5**2)
    xa = alpha * eff**2
    k = alpha * eff * (1 - 0.75*alpha**2) / (1 + xa + xa**2)
    dn = L3 * k
    base = 2*pi*alpha*(2/math.sqrt(5))*phi**2*m_p
    corr = 1 + Rs2 + 2*alpha
    return base * (1 + dn/poly_norm) * corr

poly_mu = 5*math.tan(pi/5)   # same polygon normalization as lepton_mass.py

def err_fn(x):
    return (mass_from_eff(eff_mu(x), poly_mu) - m_mu_pdg) / m_mu_pdg * 100

print(f"\n  {'x':>10}  {'eff_mu':>10}  {'m_mu (MeV)':>12}  {'err %':>10}")
xs_report = [0.30, 0.40, 0.45, 0.50, 0.55, 0.58, 0.60, x0, 0.64, 0.66, 0.70, 0.75, 0.80, 0.90]
for x in xs_report:
    eff = eff_mu(x)
    m = mass_from_eff(eff, poly_mu)
    err = err_fn(x)
    marker = "  <-- 1/phi (documented, used in mass formula)" if abs(x - x0) < 1e-6 else ""
    print(f"  {x:>10.6f}  {eff:>10.6f}  {m:>12.4f}  {err:>+10.4f}{marker}")

def absf(x):
    return abs(err_fn(x))

# golden-section search for the x minimizing |err(x)|; bounded away from the
# x=1 singularity (cos_eq(1)=0 -> eff_mu undefined there).
lo, hi = 0.30, 0.95
gr = (math.sqrt(5) - 1) / 2
c = hi - gr*(hi - lo)
d = lo + gr*(hi - lo)
for _ in range(200):
    if absf(c) < absf(d):
        hi = d
    else:
        lo = c
    c = hi - gr*(hi - lo)
    d = lo + gr*(hi - lo)
    if hi - lo < 1e-12:
        break
x_bestfit = (lo + hi) / 2
print(f"\n  Best-fit x (minimizes |err(x)|, search range [0.30,0.95]): x_bestfit = {x_bestfit:.8f}")
print(f"    err(x_bestfit) = {err_fn(x_bestfit):+.6f}%")
print(f"    1/phi          = {x0:.8f}   (x_bestfit - 1/phi = {x_bestfit - x0:+.8f})")

check("BR2: best-fit x does NOT coincide with 1/phi (mass-fit is not the justification)",
      abs(x_bestfit - x0) > 1e-4,
      f"|x_bestfit - 1/phi| = {abs(x_bestfit - x0):.6f} (nonzero, as expected)")


print()
print(SEP)
print("BR3 -- IS x=1/phi THE EQUAL-EDGE-LENGTH CONDITION? (independent of mass fit)")
print(SEP)

def edge_apex_eq(x):
    """Apex-to-equatorial-vertex edge length, r_e=1."""
    return math.sqrt(1 + x*x)

def edge_eq_eq():
    """Adjacent equatorial-vertex-to-vertex edge length (regular pentagon, r_e=1).
    Independent of x -- the equatorial ring's own radius is fixed at r_e=1."""
    return 2 * math.sin(pi/5)

eq_eq = edge_eq_eq()
print(f"\n  edge_eq_eq (pentagon side, r_e=1)  = {eq_eq:.8f}   (constant, does not depend on x)")
print(f"  edge_apex_eq(x0=1/phi)             = {edge_apex_eq(x0):.8f}")

check("BR3a: edge_apex_eq(1/phi) == edge_eq_eq (all-edges-equal, the J13 defining condition)",
      abs(edge_apex_eq(x0) - eq_eq) < 1e-10,
      f"apex-eq edge = {edge_apex_eq(x0):.8f}, eq-eq edge = {eq_eq:.8f}")

# Solve the equal-edge condition directly (no reference to phi anywhere) and
# check it independently reproduces 1/phi:
x_from_equal_edges = math.sqrt(eq_eq**2 - 1)
print(f"\n  Solving edge_apex_eq(x) = edge_eq_eq for x directly (no phi in the setup):")
print(f"    x_from_equal_edges = sqrt(4*sin^2(36deg) - 1) = {x_from_equal_edges:.8f}")
print(f"    1/phi                                          = {x0:.8f}")

check("BR3b: solving the equal-edge condition directly (no phi input) reproduces 1/phi",
      abs(x_from_equal_edges - x0) < 1e-10,
      f"diff = {x_from_equal_edges - x0:+.2e}")


print()
print(SEP)
print("INTERPRETATION")
print(SEP)
print("""
  BR1: the generalized formulas reduce exactly to the documented bipyramid
       deflection values at x=1/phi -- the x-parameterization is correct.
  BR2: x=1/phi is NOT secretly the best fit to the measured muon mass --
       the best fit is a nearby but DIFFERENT x. Rules out a circular
       "chose 1/phi because it matches the mass" justification.
  BR3: x=1/phi is exactly and ONLY the equal-edge-length condition for
       this 2-apex/1-ring path topology -- an independent, non-mass-fit
       geometric reason, matching the medium's own uniform-edge-length
       convention (every icosahedral edge = L_J elsewhere in the
       framework) carried over as the one surviving constraint once the
       free muon is no longer pinned to exact lattice vertex positions.

  STILL OPEN (see notes/koide_investigation_history.txt Part 1): WHY the
  free muon's path relaxes to this SIMPLER 2-apex/1-ring (4-vertex)
  topology at all, rather than the icosahedron's own 2-apex/2-ring
  (6-vertex, antiprism-twisted) topology used by the structural G32 mode.
  This script assumes the bipyramid topology and only tests the height
  ratio within it.
""")

print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
for name, status, detail in results:
    print(f"  {status}: {name}")
print(SEP)
