"""
mesh_grind_intervening_cell_check.py
=====================================
Stage 2 of the "chirality -> mesh/grind -> pressure sign" derivation plan
(see /memories/repo/mesh-grind-chirality-derivation.md and Stage 1,
analysis/nuclear/mesh_grind_contact_kinematics.py).

GAP THIS ADDRESSES:
  Stage 1 showed two rotating bodies in DIRECT contact either slide
  (same chirality, grind) or roll without slipping (opposite chirality,
  mesh). But every doc actually describes cells BETWEEN two chirality
  sources being dragged, not the sources touching directly. This script
  refines Stage 1 to the physically correct picture: one or more
  intervening Jobson cells sandwiched between two rotating boundary
  regions (the Zone 2/3 shells, locally flattened -- a small contact
  patch on a much larger curved rotating shell looks, to leading order,
  like a flat surface translating at fixed tangential speed v=Rs*c, the
  same v_surface already used in Stage 1 and elsewhere in this framework).

PHYSICAL JUSTIFICATION FOR THE MODEL (not a new assumption -- already
established elsewhere in this framework, doc_nucleus.txt / MG-J9):
  cells are TRANSLATIONALLY jammed (locked in place, cannot rearrange/flow
  past neighbors) but retain zero-cost ROTATIONAL freedom (spin in place).
  This licenses exactly the model built here: each intervening grain's
  CENTER is fixed (v_c = 0, confined by the surrounding dense packing);
  only its spin omega is a free variable. This is a real physical
  constraint, not a simplifying convenience -- it is what the framework
  already says about this medium.

METHOD:
  N grains (radius r, centers fixed) stacked between boundary 1 (velocity
  v1) and boundary 2 (velocity v2). Adjacent bodies contact via ordinary
  rigid-body kinematics (v = v_c + omega x r); each contact gives ONE
  linear equation in the grains' spins. N grains touching 2 boundaries
  give N+1 contacts but only N free spins -- generically OVER-determined
  by exactly one equation. Built and solved as a real (N+1)xN linear
  least-squares system (numpy), not hand-derived -- the residual at each
  contact IS the forced (unavoidable) slip velocity there.

WHAT THIS ESTABLISHES:
  For N=1 (a single intervening cell): opposite chirality (v2=-v1) admits
  an EXACT zero-residual solution (perfect mesh, no forced slip anywhere).
  Same chirality (v2=v1) admits NO exact solution -- forced residual slip
  of exactly v_surface at BOTH contacts (total slip budget 2*v_surface,
  conserved from Stage 1's single-contact 2*v_surface, now split across
  2 contacts instead of concentrated at 1).

WHAT THIS SCRIPT ALSO SURFACES (an honest complication, not swept under
the rug): the mesh/grind condition's dependence on chain length N is
tested directly (not assumed) in Section 2 -- see the printed verdict for
whether "opposite chirality meshes" continues to hold, or flips, as N
grows (a real gear-train fact: idler gears can reverse the required sense).

WHAT REMAINS OPEN (Stage 2 does NOT close this): converting a nonzero
forced-slip velocity into an actual PRESSURE VALUE requires a friction/
contact-stiffness law for the Jobson-cell medium that does not yet exist
in this framework. Stage 2 establishes WHERE dissipation/stress must occur
(mesh: nowhere: grind: at every contact) -- not yet HOW MUCH pressure that
produces.

Run: python analysis/nuclear/mesh_grind_intervening_cell_check.py
Reference: docs/series1/doc_electron.txt Section 4-5, doc_orbit_pressure.txt
  Section 1.2/1.2a, doc_nucleus.txt Section 5.4 (N2-1)
"""

import sys, os, math
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'higgs'))

pi = math.pi
Rs = math.sqrt(5) / (4 * pi)
c = 2.99792458e8    # m/s, exact by SI definition
v_surface = Rs * c  # Zone 3 co-rotation speed [Stage 1; entanglement_phase_lock.py]
r = 1.0             # grain radius -- Stage 1's MGK4 already confirmed R-independence
TOL = 1e-6 * v_surface   # "exact zero" tolerance, RELATIVE to v_surface (~5.3e7 m/s) --
                          # an absolute 1e-9 threshold is meaningless at this magnitude;
                          # real floating-point noise here is ~1e-8 m/s (~16 sig figs).

SEP = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

def build_system(N, r):
    """
    N grains stacked between two boundaries. Rows = contacts (N+1 total),
    columns = grain spins omega_1..omega_N (N unknowns). Grain i's
    bottom-facing contact velocity = +omega_i*r, top-facing = -omega_i*r
    (v = v_c + omega x r, v_c=0, r_bottom=(0,-r), r_top=(0,r)). Boundaries
    are flat/non-rotating: contribute only their own translation velocity.
    """
    A = np.zeros((N + 1, N))
    A[0, 0] = r                      # contact 0: boundary1-grain1 -> omega_1*r = v1
    for i in range(1, N):            # contacts 1..N-1: grain_i - grain_(i+1)
        A[i, i - 1] = -r             # grain i's top-facing velocity term
        A[i, i] = -r                 # grain (i+1)'s bottom-facing velocity term
    A[N, N - 1] = -r                 # contact N: grain_N-boundary2 -> -omega_N*r = v2
    return A

def solve_and_residual(N, v1, v2):
    """Returns (omega_solution, residual_vector) for N intervening grains."""
    if N == 0:
        # No intervening grain: direct boundary-boundary contact, 1 constraint, 0 unknowns.
        return np.array([]), np.array([v1 - v2])
    A = build_system(N, r)
    b = np.zeros(N + 1)
    b[0] = v1
    b[N] = v2
    omega, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    residual = A @ omega - b
    return omega, residual

# ── Section 1: single intervening grain (N=1), the direct Stage 1 refinement ─
print(SEP)
print("SECTION 1: ONE INTERVENING CELL BETWEEN TWO CHIRALITY SOURCES")
print(SEP2)

omega_mesh, resid_mesh = solve_and_residual(1, v_surface, -v_surface)   # opposite chirality
omega_grind, resid_grind = solve_and_residual(1, v_surface, v_surface)  # same chirality

print(f"  v_surface = Rs*c = {v_surface:.6e} m/s = {Rs:.4f}*c")
print()
print(f"  OPPOSITE chirality (v1={v_surface:.3e}, v2={-v_surface:.3e}):")
print(f"    grain spin omega = {omega_mesh[0]:.6e}")
print(f"    residual slip at each contact: {resid_mesh}")
print()
print(f"  SAME chirality (v1=v2={v_surface:.3e}):")
print(f"    grain spin omega = {omega_grind[0]:.6e}  (least-squares optimum)")
print(f"    residual slip at each contact: {resid_grind}")

check("MGC1: opposite chirality (N=1) admits an EXACT zero-slip solution (perfect mesh)",
      np.max(np.abs(resid_mesh)) < TOL,
      f"max residual = {np.max(np.abs(resid_mesh)):.3e} m/s  (tol={TOL:.3e})")

check("MGC2: same chirality (N=1) has NO exact solution -- forced residual = v_surface at BOTH contacts",
      np.allclose(np.abs(resid_grind), v_surface, rtol=1e-9),
      f"residuals = {resid_grind}  (expect +/-{v_surface:.4e} at each)")

total_grind_N1 = np.sum(np.abs(resid_grind))
check("MGC3: total forced-slip budget conserved: N=1's 2 contacts sum to 2*v_surface (Stage 1's single-contact total)",
      abs(total_grind_N1 - 2 * v_surface) < 1e-6 * v_surface,
      f"sum|residual| = {total_grind_N1:.6e}  2*v_surface = {2*v_surface:.6e}")

# ── Section 2: does the mesh/grind condition depend on chain length N? ───────
print()
print(SEP)
print("SECTION 2: DOES ADDING MORE INTERVENING CELLS CHANGE WHICH CHIRALITY MESHES?")
print(SEP2)
print(f"  Testing N=0..6 directly (not assumed) -- for each N, compute the forced")
print(f"  residual slip under 'same absolute chirality' (v1=v2) and 'opposite")
print(f"  absolute chirality' (v1=-v2), report which one is the zero-residual (mesh) case.")
print()

pattern = []
for N in range(0, 7):
    _, r_same = solve_and_residual(N, v_surface, v_surface)
    _, r_opp = solve_and_residual(N, v_surface, -v_surface)
    same_is_mesh = np.max(np.abs(r_same)) < TOL
    opp_is_mesh = np.max(np.abs(r_opp)) < TOL
    verdict = "SAME meshes" if same_is_mesh else ("OPPOSITE meshes" if opp_is_mesh else "NEITHER exact")
    pattern.append(verdict)
    print(f"  N={N}:  max|resid|(same)={np.max(np.abs(r_same)):.3e}   "
          f"max|resid|(opposite)={np.max(np.abs(r_opp)):.3e}   -> {verdict}")

print()
print(f"  Pattern across N=0..6: {pattern}")

odd_N_mesh_opposite = all(pattern[N] == "OPPOSITE meshes" for N in range(0, 7) if N % 2 == 1)
even_N_mesh_same = all(pattern[N] == "SAME meshes" for N in range(0, 7) if N % 2 == 0 and N > 0)

check("MGC4: for ODD N (1,3,5), OPPOSITE chirality is the exact zero-slip (mesh) case",
      odd_N_mesh_opposite,
      f"odd-N verdicts: {[pattern[N] for N in range(0,7) if N%2==1]}")

check("MGC5: for EVEN N>=2 (2,4,6), SAME chirality becomes the exact zero-slip (mesh) case -- a PARITY FLIP",
      even_N_mesh_same,
      f"even-N (>0) verdicts: {[pattern[N] for N in range(0,7) if N%2==0 and N>0]}")

# ── Section 3: what this means physically ─────────────────────────────────────
print()
print(SEP)
print("SECTION 3: PHYSICAL READING AND HONEST SCOPE")
print(SEP2)
print("""
  CONFIRMED, not assumed: with an ODD number of intervening cells (matching
  the single-cell N=1 case doc_electron.txt/doc_orbit_pressure.txt/
  doc_nucleus.txt actually describe), opposite chirality gives an EXACT,
  zero-residual rolling solution (mesh) and same chirality gives a forced,
  nonzero residual slip at every contact (grind) -- exactly the established
  claim, now derived rather than asserted.

  HONEST COMPLICATION (found here, not previously flagged anywhere in this
  framework): this is a real gear-train parity effect -- an EVEN number of
  intervening cells flips which chirality meshes. This means the simple
  'same chirality always grinds, opposite always meshes' rule implicitly
  assumes an ODD contact-chain length between the two sources. Whether the
  real nucleon-nucleon Zone 2/3 contact zone has an odd or even effective
  cell count is NOT established anywhere in this repo -- Stage 1's own
  r_grind/r_mid_cell = 52.3 (mesh_grind_contact_kinematics.py MGK7) is a
  LENGTH ratio, not a validated discrete contact count, and 52 (if taken
  literally) would be EVEN. This is a genuine new open question this
  derivation surfaces, not a result that was already known.

  STILL OPEN (unchanged from Stage 1): converting a nonzero forced-slip
  velocity into an actual pressure/energy value needs a friction/contact-
  stiffness law for this medium that does not yet exist here.
""")

check("MGC6: single-cell case (N=1, matching the docs' own 'cells between them' language) "
      "is unambiguous -- odd-N result, no parity ambiguity there",
      pattern[1] == "OPPOSITE meshes",
      f"N=1 verdict: {pattern[1]}")

# ── Summary ────────────────────────────────────────────────────────────────────
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
