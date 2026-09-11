"""
mesh_grind_2d_lattice_check.py
===============================
Direct follow-up to Stage 2 (analysis/nuclear/mesh_grind_intervening_cell_check.py),
prompted by the author's question: Stage 2 was a 1D CHAIN of 2D disks (one
column, each grain touching only its neighbor directly above/below). That is
exactly the topology where gear-train parity effects (odd vs even chain
length) are expected. Does the parity effect survive once there is more than
one column -- i.e. multiple parallel contact paths between the same two
boundaries, sharing a single coherent boundary velocity each side (matching
the physical picture: the whole Zone 2/3 shell surface moves at one coherent
tangential speed, not alternating point to point)?

WHAT STAGE 2 ACTUALLY WAS (clarification): a single column, N grains stacked
vertically between boundary 1 (bottom) and boundary 2 (top). M (number of
side-by-side columns) was implicitly 1 throughout. This script generalizes
to an N (rows, stacking direction) x M (columns, side-by-side direction)
grid of grains, all sharing the SAME two boundaries.

METHOD: same rigid-body kinematics as Stage 2 (v = v_c + omega x r, grain
centers fixed -- translationally jammed, only spin free), now with TWO
contact types instead of one:
  VERTICAL   (within a column): grain(i,j)-grain(i+1,j) -- same rule as Stage 2.
  HORIZONTAL (within a row):    grain(i,j)-grain(i,j+1) -- a new constraint,
    derived the same way: touching disks side-by-side (in the x-direction)
    require omega(i,j) = -omega(i,j+1) for zero relative slip at their
    shared contact point (the SAME kind of alternation as the vertical case,
    now sideways).
Built as a real (constraints)x(N*M) linear least-squares system (numpy),
not hand-derived -- Section 1 below shows why hand-derivation is untrustworthy
here (a naive "checkerboard" ansatz looks like it should work but does NOT
satisfy the boundary condition once M>1; only the actual linear solve settles
what really happens).

Run: python analysis/nuclear/mesh_grind_2d_lattice_check.py
Reference: analysis/nuclear/mesh_grind_intervening_cell_check.py (Stage 2),
  /memories/repo/mesh-grind-chirality-derivation.md
"""

import sys, os, math
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

pi = math.pi
Rs = math.sqrt(5) / (4 * pi)
c = 2.99792458e8
v_surface = Rs * c
r = 1.0
TOL = 1e-6 * v_surface

SEP = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

def idx(i, j, M):
    return i * M + j

def build_2d_system(N, M, r):
    """
    N rows (stacking direction, boundary1->boundary2), M columns (side by
    side). Unknowns omega[i,j], flattened row-major. Returns dense matrix A
    (rows = contacts, cols = N*M unknowns) and a template for the RHS
    (list tagged 'v1', 'v2', or 0.0 per contact row).
    """
    n_unk = N * M
    rows = []
    rhs_tags = []

    # Vertical contacts (per column j): boundary1 -> (0,j) -> (1,j) -> ... -> (N-1,j) -> boundary2
    for j in range(M):
        row = np.zeros(n_unk); row[idx(0, j, M)] = r
        rows.append(row); rhs_tags.append('v1')
        for i in range(N - 1):
            row = np.zeros(n_unk)
            row[idx(i, j, M)] = -r
            row[idx(i + 1, j, M)] = -r
            rows.append(row); rhs_tags.append(0.0)
        row = np.zeros(n_unk); row[idx(N - 1, j, M)] = -r
        rows.append(row); rhs_tags.append('v2')

    # Horizontal contacts (per row i): (i,j) -- (i,j+1)
    for i in range(N):
        for j in range(M - 1):
            row = np.zeros(n_unk)
            row[idx(i, j, M)] = r
            row[idx(i, j + 1, M)] = r
            rows.append(row); rhs_tags.append(0.0)

    A = np.array(rows)
    return A, rhs_tags

def solve_2d(N, M, v1, v2):
    A, tags = build_2d_system(N, M, r)
    b = np.array([v1 if t == 'v1' else (v2 if t == 'v2' else 0.0) for t in tags])
    omega, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    residual = A @ omega - b
    return omega, residual

# ── Section 1: why the "checkerboard" hand-derivation cannot be trusted ──────
print(SEP)
print("SECTION 1: WHY THIS MUST BE COMPUTED, NOT HAND-DERIVED")
print(SEP2)
print("""
  Hand-reasoning suggests a checkerboard pattern omega[i,j] = (-1)^(i+j)*omega0
  satisfies every vertical AND horizontal contact identically. But row 0
  touches a SINGLE boundary 1 at a SINGLE velocity v1 -- while the checkerboard
  forces omega[0,j] to alternate sign with j. A sign-alternating quantity
  cannot equal one constant v1/r for every column j at once (M>1) except the
  trivial omega0=0 case. So the checkerboard guess already looks broken for
  M>1 -- but that is still hand-reasoning. The actual answer requires solving
  the real over-determined system and reading off the residual, which is what
  the rest of this script does.
""")

# ── Section 2: does a single row (N=1) still mesh exactly as M grows? ────────
print(SEP)
print("SECTION 2: N=1 (ONE LAYER, THE PHYSICALLY MOTIVATED CASE) -- VARYING M")
print(SEP2)

for M in range(1, 6):
    _, resid_opp = solve_2d(1, M, v_surface, -v_surface)
    _, resid_same = solve_2d(1, M, v_surface, v_surface)
    print(f"  M={M}:  max|resid|(opposite chirality)={np.max(np.abs(resid_opp)):.4e}   "
          f"max|resid|(same chirality)={np.max(np.abs(resid_same)):.4e}")

_, resid_opp_M1 = solve_2d(1, 1, v_surface, -v_surface)
_, resid_opp_M3 = solve_2d(1, 3, v_surface, -v_surface)
check("2DL1: N=1,M=1 (Stage 2's own case) still meshes exactly under opposite chirality",
      np.max(np.abs(resid_opp_M1)) < TOL,
      f"max residual = {np.max(np.abs(resid_opp_M1)):.3e}")

check("2DL2: N=1,M=3 (three parallel columns, single layer) -- opposite chirality mesh survives adding columns",
      np.max(np.abs(resid_opp_M3)) < TOL,
      f"max residual = {np.max(np.abs(resid_opp_M3)):.3e}")

# ── Section 3: multi-row (N>1) with multiple columns (M>1) -- does the parity survive? ─
print()
print(SEP)
print("SECTION 3: MULTI-ROW x MULTI-COLUMN GRID -- DOES THE ODD/EVEN-N PARITY SURVIVE?")
print(SEP2)

pattern2d = {}
for N in range(1, 6):
    row_report = []
    for M in range(1, 5):
        _, r_same = solve_2d(N, M, v_surface, v_surface)
        _, r_opp = solve_2d(N, M, v_surface, -v_surface)
        same_ok = np.max(np.abs(r_same)) < TOL
        opp_ok = np.max(np.abs(r_opp)) < TOL
        verdict = "SAME" if same_ok else ("OPP" if opp_ok else "NEITHER")
        row_report.append(verdict)
        pattern2d[(N, M)] = verdict
    print(f"  N={N} (rows): M=1..4 verdicts -> {row_report}")

print()
print("  Reading the table: does the verdict in each row stay constant across M")
print("  (parity survives, set only by N, as in the 1D chain), or does it change")
print("  once M>1 (parity is an artifact of the 1D topology)?")

odd_N_stable = all(pattern2d[(N, 1)] == pattern2d[(N, M)] for N in (1, 3, 5) for M in (1, 2, 3, 4))
even_N_stable = all(pattern2d[(N, 1)] == pattern2d[(N, M)] for N in (2, 4) for M in (1, 2, 3, 4))

check("2DL3: for ODD N, the M=1 verdict (OPP meshes) is unchanged for all tested M (parity survives sideways extension)",
      odd_N_stable,
      f"odd-N rows: { {N: [pattern2d[(N,M)] for M in (1,2,3,4)] for N in (1,3,5)} }")

check("2DL4: for EVEN N, the M=1 verdict (SAME meshes) is unchanged for all tested M (parity survives sideways extension)",
      even_N_stable,
      f"even-N rows: { {N: [pattern2d[(N,M)] for M in (1,2,3,4)] for N in (2,4)} }")

# ── Summary ────────────────────────────────────────────────────────────────────
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
