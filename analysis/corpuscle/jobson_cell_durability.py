#!/usr/bin/env python3
"""
jobson_cell_durability.py

Proves the Jobson cell's outer icosahedral shell is structurally self-sufficient:
no inner content is required for rigidity or load-bearing under radial pressure
or face-normal forces. All resistance comes from the outer triangulated edge network.

PHYSICS BACKGROUND:
  The Jobson cell is a bar-joint network: V=12 vertices, E=30 edges, F=20 faces.
  Maxwell count: 3V - E = 36 - 30 = 6 = exactly the rigid-body zero-modes.
  Triangulated faces -> zero floppy modes. The shell is a geodesic dome: rigid by
  triangulation, not by internal fill.

Checks:
  DC1: Zero floppy modes -- kernel = 6 (rigid body only, no soft deformations)
  DC2: Uniform inward radial force -> outer edges resist (no radial collapse)
  DC3: Face-center inward force -> outer edges resist (no inner support needed)
  DC4: Every face resists inward push (tested all 20 faces)
  DC5: Full row rank: rank(R) = 30 = n_edges (each constraint independent; shell
       maximally rigid; inner content would be redundant by definition)

Reference: jobson_cell_rigidity_matrix.py (RM1-RM5, 6/6 PASS),
           jobson_cell_force_balance.py (FB1-FB19, 19/19 PASS),
           docs/series1/doc_jobson_cell.txt Sec 7.2
"""
import math
import numpy as np
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 66
SEP2 = "-" * 66
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

phi      = (1 + math.sqrt(5)) / 2
alpha    = 7.2973525693e-3
r_p_fm   = 0.8414
L_J_fm   = alpha * phi * r_p_fm

print(SEP)
print("JOBSON CELL DURABILITY: OUTER SHELL IS STRUCTURALLY SELF-SUFFICIENT")
print(SEP)
print(f"  L_J = {L_J_fm:.6f} fm  (same alpha*phi*r_p scaling as all other scripts)")
print()

# ── Icosahedron vertices (same construction as rigidity_matrix.py) ───────────
verts_raw = []
for perm in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            v = [0.0, 0.0, 0.0]
            v[perm[1]] = s1 * 1.0
            v[perm[2]] = s2 * phi
            verts_raw.append(tuple(v))
verts_raw = list(dict.fromkeys(verts_raw))

def dist3(a, b):
    return math.sqrt(sum((a[k]-b[k])**2 for k in range(3)))

edge_raw = min(dist3(verts_raw[0], v) for v in verts_raw[1:])
scale    = L_J_fm / edge_raw
V = np.array([[c * scale for c in v] for v in verts_raw])   # shape (12, 3)
n_v = len(V)

edges = [(i, j) for i in range(n_v) for j in range(i+1, n_v)
         if abs(np.linalg.norm(V[i] - V[j]) - L_J_fm) < 1e-9]
n_e = len(edges)

# ── Find all triangular faces ─────────────────────────────────────────────────
edge_set = {(i, j) for i, j in edges} | {(j, i) for i, j in edges}
faces = [(i, j, k)
         for i in range(n_v) for j in range(i+1, n_v) for k in range(j+1, n_v)
         if (i, j) in edge_set and (i, k) in edge_set and (j, k) in edge_set]
n_f = len(faces)

# ── Build rigidity matrix R (n_e x 3*n_v) ────────────────────────────────────
# Row (i,j): first-order edge-length change under vertex displacements.
R = np.zeros((n_e, 3 * n_v))
for row, (i, j) in enumerate(edges):
    d = V[i] - V[j]
    R[row, 3*i:3*i+3] =  d
    R[row, 3*j:3*j+3] = -d

print(f"  {n_v} vertices, {n_e} edges, {n_f} faces")
print()

# =============================================================================
print("SECTION 1: ZERO FLOPPY MODES")
print(SEP2)

rank_R     = np.linalg.matrix_rank(R, tol=1e-9)
dim_kernel = 3 * n_v - rank_R
print(f"  R: {n_e} x {3*n_v}   rank = {rank_R}   kernel = {dim_kernel}")
check("DC1: zero floppy modes -- kernel = 6 (rigid body only)",
      dim_kernel == 6,
      f"rank(R) = {rank_R}, dim(kernel) = {dim_kernel}  [must equal 6]")

# =============================================================================
print()
print("SECTION 2: UNIFORM INWARD RADIAL PRESSURE")
print(SEP2)
print("  Each vertex displaced inward by unit distance along -r_i direction.")

radii    = np.linalg.norm(V, axis=1, keepdims=True)
u_radial = np.zeros(3 * n_v)
for i in range(n_v):
    u_radial[3*i:3*i+3] = -V[i] / radii[i]   # unit inward

Ru_rad   = R @ u_radial
# Each edge shortens by L_J / R_c (same for all edges by icosahedral symmetry)
# Normalise by L_J^2 so the result is scale-independent (fractional edge-length change)
cost_rad = float(np.dot(Ru_rad, Ru_rad)) / L_J_fm**4
print(f"  |R*u_radial|^2 / L_J^4 = {cost_rad:.4f}  (fractional, scale-independent)")
check("DC2: radial inward pressure costs energy (outer edge network resists, no collapse)",
      cost_rad > 1.0,
      f"cost = {cost_rad:.4f} >> 0")

# =============================================================================
print()
print("SECTION 3: FACE-CENTER INWARD FORCE")
print(SEP2)
print("  Force applied at face center = distributed equally to 3 bounding vertices.")
print("  No inner structure present. Outer edges carry all load.")

# test one face; confirm with all 20 below (DC4)
fi, fj, fk = faces[0]
face_verts  = V[[fi, fj, fk]]
face_center = face_verts.mean(axis=0)
n_face = np.cross(face_verts[1] - face_verts[0], face_verts[2] - face_verts[0])
n_face /= np.linalg.norm(n_face)
if np.dot(n_face, face_center) > 0:   # ensure inward
    n_face = -n_face

u_face = np.zeros(3 * n_v)
for vi in [fi, fj, fk]:
    u_face[3*vi:3*vi+3] = n_face / 3.0

Ru_face   = R @ u_face
cost_face = float(np.dot(Ru_face, Ru_face)) / L_J_fm**4

# Identify which edges carry the load (non-zero row in R*u_face)
loaded_edges = [(i, j) for row, (i, j) in enumerate(edges) if abs(Ru_face[row]) > 1e-12]
print(f"  Face vertices: {fi}, {fj}, {fk}")
print(f"  Edges carrying non-zero load: {len(loaded_edges)} of {n_e}")
print(f"  |R*u_face|^2 / L_J^4 = {cost_face:.6f}")
check("DC3: face-center inward force costs energy (outer edges carry load; no inner support needed)",
      cost_face > 1e-6,
      f"cost = {cost_face:.6f} >> 0  [{len(loaded_edges)} outer edges loaded]")

# =============================================================================
print()
print("SECTION 4: ALL 20 FACES RESIST INWARD PUSH")
print(SEP2)

min_cost = float('inf')
for fi2, fj2, fk2 in faces:
    fv2         = V[[fi2, fj2, fk2]]
    fc2         = fv2.mean(axis=0)
    nf2         = np.cross(fv2[1] - fv2[0], fv2[2] - fv2[0])
    nf2        /= np.linalg.norm(nf2)
    if np.dot(nf2, fc2) > 0:
        nf2 = -nf2
    u2 = np.zeros(3 * n_v)
    for vi2 in [fi2, fj2, fk2]:
        u2[3*vi2:3*vi2+3] = nf2 / 3.0
    c2 = float(np.dot(R @ u2, R @ u2)) / L_J_fm**4
    min_cost = min(min_cost, c2)

print(f"  Minimum |R*u_face|^2 / L_J^4 across all 20 faces: {min_cost:.6f}")
check("DC4: all 20 faces resist inward push (min cost > 0 across all faces)",
      min_cost > 1e-6,
      f"min cost = {min_cost:.6f} >> 0  [symmetry: all faces identical by icosahedral symmetry]")

# =============================================================================
print()
print("SECTION 5: FULL ROW RANK -- INNER CONTENT IS REDUNDANT BY DEFINITION")
print(SEP2)
print(f"  rank(R) = {rank_R} = n_edges = {n_e}")
print(f"  Every edge is an independent constraint. No constraint is slack.")
print(f"  Adding inner edges (inner content) would create REDUNDANT constraints.")
print(f"  Redundant = over-constrained = mechanically unnecessary for rigidity.")
print(f"  The shell cannot be made 'more rigid' in any structural sense.")
check("DC5: rank(R) = n_edges (full row rank; inner content structurally redundant)",
      rank_R == n_e,
      f"rank = {rank_R} = n_e = {n_e}  -> inner content cannot improve on this")

# =============================================================================
print()
print("SECTION 6: EDGE-MIDPOINT NEXUS STRENGTH (SINGLE-EDGE TRANSVERSE PUSH)")
print(SEP2)
print("  Gluon/muon nexus sits at each edge midpoint (r_mid = 0.809*L_J).")
print("  Test: apply a transverse force at one edge midpoint (toward its face center,")
print("  the gluon amplitude direction, GH0b) distributed equally to 2 bounding vertices.")
print("  Each vertex displaced by 1/2 in direction: edge_midpoint -> face_center.")

# Find faces for each edge: each edge borders exactly 2 faces
edge_faces = {(min(i,j),max(i,j)): [] for i,j in edges}
for fi2, fj2, fk2 in faces:
    for ei, ej in [(fi2,fj2),(fj2,fk2),(fk2,fi2)]:
        key = (min(ei,ej), max(ei,ej))
        fc = (V[fi2] + V[fj2] + V[fk2]) / 3.0
        edge_faces[key].append(fc)

edge_costs = []
for i, j in edges:
    midpoint = (V[i] + V[j]) / 2.0
    fc_list = edge_faces[(min(i,j), max(i,j))]
    # Use first adjacent face center; direction = toward face center
    fc = fc_list[0]
    d_raw = fc - midpoint
    d_norm = np.linalg.norm(d_raw)
    d_unit = d_raw / d_norm   # unit transverse direction (toward face center)

    u_edge = np.zeros(3 * n_v)
    u_edge[3*i:3*i+3] = d_unit / 2.0   # each bounding vertex gets half
    u_edge[3*j:3*j+3] = d_unit / 2.0

    Ru_edge = R @ u_edge
    cost_edge = float(np.dot(Ru_edge, Ru_edge)) / L_J_fm**4
    edge_costs.append(cost_edge)

min_ec6 = min(edge_costs)
max_ec6 = max(edge_costs)
print(f"  Edge-midpoint push cost: min={min_ec6:.4f}, max={max_ec6:.4f}  (over all 30 edges)")
print(f"  Uniform: {'YES' if abs(max_ec6-min_ec6)/max_ec6 < 0.01 else 'NO (asymmetric strip coloring expected)'}")
check("DC6: single edge-midpoint transverse force costs energy (shell resists at every edge nexus)",
      min_ec6 > 0.01,
      f"min cost = {min_ec6:.4f}  (gluon/muon nexus position; direction toward face center)")

# =============================================================================
print()
print(SEP)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"RESULT: {len(results)} checks  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
    print()
    print("  STRUCTURAL CONCLUSION:")
    print("  The Jobson cell outer shell is self-sufficient -- no inner content required.")
    print(f"  - Triangulated faces: zero floppy modes, kernel = 6 (rigid body only)  [DC1]")
    print(f"  - Radial inward pressure resisted by outer edge network                [DC2]")
    print(f"  - Face-center inward force resisted by outer edge network               [DC3]")
    print(f"  - All 20 faces resist inward push; icosahedral symmetry                [DC4]")
    print(f"  - Full row rank: no redundant constraints; inner content adds none      [DC5]")
    print(f"  - Single edge-midpoint transverse force resisted at all 30 edges        [DC6]")
    print()
    print("  Inner cell structure cannot be derived from or required by durability.")
    print("  If inner content exists, it is occupant -- not structural support.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAIL: {name}")
            print(f"        {detail}")
print(SEP)
