"""
mesh_grind_pybullet_icosahedron_validation.py
================================================
MIGRATION START (2026-09-16): first script in the PyBullet port of the mesh/grind
chirality-derivation line. Purpose: validate that a PyBullet convex-hull collision
shape built from the SAME 12-vertex Jobson-cell icosahedron template already used
throughout the numpy engine (mesh_grind_single_driver_parallel_sweep.py, cell_
rotation_propagation.py) reproduces the SAME three independently-known touching
distances already cross-validated there (V2-0e in mesh_grind_cluster_relaxation_v2.py:
"SAT gap is exactly 0 at each of 3 independently-known touching distances... at ALL
three") -- BEFORE trusting anything built on top of this new engine.

WHY MIGRATE AT ALL (author, 2026-09-15/16): the hand-rolled numpy SAT + Gauss-Seidel
rigid-overlap engine has a persistent, never-fixed ~2.05*r_in penetration bug in the
full two-cluster dynamical loop (root-cause diagnostic isolated the corrector and the
static geometry as NOT the cause -- see mesh_grind_penetration_root_cause_diagnostic.py
and /memories/repo/mesh-grind-chirality-derivation.md). PyBullet is a mature, widely-
used rigid-body solver with real contact resolution (LCP/PGS or Bullet's own solver) --
the hope is this class of bug simply does not occur with a real physics engine, and
that a well-modeled rigid Zone-2 driver body (a genuine compound rigid body, not the
"locked group / shared-omega-only" approximation flagged as physically incoherent last
session -- "group 5 is impossible -- how do you spin 5 icosahedra together as one
blob") becomes straightforward once PyBullet handles contacts/constraints natively.

CRITICAL DESIGN DECISION -- UNIT RESCALING (read before extending this script):
  Jobson-cell length scales are ~1e-17 m (L_J ~= 9.94e-18 m). The standard pybullet
  wheel is a SINGLE-PRECISION (float32) build of Bullet -- feeding it raw femtometre
  numbers directly would immediately underflow all meaningful precision (float32 has
  ~7 significant decimal digits; contact margins, solver tolerances, and the default
  collision margin (~0.04 units) are all tuned assuming O(1)-scale geometry). FIX:
  every length fed to PyBullet is rescaled by SIM_SCALE = 1.0 / L_J, so "1 PyBullet
  unit = 1 L_J" -- shapes, positions, and velocities are all defined in this rescaled
  system; any measured PyBullet distance is converted back to physical fm by
  multiplying by L_J (dividing by SIM_SCALE) before comparing to established constants
  or printing. This is a NEW shortcut/consideration that did not exist in the numpy
  engine (which worked directly in raw fm since numpy float64 has no such problem) --
  flagged explicitly, not silently baked in.

EXPECTED, HONEST CAVEAT (do not expect exact-to-1e-9 agreement like the numpy SAT
engine achieved): PyBullet's GJK/EPA narrow-phase collision applies a default
COLLISION MARGIN around every convex shape (~0.04 units, i.e. ~0.04*L_J here) to keep
the solver numerically robust -- this shifts the effective "touching" separation
measured via getClosestPoints outward by roughly that margin, and there is a small
inherent numerical fuzz beyond that from the reduced float32 precision. The checks
below verify the SIGN and rough MAGNITUDE match the established face/edge/vertex
touching distances to a few percent of L_J, not exact agreement -- a materially
different (looser) tolerance than the numpy engine's own 1e-16-level self-consistency
check, and that difference is the honest cost of switching to a real physics engine.

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_icosahedron_validation.py
"""

import math
import numpy as np
import pybullet as p

pi = math.pi
alpha = 7.2973525693e-3
hbar_c = 197.3269804
m_p = 938.272
lambda_p = hbar_c / m_p
r_p_fm = 4 * lambda_p
phi = (1 + math.sqrt(5)) / 2
L_J = alpha * phi * r_p_fm            # fm
R_c = L_J * math.sqrt(1 + phi**2) / 2  # circumradius, fm
r_in = L_J * phi**2 / (2 * math.sqrt(3))  # inradius (face), fm
r_mid = L_J * phi / 2                  # midradius (edge), fm

SIM_SCALE = 1.0 / L_J  # "1 PyBullet unit = 1 L_J" -- see docstring above


def norm(v):
    v = np.array(v, dtype=float)
    nv = np.linalg.norm(v)
    return v / nv if nv > 1e-15 else v


# same raw vertex template as mesh_grind_single_driver_parallel_sweep.py
scale_fm = L_J / 2.0
raw = []
for s1 in (+1, -1):
    for s2 in (+1, -1):
        raw.append((0, s1, s2 * phi))
        raw.append((s1, s2 * phi, 0))
        raw.append((s1 * phi, 0, s2))
verts_fm = [scale_fm * np.array(v, dtype=float) for v in raw]  # 12 verts, physical fm
verts_sim = [v * SIM_SCALE for v in verts_fm]                  # same verts, PyBullet units
n_t = len(verts_fm)

# rebuild edges/faces (same method as the numpy engine) to get real face-centroid and
# edge-midpoint DIRECTIONS -- not assumed, derived from the same vertex template.
pair_dists = []
for i in range(n_t):
    for j in range(i + 1, n_t):
        d = np.linalg.norm(verts_fm[i] - verts_fm[j])
        pair_dists.append((d, i, j))
pair_dists.sort(key=lambda x: x[0])
edge_len = pair_dists[0][0]
edges = [(i, j) for d, i, j in pair_dists if abs(d - edge_len) < 1e-9]
edge_set = set(edges) | set((j, i) for (i, j) in edges)

faces = []
for i in range(n_t):
    for j in range(i + 1, n_t):
        if (i, j) not in edge_set:
            continue
        for k in range(j + 1, n_t):
            if (i, k) in edge_set and (j, k) in edge_set:
                faces.append((i, j, k))

face_centroid_dir = norm(sum(verts_fm[i] for i in faces[0]))
edge_midpoint_dir = norm(verts_fm[edges[0][0]] + verts_fm[edges[0][1]])
vertex_dir = norm(verts_fm[0])


def make_icosahedron_body(client, position_sim, orientation=(0, 0, 0, 1), mass=0.0):
    """One rigid body with a convex-hull collision shape from the 12-vertex template
    (vertices-only -> PyBullet computes the convex hull internally)."""
    col = p.createCollisionShape(p.GEOM_MESH, vertices=verts_sim, physicsClientId=client)
    body = p.createMultiBody(baseMass=mass, baseCollisionShapeIndex=col,
                              basePosition=list(position_sim), baseOrientation=list(orientation),
                              physicsClientId=client)
    return body


def check_separation(client, label, direction_dir, expected_fm, tol_frac=0.08):
    """Places body B along `direction_dir` from body A (both identity orientation, both
    at the origin's local frame) at the EXPECTED touching separation (converted to sim
    units), then reports PyBullet's own measured closest-point gap there, plus at +/-10%
    of that separation, to confirm sign convention (positive=separated, negative=
    penetrating) and rough magnitude -- not exact agreement, see docstring caveat."""
    bodyA = make_icosahedron_body(client, [0, 0, 0])
    results = []
    for frac, tag in [(0.90, "-10%"), (1.00, "exact"), (1.10, "+10%")]:
        sep_fm = expected_fm * frac
        posB_sim = np.array(direction_dir) * sep_fm * SIM_SCALE
        bodyB = make_icosahedron_body(client, posB_sim)
        p.performCollisionDetection(physicsClientId=client)
        pts = p.getClosestPoints(bodyA, bodyB, distance=10.0, physicsClientId=client)
        gap_sim = min(pt[8] for pt in pts) if pts else float("nan")  # index 8 = contactDistance
        gap_fm = gap_sim / SIM_SCALE
        results.append((tag, sep_fm, gap_fm))
        p.removeBody(bodyB, physicsClientId=client)
    p.removeBody(bodyA, physicsClientId=client)

    print(f"\n{label}: expected touching separation = {expected_fm:.6f} fm")
    for tag, sep_fm, gap_fm in results:
        print(f"  sep={sep_fm:.6f} fm ({tag:>5s}) -> PyBullet closest-point gap = {gap_fm:+.6f} fm")

    exact_gap_fm = results[1][2]
    margin_frac = abs(exact_gap_fm) / expected_fm
    ok_sign_order = results[0][2] < results[1][2] < results[2][2]
    ok_margin = margin_frac < tol_frac
    passed = ok_sign_order and ok_margin
    print(f"  gap at exact separation = {exact_gap_fm:+.6f} fm "
          f"({margin_frac*100:.2f}% of touching distance -- collision-margin + "
          f"float32 fuzz, expected nonzero, tol={tol_frac*100:.0f}%)")
    print(f"  monotonic gap (-10% < exact < +10%): {ok_sign_order}")
    print(f"  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed


def main():
    client = p.connect(p.DIRECT)
    p.setPhysicsEngineParameter(numSolverIterations=50, physicsClientId=client)

    all_pass = True
    all_pass &= check_separation(client, "V1: face-to-face (2*r_in)", face_centroid_dir, 2 * r_in)
    all_pass &= check_separation(client, "V2: edge-to-edge (2*r_mid)", edge_midpoint_dir, 2 * r_mid)
    all_pass &= check_separation(client, "V3: vertex-to-vertex (2*R_c)", vertex_dir, 2 * R_c)

    print(f"\n{'='*60}")
    print(f"OVERALL: {'3/3 PASS' if all_pass else 'SOME FAILED'}")
    print(f"{'='*60}")
    p.disconnect(client)
    return all_pass


if __name__ == "__main__":
    main()
