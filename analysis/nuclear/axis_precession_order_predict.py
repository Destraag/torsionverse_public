"""
axis_precession_order_predict.py
==========================================
Author's question (2026-09-11), asked mid-Phase-C: given that the current
cluster scripts do NOT actually rotate a cell's orientation (Omega is used
only as an instantaneous kinematic quantity for contact-point velocity,
never integrated into an accumulated rotation -- confirmed by grep, no
rotation matrix/quaternion anywhere in mesh_grind_cluster_spin_test.py or
mesh_grind_two_cluster_test.py), what ORDER of the 6 discrete C5 axes
would a real, continuously-precessing rigid body naturally visit?

METHOD: compute the exact pairwise angles between all 6 C5 axes (derived
from the same icosahedron vertex template used throughout this whole
mesh/grind effort). For a regular icosahedron, vertex-pair angles take
only 2 possible values: ~63.43 deg (edge-adjacent) or ~116.57 deg
(antipode-of-edge-adjacent) -- verified computationally below, not assumed.
Then construct a nearest-neighbor tour (always step to the closest
not-yet-visited axis) as a REASONED PREDICTION for what a smoothly-
precessing spin axis would trace out -- explicitly flagged as a
hypothesis about a smooth path through orientation space, NOT a rigorous
derivation of real gyroscopic/torque-driven precession dynamics (which
would require modeling the actual constraint torques from neighboring
cells -- not attempted here).

Run: python analysis/nuclear/axis_precession_order_predict.py
Reference: /memories/repo/mesh-grind-chirality-derivation.md,
  analysis/nuclear/mesh_grind_discrete_cell_axes.py (where the 6 C5 axes were first enumerated)
"""

import sys, math
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

phi = (1 + math.sqrt(5)) / 2
SEP = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

def norm(v):
    v = np.array(v, dtype=float)
    return v / np.linalg.norm(v)

raw = []
for s1 in (+1, -1):
    for s2 in (+1, -1):
        raw.append((0, s1, s2 * phi))
        raw.append((s1, s2 * phi, 0))
        raw.append((s1 * phi, 0, s2))
verts = [norm(v) for v in raw]

seen, c5_axes = [], []
for v in verts:
    if not any(np.allclose(v, a, atol=1e-9) or np.allclose(v, -a, atol=1e-9) for a in seen):
        seen.append(v)
        c5_axes.append(v)

print(SEP)
print("SECTION 1: EXACT PAIRWISE ANGLES BETWEEN THE 6 C5 AXES")
print(SEP2)
n = len(c5_axes)
angle_matrix = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        d = np.clip(np.dot(c5_axes[i], c5_axes[j]), -1, 1)
        angle_matrix[i, j] = math.degrees(math.acos(d))

print("  Angle matrix (degrees), rows/cols = axis 0-5:")
print("        " + "".join(f"{j:>8}" for j in range(n)))
for i in range(n):
    print(f"    {i:>3} " + "".join(f"{angle_matrix[i,j]:>8.2f}" for j in range(n)))

off_diag = angle_matrix[~np.eye(n, dtype=bool)]
distinct_angles = sorted(set(round(a, 1) for a in off_diag))
print(f"\n  Distinct off-diagonal angles found: {distinct_angles}")

check("AP-1: exactly 2 distinct non-zero angles exist between C5 axis pairs "
      "(the two characteristic icosahedral vertex angles: ~63.43 and ~116.57 deg)",
      len(distinct_angles) == 2,
      f"found {len(distinct_angles)} distinct angles: {distinct_angles}")

expected = [round(math.degrees(math.acos(1 / math.sqrt(5))), 1), round(180 - math.degrees(math.acos(1 / math.sqrt(5))), 1)]
check("AP-2: the 2 angles match the exact algebraic values arccos(1/sqrt(5)) and its supplement "
      "(the standard icosahedron vertex central angle, ~63.43/116.57 deg)",
      distinct_angles == sorted(expected),
      f"found {distinct_angles}, expected {sorted(expected)}")

print()
print(SEP)
print("SECTION 2: NEAREST-NEIGHBOR TOUR -- A REASONED PREDICTION, NOT A DERIVED TRAJECTORY")
print(SEP2)
print("""
  HONEST FRAMING: this is a PREDICTION about what a smoothly-varying (continuously
  precessing) spin axis might visit, based on always moving to the geometrically
  CLOSEST not-yet-visited axis -- a plausible hypothesis for "smooth", not a
  rigorously derived gyroscopic/torque trajectory (that would require modeling
  the actual constraint torques from neighboring cells, not attempted here or
  anywhere in this repo yet).
""")

visited = [0]
remaining = set(range(1, n))
while remaining:
    last = visited[-1]
    nxt = min(remaining, key=lambda j: angle_matrix[last, j])
    visited.append(nxt)
    remaining.remove(nxt)

print(f"  Predicted order (nearest-neighbor tour, starting at axis 0): {visited}")
step_angles = [angle_matrix[visited[i], visited[i+1]] for i in range(len(visited)-1)]
print(f"  Step angles along this path: {[f'{a:.2f}' for a in step_angles]} deg")

check("AP-3: the nearest-neighbor tour only ever takes the SMALLER of the 2 possible angles "
      "(~63.43 deg) for every step -- a genuinely 'smooth' path exists using only edge-adjacent "
      "steps, not the larger 116.57 deg jumps",
      all(abs(a - min(distinct_angles)) < 0.5 for a in step_angles),
      f"step angles: {[round(a,2) for a in step_angles]}, smallest available = {min(distinct_angles)}")

print(f"""
  PREDICTION FOR THE AUTHOR'S QUESTION: if mesh_grind_cluster_spin_test.py /
  mesh_grind_two_cluster_test.py are upgraded to track REAL orientation (a rotation
  matrix/quaternion integrated from Omega, not just the current kinematic-only
  Omega), a single continuously-driven spin about ONE fixed lab-frame axis would
  already sweep the OTHER 5 axes through a cone around it automatically -- no
  manual "order" would be needed at all. If instead the intent is to manually
  step through a SEQUENCE of distinct axis choices (as the current scripts do,
  resetting between each), the order found above ({visited}) is the smoothest
  available path (every step the minimum 63.43 deg, never a 116.57 deg jump) --
  this is offered as a reasoned hypothesis for "the order that would look most
  like continuous precession", not a proven physical trajectory.
""")

print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
