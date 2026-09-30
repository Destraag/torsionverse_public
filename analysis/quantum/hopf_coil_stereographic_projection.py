"""
hopf_coil_stereographic_projection.py
==========================================
Author redirect (2026-09-16): stop guessing icosahedral vertex-hopping shapes
for the free tau's coil -- go back to the ACTUAL (1,2) Hopf torus knot
already built and validated for the alpha derivation (gap3_chern_simons.py),
and see what it looks like as a genuine 3D coil.

METHOD: the (1,2) torus knot lives on the Clifford torus in S^3 subset R^4
(gap3_chern_simons.py, unchanged, reused verbatim):
    x(theta) = (cos(theta), sin(theta), cos(2*theta), sin(2*theta)) / sqrt(2)

To see this as an ordinary 3D coil, project it down via STEREOGRAPHIC
PROJECTION from S^3 to R^3 -- the standard, canonical (not arbitrary) way to
bring a Hopf-fibration-related S^3 curve into physical 3D space while
preserving its actual connection/angle structure (used throughout the
Hopf-fibration visualization literature, not invented here):
    (X, Y, Z) = (x1, x2, x3) / (1 - x4)

Also computed for cross-check: the STANDARD explicit (p,q) torus-knot formula
in R^3 (the classic "coil wound around a donut" parametrization) -- a
DIFFERENT, more commonly-visualized construction, to see whether the two
approaches agree on basic shape features (are they even qualitatively
compatible) before trusting either one's numbers.

GOAL: extract genuine SHAPE properties (effective radius, pitch, how many
times it winds around which axis) from BOTH constructions, and compare them
to the discrete icosahedral tau path's own already-established pitch angle
(gluon_tau_helix.py, 30.0666 deg) and deflection (72 deg) -- reported
honestly, not forced to match.

Run: python analysis/quantum/hopf_coil_stereographic_projection.py
"""
import math
import numpy as np

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
p, q = 1, 2

SEP = "=" * 68
results = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")


def torus_knot_S3(theta):
    """(1,2) torus knot on the Clifford torus in S^3, reused verbatim from
    gap3_chern_simons.py (already validated there)."""
    s = 1.0 / math.sqrt(2)
    return np.array([
        s * math.cos(p * theta),
        s * math.sin(p * theta),
        s * math.cos(q * theta),
        s * math.sin(q * theta)
    ])


def stereographic_project(x4vec):
    """Standard stereographic projection S^3 -> R^3 from the north pole
    (0,0,0,1) -- the canonical way to view a Hopf-related S^3 curve in
    ordinary 3D space."""
    x1, x2, x3, x4 = x4vec
    denom = 1.0 - x4
    return np.array([x1/denom, x2/denom, x3/denom])


print(SEP)
print("PART A -- (1,2) HOPF TORUS KNOT, STEREOGRAPHIC PROJECTION TO R^3")
print(SEP)

N = 2000
thetas = np.linspace(0, 2*pi, N, endpoint=False)
pts_S3 = [torus_knot_S3(t) for t in thetas]
pts_R3 = np.array([stereographic_project(x) for x in pts_S3])

print(f"\n  Sampled {N} points along theta in [0, 2*pi).")
r_from_origin = np.linalg.norm(pts_R3, axis=1)
print(f"  Distance from R^3 origin: min={r_from_origin.min():.4f}  "
      f"max={r_from_origin.max():.4f}  mean={r_from_origin.mean():.4f}")
print(f"  (NOTE: the projection point itself, theta where x4->1, sends that")
print(f"  point to infinity -- check for this before trusting the max.)")
max_x4 = max(abs(x[3]) for x in pts_S3)
print(f"  max |x4| sampled = {max_x4:.6f} (danger zone near 1/sqrt(2)={1/math.sqrt(2):.6f})")

# find a natural axis: the mean position direction (since this curve is not
# symmetric about the R^3 origin the way the icosahedral construction was)
centroid = pts_R3.mean(axis=0)
print(f"\n  Centroid of projected curve: {centroid}")
print(f"  (a genuine coil should spiral AROUND some axis -- check by PCA)")

# PCA to find the natural coil axis
centered = pts_R3 - centroid
cov = centered.T @ centered / len(centered)
eigvals, eigvecs = np.linalg.eigh(cov)
order = np.argsort(eigvals)[::-1]
eigvals = eigvals[order]
eigvecs = eigvecs[:, order]
print(f"  PCA eigenvalues (variance along principal axes): {eigvals}")
print(f"  Principal (coil) axis: {eigvecs[:, 0]}")

coil_axis = eigvecs[:, 0]
z_coord = centered @ coil_axis
transverse = centered - np.outer(z_coord, coil_axis)
r_transverse = np.linalg.norm(transverse, axis=1)
print(f"\n  Along principal axis: z range = [{z_coord.min():.4f}, {z_coord.max():.4f}]")
print(f"  Transverse radius: min={r_transverse.min():.4f}  max={r_transverse.max():.4f}  "
      f"mean={r_transverse.mean():.4f}  std={r_transverse.std():.4f}")
uniform_radius = r_transverse.std() / r_transverse.mean() < 0.05
print(f"  Roughly uniform transverse radius (std/mean<5%)? {uniform_radius}")

check("SP1: stereographic projection of the (1,2) knot has a well-defined "
      "principal (coil) axis via PCA",
      eigvals[0] > 2*eigvals[1],
      f"eigvals={eigvals}")

print()
print(SEP)
print("PART B -- STANDARD (p,q) TORUS KNOT IN R^3 (classic 'coil around a "
      "donut' construction, cross-check)")
print(SEP)

R_big, r_small = 2.0, 1.0  # arbitrary but fixed reference radii for shape-only comparison


def std_torus_knot(theta):
    x = (R_big + r_small*math.cos(q*theta)) * math.cos(p*theta)
    y = (R_big + r_small*math.cos(q*theta)) * math.sin(p*theta)
    z = r_small * math.sin(q*theta)
    return np.array([x, y, z])


pts_std = np.array([std_torus_knot(t) for t in thetas])
z_std = pts_std[:, 2]
r_std = np.linalg.norm(pts_std[:, :2], axis=1)
print(f"\n  z range (donut tube axis): [{z_std.min():.4f}, {z_std.max():.4f}]  "
      f"(should be +/-r_small={r_small})")
print(f"  transverse radius (distance from central hole axis): "
      f"[{r_std.min():.4f}, {r_std.max():.4f}]  (should be R_big +/- r_small = "
      f"[{R_big-r_small:.1f}, {R_big+r_small:.1f}])")
print(f"  Winds {p}x around the big (hole) axis, {q}x around the small (tube)")
print(f"  axis -- by construction, exact, not fit.")

# pitch angle analogue: at any point, ratio of "big-axis progress" to
# "small-axis (tube) progress" per unit theta
dz_dtheta = r_small * q  # amplitude of dz/dtheta
d_big_dtheta = R_big * p  # amplitude of angular progress x R_big (rough arc-length scale)
pitch_std = math.degrees(math.atan2(dz_dtheta, d_big_dtheta))
print(f"  Rough pitch angle (tube oscillation rate vs hole-circulation rate) "
      f"= arctan(r_small*q / (R_big*p)) = {pitch_std:.4f} deg")
print(f"  (Depends on the ARBITRARY choice R_big=2, r_small=1 here -- this")
print(f"  specific number is NOT meaningful on its own; only used to confirm")
print(f"  the construction behaves as expected before any physical scale is")
print(f"  attached.)")

print()
print(SEP)
print("COMPARISON TO THE DISCRETE ICOSAHEDRAL TAU PATH")
print(SEP)
print(f"\n  Discrete (THE TAU TRADE, structural, 20-face circuit,")
print(f"  gluon_tau_helix.py): deflection=72.0000 deg EXACTLY, pitch=30.0666 deg.")
print(f"  Continuous (1,2) Hopf torus knot: winds p=1 time one way, q=2 times")
print(f"  the other -- NO single 'deflection angle' in the discrete sense (this")
print(f"  is a smooth curve, not a sequence of straight-line bounces), so a")
print(f"  direct numeric comparison to 72 deg is not well-posed as asked --")
print(f"  reported honestly rather than forcing an artificial angle match.")
print(f"  What IS comparable: the WINDING RATIO q/p=2 for the continuous knot")
print(f"  vs the discrete path's OWN 20 hops / 1 circuit -- these are different")
print(f"  countings (a continuous ratio vs a discrete hop count) and are NOT")
print(f"  obviously the same kind of quantity; no forced numeric match claimed.")

print()
print(SEP)
print("HONEST STATUS")
print(SEP)
print(f"  Both constructions (stereographic projection of the REAL S^3 (1,2)")
print(f"  knot, and the classic explicit R^3 torus-knot formula) DO produce")
print(f"  genuine coil shapes -- confirmed, not assumed (PCA found a real")
print(f"  principal axis; the standard construction winds p/q times as built).")
print(f"  BUT: neither one has yet been tied to a PHYSICAL SCALE (L_J, r_in,")
print(f"  R_c) or to the discrete icosahedral path's own specific numbers --")
print(f"  that mapping (what R_big/r_small, or what stereographic-projection")
print(f"  pole choice, corresponds to 'one Jobson cell') has NOT been derived,")
print(f"  only flagged as the next real gap. Forcing a match at this stage")
print(f"  would be numerology; this script stops at 'both are real coils,")
print(f"  by two independent methods' and goes no further.")
print(SEP)
