"""
mesh_grind_pybullet_zone2_driver_test.py
===========================================
SECOND script in the PyBullet migration (see mesh_grind_pybullet_icosahedron_
validation.py for the geometry-validation step this builds on -- 3/3 PASS, PyBullet's
convex-hull collision shape reproduces the established face/edge/vertex touching
distances to <0.2% of L_J). This script builds the actual target the migration was
started for: A GENUINELY RIGID, WELDED ZONE-2 DRIVER BODY, replacing last session's
"locked group / shared-omega-only" approximation that the author rejected as
physically incoherent ("group 5 is impossible -- how do you spin 5 icosahedra
together as one blob").

ZONE-2 DRIVER: a real PyBullet COMPOUND rigid body -- one base cell (center) plus 12
shell cells attached via JOINT_FIXED links (0 relative DOF, so the whole 13-cell
cluster moves and rotates as a single rigid object by construction, not an
approximation). 13 = 1 + kissing number (12, the max spheres touching one central
sphere in 3D) -- reusing the SAME 12-vertex icosahedral direction set already
established throughout this repo (cog_geometry.py, the C5-axis template) as the
principled shell-placement directions, per last session's author-endorsed reasoning.
Shell radius = 2*r_in (the established face-touching distance, V1 in the validation
script) -- a "packed shell" choice, NOT a claim of exact space-filling tiling (real
icosahedra do not tile 3D space by simple repetition -- that broader packing
question remains open repo-wide, tracked internally).

DRIVER KINEMATICS (author-mandated convention, unchanged from the numpy engine):
ROTATION IS PRESCRIBED, POSITION IS FULLY FREE. Every step, the driver's angular
velocity is reset to a fixed Omega (kinematic driving), while its linear velocity is
read back from the solver and resupplied unchanged -- i.e. NOTHING holds its position
fixed; it responds to ordinary contact forces exactly like any other dynamic body.

MEDIUM + CONFINEMENT: N_FREE ordinary dynamic icosahedra seeded in a ball around the
driver (small initial overlaps expected and left for PyBullet's OWN solver to resolve
-- this is the actual point of the migration, letting a mature contact solver do what
the numpy engine's hand-rolled Gauss-Seidel correction could not: converge to ~0
one-cluster penetration inside a live dynamical loop with continuous forcing). A
constant-magnitude inward force (-P_conf * direction-to-center) stands in for "one
common external bounding pressure" (the corrected architecture from last session --
NOT two separate clusters, one shared confined medium with an embedded driver).

SCOPE FOR THIS FIRST RUN (deliberately NOT yet the full N~830 kissing-number-scaled
real-Zone-2-radius target from last session's driver-group critique): N_FREE=100,
matching the scale already used/validated in the numpy engine's own incremental
history (N=24 -> N=150) -- prove the compound-body + solver combination is stable and
low-penetration FIRST, then scale N up as a separate, explicit next step.

MEASUREMENTS: (1) max penetration depth across the WHOLE run (via getContactPoints'
contactDistance, all pairs) -- the direct test of whether PyBullet's real solver
avoids the numpy engine's persistent ~2.05*r_in bug; (2) spread (mean free-cell
displacement from start); (3) local spacing change near the driver (thinning/
expansion test, same diagnostic the numpy sweep used).

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_zone2_driver_test.py
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
L_J = alpha * phi * r_p_fm
R_c = L_J * math.sqrt(1 + phi**2) / 2
r_in = L_J * phi**2 / (2 * math.sqrt(3))
r_mid = L_J * phi / 2

SIM_SCALE = 1.0 / L_J  # "1 PyBullet unit = 1 L_J" -- see validation script docstring

N_FREE = 100
N_STEPS = 400
CELL_MASS = 1.0          # placeholder, arbitrary units -- real mass scale not yet derived
# P_CONF: the original 0.05 placeholder was checked against the sim's own contact-
# stiffness force scale (jobson_cell_natural_force_pressure_scale.py) and found almost
# certainly too weak to be a meaningful background confinement relative to contacts
# (contactStiffness=1e6 dwarfs it by many orders of magnitude) -- raised to a still-
# moderate, not-yet-fully-pinned-down value, clearly flagged as informed-but-provisional.
P_CONF = 2.0
OMEGA_MAG = 0.3          # placeholder driver angular speed, rad per sim step-ish (arbitrary units)
# LINEAR/ANGULAR DAMPING: added after the first stiffened-solver run ejected cells
# violently (spacing +1076%, 56 r_in mean displacement) -- the likely cause is a stiff
# contact solver converting corrective overlap into velocity with nothing to dissipate
# it. Real materials have damping; the numpy engine had none either (a standing,
# unaddressed gap there too) -- added explicitly here rather than assumed away.
LINEAR_DAMPING = 0.4
ANGULAR_DAMPING = 0.4
SEED_RADIUS_FM = 10 * r_in
rng = np.random.RandomState(42)


def norm(v):
    v = np.array(v, dtype=float)
    nv = np.linalg.norm(v)
    return v / nv if nv > 1e-15 else v


scale_fm = L_J / 2.0
raw = []
for s1 in (+1, -1):
    for s2 in (+1, -1):
        raw.append((0, s1, s2 * phi))
        raw.append((s1, s2 * phi, 0))
        raw.append((s1 * phi, 0, s2))
verts_fm = [scale_fm * np.array(v, dtype=float) for v in raw]
verts_sim = [v * SIM_SCALE for v in verts_fm]

# 12 distinct vertex DIRECTIONS (already deduplicated -- raw has each direction once,
# no antipodal duplicates by construction of the (0,+-1,+-phi)-type template)
shell_directions = [norm(v) for v in verts_fm]
assert len(shell_directions) == 12


def make_ico_shape(client):
    return p.createCollisionShape(p.GEOM_MESH, vertices=verts_sim, physicsClientId=client)


def build_zone2_driver(client, ico_shape, position_sim):
    """ONE compound rigid body: base (center cell) + 12 shell cells on JOINT_FIXED
    links -- a genuine rigid weld (0 relative DOF), not a shared-omega approximation."""
    n_shell = 12
    link_masses = [CELL_MASS] * n_shell
    link_collision_shapes = [ico_shape] * n_shell
    link_visual_shapes = [-1] * n_shell
    link_positions = [list(np.array(d) * (2 * r_in * SIM_SCALE)) for d in shell_directions]
    link_orientations = [[0, 0, 0, 1]] * n_shell
    link_inertial_positions = [[0, 0, 0]] * n_shell
    link_inertial_orientations = [[0, 0, 0, 1]] * n_shell
    link_parent_indices = [0] * n_shell
    link_joint_types = [p.JOINT_FIXED] * n_shell
    link_joint_axis = [[0, 0, 0]] * n_shell

    body = p.createMultiBody(
        baseMass=CELL_MASS,
        baseCollisionShapeIndex=ico_shape,
        baseVisualShapeIndex=-1,
        basePosition=list(position_sim),
        baseOrientation=[0, 0, 0, 1],
        linkMasses=link_masses,
        linkCollisionShapeIndices=link_collision_shapes,
        linkVisualShapeIndices=link_visual_shapes,
        linkPositions=link_positions,
        linkOrientations=link_orientations,
        linkInertialFramePositions=link_inertial_positions,
        linkInertialFrameOrientations=link_inertial_orientations,
        linkParentIndices=link_parent_indices,
        linkJointTypes=link_joint_types,
        linkJointAxis=link_joint_axis,
        physicsClientId=client,
    )
    return body


def sample_uniform_ball(n, r_max, rng):
    u = rng.random(n)
    r = r_max * u ** (1.0 / 3.0)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return r[:, None] * d


def main():
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, 0, physicsClientId=client)
    # steady-state penetration under a SUSTAINED confinement force is expected with
    # Bullet's default soft/compliant contact solver (a continuous inward force will
    # always produce some equilibrium overlap unless contacts are made close to rigid)
    # -- tighten solver/contact parameters rather than accept the default softness,
    # since a real rigid Zone-2 boundary should not overlap under an ordinary applied
    # pressure (a genuine yield/jamming threshold is a separate, not-yet-modeled
    # question -- open_items.txt N-12 -- not what's being tested here).
    p.setPhysicsEngineParameter(numSolverIterations=100, numSubSteps=4,
                                 fixedTimeStep=1.0 / 60.0, erp=0.9, contactERP=0.9,
                                 globalCFM=1e-7, contactSlop=0.0001,
                                 physicsClientId=client)

    ico_shape = make_ico_shape(client)
    driver_id = build_zone2_driver(client, ico_shape, [0, 0, 0])

    # driver's real max radial footprint: shell-center offset (2*r_in) + a shell cell's
    # own circumradius (R_c) + a small safety margin -- properly RESAMPLE (not just
    # nudge) any free cell that lands inside this, so nothing starts pre-overlapped
    driver_footprint_fm = 2 * r_in + R_c + r_in
    free_start_fm = sample_uniform_ball(N_FREE, SEED_RADIUS_FM, rng)
    dist_from_origin = np.linalg.norm(free_start_fm, axis=1)
    too_close = dist_from_origin < driver_footprint_fm
    n_resampled = 0
    while np.any(too_close):
        n_bad = int(too_close.sum())
        n_resampled += n_bad
        free_start_fm[too_close] = sample_uniform_ball(n_bad, SEED_RADIUS_FM, rng)
        dist_from_origin = np.linalg.norm(free_start_fm, axis=1)
        too_close = dist_from_origin < driver_footprint_fm
    print(f"resampled {n_resampled} initial free-cell positions clear of the driver footprint "
          f"({driver_footprint_fm:.6f} fm)")

    # continuous collision detection: without this, a shape moving more than its own
    # size in one timestep can tunnel through a neighbor before discrete detection
    # sees it -- cheap insurance against exactly the kind of penetration this test is
    # meant to diagnose, not assumed unnecessary
    ccd_radius_sim = r_in * SIM_SCALE * 0.5
    for bid in [driver_id]:
        p.changeDynamics(bid, -1, ccdSweptSphereRadius=ccd_radius_sim,
                          linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                          physicsClientId=client)
        for link in range(12):
            p.changeDynamics(bid, link, ccdSweptSphereRadius=ccd_radius_sim,
                              linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                              physicsClientId=client)

    free_ids = []
    for i in range(N_FREE):
        pos_sim = free_start_fm[i] * SIM_SCALE
        body = p.createMultiBody(baseMass=CELL_MASS, baseCollisionShapeIndex=ico_shape,
                                  baseVisualShapeIndex=-1, basePosition=list(pos_sim),
                                  baseOrientation=[0, 0, 0, 1], physicsClientId=client)
        p.changeDynamics(body, -1, ccdSweptSphereRadius=ccd_radius_sim,
                          linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                          physicsClientId=client)
        free_ids.append(body)

    omega_vec = OMEGA_MAG * norm(np.array([1.0, 1.0, 1.0]))  # placeholder fixed axis

    start_pos_sim = {bid: np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0])
                      for bid in free_ids}

    # cells originally near the driver (within its shell radius + 2 cell widths) -- for
    # the local spacing/thinning diagnostic, same idea as the numpy sweep's readout
    near_ids = [bid for bid in free_ids
                if np.linalg.norm(start_pos_sim[bid]) / SIM_SCALE < (2 * r_in + 6 * r_in)]

    def nn_spacing(ids):
        if len(ids) < 2:
            return float("nan")
        pts = np.array([p.getBasePositionAndOrientation(bid, physicsClientId=client)[0] for bid in ids])
        d = np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=2)
        np.fill_diagonal(d, np.inf)
        return float(d.min(axis=1).mean()) / SIM_SCALE  # back to fm

    nn_spacing_before = nn_spacing(near_ids)

    max_pen_fm_overall = 0.0
    for step in range(N_STEPS):
        for bid in free_ids:
            pos = np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0])
            r = np.linalg.norm(pos)
            if r > 1e-9:
                force_sim = -P_CONF * (pos / r)
                p.applyExternalForce(bid, -1, list(force_sim), list(pos), p.WORLD_FRAME,
                                      physicsClientId=client)

        lin_vel, _ = p.getBaseVelocity(driver_id, physicsClientId=client)
        p.resetBaseVelocity(driver_id, linearVelocity=lin_vel, angularVelocity=list(omega_vec),
                             physicsClientId=client)

        p.stepSimulation(physicsClientId=client)

        contacts = p.getContactPoints(physicsClientId=client)
        if contacts:
            step_min = min(c[8] for c in contacts)  # index 8 = contactDistance
            if step_min < max_pen_fm_overall:
                max_pen_fm_overall = step_min

    max_pen_fm_overall = -max_pen_fm_overall / SIM_SCALE  # positive = penetration depth, in fm

    end_pos_sim = {bid: np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0])
                    for bid in free_ids}
    disp_fm = np.array([np.linalg.norm(end_pos_sim[bid] - start_pos_sim[bid]) / SIM_SCALE
                         for bid in free_ids])
    spread_mean_fm = float(disp_fm.mean())

    nn_spacing_after = nn_spacing(near_ids)
    thinning_pct = 100.0 * (nn_spacing_after - nn_spacing_before) / nn_spacing_before

    print(f"N_FREE={N_FREE}, N_STEPS={N_STEPS}, driver=13-cell rigid compound (1+12 kissing shell)")
    print(f"max penetration over whole run: {max_pen_fm_overall:.6f} fm "
          f"({100*max_pen_fm_overall/r_in:.2f}% of r_in)")
    print(f"spread_mean (free-cell displacement): {spread_mean_fm:.6f} fm "
          f"({spread_mean_fm/r_in:.3f} r_in)")
    print(f"near-driver NN spacing: {nn_spacing_before:.6f} -> {nn_spacing_after:.6f} fm "
          f"({thinning_pct:+.2f}%)")
    print(f"cells tracked near driver: {len(near_ids)}")

    p.disconnect(client)


if __name__ == "__main__":
    main()
