"""
mesh_grind_pybullet_two_driver_chirality_test.py
================================================
THIRD script in the PyBullet migration -- the actual chirality target ("model
Zone2 + observe how spin affects cells to see how chirality can repel same or
attract opposite", author, 2026-09-20). Builds DIRECTLY on the validated
single-driver engine (mesh_grind_pybullet_icosahedron_validation.py's convex
hull, mesh_grind_pybullet_zone2_driver_test.py's rigid welded 13-cell compound
+ damping(0.4/0.4) + grounded P_CONF=2.0 + CCD + resampling -- 15.64% r_in max
penetration, stable, the best-known PyBullet config).

ARCHITECTURE CORRECTION (author, 2026-09-20): this is NOT "two clusters" --
it is ONE shared pool of free cells with TWO embedded Zone-2 drivers, matching
the numpy engine's own "UNIFIED CONFINEMENT" fix (mesh-grind-chirality-
derivation.md, 2026-09-13): every free cell is confined toward the SHARED
midpoint of both driver centers, not toward whichever driver it started near
-- removing the per-driver "homing bias" that the numpy engine found could by
itself manufacture a spurious repel-looking result.

DRIVING RATE: OMEGA_MAG derived in proton_spin_omega_scale.py (same folder),
NOT a fresh guess. The proton's real measured total spin (1/2 hbar) implies a
bulk Zone-2 rotation rate (2.82e17 rad/s) 17 orders of magnitude beyond what
any discrete-time contact solver can run directly (same failure mode already
hit once before in the numpy engine). So: drive at omega_sim_safe (the
largest rate the "1 degree per physics step" numerical-safety cap allows,
already-accepted convention for this derivation line) -- the real number sets
the ceiling/justification, not a literal injected value. Axis: 3-axis
(Cartesian x/y/z, equal-weighted), the author's last explicit priority.

APPROACH FORCE: a constant attractive force between the two driver centers
(magnitude = P_CONF, same order as the per-cell confinement) stands in for
"the natural pressure-well attraction that draws nucleons together" (numpy
engine's own established mechanism, chirality-BLIND by construction) --
matches the SIMPLEST/first version of that line's own approach-force choice,
not yet the later inverse-square refinement (a natural next step if this
first pass is informative).

CONDITIONS: CONTROL (both drivers omega=0), SAME (both +omega0*axis), HETERO
(+omega0*axis / -omega0*axis) -- 3 conditions, same free-cell seed reused
across all three for a fair comparison (ORIENT_SEED convention, numpy engine).

DIAGNOSTICS: separation(t) between driver centers (the primary repel/attract
observable), max penetration (engine-health check, same as the single-driver
script), mixed_fraction (fraction of free cells that end up closer to the
OTHER driver than the one they started nearest -- numpy engine's own PC-7b
diagnostic, directly meaningful now that confinement is unified/shared).

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_two_driver_chirality_test.py
"""

import math
import os
import time
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

SIM_SCALE = 1.0 / L_J

N_STEPS = 800
CELL_MASS = 1.0
P_CONF = 2.0                 # grounded per jobson_cell_natural_force_pressure_scale.py
APPROACH_MAG = P_CONF        # same order as confinement -- constant approach force placeholder
# omega_sim_safe from proton_spin_omega_scale.py STEP 6 (1 deg/step cap, fixedTimeStep=1/60):
OMEGA_MAG = math.radians(1.0) / (1.0 / 60.0)   # = 1.047198 rad/(sim step-ish)
# real seconds per sim step, from proton_spin_omega_scale.py STEP 7 (radians(1deg)/omega_bulk) --
# keep this metric visible in every run, per author request 2026-09-20.
REAL_S_PER_STEP = 6.179739e-20
MAX_PLATEAU_WALL_SECONDS = 30.0   # author, 2026-09-20: stop once trend plateaus this long (wall-clock)
LINEAR_DAMPING = 0.4
ANGULAR_DAMPING = 0.4
# LATERAL FRICTION (2026-09-20, author: "we expect a perfect energy transfer... should
# probably be close to 100"): NEVER explicitly set before now -- PyBullet silently used
# its own default (verified below, printed at runtime, not assumed). PyBullet models
# COULOMB friction: tangential contact force is capped at mu*normal_force: mu=100 does NOT
# mean "100% energy transfer" literally -- it means the contact can sustain up to a huge
# tangential force before slipping, i.e. a firm/near-rigid no-slip grip (what the framework
# wants: RP1-RP4's zero-attenuation handoff, AUTHOR CORRECTION #4's lossless-coupling fix).
# Whether contact is actually lossless depends on whether it achieves zero relative sliding
# velocity (pure rolling/rotation), not on mu alone -- a high mu makes that MORE likely,
# it doesn't guarantee it. Adopting mu=100 per the author's instinct; comparing against the
# unset-default run to see if it changes anything.
LATERAL_FRICTION = 100.0
# CONFINEMENT RADIUS GATE (2026-09-20, author: "a specific pressure on each cell rather
# than only cells that go outside a specified radius?"): confirmed the prior version DID
# apply P_CONF to every cell everywhere -- fixed to a real boundary/wall model: zero force
# on any cell inside the pool's own seeding envelope, only cells that stray OUTSIDE it feel
# a constant inward pull -- matches "one common EXTERNAL bounding pressure" (a container
# wall) rather than a uniform interior body-force (which artificially compresses cells
# already near the drivers, not just strays at the edge).
#
# GEOMETRY REDERIVED (2026-09-20, author: "set them to a meaningful distance -- also make
# sure the wall radius is appropriate to the total number of cells"): the previous
# SEPARATION_FM=2.2x footprint was barely past the drivers-touch distance (2x footprint),
# so HETERO's "touch" test needed almost no travel. Rederived from TOUCH_SEP_FM (drivers
# just touching) so SAME (grow to wall) and HETERO (shrink to touch) have EQUAL travel
# distance (=TOUCH_SEP_FM each): start at 2x touch, wall at 3x touch.
DRIVER_FOOTPRINT_FM = 2 * r_in + R_c + r_in
TOUCH_SEP_FM = 2 * DRIVER_FOOTPRINT_FM
SEPARATION_FM = 2.0 * TOUCH_SEP_FM
# CONFINE_RADIUS widened (2026-09-20, author: "containment range... sufficiently far enough
# back not to be the deciding factor as the distance increases"): was 1.5*TOUCH_SEP_FM (wall
# diameter=3x touch) -- the sweep's own largest tested start (mult=2.9) sat almost ON the
# wall (driver-to-center distance 1.45x touch vs wall radius 1.5x touch, 96.7% of the way
# there). Tried 6x touch first -- but a first full decoupling (seed pool grown to match, no
# outer margin shell) made EVERY cell start inside the wall, so the confinement force never
# triggered for ANY cell -- confirmed by direct test: coord collapsed to 0.04 (no consolidating
# force at all, not just dilution). The wall-triggered squeeze on an outer margin shell is the
# medium's ONLY consolidating mechanism (no gravity analog) -- removing it is worse than being
# close to it. Settled on 2x touch (driver@mult=2.9 is 1.45/2.0=72.5% of the way to the wall --
# a real, meaningful improvement over 96.7% -- without discarding the margin-shell mechanism.
CONFINE_RADIUS_FM = 2.0 * TOUCH_SEP_FM
SEED_MARGIN_FM = 3 * r_in   # outer shell beyond the wall -- generates the "universal pressure"
                            # (shrunk from 10*r_in: that made N_FREE=3795, too slow for the
                            # unvectorized per-cell Python loop -- iterate-faster > perfect density)
SEED_RADIUS_FM = CONFINE_RADIUS_FM + SEED_MARGIN_FM

# N_FREE sized so the INTERIOR (r<CONFINE_RADIUS_FM) starts reasonably filled -- NOT seeded at
# full 0.64 random-close-packing (a naive rejection sampler stalls badly approaching that
# density -- confirmed: N_FREE=1261 at 0.64 hung for minutes just placing cells). Use a lower
# SEEDING target (0.20) and let the confinement force + contacts consolidate it further during
# the run, same as the earlier N=400 runs (seeded sparse, ended at coord~6-7 by step ~200-400).
_PHI_RCP = 0.20
_percell_vol_fm3 = (4.0 / 3.0) * pi * R_c**3
_interior_vol_fm3 = (4.0 / 3.0) * pi * CONFINE_RADIUS_FM**3
_n_interior_needed = _interior_vol_fm3 * _PHI_RCP / _percell_vol_fm3
_interior_volume_fraction = (CONFINE_RADIUS_FM / SEED_RADIUS_FM) ** 3
N_FREE = int(math.ceil(_n_interior_needed / _interior_volume_fraction))
ORIENT_SEED = 456                            # matches numpy engine's own fair-comparison convention


def n_free_for_confine(confine_radius_fm):
    """Same N_FREE derivation as above, parameterized by confine_radius_fm -- needed because
    the wall now scales with start_sep per-run (2026-09-21 wall-scaling fix: driver-to-wall
    margin was growing from 27.5% to 72.5% across the sweep with a FIXED wall, a real confound
    on the distance trend -- see two_driver_wall_geometry_tuning.py's follow-up section)."""
    seed_radius_fm = confine_radius_fm + SEED_MARGIN_FM
    interior_vol_fm3 = (4.0 / 3.0) * pi * confine_radius_fm**3
    interior_volume_fraction = (confine_radius_fm / seed_radius_fm) ** 3
    n_needed = interior_vol_fm3 * _PHI_RCP / _percell_vol_fm3
    return int(math.ceil(n_needed / interior_volume_fraction))


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
shell_directions = [norm(v) for v in verts_fm]
assert len(shell_directions) == 12

AXIS_3 = norm(np.array([1.0, 1.0, 1.0]))   # 3-axis (Cartesian) convention, author's priority


def make_ico_shape(client):
    return p.createCollisionShape(p.GEOM_MESH, vertices=verts_sim, physicsClientId=client)


def build_zone2_driver(client, ico_shape, position_sim):
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
        baseMass=CELL_MASS, baseCollisionShapeIndex=ico_shape, baseVisualShapeIndex=-1,
        basePosition=list(position_sim), baseOrientation=[0, 0, 0, 1],
        linkMasses=link_masses, linkCollisionShapeIndices=link_collision_shapes,
        linkVisualShapeIndices=link_visual_shapes, linkPositions=link_positions,
        linkOrientations=link_orientations, linkInertialFramePositions=link_inertial_positions,
        linkInertialFrameOrientations=link_inertial_orientations,
        linkParentIndices=link_parent_indices, linkJointTypes=link_joint_types,
        linkJointAxis=link_joint_axis, physicsClientId=client,
    )
    return body


def sample_uniform_ball(n, r_max, rng, center=np.zeros(3)):
    u = rng.random(n)
    r = r_max * u ** (1.0 / 3.0)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    return center + r[:, None] * d


def build_grid_seed(seed_radius_fm, posA_fm, posB_fm, driver_footprint_fm, rng, max_n,
                     jitter_frac=0.15):
    """Deterministic simple-cubic grid seeding, spacing=1.05*2*R_c (guarantees no overlap
    regardless of relative orientation) -- REPLACES random rejection-sampling, which stalls
    badly approaching a high packing-fraction target (confirmed: N=1261 at phi_rcp=0.64 hung
    for minutes; this is O(N) grid generation + vectorized filtering, no rejection loop at
    all). Small random jitter breaks exact-lattice symmetry (avoids degenerate-SAT artifacts,
    same lesson as the numpy engine's own random-orientation-seed precedent). SUBSAMPLED down
    to max_n if the raw grid exceeds it -- the per-step confinement loop is an unvectorized
    Python loop (PyBullet has no batch position/force query for independent bodies), so keeping
    n_free bounded is what keeps a multi-run SWEEP tractable (confirmed: an unsubsampled dense
    grid, ~1000+ cells, made even ONE 800-step run impractically slow)."""
    spacing_fm = 1.05 * 2 * R_c
    n_side = int(math.ceil(seed_radius_fm / spacing_fm))
    axis_1d = np.arange(-n_side, n_side + 1) * spacing_fm
    gx, gy, gz = np.meshgrid(axis_1d, axis_1d, axis_1d, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    pts = pts + rng.uniform(-jitter_frac * spacing_fm, jitter_frac * spacing_fm, size=pts.shape)
    r = np.linalg.norm(pts, axis=1)
    pts = pts[r <= seed_radius_fm]
    dA = np.linalg.norm(pts - posA_fm, axis=1)
    dB = np.linalg.norm(pts - posB_fm, axis=1)
    pts = pts[(dA >= driver_footprint_fm) & (dB >= driver_footprint_fm)]
    if len(pts) > max_n:
        idx = rng.choice(len(pts), size=max_n, replace=False)
        pts = pts[idx]
    return pts


def run_condition(label, omega_A, omega_B, stop_mode="fixed", start_sep_fm=None, n_free_target=N_FREE,
                   approach_mag=None, confine_radius_fm=None):
    approach_mag = APPROACH_MAG if approach_mag is None else approach_mag
    confine_radius_fm = CONFINE_RADIUS_FM if confine_radius_fm is None else confine_radius_fm
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, 0, physicsClientId=client)
    p.setPhysicsEngineParameter(numSolverIterations=100, numSubSteps=4,
                                 fixedTimeStep=1.0 / 60.0, erp=0.9, contactERP=0.9,
                                 globalCFM=1e-7, contactSlop=0.0001, physicsClientId=client)

    sep_fm = start_sep_fm if start_sep_fm is not None else SEPARATION_FM
    ico_shape = make_ico_shape(client)
    posA_fm = np.array([-sep_fm / 2, 0.0, 0.0])
    posB_fm = np.array([+sep_fm / 2, 0.0, 0.0])
    driver_A = build_zone2_driver(client, ico_shape, list(posA_fm * SIM_SCALE))
    driver_B = build_zone2_driver(client, ico_shape, list(posB_fm * SIM_SCALE))

    rng = np.random.RandomState(ORIENT_SEED)
    seed_radius_fm = confine_radius_fm + SEED_MARGIN_FM
    free_start_fm = build_grid_seed(seed_radius_fm, posA_fm, posB_fm, DRIVER_FOOTPRINT_FM, rng,
                                     max_n=n_free_target)
    n_free = len(free_start_fm)

    ccd_radius_sim = r_in * SIM_SCALE * 0.5
    for bid in (driver_A, driver_B):
        p.changeDynamics(bid, -1, ccdSweptSphereRadius=ccd_radius_sim,
                          linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                          lateralFriction=LATERAL_FRICTION, physicsClientId=client)
        for link in range(12):
            p.changeDynamics(bid, link, ccdSweptSphereRadius=ccd_radius_sim,
                              linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                              lateralFriction=LATERAL_FRICTION, physicsClientId=client)

    free_ids = []
    for i in range(n_free):
        pos_sim = free_start_fm[i] * SIM_SCALE
        body = p.createMultiBody(baseMass=CELL_MASS, baseCollisionShapeIndex=ico_shape,
                                  baseVisualShapeIndex=-1, basePosition=list(pos_sim),
                                  baseOrientation=[0, 0, 0, 1], physicsClientId=client)
        if i == 0 and label == "CONTROL":
            default_friction = p.getDynamicsInfo(body, -1, physicsClientId=client)[1]
            print(f"[DIAG] PyBullet's default lateralFriction (before we override it) = "
                  f"{default_friction} -- now explicitly set to {LATERAL_FRICTION}")
        p.changeDynamics(body, -1, ccdSweptSphereRadius=ccd_radius_sim,
                          linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                          lateralFriction=LATERAL_FRICTION, physicsClientId=client)
        free_ids.append(body)

    # DIAG: was any free cell seeded overlapping ANOTHER free cell (no mutual-exclusion check
    # exists in the resampling above, only vs the two drivers) -- a real candidate cause of
    # early, seed-driven penetration that would look identical across conditions (same seed).
    pair_d = np.linalg.norm(free_start_fm[:, None, :] - free_start_fm[None, :, :], axis=2)
    np.fill_diagonal(pair_d, np.inf)
    min_seed_gap_fm = float(pair_d.min())

    home_is_A = {}
    for i, bid in enumerate(free_ids):
        home_is_A[bid] = np.linalg.norm(free_start_fm[i] - posA_fm) < np.linalg.norm(free_start_fm[i] - posB_fm)

    sep_start_fm = np.linalg.norm(posB_fm - posA_fm)
    pen_history_fm = []   # per-step worst (most negative) contactDistance, for trend diagnosis
    contact_samples = []  # (step, mismatch_mag, normal_force) at driver-driver contacts only --
                           # empirical check of the mismatch->force constitutive law (2026-09-21)
    omega_free_start = {bid: np.array(p.getBaseVelocity(bid, physicsClientId=client)[1]) for bid in free_ids}

    # ADAPTIVE STOP (author, 2026-09-20): stop_mode="repel" (SAME) runs until a driver reaches
    # the confinement wall OR separation plateaus (no new max) for >30 wall-clock seconds.
    # stop_mode="attract" (HETERO) runs until the drivers touch OR separation plateaus (no new
    # min) for >30 wall-clock seconds. stop_mode="fixed" (CONTROL) keeps the old N_STEPS loop.
    best_sep = sep_start_fm
    last_progress_wall = time.time()
    stop_reason = "max_steps"
    step = 0
    sep_history_fm = [sep_start_fm]   # full time series, not just start/end -- author asked
    posA_history_fm = [posA_fm.copy()]   # FULL 3D trajectory (author, 2026-09-20: "how did the
    posB_history_fm = [posB_fm.copy()]   # detailed positional logging suggest the maintained
                                          # orientation? or did we only track distance?" -- we
                                          # only tracked scalar separation before; a spiral/orbit
                                          # would be invisible to that metric. Track full (x,y,z).
    # FIXED CONFINEMENT CENTER (2026-09-20, author correction): was recomputed every step as
    # the LIVE midpoint of the two drivers' CURRENT positions -- a self-referential feedback
    # (cells get pulled toward wherever the drivers currently average out to, then push the
    # drivers back toward that same point via contact, artificially coupling their motion
    # through the medium). Anchored ONCE here, from the INITIAL positions, in the lab frame --
    # doesn't respond to what the drivers do, same as any ordinary bounded-container setup.
    fixed_center = 0.5 * (posA_fm + posB_fm)
    while True:
        cA = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
        cB = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE

        for bid in free_ids:
            pos = np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0]) / SIM_SCALE
            rel = pos - fixed_center
            r = np.linalg.norm(rel)
            if r > confine_radius_fm:
                force_sim = -P_CONF * (rel / r)
                p.applyExternalForce(bid, -1, list(force_sim), list(pos * SIM_SCALE), p.WORLD_FRAME,
                                      physicsClientId=client)

        d_ab = cB - cA
        d_norm = np.linalg.norm(d_ab)
        if d_norm > 1e-9:
            approach_dir = d_ab / d_norm
            # INVERSE-SQUARE (2026-09-20, quick-win alignment fix before the distance sweep):
            # was a flat constant -- the established mechanism (numpy engine's own "AUTHOR
            # CORRECTION", mesh-grind-chirality-derivation.md) is Coulomb-like, F~1/r^2, matching
            # V=-alpha_em*hbar*c/r. A flat force would confound the very distance-dependence this
            # sweep measures (relatively too weak close in, too strong far out). Anchored so the
            # magnitude still equals APPROACH_MAG exactly AT the touch distance, for continuity.
            force_mag = approach_mag * (TOUCH_SEP_FM / d_norm) ** 2
            p.applyExternalForce(driver_A, -1, list(force_mag * approach_dir), list(cA * SIM_SCALE),
                                  p.WORLD_FRAME, physicsClientId=client)
            p.applyExternalForce(driver_B, -1, list(-force_mag * approach_dir), list(cB * SIM_SCALE),
                                  p.WORLD_FRAME, physicsClientId=client)

        for bid, omega_vec in ((driver_A, omega_A), (driver_B, omega_B)):
            lin_vel, _ = p.getBaseVelocity(bid, physicsClientId=client)
            p.resetBaseVelocity(bid, linearVelocity=lin_vel, angularVelocity=list(omega_vec),
                                 physicsClientId=client)

        p.stepSimulation(physicsClientId=client)
        step += 1

        cA_now = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
        cB_now = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE
        posA_history_fm.append(cA_now)
        posB_history_fm.append(cB_now)

        contacts = p.getContactPoints(physicsClientId=client)
        step_min = min((c[8] for c in contacts), default=0.0)
        pen_history_fm.append(-step_min / SIM_SCALE)   # positive = penetration depth, in fm

        # CONTACT-FORCE LOGGING (2026-09-21): driver-driver contacts only (not free-cell
        # contacts) -- pair the REAL PyBullet normal force (index 9) with the analytic
        # kinematic mismatch (same cross-product formula as mesh_grind_patch_integration.py)
        # at that same contact point, using the actual omega_A/omega_B and current centers.
        driver_pair = {driver_A, driver_B}
        for c in contacts:
            if {c[1], c[2]} == driver_pair:
                contact_pos_fm = np.array(c[5]) / SIM_SCALE   # positionOnA, world frame
                r_A = contact_pos_fm - cA_now
                r_B = contact_pos_fm - cB_now
                mismatch_vec = np.cross(omega_A, r_A) - np.cross(omega_B, r_B)
                contact_samples.append((step, float(np.linalg.norm(mismatch_vec)), float(c[9])))

        cA2 = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
        cB2 = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE
        sep_now = np.linalg.norm(cB2 - cA2)
        sep_history_fm.append(sep_now)

        if stop_mode == "fixed":
            if step >= N_STEPS:
                break
            continue

        if stop_mode == "repel":
            if sep_now > best_sep + 1e-9:
                best_sep = sep_now
                last_progress_wall = time.time()
            if sep_now / 2 >= confine_radius_fm:
                stop_reason = "wall_touch"
                break
        elif stop_mode == "attract":
            if sep_now < best_sep - 1e-9:
                best_sep = sep_now
                last_progress_wall = time.time()
            if sep_now <= 2 * DRIVER_FOOTPRINT_FM:
                stop_reason = "drivers_touch"
                break
        if time.time() - last_progress_wall > MAX_PLATEAU_WALL_SECONDS:
            stop_reason = "plateau"
            break

    n_steps_run = step
    pen_history_fm = np.array(pen_history_fm)
    max_pen_fm_overall = float(pen_history_fm.max())
    step_of_max = int(pen_history_fm.argmax())

    cA_end = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
    cB_end = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE
    sep_end_fm = np.linalg.norm(cB_end - cA_end)
    driver_A_disp_fm = float(np.linalg.norm(cA_end - posA_fm))
    driver_B_disp_fm = float(np.linalg.norm(cB_end - posB_fm))

    # ORBITAL/SPIRAL CHECK (author, 2026-09-20): scalar separation is blind to tangential motion
    # -- two drivers spiraling around each other could show the same net delta_sep as a straight
    # approach. Measure max perpendicular deviation from the INITIAL A-B axis for each driver.
    posA_history_fm = np.array(posA_history_fm)
    posB_history_fm = np.array(posB_history_fm)
    init_axis = norm(posB_fm - posA_fm)
    initA_mid = 0.5 * (posA_fm + posB_fm)
    relA = posA_history_fm - initA_mid
    relB = posB_history_fm - initA_mid
    perpA = relA - np.outer(relA @ init_axis, init_axis)
    perpB = relB - np.outer(relB @ init_axis, init_axis)
    max_perp_fm = float(max(np.linalg.norm(perpA, axis=1).max(), np.linalg.norm(perpB, axis=1).max()))

    n_crossed = 0
    free_disp_fm = np.zeros(n_free)
    omega_free_end_mag = np.zeros(n_free)
    free_end_fm = np.zeros((n_free, 3))
    for i, bid in enumerate(free_ids):
        pos_end = np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0]) / SIM_SCALE
        free_end_fm[i] = pos_end
        now_closer_to_A = np.linalg.norm(pos_end - cA_end) < np.linalg.norm(pos_end - cB_end)
        if now_closer_to_A != home_is_A[bid]:
            n_crossed += 1
        free_disp_fm[i] = np.linalg.norm(pos_end - free_start_fm[i])
        omega_free_end_mag[i] = np.linalg.norm(p.getBaseVelocity(bid, physicsClientId=client)[1])
    mixed_fraction = n_crossed / n_free

    # MAXWELL-CRITICALITY / ISOSTATIC-JAMMING CHECK (author, 2026-09-20): reuses the SAME
    # coordination-number convention already established in mesh_grind_real_rotation_spin_test.py
    # (dist < 1.05*(2*R_c) counts as a touching neighbor) -- target is the isostatic ~6
    # (random close packing of 3D spheres, Z_iso=2*d=6) that mesh_grind_cluster_relaxation_v2.py
    # already found (coord mean=6.63) for a single relaxed jam. If this pool is NOT near that
    # coordination number, it isn't genuinely "Maxwell-critical jammed" (translationally locked,
    # rotationally free) the way Zone 2 is supposed to be -- it's just a loose gas of cells
    # under a weak confinement force, a real candidate explanation for a muddled chirality signal.
    diff_end = free_end_fm[:, None, :] - free_end_fm[None, :, :]
    dist_end = np.linalg.norm(diff_end, axis=2)
    np.fill_diagonal(dist_end, np.inf)
    coord_end = np.sum(dist_end < 1.05 * (2 * R_c), axis=1)

    p.disconnect(client)

    sep_history_fm = np.array(sep_history_fm)
    sep_diffs = np.diff(sep_history_fm)
    non_monotonic_frac = float(np.mean(np.sign(sep_diffs) != np.sign(sep_end_fm - sep_start_fm))) if len(sep_diffs) else 0.0

    # CONSTITUTIVE-LAW CHECK (2026-09-21): through-origin least-squares fit normal_force =
    # k_fit * mismatch, over every driver-driver contact sample this run -- a first empirical
    # read on the still-missing mismatch->force conversion (mesh_grind_stiffness_close.py's
    # open question), directly from real PyBullet contact forces rather than an assumed law.
    contact_samples_arr = np.array(contact_samples) if contact_samples else np.zeros((0, 3))
    if len(contact_samples_arr):
        mismatch_vals = contact_samples_arr[:, 1]
        force_vals = contact_samples_arr[:, 2]
        denom = float(np.sum(mismatch_vals**2))
        k_fit = float(np.sum(mismatch_vals * force_vals) / denom) if denom > 1e-300 else float("nan")
    else:
        k_fit = float("nan")

    # TRAJECTORY PERSISTENCE (2026-09-21, author: raw per-frame data was in-memory only
    # before): full position/separation/penetration histories + contact samples, not
    # committed to the repo (matches /tmp/sweep_full_v2.log precedent for run artifacts).
    out_dir = "/tmp/mesh_grind_hetero_v3_data"
    os.makedirs(out_dir, exist_ok=True)
    mult_tag = sep_start_fm / TOUCH_SEP_FM
    np.savez_compressed(f"{out_dir}/{label}_mult{mult_tag:.1f}.npz",
                         posA_history_fm=posA_history_fm, posB_history_fm=posB_history_fm,
                         sep_history_fm=sep_history_fm, pen_history_fm=pen_history_fm,
                         contact_samples=contact_samples_arr)

    print(f"[{label}] n_free={n_free} stop={stop_reason} steps={n_steps_run} "
          f"real_t={n_steps_run*REAL_S_PER_STEP:.3e}s  "
          f"delta_sep={sep_end_fm - sep_start_fm:+.6f}fm  sep_range=[{sep_history_fm.min():.4f},"
          f"{sep_history_fm.max():.4f}]  non_monotonic_frac={non_monotonic_frac:.2f}  "
          f"max_perp={max_perp_fm:.4f}fm  "
          f"max_pen={100*max_pen_fm_overall/r_in:.2f}%r_in "
          f"(step{step_of_max})  mixed={mixed_fraction:.3f}  coord={coord_end.mean():.2f}  "
          f"disp_mean={free_disp_fm.mean():.4f}  omega_end={omega_free_end_mag.mean():.3f}  "
          f"driverAB_disp=({driver_A_disp_fm:.4f},{driver_B_disp_fm:.4f})  "
          f"n_contact_samples={len(contact_samples_arr)}  k_fit={k_fit:.4e}")

    return {"label": label, "sep_start": sep_start_fm, "sep_end": sep_end_fm,
            "delta_sep": sep_end_fm - sep_start_fm, "max_pen_pct": 100 * max_pen_fm_overall / r_in,
            "mixed_fraction": mixed_fraction, "step_of_max": step_of_max,
            "sep_history": sep_history_fm, "max_perp_fm": max_perp_fm,
            "coord_mean": float(coord_end.mean()), "n_steps_run": n_steps_run,
            "stop_reason": stop_reason, "n_free": n_free, "k_fit": k_fit,
            "n_contact_samples": len(contact_samples_arr)}


def main():
    # ORIENTATIONS expanded (author, 2026-09-20: "only 3 orientations doesn't quite seem
    # ample") -- added perpZ (the OTHER perpendicular direction, should match perp/Y by
    # symmetry -- a sanity check) and perp45 (an intermediate angle, not axis-aligned at all).
    ORIENTATIONS = {
        "coaxial": norm(np.array([1.0, 0.0, 0.0])),   # parallel to the driver separation axis
        "perp": norm(np.array([0.0, 1.0, 0.0])),      # perpendicular to it (Y)
        "perpZ": norm(np.array([0.0, 0.0, 1.0])),     # perpendicular to it (Z) -- symmetry check vs perp
        "perp45": norm(np.array([1.0, 1.0, 0.0])),    # intermediate, 45deg off-axis, not symmetry-special
        "diag": AXIS_3,                                # (1,1,1)/sqrt(3), prior default
    }
    print("=" * 78)
    print("CONFIG (2026-09-21) -- FULL re-run (CONTROL+SAME+HETERO), 3 changes vs the 209-run")
    print("2026-09-20 sweep (sweep_full_v2.log):")
    print(f"  (a) WALL SCALING: CONFINE_RADIUS_FM now = start_sep_fm exactly (was FIXED at"
          f" {CONFINE_RADIUS_FM:.4f}fm=2.0x touch for every mult) -- driver-to-wall margin is"
          f" now a CONSTANT 50% at every tested distance (was 27.5%->72.5%, a real confound on"
          f" the distance trend) -- see two_driver_wall_geometry_tuning.py's follow-up section.")
    print(f"  (b) TRAJECTORY PERSISTENCE: full posA/posB/sep/penetration histories now saved to"
          f" /tmp/mesh_grind_hetero_v3_data/<label>_mult<mult>.npz (previously in-memory only).")
    print(f"  (c) CONTACT-FORCE LOGGING: at every driver-driver contact, the REAL PyBullet normal"
          f" force is logged alongside the analytic kinematic mismatch (same formula as"
          f" mesh_grind_patch_integration.py) at that point -- empirical check of the still-"
          f" missing mismatch->force law, k_fit printed per run (through-origin F=k*mismatch fit;"
          f" CONTROL has omega=0 so its mismatch is always 0 and k_fit is expectedly nan).")
    print(f"  DRAW (approach force): inverse-square, magnitude={APPROACH_MAG} anchored AT touch"
          f" distance={TOUCH_SEP_FM:.4f}fm, chirality-BLIND. CHIRALITY: OMEGA_MAG={OMEGA_MAG:.4f}"
          f" (sim-safety ceiling, proton_spin_omega_scale.py). Orientations tested (axis vectors):")
    for name, axis in ORIENTATIONS.items():
        print(f"     {name:7s} = {np.array2string(axis, precision=4)}  "
              f"(SAME uses +{OMEGA_MAG:.3f}*axis both drivers, HETERO uses "
              f"+{OMEGA_MAG:.3f}*axis / -{OMEGA_MAG:.3f}*axis)")
    print(f"  friction={LATERAL_FRICTION}")
    print("=" * 78)

    zero = np.zeros(3)

    # SWEEP: CONTROL+SAME+HETERO restored (author, 2026-09-21: "the overnight full run should
    # probably turn the other checks back on which is fine for a true long overnight run") --
    # wall now scales with start_sep for ALL THREE conditions (confine_mult=mult, constant 50%
    # margin -- see two_driver_wall_geometry_tuning.py follow-up section for the full
    # N_FREE-vs-mult table, range 190 (mult=1.1) to 2133 (mult=2.9), both tractable per that
    # script's own precedent). CONTROL still runs ONCE per mult, shared across orientations.
    SWEEP_MULTIPLIERS = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                         2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9]
    sweep_results = []
    for mult in SWEEP_MULTIPLIERS:
        start_sep = mult * TOUCH_SEP_FM
        confine_radius_fm = start_sep   # (a) wall scaling: constant 50% margin at every mult
        n_free_target = n_free_for_confine(confine_radius_fm)
        control = run_condition("CONTROL", zero, zero, stop_mode="fixed", start_sep_fm=start_sep,
                                 confine_radius_fm=confine_radius_fm, n_free_target=n_free_target)
        for orient_name, axis in ORIENTATIONS.items():
            omega_vec = OMEGA_MAG * axis
            same = run_condition(f"SAME-{orient_name}", omega_vec, omega_vec, stop_mode="repel",
                                  start_sep_fm=start_sep, confine_radius_fm=confine_radius_fm,
                                  n_free_target=n_free_target)
            hetero = run_condition(f"HETERO-{orient_name}", omega_vec, -omega_vec, stop_mode="attract",
                                    start_sep_fm=start_sep, confine_radius_fm=confine_radius_fm,
                                    n_free_target=n_free_target)
            same_vs_ctrl = same["delta_sep"] - control["delta_sep"]
            hetero_vs_ctrl = hetero["delta_sep"] - control["delta_sep"]
            print(f"[mult={mult} orient={orient_name}] same>hetero={same['delta_sep'] > hetero['delta_sep']}  "
                  f"same_vs_ctrl={same_vs_ctrl:+.4f}  hetero_vs_ctrl={hetero_vs_ctrl:+.4f}  "
                  f"n_free={n_free_target}  k_fit_same={same['k_fit']:.4e}  k_fit_hetero={hetero['k_fit']:.4e}")
            sweep_results.append({"mult": mult, "orient": orient_name, "same_vs_ctrl": same_vs_ctrl,
                                   "hetero_vs_ctrl": hetero_vs_ctrl,
                                   "same_gt_hetero": same["delta_sep"] > hetero["delta_sep"],
                                   "k_fit_same": same["k_fit"], "k_fit_hetero": hetero["k_fit"]})

    print("---")
    for r in sweep_results:
        print(f"start={r['mult']}x orient={r['orient']:7s}: same>hetero={r['same_gt_hetero']}  "
              f"same_vs_ctrl={r['same_vs_ctrl']:+.4f}  hetero_vs_ctrl={r['hetero_vs_ctrl']:+.4f}  "
              f"k_fit_same={r['k_fit_same']:.4e}  k_fit_hetero={r['k_fit_hetero']:.4e}")


def quick_test_no_approach_force():
    """One-off (author, 2026-09-20): "not convinced that force shouldn't actually exist" --
    HETERO at mult=2.5*TOUCH_SEP_FM (the strongest recorded draw effect so far,
    hetero_vs_ctrl=-0.0333 WITH the approach force) but with approach_mag=0.0 -- does
    chirality/rotation ALONE produce any net attraction, or was the earlier effect riding
    entirely on the pre-imposed, chirality-blind approach force?
    RESULT (run 1): HETERO-noapproach moved together on its own, delta_sep=-0.010966fm,
    stop=plateau (never reached touch), steps=2190. REAL but needs the obvious control:
    does a non-spinning (omega=0) pair ALSO drift this much from generic settling/jostling
    with no approach force, at the same distance? Added CONTROL-noapproach (same stop_mode
    logic, same start_sep, omega=0) for a fair, apples-to-apples comparison.
    """
    omega_vec = OMEGA_MAG * AXIS_3
    start_sep = 2.5 * TOUCH_SEP_FM
    zero = np.zeros(3)
    control = run_condition("CONTROL-noapproach", zero, zero, stop_mode="attract",
                             start_sep_fm=start_sep, approach_mag=0.0)
    hetero = run_condition("HETERO-noapproach", omega_vec, -omega_vec, stop_mode="attract",
                            start_sep_fm=start_sep, approach_mag=0.0)
    print(f"\nNO APPROACH FORCE (start={start_sep:.4f}fm):")
    print(f"  CONTROL delta_sep={control['delta_sep']:+.6f}fm (stop={control['stop_reason']})")
    print(f"  HETERO  delta_sep={hetero['delta_sep']:+.6f}fm (stop={hetero['stop_reason']})")
    print(f"  hetero_vs_ctrl={hetero['delta_sep'] - control['delta_sep']:+.6f}fm -- "
          f"{'chirality itself pulls beyond generic settling' if hetero['delta_sep'] < control['delta_sep'] - 1e-4 else 'no clear chirality-specific effect vs generic settling'}")


def quick_sanity_check_new_geometry():
    """One-off (author, 2026-09-20): verify the widened CONFINE_RADIUS_FM (1.5x->6x touch)
    and new max_perp_fm diagnostic both work before committing to the full (now 5-orientation)
    sweep. Single fast CONTROL run only."""
    r = run_condition("CONTROL", np.zeros(3), np.zeros(3), stop_mode="fixed",
                       start_sep_fm=2.5 * TOUCH_SEP_FM)
    print(f"\nSANITY: coord={r['coord_mean']:.2f} (target ~6, was ~4 before widening -- watch for "
          f"further dilution) max_perp_fm={r['max_perp_fm']:.4f} n_free={r['n_free']}")


def quick_sanity_check_wall_scaling():
    """One-off (author, 2026-09-21): before committing to the full 95-run HETERO-only overnight
    sweep, verify the 2026-09-21 wall-scaling fix (confine_mult=mult) + new contact-force/
    mismatch logging + trajectory persistence all work, at the WORST-CASE mult=2.9 (largest
    N_FREE=2133 per two_driver_wall_geometry_tuning.py's follow-up table -- if any mult is
    going to crash, run slowly, or produce empty/garbage contact data, this is the one)."""
    mult = 2.9
    start_sep = mult * TOUCH_SEP_FM
    confine_radius_fm = start_sep
    n_free_target = n_free_for_confine(confine_radius_fm)
    omega_vec = OMEGA_MAG * AXIS_3
    t0 = time.time()
    r = run_condition("HETERO-diag", omega_vec, -omega_vec, stop_mode="attract",
                       start_sep_fm=start_sep, confine_radius_fm=confine_radius_fm,
                       n_free_target=n_free_target)
    wall_s = time.time() - t0
    print(f"\nSANITY (mult={mult}, worst-case N_FREE={n_free_target}): wall_clock={wall_s:.1f}s  "
          f"steps={r['n_steps_run']}  stop={r['stop_reason']}  delta_sep={r['delta_sep']:+.4f}fm  "
          f"n_contact_samples={r['n_contact_samples']}  k_fit={r['k_fit']:.4e}  "
          f"est_full_sweep_min={95*wall_s/60:.1f} (upper bound, assumes every run this slow)")


if __name__ == "__main__":
    main()
