"""
mesh_grind_pybullet_two_driver_contact_gated_test.py
======================================================
FIFTH script in the PyBullet two-driver chirality line. Fixes two real
issues found in the extended-range "v4" script (mesh_grind_pybullet_two_
driver_extended_range_test.py, sweep_extended_range_v1.log) before running
it again, per author instruction (2026-09-23) -- does NOT overwrite v4,
which already produced a real, already-logged result; kept as-is.

FIX 1 -- STOP CONDITION WAS DISTANCE-BASED, NOT CONTACT-BASED (why
n_contact_samples=0 in EVERY v4 row, even "drivers_touch" ones): v4's
"attract" stop condition was `sep_now <= 2*DRIVER_FOOTPRINT_FM` -- a
CENTER-TO-CENTER heuristic based on the driver's own nominal geometric
footprint, not on whether PyBullet's real collision solver actually
registered a driver-driver contact point. Since DRIVER_FOOTPRINT_FM is an
approximate bound (12 shell cells around a center, at a FIXED, never-
randomized relative orientation -- see FIX 2), the true rigid-body surfaces
are often still not touching when this heuristic fires, especially at
larger separations where more steps elapse before the threshold is crossed.
The v3 script (sweep_full_v3.log) DID capture real, substantial contact
counts (15 to 13,979 samples per condition, in ~21 of 209 rows) at its own
shorter-range scale -- confirming the underlying contact-detection code
works fine; this is specifically a v4 stop-condition scale/tuning issue,
not a repo-wide "contact detection never works" problem.
FIX: "attract" mode's PRIMARY stop trigger is now real accumulated driver-
driver contact samples (>= N_CONTACT_MIN, a small buffer against a single
noisy/transient contact), not center-to-center distance. The old distance
check is kept ONLY as a much TIGHTER safety-net fallback (half of v4's own
threshold) in case real contact never registers even at very close
approach -- itself a genuinely informative outcome (would mean free cells
stay wedged between the drivers, preventing direct contact) rather than a
silent stop-condition artifact.

FIX 2 -- DRIVER ORIENTATION WAS NEVER RANDOMIZED (both v3 and v4 built
every driver at baseOrientation=[0,0,0,1], identity, always): this is the
SAME class of design flaw already caught once before and paused over in
the numpy-engine Phase C investigation (mesh-grind-chirality-derivation.md,
2026-09-12: "the Phase C test held ONE C5 axis fixed... this static setup
may not represent the real mechanism"). With identity orientation always,
the specific shell-cell-to-shell-cell facing alignment at the moment of
closest approach is fully deterministic (fixed start + fixed spin axis +
elapsed time) -- v4's clean same>hetero=True in 7/7 (with approach) could,
as tested, be an artifact of this ONE relative alignment, not a general
result across the many possible relative orientations two real rotating
Zone-2 compounds could have.
FIX: each driver gets a genuinely random starting orientation (uniform
quaternion, via 4D-Gaussian-then-normalize, the standard Haar-uniform-on-
SO(3) sampling method), drawn from a seed that is SHARED across CONTROL/
SAME/HETERO at the same mult (fair same-configuration comparison) but
DIFFERENT across mult values (so the sweep as a whole samples multiple
relative orientations, not just one).

SCOPE NARROWING (author, 2026-09-23): "we already proved this with the
previous script [v3] -- this one was about getting further distances and
testing approach off... a new run could probably disable approach off."
v3 already established the base chirality-gates-approach finding at short
range; approach=without was already shown negligible in v4 (deltas ~1e-4
to 1e-3, never reaching drivers_touch, no clear same>hetero majority).
approach=without is DROPPED here to spend the compute budget on validating
the with-approach finding at extended range, with both fixes in place.

VALIDATION-FIRST WORKFLOW (author instruction: "run a short test to
validate" before any full run): --sanity runs ONE small-N, single-mult
condition set (CONTROL/SAME/HETERO at mult=2.0) to confirm (a) contact
samples now accumulate for at least the SAME condition (same-chirality
grinding contact is the easier case to trigger first) and (b) driver
orientation actually differs between mult values, before committing to
the full mult sweep.

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_two_driver_contact_gated_test.py --sanity
     python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_two_driver_contact_gated_test.py
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
P_CONF = 2.0
APPROACH_MAG = P_CONF
OMEGA_MAG = math.radians(1.0) / (1.0 / 60.0)
REAL_S_PER_STEP = 6.179739e-20
MAX_PLATEAU_WALL_SECONDS = 30.0
MAX_STEPS_ADAPTIVE = 60000
LINEAR_DAMPING = 0.4
ANGULAR_DAMPING = 0.4
LATERAL_FRICTION = 100.0
N_CONTACT_MIN = 3   # FIX 1: min accumulated driver-driver contact samples to declare touch

DRIVER_FOOTPRINT_FM = 2 * r_in + R_c + r_in
TOUCH_SEP_FM = 2 * DRIVER_FOOTPRINT_FM
B_FIXED_FM = 2.0 * TOUCH_SEP_FM
SEED_MARGIN_FM = 3 * r_in

_PHI_RCP = 0.20
_percell_vol_fm3 = (4.0 / 3.0) * pi * R_c**3
ORIENT_SEED = 456


def ellipsoid_semi_axes(mult):
    a = max(B_FIXED_FM, mult * TOUCH_SEP_FM)
    b = B_FIXED_FM * (a / B_FIXED_FM) ** 0.4
    return a, b


def n_free_for_ellipsoid(a_confine, b_confine):
    a_seed = a_confine + SEED_MARGIN_FM
    b_seed = b_confine + SEED_MARGIN_FM
    interior_vol_fm3 = (4.0 / 3.0) * pi * a_confine * b_confine * b_confine
    seed_vol_fm3 = (4.0 / 3.0) * pi * a_seed * b_seed * b_seed
    interior_volume_fraction = interior_vol_fm3 / seed_vol_fm3
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

AXIS_COAXIAL = norm(np.array([1.0, 0.0, 0.0]))


def make_ico_shape(client):
    return p.createCollisionShape(p.GEOM_MESH, vertices=verts_sim, physicsClientId=client)


def random_quaternion(rng):
    """FIX 2: uniform-on-SO(3) (Haar) random orientation -- normalize a 4D Gaussian vector."""
    q = rng.normal(size=4)
    q = q / np.linalg.norm(q)
    return list(q)


def build_zone2_driver(client, ico_shape, position_sim, orientation):
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
        basePosition=list(position_sim), baseOrientation=list(orientation),
        linkMasses=link_masses, linkCollisionShapeIndices=link_collision_shapes,
        linkVisualShapeIndices=link_visual_shapes, linkPositions=link_positions,
        linkOrientations=link_orientations, linkInertialFramePositions=link_inertial_positions,
        linkInertialFrameOrientations=link_inertial_orientations,
        linkParentIndices=link_parent_indices, linkJointTypes=link_joint_types,
        linkJointAxis=link_joint_axis, physicsClientId=client,
    )
    return body


def build_grid_seed_ellipsoid(a_seed_fm, b_seed_fm, posA_fm, posB_fm, driver_footprint_fm,
                               rng, max_n, jitter_frac=0.15):
    spacing_fm = 1.05 * 2 * R_c
    n_side_a = int(math.ceil(a_seed_fm / spacing_fm))
    n_side_b = int(math.ceil(b_seed_fm / spacing_fm))
    axis_x = np.arange(-n_side_a, n_side_a + 1) * spacing_fm
    axis_yz = np.arange(-n_side_b, n_side_b + 1) * spacing_fm
    gx, gy, gz = np.meshgrid(axis_x, axis_yz, axis_yz, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    pts = pts + rng.uniform(-jitter_frac * spacing_fm, jitter_frac * spacing_fm, size=pts.shape)
    normalized_r = np.sqrt((pts[:, 0] / a_seed_fm) ** 2 + (pts[:, 1] / b_seed_fm) ** 2
                            + (pts[:, 2] / b_seed_fm) ** 2)
    pts = pts[normalized_r <= 1.0]
    dA = np.linalg.norm(pts - posA_fm, axis=1)
    dB = np.linalg.norm(pts - posB_fm, axis=1)
    pts = pts[(dA >= driver_footprint_fm) & (dB >= driver_footprint_fm)]
    if len(pts) > max_n:
        idx = rng.choice(len(pts), size=max_n, replace=False)
        pts = pts[idx]
    return pts


def run_condition(label, omega_A, omega_B, stop_mode="fixed", start_sep_fm=None,
                   n_free_target=None, approach_mag=None, confine_a_fm=None, confine_b_fm=None,
                   condition_seed=ORIENT_SEED):
    approach_mag = APPROACH_MAG if approach_mag is None else approach_mag
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, 0, physicsClientId=client)
    p.setPhysicsEngineParameter(numSolverIterations=100, numSubSteps=4,
                                 fixedTimeStep=1.0 / 60.0, erp=0.9, contactERP=0.9,
                                 globalCFM=1e-7, contactSlop=0.0001, physicsClientId=client)

    sep_fm = start_sep_fm
    ico_shape = make_ico_shape(client)
    posA_fm = np.array([-sep_fm / 2, 0.0, 0.0])
    posB_fm = np.array([+sep_fm / 2, 0.0, 0.0])

    # FIX 2: shared rng per (mult, approach) -- same starting orientations across
    # CONTROL/SAME/HETERO at one configuration (fair comparison), different across mult.
    rng = np.random.RandomState(condition_seed)
    orient_A = random_quaternion(rng)
    orient_B = random_quaternion(rng)

    driver_A = build_zone2_driver(client, ico_shape, list(posA_fm * SIM_SCALE), orient_A)
    driver_B = build_zone2_driver(client, ico_shape, list(posB_fm * SIM_SCALE), orient_B)

    a_seed_fm = confine_a_fm + SEED_MARGIN_FM
    b_seed_fm = confine_b_fm + SEED_MARGIN_FM
    free_start_fm = build_grid_seed_ellipsoid(a_seed_fm, b_seed_fm, posA_fm, posB_fm,
                                               DRIVER_FOOTPRINT_FM, rng, max_n=n_free_target)
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
        p.changeDynamics(body, -1, ccdSweptSphereRadius=ccd_radius_sim,
                          linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                          lateralFriction=LATERAL_FRICTION, physicsClientId=client)
        free_ids.append(body)

    home_is_A = {}
    for i, bid in enumerate(free_ids):
        home_is_A[bid] = np.linalg.norm(free_start_fm[i] - posA_fm) < np.linalg.norm(free_start_fm[i] - posB_fm)

    sep_start_fm = np.linalg.norm(posB_fm - posA_fm)
    pen_history_fm = []
    contact_samples = []

    best_sep = sep_start_fm
    last_progress_wall = time.time()
    stop_reason = "max_steps"
    step = 0
    sep_history_fm = [sep_start_fm]
    posA_history_fm = [posA_fm.copy()]
    posB_history_fm = [posB_fm.copy()]
    fixed_center = 0.5 * (posA_fm + posB_fm)

    while True:
        cA = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
        cB = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE

        for bid in free_ids:
            pos = np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0]) / SIM_SCALE
            rel = pos - fixed_center
            normalized_r = math.sqrt((rel[0] / confine_a_fm) ** 2 + (rel[1] / confine_b_fm) ** 2
                                      + (rel[2] / confine_b_fm) ** 2)
            if normalized_r > 1.0:
                r_mag = np.linalg.norm(rel)
                if r_mag > 1e-12:
                    force_sim = -P_CONF * (rel / r_mag)
                    p.applyExternalForce(bid, -1, list(force_sim), list(pos * SIM_SCALE), p.WORLD_FRAME,
                                          physicsClientId=client)

        d_ab = cB - cA
        d_norm = np.linalg.norm(d_ab)
        if d_norm > 1e-9:
            approach_dir = d_ab / d_norm
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
        pen_history_fm.append(-step_min / SIM_SCALE)

        driver_pair = {driver_A, driver_B}
        for c in contacts:
            if {c[1], c[2]} == driver_pair:
                contact_pos_fm = np.array(c[5]) / SIM_SCALE
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
            if sep_now / 2 >= confine_a_fm:
                stop_reason = "wall_touch"
                break
        elif stop_mode == "attract":
            if sep_now < best_sep - 1e-9:
                best_sep = sep_now
                last_progress_wall = time.time()
            # FIX 1: gate primarily on real accumulated contact, not raw distance.
            if len(contact_samples) >= N_CONTACT_MIN:
                stop_reason = "drivers_touch"
                break
            # Tight distance fallback (half v4's threshold) -- only fires if the
            # drivers get very close WITHOUT ever registering real contact, which
            # would itself be a genuine finding (free cells wedged between them).
            if sep_now <= DRIVER_FOOTPRINT_FM:
                stop_reason = "close_no_contact"
                break
        if time.time() - last_progress_wall > MAX_PLATEAU_WALL_SECONDS:
            stop_reason = "plateau"
            break
        if step >= MAX_STEPS_ADAPTIVE:
            stop_reason = "max_steps_adaptive_cap"
            break

    n_steps_run = step
    pen_history_fm = np.array(pen_history_fm)
    max_pen_fm_overall = float(pen_history_fm.max())
    step_of_max = int(pen_history_fm.argmax())

    cA_end = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
    cB_end = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE
    sep_end_fm = np.linalg.norm(cB_end - cA_end)

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
    free_end_fm = np.zeros((n_free, 3))
    for i, bid in enumerate(free_ids):
        pos_end = np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0]) / SIM_SCALE
        free_end_fm[i] = pos_end
        now_closer_to_A = np.linalg.norm(pos_end - cA_end) < np.linalg.norm(pos_end - cB_end)
        if now_closer_to_A != home_is_A[bid]:
            n_crossed += 1
    mixed_fraction = n_crossed / n_free if n_free else 0.0

    diff_end = free_end_fm[:, None, :] - free_end_fm[None, :, :]
    dist_end = np.linalg.norm(diff_end, axis=2)
    np.fill_diagonal(dist_end, np.inf)
    coord_end = np.sum(dist_end < 1.05 * (2 * R_c), axis=1)

    p.disconnect(client)

    sep_history_fm = np.array(sep_history_fm)

    contact_samples_arr = np.array(contact_samples) if contact_samples else np.zeros((0, 3))
    if len(contact_samples_arr):
        mismatch_vals = contact_samples_arr[:, 1]
        force_vals = contact_samples_arr[:, 2]
        denom = float(np.sum(mismatch_vals**2))
        k_fit = float(np.sum(mismatch_vals * force_vals) / denom) if denom > 1e-300 else float("nan")
    else:
        k_fit = float("nan")

    out_dir = "/tmp/mesh_grind_contact_gated_data"
    os.makedirs(out_dir, exist_ok=True)
    mult_tag = sep_start_fm / TOUCH_SEP_FM
    np.savez_compressed(f"{out_dir}/{label}_mult{mult_tag:.1f}_seed{condition_seed}.npz",
                         posA_history_fm=posA_history_fm, posB_history_fm=posB_history_fm,
                         sep_history_fm=sep_history_fm, pen_history_fm=pen_history_fm,
                         contact_samples=contact_samples_arr, orient_A=orient_A, orient_B=orient_B)

    print(f"[{label}] n_free={n_free} stop={stop_reason} steps={n_steps_run} "
          f"delta_sep={sep_end_fm - sep_start_fm:+.6f}fm  sep_range=[{sep_history_fm.min():.4f},"
          f"{sep_history_fm.max():.4f}]  max_perp={max_perp_fm:.4f}fm  "
          f"max_pen={100*max_pen_fm_overall/r_in:.2f}%r_in (step{step_of_max})  "
          f"mixed={mixed_fraction:.3f}  coord={coord_end.mean() if n_free else 0:.2f}  "
          f"n_contact_samples={len(contact_samples_arr)}  k_fit={k_fit:.4e}")

    return {"label": label, "sep_start": sep_start_fm, "sep_end": sep_end_fm,
            "delta_sep": sep_end_fm - sep_start_fm, "max_pen_pct": 100 * max_pen_fm_overall / r_in,
            "mixed_fraction": mixed_fraction, "n_steps_run": n_steps_run,
            "stop_reason": stop_reason, "n_free": n_free, "k_fit": k_fit,
            "n_contact_samples": len(contact_samples_arr), "max_perp_fm": max_perp_fm}


def main():
    print("=" * 78)
    print("CONFIG (2026-09-23) -- contact-gated stop + randomized driver orientation,")
    print("approach=without DROPPED (already shown negligible in v4). Fixes applied:")
    print(f"  FIX 1: 'attract' stop gates on n_contact_samples>={N_CONTACT_MIN}, not raw distance.")
    print("  FIX 2: driver orientation is a random (Haar-uniform) quaternion, shared across")
    print("         CONTROL/SAME/HETERO at one mult, different across mult values.")
    print("=" * 78)

    zero = np.zeros(3)
    omega_vec = OMEGA_MAG * AXIS_COAXIAL

    SWEEP_MULTIPLIERS = [1.5, 2.0, 2.9, 4.0, 6.0, 8.0, 10.0]
    sweep_results = []
    for mi, mult in enumerate(SWEEP_MULTIPLIERS):
        start_sep = mult * TOUCH_SEP_FM
        a_confine, b_confine = ellipsoid_semi_axes(mult)
        n_free_target = n_free_for_ellipsoid(a_confine, b_confine)
        condition_seed = ORIENT_SEED + mi * 100
        control = run_condition("CONTROL", zero, zero, stop_mode="fixed", start_sep_fm=start_sep,
                                 n_free_target=n_free_target, approach_mag=APPROACH_MAG,
                                 confine_a_fm=a_confine, confine_b_fm=b_confine, condition_seed=condition_seed)
        same = run_condition("SAME-coaxial", omega_vec, omega_vec, stop_mode="repel",
                              start_sep_fm=start_sep, n_free_target=n_free_target,
                              approach_mag=APPROACH_MAG, confine_a_fm=a_confine, confine_b_fm=b_confine,
                              condition_seed=condition_seed)
        hetero = run_condition("HETERO-coaxial", omega_vec, -omega_vec, stop_mode="attract",
                                start_sep_fm=start_sep, n_free_target=n_free_target,
                                approach_mag=APPROACH_MAG, confine_a_fm=a_confine, confine_b_fm=b_confine,
                                condition_seed=condition_seed)
        same_vs_ctrl = same["delta_sep"] - control["delta_sep"]
        hetero_vs_ctrl = hetero["delta_sep"] - control["delta_sep"]
        print(f"[mult={mult}] same>hetero={same['delta_sep'] > hetero['delta_sep']}  "
              f"same_vs_ctrl={same_vs_ctrl:+.4f}  hetero_vs_ctrl={hetero_vs_ctrl:+.4f}  "
              f"n_free={n_free_target}  hetero_stop={hetero['stop_reason']}  "
              f"hetero_contacts={hetero['n_contact_samples']}  same_contacts={same['n_contact_samples']}")
        sweep_results.append({"mult": mult, "same_vs_ctrl": same_vs_ctrl,
                               "hetero_vs_ctrl": hetero_vs_ctrl,
                               "same_gt_hetero": same["delta_sep"] > hetero["delta_sep"],
                               "hetero_stop": hetero["stop_reason"], "n_free": n_free_target,
                               "hetero_contacts": hetero["n_contact_samples"],
                               "same_contacts": same["n_contact_samples"]})

    print("---")
    for r in sweep_results:
        print(f"start={r['mult']}x: same>hetero={r['same_gt_hetero']}  "
              f"same_vs_ctrl={r['same_vs_ctrl']:+.4f}  hetero_vs_ctrl={r['hetero_vs_ctrl']:+.4f}  "
              f"hetero_stop={r['hetero_stop']:22s}  n_free={r['n_free']}  "
              f"contacts(same/hetero)={r['same_contacts']}/{r['hetero_contacts']}")


def quick_sanity_check(mult=2.0):
    """Short validation run (author instruction, before any full sweep): confirms
    (a) contact samples now accumulate for at least the SAME condition, (b) driver
    orientations actually differ, (c) nothing crashes -- on ONE small mult value."""
    start_sep = mult * TOUCH_SEP_FM
    a_confine, b_confine = ellipsoid_semi_axes(mult)
    n_free_target = n_free_for_ellipsoid(a_confine, b_confine)
    omega_vec = OMEGA_MAG * AXIS_COAXIAL
    zero = np.zeros(3)
    condition_seed = ORIENT_SEED

    print(f"SANITY CHECK: mult={mult}, N_FREE={n_free_target}")
    for label, oA, oB, mode in [("CONTROL", zero, zero, "fixed"),
                                 ("SAME-coaxial", omega_vec, omega_vec, "repel"),
                                 ("HETERO-coaxial", omega_vec, -omega_vec, "attract")]:
        t0 = time.time()
        r = run_condition(label, oA, oB, stop_mode=mode, start_sep_fm=start_sep,
                           n_free_target=n_free_target, approach_mag=APPROACH_MAG,
                           confine_a_fm=a_confine, confine_b_fm=b_confine, condition_seed=condition_seed)
        wall_s = time.time() - t0
        print(f"  {label}: wall={wall_s:.1f}s steps={r['n_steps_run']} stop={r['stop_reason']} "
              f"n_contact_samples={r['n_contact_samples']} k_fit={r['k_fit']:.4e}")

    # Confirm orientation actually varies across mult (FIX 2 sanity)
    rng1 = np.random.RandomState(ORIENT_SEED)
    q1 = random_quaternion(rng1)
    rng2 = np.random.RandomState(ORIENT_SEED + 100)
    q2 = random_quaternion(rng2)
    print(f"  orientation check: mult-index-0 seed quaternion={np.round(q1,4)}, "
          f"mult-index-1 seed quaternion={np.round(q2,4)}  (should differ)")


if __name__ == "__main__":
    import sys
    if "--sanity" in sys.argv:
        args = [a for a in sys.argv[1:] if a != "--sanity"]
        mult_arg = float(args[0]) if args else 2.0
        quick_sanity_check(mult_arg)
    else:
        main()
