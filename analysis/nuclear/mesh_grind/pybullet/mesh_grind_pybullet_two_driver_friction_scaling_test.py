"""
mesh_grind_pybullet_two_driver_friction_scaling_test.py
==========================================================
EIGHTH script in the PyBullet two-driver chirality line. Follows directly
from mesh_grind_pybullet_two_driver_pressure_scaling_test.py's own finding
(2026-09-23): pure ambient confinement pressure alone (no explicit
approach force) showed the correct same>hetero direction at mult=2.0
under 10-50x baseline pressure, but never reached real contact, and did
not generalize to mult=3.0/7.0 -- a weaker, secondary signal than the
orientation-sweep sibling script's much cleaner 0/20 vs 19/20 result.

HYPOTHESIS TESTED HERE (2026-09-28): the gap is not pressure magnitude
ALONE but RIGIDITY REPRESENTATION. sandbox/corpuscle/
cell_rigidity_vs_macroscopic_stress_check.py already established that a
real Jobson cell is rigid to ~30 orders of magnitude beyond floating-
point precision under any medium-scale stress -- effectively perfectly
rigid, zero-slip. PyBullet's LATERAL_FRICTION=100.0 (used unchanged in
every prior script in this line) is a large but FINITE Coulomb-friction
coefficient, not a representation of true near-infinite rigidity.

FIRST ATTEMPT (1D, friction only, P_CONF held fixed at 100.0 -- the best-
tested level from the pressure-scaling sweep) was run and found
UNINFORMATIVE: a --sanity check showed IDENTICAL results at friction=100
vs friction=10000 (delta_sep matched to 4 decimal places), because BOTH
runs never reached real contact (n_contact_samples=0, stop=plateau) at
any friction level -- friction only acts at contact, so if contact never
happens, friction cannot matter regardless of its value. The 1D design
was testing the wrong thing: it couldn't distinguish "friction is fine"
from "friction is irrelevant because pressure alone never gets close
enough to touch."

CORRECTED DESIGN (2D, this version): scan FRICTION and P_CONF TOGETHER,
both pushed well past their previously-tested maxima, restricted to
mult=2.0 only (the one separation that showed ANY partial signal before;
farther separations are even less likely to help and are not worth the
compute until mult=2.0 itself is resolved). FRICTION_SCALES = [1x, 10x,
100x] the LATERAL_FRICTION baseline (100, 1000, 10000) CROSSED WITH
P_CONF_SCALES = [50x, 200x, 1000x] the P_CONF baseline (100, 400, 2000) --
9 (friction, pressure) combinations x 2 conditions (same/hetero) = 18
runs, approach_mag=0.0 throughout (pure ambient pressure, no explicit
force). Reports whether ANY combination reaches real contact
(n_contact_samples>=3, stop=drivers_touch) for HETERO, and whether SAME
ever does (it should not, per every prior result in this line).

Reuses ALL geometry/physics/contact-gating code from
mesh_grind_pybullet_two_driver_pressure_scaling_test.py verbatim, with
ONE deliberate change: run_condition() now takes an explicit friction
parameter (default = the module LATERAL_FRICTION constant, so the
function is unchanged for any caller that does not pass it), used in the
3 changeDynamics() calls instead of always reading the module-level
constant directly -- the same pattern p_conf already used for pressure.

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_two_driver_friction_scaling_test.py --sanity
     python analysis/nuclear/mesh_grind/pybullet/mesh_grind_pybullet_two_driver_friction_scaling_test.py
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
N_CONTACT_MIN = 3   # min accumulated driver-driver contact samples to declare touch

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
    """Uniform-on-SO(3) (Haar) random orientation -- normalize a 4D Gaussian vector."""
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
                   condition_seed=ORIENT_SEED, p_conf=None, friction=None):
    approach_mag = APPROACH_MAG if approach_mag is None else approach_mag
    p_conf = P_CONF if p_conf is None else p_conf
    friction = LATERAL_FRICTION if friction is None else friction
    client = p.connect(p.DIRECT)
    p.setGravity(0, 0, 0, physicsClientId=client)
    p.setPhysicsEngineParameter(numSolverIterations=100, numSubSteps=4,
                                 fixedTimeStep=1.0 / 60.0, erp=0.9, contactERP=0.9,
                                 globalCFM=1e-7, contactSlop=0.0001, physicsClientId=client)

    sep_fm = start_sep_fm
    ico_shape = make_ico_shape(client)
    posA_fm = np.array([-sep_fm / 2, 0.0, 0.0])
    posB_fm = np.array([+sep_fm / 2, 0.0, 0.0])

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
                          lateralFriction=friction, physicsClientId=client)
        for link in range(12):
            p.changeDynamics(bid, link, ccdSweptSphereRadius=ccd_radius_sim,
                              linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                              lateralFriction=friction, physicsClientId=client)

    free_ids = []
    for i in range(n_free):
        pos_sim = free_start_fm[i] * SIM_SCALE
        body = p.createMultiBody(baseMass=CELL_MASS, baseCollisionShapeIndex=ico_shape,
                                  baseVisualShapeIndex=-1, basePosition=list(pos_sim),
                                  baseOrientation=[0, 0, 0, 1], physicsClientId=client)
        p.changeDynamics(body, -1, ccdSweptSphereRadius=ccd_radius_sim,
                          linearDamping=LINEAR_DAMPING, angularDamping=ANGULAR_DAMPING,
                          lateralFriction=friction, physicsClientId=client)
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
                    force_sim = -p_conf * (rel / r_mag)
                    p.applyExternalForce(bid, -1, list(force_sim), list(pos * SIM_SCALE), p.WORLD_FRAME,
                                          physicsClientId=client)

        d_ab = cB - cA
        d_norm = np.linalg.norm(d_ab)
        if d_norm > 1e-9 and approach_mag != 0.0:
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
            if len(contact_samples) >= N_CONTACT_MIN:
                stop_reason = "drivers_touch"
                break
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

    cA_end = np.array(p.getBasePositionAndOrientation(driver_A, physicsClientId=client)[0]) / SIM_SCALE
    cB_end = np.array(p.getBasePositionAndOrientation(driver_B, physicsClientId=client)[0]) / SIM_SCALE
    sep_end_fm = np.linalg.norm(cB_end - cA_end)

    n_crossed = 0
    free_end_fm = np.zeros((n_free, 3))
    for i, bid in enumerate(free_ids):
        pos_end = np.array(p.getBasePositionAndOrientation(bid, physicsClientId=client)[0]) / SIM_SCALE
        free_end_fm[i] = pos_end
        now_closer_to_A = np.linalg.norm(pos_end - cA_end) < np.linalg.norm(pos_end - cB_end)
        if now_closer_to_A != home_is_A[bid]:
            n_crossed += 1
    mixed_fraction = n_crossed / n_free if n_free else 0.0

    p.disconnect(client)

    sep_history_fm = np.array(sep_history_fm)

    contact_samples_arr = np.array(contact_samples) if contact_samples else np.zeros((0, 3))

    out_dir = "/tmp/mesh_grind_friction_scaling_data"
    os.makedirs(out_dir, exist_ok=True)
    mult_tag = sep_start_fm / TOUCH_SEP_FM
    np.savez_compressed(f"{out_dir}/{label}_mult{mult_tag:.1f}_fric{friction:.1f}_seed{condition_seed}.npz",
                         sep_history_fm=sep_history_fm, pen_history_fm=pen_history_fm,
                         contact_samples=contact_samples_arr, orient_A=orient_A, orient_B=orient_B)

    return {"label": label, "sep_start": sep_start_fm, "sep_end": sep_end_fm,
            "delta_sep": sep_end_fm - sep_start_fm, "max_pen_pct": 100 * max_pen_fm_overall / r_in,
            "mixed_fraction": mixed_fraction, "n_steps_run": n_steps_run,
            "stop_reason": stop_reason, "n_free": n_free,
            "n_contact_samples": len(contact_samples_arr)}


def main():
    print("=" * 78)
    print("FRICTION x PRESSURE 2D SCAN (2026-09-28) -- same/hetero only (no control),")
    print("mult=2.0 ONLY (the one separation showing partial signal before). Explicit")
    print("approach force OFF (approach_mag=0). FRICTION_SCALES=[1x,10x,100x] baseline")
    print("100 (->100,1000,10000) CROSSED WITH P_CONF_SCALES=[50x,200x,1000x] baseline")
    print("2.0 (->100,400,2000) -- both pushed well past their prior tested maxima.")
    print("Tests whether ANY combination reaches real contact (the 1D friction-only")
    print("version never did, at any friction level, because P_CONF was fixed at a")
    print("level that itself never produced contact -- uninformative by construction).")
    print("=" * 78)

    omega_vec = OMEGA_MAG * AXIS_COAXIAL
    mult = 2.0
    FRICTION_SCALES = [1.0, 10.0, 100.0]     # x100 -> 100,1000,10000
    P_CONF_SCALES = [50.0, 200.0, 1000.0]    # x2.0 -> 100,400,2000

    start_sep = mult * TOUCH_SEP_FM
    a_confine, b_confine = ellipsoid_semi_axes(mult)
    n_free_target = n_free_for_ellipsoid(a_confine, b_confine)
    print(f"\n--- mult={mult}  n_free={n_free_target} ---")

    all_results = []
    for fj, f_scale in enumerate(FRICTION_SCALES):
        friction_value = LATERAL_FRICTION * f_scale
        for pj, p_scale in enumerate(P_CONF_SCALES):
            p_conf_value = P_CONF * p_scale
            condition_seed = ORIENT_SEED + fj * 10000 + pj * 137
            same = run_condition("SAME-coaxial", omega_vec, omega_vec, stop_mode="repel",
                                  start_sep_fm=start_sep, n_free_target=n_free_target,
                                  approach_mag=0.0, confine_a_fm=a_confine, confine_b_fm=b_confine,
                                  condition_seed=condition_seed, p_conf=p_conf_value,
                                  friction=friction_value)
            hetero = run_condition("HETERO-coaxial", omega_vec, -omega_vec, stop_mode="attract",
                                    start_sep_fm=start_sep, n_free_target=n_free_target,
                                    approach_mag=0.0, confine_a_fm=a_confine, confine_b_fm=b_confine,
                                    condition_seed=condition_seed, p_conf=p_conf_value,
                                    friction=friction_value)
            result = {"friction": friction_value, "p_conf": p_conf_value,
                      "same_delta": same["delta_sep"], "hetero_delta": hetero["delta_sep"],
                      "same_stop": same["stop_reason"], "hetero_stop": hetero["stop_reason"],
                      "same_contacts": same["n_contact_samples"],
                      "hetero_contacts": hetero["n_contact_samples"],
                      "same_gt_hetero": same["delta_sep"] > hetero["delta_sep"]}
            all_results.append(result)
            print(f"  friction={friction_value:.1f} p_conf={p_conf_value:.1f}: "
                  f"same_delta={same['delta_sep']:+.4f}(stop={same['stop_reason']}, "
                  f"contacts={same['n_contact_samples']})  "
                  f"hetero_delta={hetero['delta_sep']:+.4f}(stop={hetero['stop_reason']}, "
                  f"contacts={hetero['n_contact_samples']})  "
                  f"same>hetero={result['same_gt_hetero']}")

    print("\n" + "=" * 78)
    print("FINAL: friction x pressure 2D scan summary (mult=2.0, approach force OFF)")
    print("=" * 78)
    any_hetero_contact = False
    any_same_contact = False
    for r in all_results:
        print(f"friction={r['friction']:.1f} p_conf={r['p_conf']:.1f}: "
              f"same>hetero={r['same_gt_hetero']}  same_delta={r['same_delta']:+.4f}  "
              f"hetero_delta={r['hetero_delta']:+.4f}  hetero_stop={r['hetero_stop']}  "
              f"hetero_contacts={r['hetero_contacts']}  same_contacts={r['same_contacts']}")
        if r["hetero_stop"] == "drivers_touch":
            any_hetero_contact = True
        if r["same_stop"] == "drivers_touch":
            any_same_contact = True
    print(f"\nAny HETERO combination reached real contact? {any_hetero_contact}")
    print(f"Any SAME combination reached real contact (should be False)? {any_same_contact}")


def quick_sanity_check(mult=2.0):
    """Fast smoke test: confirms (a) friction AND p_conf parameters are both
    genuinely applied (no silent fallback to module defaults), (b) the LOW
    corner (friction=100,p_conf=100 -- matches the already-known plateau/no-
    contact result) differs from the HIGH corner (friction=10000,p_conf=2000
    -- the new extreme), before committing to the full 9-combination sweep."""
    start_sep = mult * TOUCH_SEP_FM
    a_confine, b_confine = ellipsoid_semi_axes(mult)
    n_free_target = n_free_for_ellipsoid(a_confine, b_confine)
    omega_vec = OMEGA_MAG * AXIS_COAXIAL

    print(f"SANITY CHECK: mult={mult}, N_FREE={n_free_target}, approach_mag=0.0")
    for friction_value, p_conf_value in [(LATERAL_FRICTION, P_CONF * 50.0),
                                          (LATERAL_FRICTION * 100.0, P_CONF * 1000.0)]:
        condition_seed = ORIENT_SEED
        for label, oA, oB, mode in [("SAME-coaxial", omega_vec, omega_vec, "repel"),
                                     ("HETERO-coaxial", omega_vec, -omega_vec, "attract")]:
            t0 = time.time()
            r = run_condition(label, oA, oB, stop_mode=mode, start_sep_fm=start_sep,
                               n_free_target=n_free_target, approach_mag=0.0,
                               confine_a_fm=a_confine, confine_b_fm=b_confine,
                               condition_seed=condition_seed, p_conf=p_conf_value,
                               friction=friction_value)
            wall_s = time.time() - t0
            print(f"  friction={friction_value:.1f} p_conf={p_conf_value:.1f} {label}: "
                  f"wall={wall_s:.1f}s steps={r['n_steps_run']} stop={r['stop_reason']} "
                  f"delta_sep={r['delta_sep']:+.4f} n_contact_samples={r['n_contact_samples']}")


if __name__ == "__main__":
    import sys
    if "--sanity" in sys.argv:
        args = [a for a in sys.argv[1:] if a != "--sanity"]
        mult_arg = float(args[0]) if args else 2.0
        quick_sanity_check(mult_arg)
    else:
        main()
