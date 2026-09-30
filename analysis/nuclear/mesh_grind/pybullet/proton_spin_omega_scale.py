"""
proton_spin_omega_scale.py
================================================
Author's direction (2026-09-20): "axis and spin should be based on the total
measured spin of the proton scaled to our model" -- replaces the PyBullet
Zone-2 driver test's OMEGA_MAG=0.3 / axis=[1,1,1] (both explicit placeholders,
"arbitrary units" per mesh_grind_pybullet_zone2_driver_test.py's own comments)
with a value chained from real, established quantities. This is ONLY about
the driver's spin-rate magnitude -- NOT an attempt to resolve the proton
spin puzzle itself (docs/open_items.txt N-15: quark spin alone is only
~25-30% of the proton's total 1/2 hbar spin) -- that stays explicitly parked
per the author's own framing this session. The quantity used here (the
proton's TOTAL spin, 1/2 hbar) is the uncontroversial, textbook part, not
the disputed internal quark/gluon decomposition.

CHAIN (every step real, no inline guessing -- ABSOLUTE RULE compliant):
  1. m_cell = E_cell/c^2 -- E=mc^2, the SAME convention already used for
     neutron/proton mass-from-volume (claims_audit.txt M8e). This framework
     has no OTHER independent cell-mass convention (already flagged in
     jobson_cell_natural_force_pressure_scale.py, not invented here).
  2. M_zone2 = N_zone2 * m_cell. N_zone2 = 2.846e4 is REUSED directly from
     zone2_cell_count_estimate.py's own printed output (literal Zone 2,
     r~lambda_p, TRUE icosahedron volume + phi_rcp=0.64 realistic packing)
     -- not recomputed by hand, taken from that already-committed script.
  3. I_zone2 = (2/3)*M_zone2*lambda_p^2 -- THIN SPHERICAL SHELL moment of
     inertia, because doc_nucleus.txt Section 1.1 explicitly describes
     Zone 2 as "a thin transition" at r~lambda_p, not a filled ball
     (which would need the (2/5)*M*R^2 solid-sphere formula instead).
  4. omega_bulk = (hbar/2) / I_zone2 -- the real angular rate a rigid
     Zone-2 shell of THIS mass, at THIS radius, must have to carry the
     proton's actual measured total spin (1/2 hbar, matching
     open_items.txt N-15's own wording).
  5. CROSS-CHECK (reported plainly, not interpreted -- N-15 stays parked):
     compare omega_bulk*lambda_p (the bulk equatorial speed this implies)
     against the framework's OWN already-established v_surface = Rs*c
     local co-rotation speed (used throughout the mesh-grind derivation,
     mesh-grind-chirality-derivation.md), evaluated at the SAME radius
     lambda_p (apples-to-apples -- NOT r_mid, the much-smaller single-cell
     contact scale used elsewhere, which would compare a bulk-structure
     rate against a cell-contact rate, a different question).
  6. SIM MAPPING: this framework has NO established time non-dimensionalization
     for the PyBullet engine (only length is rescaled, SIM_SCALE=1/L_J;
     jobson_cell_natural_force_pressure_scale.py flagged the same gap for
     force/mass). Injecting a raw SI omega_bulk directly already failed
     catastrophically once before (mesh_grind_cluster_spin_test.py's own
     "BUG CAUGHT AND FIXED" entry -- ~9 orders of magnitude too large,
     simulation exploded). So: derive the LARGEST numerically-safe sim
     rate under the same "1 degree per physics step" cap already accepted
     as standing policy for this derivation line (mesh-grind-chirality-
     derivation.md SHORTCUTS LOG item 9), and report explicitly how many
     orders of magnitude short of the real omega_bulk that safe rate is --
     an honest gap, not glossed over.

Run: python analysis/nuclear/mesh_grind/pybullet/proton_spin_omega_scale.py
"""

import math

pi = math.pi
alpha = 7.2973525693e-3
hbar_c_MeVfm = 197.3269804          # MeV*fm
hbar_c_GeVfm = hbar_c_MeVfm / 1000.0
m_p_MeV = 938.272
phi = (1 + math.sqrt(5)) / 2
c_SI = 2.99792458e8                  # m/s
hbar_SI = 1.054571817e-34            # J*s, CODATA
GeV_to_J = 1.602176634e-10
MeV_to_J = GeV_to_J / 1000.0
fm_to_m = 1e-15

r_p_fm = 0.8414                                  # CODATA proton charge radius
L_J_fm = alpha * phi * r_p_fm
E_cell_GeV = 2 * pi * hbar_c_GeVfm / L_J_fm       # canonical formula, matches every doc/script
lambda_p_fm = hbar_c_MeVfm / m_p_MeV              # Zone 1/2 boundary = 0.2103 fm

# REUSED directly from zone2_cell_count_estimate.py's own printed run
# ("N cells, TRUE volume + 0.64 random-packing frac" for the ZONE 2 LITERAL row):
N_ZONE2_REAL_PACKED = 2.846e4

Rs = math.sqrt(5) / (4 * pi)                      # established elsewhere (proton_g_factor.py)

MAX_DEG_PER_STEP = 1.0                            # already-accepted numerical safety cap
FIXED_TIMESTEP_SIM = 1.0 / 60.0                   # matches mesh_grind_pybullet_zone2_driver_test.py


def check(label, value, unit=""):
    print(f"  {label} = {value:.6e} {unit}")


def main():
    print("=" * 78)
    print("PROTON TOTAL-SPIN -> ZONE-2 BULK ANGULAR RATE -> PYBULLET OMEGA_MAG")
    print("=" * 78)

    E_cell_J = E_cell_GeV * GeV_to_J
    m_cell_kg = E_cell_J / c_SI**2
    print("\nSTEP 1: cell mass from E=mc^2 (no other convention exists in this repo)")
    check("E_cell", E_cell_GeV, "GeV")
    check("m_cell", m_cell_kg, "kg")

    M_zone2_kg = N_ZONE2_REAL_PACKED * m_cell_kg
    print("\nSTEP 2: total Zone-2 mass (N reused from zone2_cell_count_estimate.py)")
    check("N_zone2 (real, packed)", N_ZONE2_REAL_PACKED)
    check("M_zone2", M_zone2_kg, "kg")

    lambda_p_m = lambda_p_fm * fm_to_m
    I_zone2 = (2.0 / 3.0) * M_zone2_kg * lambda_p_m**2   # thin spherical shell
    print("\nSTEP 3: moment of inertia (thin spherical shell at r=lambda_p)")
    check("lambda_p", lambda_p_m, "m")
    check("I_zone2 = (2/3)*M*lambda_p^2", I_zone2, "kg*m^2")

    L_target = hbar_SI / 2.0
    omega_bulk = L_target / I_zone2
    print("\nSTEP 4: real bulk angular rate carrying the proton's measured total spin")
    check("L_target = hbar/2", L_target, "J*s")
    check("omega_bulk = L_target / I_zone2", omega_bulk, "rad/s")

    v_bulk_at_lambda_p = omega_bulk * lambda_p_m
    v_surface = Rs * c_SI
    omega_local_established = v_surface / lambda_p_m
    ratio = omega_bulk / omega_local_established
    print("\nSTEP 5: cross-check against the already-established v_surface=Rs*c "
          "(same radius lambda_p, apples-to-apples -- reported plainly, N-15 stays parked)")
    check("v_bulk_at_lambda_p = omega_bulk*lambda_p", v_bulk_at_lambda_p, "m/s")
    check("v_surface = Rs*c", v_surface, "m/s")
    check("omega_local_established = v_surface/lambda_p", omega_local_established, "rad/s")
    check("ratio omega_bulk / omega_local_established", ratio)

    omega_sim_safe = math.radians(MAX_DEG_PER_STEP) / FIXED_TIMESTEP_SIM
    gap_orders_of_magnitude = math.log10(omega_bulk / omega_sim_safe)
    print("\nSTEP 6: largest numerically-safe PyBullet OMEGA_MAG "
          "(1 deg/step cap, already-accepted convention)")
    check("omega_sim_safe = radians(1deg)/fixedTimeStep", omega_sim_safe, "rad/(sim step-ish)")
    print(f"  gap vs real omega_bulk: {gap_orders_of_magnitude:.2f} orders of magnitude "
          f"(real rate CANNOT be run directly -- same failure mode already hit once before)")

    print("\nSTEP 7 (2026-09-20, author correction -- this IS a real timescale, don't dismiss it):")
    print("  1 sim step = 1 degree of driver rotation, by construction (STEP 6). If that 1 degree")
    print("  represents the SAME physical rotation the real Zone-2 shell undergoes at omega_bulk,")
    print("  then 1 sim step maps to a real duration:")
    real_s_per_step = math.radians(MAX_DEG_PER_STEP) / omega_bulk
    check("real seconds per sim step = radians(1deg)/omega_bulk", real_s_per_step, "s")
    r_p_m = r_p_fm * fm_to_m
    t_nuclear = r_p_m / c_SI
    t_protonium = 1.0e-6   # measured protonium mean lifetime (Wikipedia, Abdel-Raouf 2009)
    for n_steps in (400, 800):
        real_s = n_steps * real_s_per_step
        print(f"  N_STEPS={n_steps}: covers {real_s:.4e} real s = "
              f"{real_s/t_nuclear:.4e}x the nuclear/strong timescale (r_p/c={t_nuclear:.4e}s), "
              f"{real_s/t_protonium:.4e}x the protonium orbital lifetime ({t_protonium:.1e}s)")

    print("\n" + "=" * 78)
    print("CONCLUSION")
    print("=" * 78)
    print("The proton's measured total spin (1/2 hbar) implies a real Zone-2 bulk")
    print("rotation rate many orders of magnitude beyond what any discrete-time")
    print("contact solver can run directly (same conclusion as F_cell's own force-")
    print("scale mismatch in jobson_cell_natural_force_pressure_scale.py -- this repo")
    print("has no time/mass non-dimensionalization yet, only length via SIM_SCALE).")
    print("PROPOSAL: drive the PyBullet test at omega_sim_safe (the largest rate the")
    print("solver can resolve without tunneling/aliasing) -- i.e. the physically-real")
    print("number sets the CEILING/justification for 'run it as fast as numerically")
    print("possible' rather than supplying a literal, directly-usable magnitude.")
    print("Axis: apply this magnitude along the already-agreed 3-axis (Cartesian)")
    print("convention, sign-flipped between same/hetero drivers.")


if __name__ == "__main__":
    main()
