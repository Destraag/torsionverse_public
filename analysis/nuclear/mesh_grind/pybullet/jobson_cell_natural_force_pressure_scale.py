"""
jobson_cell_natural_force_pressure_scale.py
================================================
Author's request (2026-09-16): "the trick is trying to keep the balance as closely
rooted to real expectations as possible based on what we have already calculated" --
this replaces the PyBullet Zone-2 driver test's placeholder P_CONF (an arbitrary
confinement-force magnitude, "not yet derived" per the standing shortcut already
logged for the numpy engine) with a value derived from ALREADY-ESTABLISHED
torsionverse constants: E_cell (Jobson cell restoration energy) and L_J (cell edge
length), plus a comparison against the medium's own established bulk modulus K.

NATURAL SCALES FROM E_cell AND L_J:
  F_cell = E_cell / L_J     -- natural FORCE scale (energy per unit length) for a
                               Jobson cell being displaced by order its own size.
  P_cell = E_cell / L_J^3   -- natural PRESSURE scale (energy density) for the cell.

WHY THIS MATTERS FOR THE CONFINEMENT-FORCE TUNING PROBLEM: an "ambient bounding
pressure" that is supposed to hold a medium together WITHOUT itself compressing it
significantly must be much smaller than the medium's own restoring/elastic force
scale -- otherwise the confinement itself would be strong enough to jam/deform the
cells on its own, which is not what "one common external bounding pressure" is meant
to represent (it should be a background confinement, not a dominant force). This
script computes the RATIO P_cell/K (K = 1/eps_0, the already-established medium bulk
modulus) to see how "gentle" the cell's own natural pressure scale already is relative
to the bulk medium -- and proposes that any CONFINEMENT pressure used in a simulation
should sit well BELOW P_cell itself (a cell's own restoring pressure), not above it,
since confinement is an external, softer effect layered on top of the cell's own
already-established stiffness.

CONCRETE OUTPUT FOR THE PYBULLET SIM: converts F_cell into the sim's own SIM_SCALE
unit system (length only rescaled by 1/L_J; this script does NOT attempt a full
mass/time non-dimensionalization -- that would require a real "mass of a Jobson
cell" convention this framework does not currently define outside of E=mc^2, which
is noted but not adopted as a literal dynamical mass here) and proposes P_CONF as a
small FRACTION of the sim's own contact-stiffness force scale, rather than an
unconstrained guess.

CORRECTION (2026-09-20, author asked "do we have a torsionverse scalable pressure" --
re-audit of this script's own premise): the "contact_stiffness_sim = 1e6 (already set
in mesh_grind_pybullet_zone2_driver_test.py)" claim below is WRONG -- grepped the real
PyBullet scripts directly: NONE of them ever call changeDynamics with contactStiffness/
contactDamping. They use PyBullet's DEFAULT constraint-based (ERP/CFM) contact solver
(erp=0.9, contactERP=0.9, globalCFM=1e-7), not a literal penalty spring -- there is no
real "1e6" force-per-overlap constant anywhere in this sim. The P_CONF=1000 proposal
computed below rests on this fictional premise and was NEVER actually adopted anyway
(the real scripts use P_CONF=2.0, carried over from the single-driver test's own
"informed but provisional" placeholder, not this script's own output). SECTION 3 below
adds the real mass/time non-dimensionalization attempt this docstring says wasn't
tried, and a check that uses the sim's ACTUAL contact model instead of an imagined one.

Run: python analysis/nuclear/mesh_grind/pybullet/jobson_cell_natural_force_pressure_scale.py
"""

import math

pi = math.pi
alpha = 7.2973525693e-3
hbar_c_MeVfm = 197.3269804        # MeV*fm
hbar_c_GeVfm = hbar_c_MeVfm / 1000.0
m_p = 938.272                     # MeV
phi = (1 + math.sqrt(5)) / 2

r_p_fm = 0.8414
L_J_fm = alpha * phi * r_p_fm
r_in_fm = L_J_fm * phi**2 / (2 * math.sqrt(3))   # matches the PyBullet scripts' own formula
E_cell_GeV = 2 * pi * hbar_c_GeVfm / L_J_fm   # canonical formula, matches every doc script

# unit conversions
GeV_to_J = 1.602176634e-10
fm_to_m = 1e-15
c_SI = 2.99792458e8
eps_0 = 8.8541878128e-12   # F/m
K_medium_Pa = 1.0 / eps_0  # already-established medium bulk modulus (doc_higgs C7)

E_cell_J = E_cell_GeV * GeV_to_J
L_J_m = L_J_fm * fm_to_m

F_cell_N = E_cell_J / L_J_m              # natural force scale, Newtons
P_cell_Pa = E_cell_J / L_J_m**3           # natural pressure scale, Pa

SIM_SCALE = 1.0 / L_J_fm  # matches the PyBullet scripts' own convention (length only)


def main():
    print("=" * 78)
    print("NATURAL JOBSON-CELL FORCE/PRESSURE SCALE FROM E_cell AND L_J")
    print("=" * 78)
    print(f"E_cell = {E_cell_GeV:.4f} GeV = {E_cell_J:.4e} J")
    print(f"L_J    = {L_J_fm:.6f} fm = {L_J_m:.4e} m")
    print()
    print(f"F_cell = E_cell / L_J = {F_cell_N:.4e} N   (natural cell restoring-force scale)")
    print(f"P_cell = E_cell / L_J^3 = {P_cell_Pa:.4e} Pa   (natural cell pressure scale)")
    print()
    print(f"Medium bulk modulus K = 1/eps_0 = {K_medium_Pa:.4e} Pa  [doc_higgs C7, established]")
    print(f"P_cell / K = {P_cell_Pa/K_medium_Pa:.4e}")
    print()
    print("READING: P_cell is many orders of magnitude LARGER than K, not smaller --")
    print("P_cell is the cell's own internal (Zone-1/2-scale) restoring pressure, a")
    print("hadronic-scale quantity; K is the medium's BULK (macroscopic/EM-scale)")
    print("modulus -- these are the SAME two-regime distinction already established")
    print("in doc_torsion.txt Section 3.3 (K_fluid vs K_jammed) applied one level")
    print("further: an 'ambient CONFINEMENT pressure' standing in for whatever holds")
    print("a many-cell cluster together should be a SMALL FRACTION of P_cell (the")
    print("cell's OWN hadronic-scale restoring force), not of K (the wrong, much")
    print("softer macroscopic scale) -- using K here would make confinement")
    print("comically weak relative to a cell's actual stiffness.")
    print()
    print("=" * 78)
    print("CONCRETE PYBULLET SIM PROPOSAL")
    print("=" * 78)
    print(f"SIM_SCALE = 1/L_J = {SIM_SCALE:.4e} (fm^-1) -- length rescaling only,")
    print("matching the existing PyBullet scripts' convention.")
    print()
    print("This framework does not (yet) define a literal dynamical 'mass of a")
    print("Jobson cell' independent of E=mc^2 (m_cell = E_cell/c^2 would give a")
    print("mass, but nothing in this repo currently uses cell mass dynamically --")
    print("flagged honestly, NOT invented here). Without a matched mass/time")
    print("non-dimensionalization, F_cell cannot be converted into an exact PyBullet")
    print("force number -- but the RATIO argument above stands regardless of that:")
    print()
    print("PROPOSAL: set the sim's confinement force (P_CONF, per free cell) to a")
    print("SMALL FRACTION (order 1-5%) of the sim's OWN contact-stiffness force scale")
    print("at a typical small overlap, rather than an independent unconstrained guess")
    print("-- i.e. confinement should be a gentle background pull, clearly")
    print("subdominant to contact repulsion, consistent with the P_cell >> K reading")
    print("above (the cell's own restoring stiffness dominates; confinement is the")
    print("much weaker, purely orienting background effect).")
    print()
    print("[ORIGINAL PROPOSAL, KEPT + LABELED WRONG, not deleted, to show the error]")
    frac = 0.02
    contact_stiffness_sim = 1e6      # FICTIONAL -- grep confirms no real script sets this
    typical_overlap_sim = 0.05
    contact_force_scale_sim = contact_stiffness_sim * typical_overlap_sim
    p_conf_proposed_wrong = frac * contact_force_scale_sim
    print(f"  [WRONG PREMISE] assumed contact stiffness (sim) = {contact_stiffness_sim:.1e} "
          f"-- NOT actually configured anywhere")
    print(f"  [WRONG PREMISE] proposed P_CONF = {frac*100:.0f}% of that = {p_conf_proposed_wrong:.2e} "
          f"-- DO NOT USE, rests on a fictional parameter")

    print()
    print("=" * 78)
    print("SECTION 3 (2026-09-20): REAL MASS/TIME NON-DIMENSIONALIZATION ATTEMPT")
    print("=" * 78)
    m_cell_kg = E_cell_J / c_SI**2
    t_cell_s = math.sqrt(m_cell_kg * L_J_m / F_cell_N)   # dimensional analysis: F=m*L/t^2
    print(f"m_cell = E_cell/c^2 = {m_cell_kg:.4e} kg  (same convention as proton_spin_omega_scale.py)")
    print(f"t_cell = sqrt(m_cell*L_J/F_cell) = {t_cell_s:.4e} s")
    print("  (natural timescale for a cell's own mass to be accelerated across its own")
    print("  size by its own restoring force -- by construction F_cell=1 EXACTLY in these units)")
    print()
    fixed_timestep_sim = 1.0 / 60.0
    real_seconds_per_step = fixed_timestep_sim * t_cell_s
    print(f"If 1 sim time unit = 1 t_cell, one physics step (dt=1/60) covers "
          f"{real_seconds_per_step:.4e} real seconds --")
    print(f"covering even 1 real microsecond would need {1e-6/real_seconds_per_step:.4e} steps.")
    print("SAME conclusion as proton_spin_omega_scale.py's omega finding: there is NO fully")
    print("self-consistent 'natural' unit system usable by a discrete solver -- hadronic-scale")
    print("dynamics are astronomically fast relative to anything a ~100-1000 step run can")
    print("resolve. P_CONF cannot be a literal unit conversion of F_cell, for the same reason")
    print("omega_bulk could not be injected directly into OMEGA_MAG.")

    print()
    print("=" * 78)
    print("SECTION 4: IS THE CURRENT P_CONF=2.0 'TOO HIGH' -- CHECKED AGAINST WHAT THE SIM")
    print("ACTUALLY DOES (real F=ma, constraint-based contacts, NOT a fictional spring)")
    print("=" * 78)
    P_CONF_current = 2.0
    CELL_MASS_sim = 1.0
    accel_sim = P_CONF_current / CELL_MASS_sim
    dv_per_step = accel_sim * fixed_timestep_sim
    r_in_sim = r_in_fm * SIM_SCALE
    print(f"acceleration from P_CONF alone = P_CONF/CELL_MASS = {accel_sim:.4f} (sim units)")
    print(f"velocity gained in ONE physics step (dt=1/60)      = {dv_per_step:.4f} (sim units)")
    print(f"r_in in sim units (a cell's own inradius)          = {r_in_sim:.4f}")
    print(f"dv_per_step / r_in = {dv_per_step/r_in_sim:.4f} -- confinement alone could move a")
    print("cell roughly this fraction of its own size (in velocity-per-step terms) before any")
    print("contact resistance -- the REAL, meaningful comparison given the sim's actual")
    print("(constraint-based) contact model, not a comparison to a stiffness that doesn't exist.")


if __name__ == "__main__":
    main()
