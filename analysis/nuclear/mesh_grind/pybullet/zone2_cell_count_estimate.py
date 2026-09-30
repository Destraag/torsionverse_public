"""
zone2_cell_count_estimate.py
================================
Direct answer to the author's question (2026-09-16): "the zone 2 driver in these
tests is meant to be a full zone 2 model. how many cells are required to make a
spherical shape equivalent to our zone 2 radius?"

REUSES the already-established, canonical geometry from analysis/nuclear/
proton_structure.py (PS1-PS4, all PASS) rather than re-deriving lambda_p, r_p, L_J
from scratch. That script already computes N_disp = V_sphere/V_cell for r < r_p using
a SIMPLE CUBE approximation for the cell volume (V_cell = L_J^3, no packing-fraction
correction, i.e. an idealized "no gaps" count) -- this script:
  (1) reproduces that r_p number for a direct cross-check,
  (2) does the SAME calculation for r ~ lambda_p specifically, since "Zone 2" per
      doc_nucleus.txt's own zone definitions IS the r~lambda_p Maxwell-jamming
      BOUNDARY (a thin transition), not r_p (that is Zone 3's OUTER edge) -- the
      literal "zone 2 radius" is lambda_p, not r_p, though both are reported here
      since casual conversation this session used "zone 2" loosely for the whole
      lambda_p-to-r_p structure,
  (3) adds a REALISTIC packing correction: true icosahedron volume (not a cube
      approximation) and a real, external, literature-standard packing fraction
      (random close packing of spheres, phi_rcp = 0.64 -- NOT a torsionverse-derived
      number, cited as such, and the closest available real-physics analogue given
      no torsionverse-native disordered-icosahedra packing-fraction result exists) --
      to see how much the idealized cube-count over/under-estimates a realistic count.

IMPORTANT CORRECTION (flagged honestly): last session's "~830 cells" estimate for a
"real-scale Zone 2 shell" (mesh-grind-chirality-derivation.md, kissing-number-based
reasoning: ~13-cell shell, loosely scaled by r_p/lambda_p) was NOT a real volume
calculation -- it was a rough extrapolation. This script's actual volume-based count
is many orders of magnitude larger (see RESULTS). The ~830 figure should be treated as
superseded/wrong, not as a lower-effort target scale for the PyBullet migration.

Run: python analysis/nuclear/mesh_grind/pybullet/zone2_cell_count_estimate.py
"""

import math

pi = math.pi
alpha = 7.2973525693e-3
hbar_c = 197.3269804          # MeV*fm
m_p = 938.272                 # MeV
phi = (1 + math.sqrt(5)) / 2

lambda_p_fm = hbar_c / m_p            # Zone 1/2 boundary = 0.2103 fm
r_p_fm = 0.8414                       # fm, CODATA 2018 (measured, canonical -- matches proton_structure.py)
L_J_fm = alpha * phi * r_p_fm         # Jobson cell edge length

r_in_fm = L_J_fm * phi**2 / (2 * math.sqrt(3))     # inradius
R_c_fm = L_J_fm * math.sqrt(1 + phi**2) / 2         # circumradius

# TRUE icosahedron volume (edge length a): V = (5/12)*(3+sqrt(5))*a^3 -- exact,
# standard solid-geometry formula, not a cube approximation.
V_cell_true_fm3 = (5.0 / 12.0) * (3 + math.sqrt(5)) * L_J_fm**3
V_cell_cube_fm3 = L_J_fm**3  # proton_structure.py's own simpler convention, reproduced for cross-check

PHI_RCP = 0.64  # random close packing fraction for spheres -- EXTERNAL, standard
                # physics constant (Bernal/Scott, ~0.64), NOT torsionverse-derived;
                # no torsionverse-native disordered-icosahedra packing result exists
                # to use instead, so this is cited as an external cross-check only.


def sphere_volume(r_fm):
    return (4.0 / 3.0) * pi * r_fm**3


def report(label, r_fm):
    V_sphere = sphere_volume(r_fm)
    n_cube_idealized = V_sphere / V_cell_cube_fm3          # matches proton_structure.py's own convention
    n_true_idealized = V_sphere / V_cell_true_fm3           # true icosahedron volume, still 0-gap idealized
    n_true_packed = n_true_idealized / PHI_RCP               # true volume + realistic packing fraction
    print(f"\n{label}: R = {r_fm:.4f} fm  (R/L_J = {r_fm/L_J_fm:.2f})")
    print(f"  sphere volume                                = {V_sphere:.6e} fm^3")
    print(f"  N cells, cube approx, 0 gaps (proton_structure.py convention) = {n_cube_idealized:.4e}")
    print(f"  N cells, TRUE icosahedron volume, 0 gaps      = {n_true_idealized:.4e}")
    print(f"  N cells, TRUE volume + {PHI_RCP} random-packing frac  = {n_true_packed:.4e}")
    return n_true_packed


def main():
    print("=" * 78)
    print("ZONE 2/3 REAL-SCALE CELL-COUNT ESTIMATE")
    print("=" * 78)
    print(f"L_J = {L_J_fm:.6f} fm, r_in = {r_in_fm:.6f} fm, R_c = {R_c_fm:.6f} fm")
    print(f"lambda_p (Zone 1/2 boundary) = {lambda_p_fm:.4f} fm")
    print(f"r_p (Zone 3 outer edge)      = {r_p_fm:.4f} fm")

    n_lambda_p = report("ZONE 2 LITERAL (r ~ lambda_p, the Maxwell-jamming boundary "
                         "itself, per doc_nucleus.txt's own zone definitions)", lambda_p_fm)
    n_r_p = report("ZONE 3 OUTER EDGE (r ~ r_p, what was loosely called \"zone 2 scale\" "
                    "last session)", r_p_fm)

    print("\n" + "=" * 78)
    print("CROSS-CHECK against analysis/nuclear/proton_structure.py's own N_disp:")
    print(f"  that script's r_p cube-approx count = 2.545e+06 (reproduced above as "
          f"{sphere_volume(r_p_fm)/V_cell_cube_fm3:.3e})")
    print("\nCORRECTION vs last session's ~830-cell estimate (kissing-number "
          "extrapolation, NOT a volume calc):")
    print(f"  literal Zone 2 (lambda_p) real-packing count = {n_lambda_p:.3e} -- "
          f"~{n_lambda_p/830:.1e}x last session's rough figure")
    print(f"  Zone 3 outer edge (r_p) real-packing count   = {n_r_p:.3e} -- "
          f"~{n_r_p/830:.1e}x last session's rough figure")
    print("Both real numbers are MILLIONS of cells, not hundreds -- the ~830 estimate")
    print("is superseded/wrong and should not be used as a PyBullet scaling target.")
    print("=" * 78)


if __name__ == "__main__":
    main()
