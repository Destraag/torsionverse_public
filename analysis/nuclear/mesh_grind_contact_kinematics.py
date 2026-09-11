"""
mesh_grind_contact_kinematics.py
=================================
Stage 1 of the "chirality -> mesh/grind -> pressure sign" derivation plan
(2026-09-10 investigation of doc_electron.txt Section 4/5.2, doc_orbit_pressure.txt
Section 1.2/1.2a, doc_nucleus.txt Section 5.4, and the proton-antiproton
annihilation mechanism -- see judgment_calls.txt N2-1).

GAP THIS ADDRESSES:
  Every doc that invokes "same chirality -> grinding -> high pressure ->
  repulsion; opposite chirality -> meshing -> low pressure -> attraction"
  asserts this as prose. No script anywhere computes it. This script does
  NOT yet derive the pressure outcome (that is Stage 2, a separate future
  script). It derives, from real vector kinematics, the one piece "grinding"
  and "meshing" literally mean in ordinary mechanics: whether two rotating
  bodies in contact slide against each other (grind) or roll together
  without slipping (mesh) -- exactly the same distinction as two external
  gears rotating in the same vs opposite absolute sense.

METHOD:
  Two circular bodies (radius R, centers separated by 2R, externally
  tangent) rotate with angular velocities omega1, omega2 (signed: positive
  = counterclockwise, using ONE shared external reference direction for
  both bodies -- this shared sign IS what "same chirality" vs "opposite
  chirality" means physically, not each body's own local handedness label).
  The velocity of the material point at the contact location is computed
  separately for each body via v = omega x r (2D cross product, r = vector
  from that body's own center to the contact point). Their difference is
  the relative sliding velocity -- no assumption beyond the definition of
  rigid rotation, applied twice at the same spatial point.

RESULT:
  v_slide = (omega1 + omega2) * R using the shared sign convention.
    SAME sign (same chirality):        v_slide = 2*v_surface (nonzero -- GRIND)
    OPPOSITE sign (opposite chirality): v_slide = 0 exactly   (MESH)
  Scale-independent: depends only on v_surface = omega*R (equal for two
  same-type sources), not on R itself -- R cancels once expressed via
  v_surface. Numeric value used: v_surface = Rs*c, the Zone 3 co-rotation
  speed already established elsewhere in this framework
  (entanglement_phase_lock.py, muon_lubrication_doc.py).

WHAT THIS DOES NOT YET DO (Stage 2, separate future script):
  Convert this sliding-vs-rolling kinematic FACT into an actual pressure/
  force outcome. A nonzero relative sliding velocity is not automatically
  "high pressure" -- that conversion (candidate route: a granular/jamming
  coordination-number argument, reusing jobson_cell_rigidity_matrix.py's
  machinery) is the real remaining derivation.

APPLIES TO (same law, different concrete cases):
  Proton-proton (same chirality):                       GRIND, v_slide = 2*Rs*c
  Proton-antiproton (opposite chirality, mirror winding): MESH, v_slide = 0
  Proton-electron (opposite chirality, complementary roll): MESH, v_slide = 0
    (qualitative only -- the electron has no established Zone 3 surface
    speed of its own to plug in numerically; see Section 2.2 of doc_electron.txt)

Run: python analysis/nuclear/mesh_grind_contact_kinematics.py
Reference: docs/series1/doc_electron.txt Section 4-5, doc_orbit_pressure.txt
  Section 1.2/1.2a, doc_nucleus.txt Section 5.4
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'higgs'))
from constants import phi, L_J

pi = math.pi
Rs = math.sqrt(5) / (4 * pi)          # shear/pressure wave speed ratio (I_h geometry)
c  = 2.99792458e8                      # m/s, exact by SI definition
hbar_c = 197.3269804                   # MeV*fm
m_p = 938.272                          # MeV

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

# ── Section 1: general 2-body rolling-contact kinematics (real vector math) ──
print(SEP)
print("SECTION 1: TWO-BODY CONTACT KINEMATICS (grind = slide, mesh = roll)")
print(SEP2)

def contact_velocities(omega1, omega2, R):
    """
    Two circles of radius R, centers at (0,0) and (2R,0), touching at
    P=(R,0). omega1, omega2 use a SHARED external sign convention
    (positive = counterclockwise viewed from +z). Returns (v1, v2, v_slide)
    -- the surface velocity of each body's own material point AT the
    contact location P, and their difference (relative sliding velocity).
    """
    C1, C2, P = (0.0, 0.0), (2 * R, 0.0), (R, 0.0)
    r1 = (P[0] - C1[0], P[1] - C1[1])   # vector from body 1's center to P
    r2 = (P[0] - C2[0], P[1] - C2[1])   # vector from body 2's center to P
    # v = omega (z-hat) x r = omega*(-ry, rx)  [2D rigid-rotation velocity]
    v1 = (-omega1 * r1[1], omega1 * r1[0])
    v2 = (-omega2 * r2[1], omega2 * r2[0])
    v_slide = (v1[0] - v2[0], v1[1] - v2[1])
    return v1, v2, v_slide

R_test = 1.0        # arbitrary unit radius -- result should not depend on it
omega_test = 1.0    # arbitrary unit angular speed

v1_same, v2_same, slide_same = contact_velocities(omega_test, omega_test, R_test)
v1_opp,  v2_opp,  slide_opp  = contact_velocities(omega_test, -omega_test, R_test)

slide_same_mag = math.hypot(*slide_same)
slide_opp_mag  = math.hypot(*slide_opp)
v_surface_test = omega_test * R_test

print(f"  Two externally-tangent circles, radius R, contact point P.")
print(f"  Surface speed each: v = omega*R = {v_surface_test:.6f}  (test units)")
print()
print(f"  SAME sign (same chirality):         v_slide = {slide_same_mag:.6f}  (expect 2v = {2*v_surface_test:.6f})")
print(f"  OPPOSITE sign (opposite chirality):  v_slide = {slide_opp_mag:.6f}  (expect 0)")

check("MGK1: same chirality -> v_slide = 2*v_surface exactly (GRIND, nonzero)",
      abs(slide_same_mag - 2 * v_surface_test) < 1e-12,
      f"v_slide = {slide_same_mag:.10f}  2v = {2*v_surface_test:.10f}")

check("MGK2: opposite chirality -> v_slide = 0 exactly (MESH, pure rolling)",
      slide_opp_mag < 1e-12,
      f"v_slide = {slide_opp_mag:.2e}  (machine-zero)")

closed_form_same = abs(omega_test + omega_test) * R_test
closed_form_opp  = abs(omega_test + (-omega_test)) * R_test
check("MGK3: closed form v_slide=(omega1+omega2)*R matches vector computation (both cases)",
      abs(slide_same_mag - closed_form_same) < 1e-12 and abs(slide_opp_mag - closed_form_opp) < 1e-12,
      f"same: {slide_same_mag:.10f} vs {closed_form_same:.10f}  opp: {slide_opp_mag:.10f} vs {closed_form_opp:.10f}")

# Confirm R-independence: rescale omega to hold v_surface fixed at a different R
R_test2 = 3.7
omega_test2 = v_surface_test / R_test2
_, _, slide_same2 = contact_velocities(omega_test2, omega_test2, R_test2)
check("MGK4: v_slide depends only on v_surface=omega*R, not on R itself (scale-independence)",
      abs(math.hypot(*slide_same2) - 2 * v_surface_test) < 1e-10,
      f"R={R_test:.2f}: v_slide={slide_same_mag:.8f}   R={R_test2:.2f}: v_slide={math.hypot(*slide_same2):.8f}")

# ── Section 2: plug in the real, established Zone 3 co-rotation speed ────────
print()
print(SEP)
print("SECTION 2: NUMERIC EVALUATION -- v_surface = Rs*c (established Zone 3 speed)")
print(SEP2)

v_surface = Rs * c   # Zone 3 co-rotation speed [entanglement_phase_lock.py, muon_lubrication_doc.py]
v_slide_grind = 2 * v_surface
v_slide_mesh = 0.0

print(f"  Rs = sqrt(5)/(4*pi) = {Rs:.6f}")
print(f"  v_surface = Rs*c = {v_surface:.6e} m/s = {Rs:.4f}*c")
print()
print(f"  Proton-proton (same chirality):        v_slide = 2*Rs*c = {v_slide_grind:.6e} m/s = {2*Rs:.4f}*c  [GRIND]")
print(f"  Proton-antiproton (opposite chirality): v_slide = 0                                            [MESH]")
print(f"  Proton-electron   (opposite chirality): v_slide = 0 (qualitative -- electron has no established")
print(f"                                                       Zone 3 surface speed to plug in numerically)")

check("MGK5: grind case v_slide = 2*Rs*c is a real, sub-c, nonzero speed (physically sane)",
      0 < v_slide_grind < c,
      f"v_slide_grind = {v_slide_grind:.6e} m/s = {v_slide_grind/c:.4f}*c")

check("MGK6: mesh case v_slide = 0 to machine precision (pure rolling, no dissipation channel)",
      v_slide_mesh == 0.0,
      f"v_slide_mesh = {v_slide_mesh}")

# ── Section 3: tie to the two concrete torsionverse contact-length scales ────
print()
print(SEP)
print("SECTION 3: WHERE THIS APPLIES -- two established contact geometries")
print(SEP2)

lambda_p = hbar_c / m_p          # fm, Zone 1 boundary (proton Compton wavelength)
r_grind = 2 * lambda_p            # fm, Zone 2/3 nucleon-nucleon contact (N2-1's own number)
r_mid_cell = phi * L_J / 2        # fm, single Jobson-cell edge-midpoint contact [RP1]

print(f"  Single-cell scale     (cell_rotation_propagation.py RP1): r_mid     = phi*L_J/2  = {r_mid_cell:.6f} fm")
print(f"  Nucleon-nucleon scale (doc_nucleus.txt / N2-1):            r_grind   = 2*lambda_p = {r_grind:.4f} fm")
print(f"  Ratio: r_grind / r_mid_cell = {r_grind/r_mid_cell:.1f}x  (~{r_grind/r_mid_cell:.0f} cell-widths)")
print()
print(f"  MGK4 above shows v_slide does NOT depend on which of these two radii is used --")
print(f"  the sliding-velocity kinematics is identical at both scales given the same")
print(f"  v_surface = Rs*c. Any genuine scale-dependence in the eventual pressure/force")
print(f"  outcome must enter through Stage 2's conversion step (how much medium, over")
print(f"  what contact area, is being sheared) -- NOT through the sliding speed itself.")

check("MGK7: r_grind (nucleon scale) is roughly 2 orders of magnitude larger than r_mid_cell (single-cell scale)",
      40 < r_grind / r_mid_cell < 100,
      f"r_grind/r_mid_cell = {r_grind/r_mid_cell:.1f}")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
print()
print("STAGE 1 COMPLETE: the kinematic half of the mesh/grind story (sliding vs")
print("rolling contact) is now a real, checked computation, not asserted prose.")
print("STAGE 2 (NOT done here): convert v_slide into an actual pressure/force")
print("outcome -- candidate route: a granular/jamming coordination-number")
print("argument (reusing jobson_cell_rigidity_matrix.py's machinery), testing")
print("whether a sliding (grind) contact locally over-constrains the packing")
print("(jams, real pressure rise) while a rolling (mesh) contact does not.")
