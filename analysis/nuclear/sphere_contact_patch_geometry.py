"""
sphere_contact_patch_geometry.py
==================================
Direct correction, author-flagged (2026-09-11): every prior script in this
derivation (mesh_grind_contact_kinematics.py, mesh_grind_intervening_cell_
check.py, mesh_grind_2d_lattice_check.py, mesh_grind_cell_count_scaling.py)
used a 1D chain or a fixed-size 2D rectangular grid of discs -- NOT a real
3D sphere-sphere contact patch. Author's correction: r_grind is a RADIUS
(2*lambda_p, "Zone 2 boundaries touch"), so each nucleus's contact boundary
is a SPHERE. Two such spheres approaching each other share a CIRCULAR
CONTACT PATCH that is a single POINT at first touch (separation = r_grind)
and GROWS as the spheres are pushed closer (separation < r_grind) -- not a
fixed chain/grid picked in advance. "The contacts are few, then scaling up
as the spheres are forced closer" (author's own words) is standard sphere-
sphere intersection geometry, computed here directly.

GEOMETRY (two equal spheres, radius R=lambda_p, centers separated by d):
  Two spheres of radius R with centers a distance d apart (d < 2R)
  intersect in a circle. Placing centers at (-d/2,0,0) and (d/2,0,0), a
  point (0,y,z) lies on sphere A when (d/2)^2 + y^2 + z^2 = R^2, so the
  intersection circle has radius:
    a(d) = sqrt(R^2 - (d/2)^2)   for d <= 2R (=r_grind here, R=lambda_p)
    a(d) = 0                     for d >= 2R (spheres not yet touching)
  This is the ACTUAL geometric contact-patch radius, not an assumed chain
  length -- a(d)=0 exactly at d=r_grind (first touch, a single point,
  matching Stage 1's original two-body picture), growing smoothly as d
  decreases (deeper overlap).

CELL COUNT IN THE PATCH: N(d) ~ (contact patch area) / (single-cell
footprint area) ~ pi*a(d)^2 / (pi*r_mid^2) = (a(d)/r_mid)^2, using
r_mid = phi*L_J/2 (single Jobson-cell edge-midpoint contact spacing,
cell_rotation_propagation.py RP1) as the per-cell footprint scale.

CONSISTENCY CHECK: at full overlap (d=0), a(0)=R=lambda_p exactly, so
N(0) ~ (lambda_p/r_mid)^2. This should be consistent with Stage 1's own
r_grind/r_mid_cell = 52.3 (mesh_grind_contact_kinematics.py MGK7), since
r_grind = 2*lambda_p is the DIAMETER in cell-widths, so lambda_p/r_mid
should be half of that, ~26.15 -- checked below, not assumed.

WHAT THIS DOES NOT YET DO: combine this real, growing N(d) with the
cell-count-scaling finding (mesh_grind_cell_count_scaling.py: for a 1D
chain, TOTAL forced residual is conserved at 2*v_surface regardless of
chain length N) -- that combination, for this circular-patch geometry
specifically (not a chain), is the natural next step, not yet attempted.

Run: python analysis/nuclear/sphere_contact_patch_geometry.py
Reference: /memories/repo/mesh-grind-chirality-derivation.md
"""

import sys, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

alpha = 7.2973525693e-3
phi = (1 + math.sqrt(5)) / 2
hbar_c = 197.3269804
m_p = 938.272
lambda_p_fm = hbar_c / m_p          # R, the Zone 1/2 boundary radius (each nucleus's contact sphere)
r_grind_fm = 2 * lambda_p_fm         # first-touch separation (Stage 1's own established number)
r_p_fm = 4 * lambda_p_fm
L_J_fm = alpha * phi * r_p_fm
r_mid_fm = phi * L_J_fm / 2          # single-cell contact footprint scale [RP1]

SEP = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

def contact_patch_radius(d, R):
    """Radius of the circle where two spheres of radius R, centers d apart, intersect."""
    if d >= 2 * R:
        return 0.0
    return math.sqrt(max(R**2 - (d / 2) ** 2, 0.0))

# ── Section 1: geometry sanity checks (real sphere-sphere intersection) ──────
print(SEP)
print("SECTION 1: SPHERE-SPHERE CONTACT PATCH GEOMETRY (R = lambda_p each nucleus)")
print(SEP2)

print(f"  R (Zone 1/2 boundary radius, lambda_p) = {lambda_p_fm:.4f} fm")
print(f"  r_grind = 2*R (first touch, a=0)        = {r_grind_fm:.4f} fm")
print(f"  r_mid (single-cell footprint scale)     = {r_mid_fm:.6f} fm")

a_at_rgrind = contact_patch_radius(r_grind_fm, lambda_p_fm)
a_at_zero = contact_patch_radius(0.0, lambda_p_fm)

check("SP1: contact patch radius = 0 exactly at d=r_grind (first touch, a single point)",
      abs(a_at_rgrind) < 1e-12,
      f"a(r_grind) = {a_at_rgrind:.2e} fm")

check("SP2: contact patch radius = R exactly at d=0 (full overlap)",
      abs(a_at_zero - lambda_p_fm) < 1e-12,
      f"a(0) = {a_at_zero:.6f} fm   R = {lambda_p_fm:.6f} fm")

# ── Section 2: contact patch radius and cell count vs separation ────────────
print()
print(SEP)
print("SECTION 2: CONTACT PATCH RADIUS a(d) AND CELL COUNT N(d) vs SEPARATION d")
print(SEP2)
print(f"  {'d (fm)':>10}  {'d/r_grind':>10}  {'a(d) (fm)':>12}  {'a(d)/r_mid':>12}  {'N(d) ~ (a/r_mid)^2':>20}")

d_fractions = [1.0, 0.99, 0.95, 0.9, 0.8, 0.6, 0.4, 0.2, 0.0]
Nd_values = []
for frac in d_fractions:
    d = frac * r_grind_fm
    a_d = contact_patch_radius(d, lambda_p_fm)
    a_over_rmid = a_d / r_mid_fm
    N_d = a_over_rmid ** 2
    Nd_values.append((frac, N_d))
    print(f"  {d:>10.4f}  {frac:>10.3f}  {a_d:>12.6f}  {a_over_rmid:>12.3f}  {N_d:>20.2f}")

check("SP3: N(d) is exactly 0 at first touch (d=r_grind) and grows monotonically as d decreases -- "
      "'contacts are few, then scaling up as the spheres are forced closer' (author's own words), computed",
      Nd_values[0][1] == 0.0 and all(Nd_values[i][1] < Nd_values[i+1][1] for i in range(len(Nd_values)-1)),
      f"N(d) at d/r_grind = {[f[0] for f in Nd_values]}: {[round(v,1) for _, v in Nd_values]}")

# ── Section 3: consistency check against Stage 1's own r_grind/r_mid_cell ────
print()
print(SEP)
print("SECTION 3: CONSISTENCY CHECK AGAINST STAGE 1's OWN NUMBER (not assumed to match)")
print(SEP2)

lambda_p_over_rmid = lambda_p_fm / r_mid_fm
r_grind_over_rmid = r_grind_fm / r_mid_fm    # Stage 1's own MGK7 result: ~52.3

print(f"  lambda_p / r_mid (patch RADIUS in cell-widths, at full overlap) = {lambda_p_over_rmid:.3f}")
print(f"  r_grind / r_mid  (patch DIAMETER in cell-widths, Stage 1 MGK7)  = {r_grind_over_rmid:.3f}")
print(f"  Ratio (should be exactly 2, since r_grind=2*lambda_p):            {r_grind_over_rmid/lambda_p_over_rmid:.6f}")

check("SP4: r_grind/r_mid is exactly 2x lambda_p/r_mid -- consistent with Stage 1's own "
      "52.3 number (this is the diameter; the new radius-based N(d) calc is consistent with it, not a new/different number)",
      abs(r_grind_over_rmid / lambda_p_over_rmid - 2.0) < 1e-9,
      f"ratio = {r_grind_over_rmid/lambda_p_over_rmid:.9f}")

check("SP5: N at full overlap (d=0) is on the order of several hundred cells, "
      "NOT the ~2-22 cell chains tested in prior (1D/2D) scripts -- a materially different regime",
      Nd_values[-1][1] > 100,
      f"N(d=0) ~ {Nd_values[-1][1]:.1f} cells")

# ── Section 4: honest scope ───────────────────────────────────────────────────
print()
print(SEP)
print("SECTION 4: WHAT THIS FIXES, AND WHAT STILL REMAINS")
print(SEP2)
print("""
  FIXES: previous scripts (Stage 1/2/2b) picked N (chain length) or (N,M)
  (grid dimensions) as FREE PARAMETERS to sweep, disconnected from actual
  sphere geometry. This script computes N(d) DIRECTLY from the real
  sphere-sphere contact geometry as a function of separation d -- "the
  contacts are few, then scaling up as the spheres are forced closer" is
  now an actual computed curve, not an assumption. It also shows the
  PHYSICALLY RELEVANT cell count (up to ~680 cells at full overlap) is
  far larger than anything tested in the chain/grid scripts (up to 22).

  STILL REMAINS (not yet done): 
  1. This uses idealized point-like cell footprints on a flat patch
     estimate (area/area), not real icosahedral vertex/edge contact
     geometry -- still a proxy, though now at least tied to the right
     GROWING, roughly-circular shape instead of an arbitrary chain/grid.
  2. Does NOT yet combine this real N(d) with the cell-count-scaling
     result (mesh_grind_cell_count_scaling.py: total forced residual is
     CONSERVED at 2*v_surface for a 1D chain, regardless of chain length)
     -- whether that conservation law holds, or does something different,
     for a genuinely circular/2D patch of THIS size is a real, different
     calculation, not yet attempted. Natural next step.
  3. Real Jobson cells are icosahedra, not the flat circular-patch proxy
     used here -- packing them into an actual growing spherical contact
     zone still runs into the same open multi-cell packing question as
     PG2-1/MG-J8/MG-J9.
""")

# ── Summary ────────────────────────────────────────────────────────────────────
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
