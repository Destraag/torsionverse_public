"""
crys1_tool2_acceptance.py — [crys1] Tool 2: Cut-and-project acceptance window
for icosahedral quasicrystal and its connection to Gap 1.

[crys1] TOOL SERIES:
  Tool 1: WZW SU(2)_2 correlator  — COMPLETE (h(j=1)=1/2 marginal)
  Tool 2: Cut-and-project acceptance window  — THIS SCRIPT
  Tool 3: Phason gap  — COMPLETE (phason ratio ≠ C_geo)
  Tool 4: 6D Brillouin zone / EW spectrum  — PENDING

HYPOTHESIS (from gap1_common_factor.py, 2026-08-18):
  The electron torus knot is a path on a (1,2) torus.  In the quasicrystal
  picture the underlying medium is icosahedral: grains are vertices of a
  Z^6-projected icosahedral quasicrystal.  The PHYSICAL SPACE projection
  gives n_EM.  But the PERPENDICULAR SPACE component (phason winding) is
  ignored in the current model.

  If the perpendicular-space component of the winding contributes:
    delta_n_perp = n_exact * |v_perp|/|v_par|
  then the acceptance window |W| determines |v_perp| and hence delta_n.

PARTS:
  I   — Z^6 icosahedral embedding: projection matrices E_par, E_perp
  II  — Standard acceptance window (triacontahedron): size and volume
  III — Vertex density: how |W| determines grain density and N_lock
  IV  — Perpendicular-space winding: phason contribution to n
  V   — Required |W| shift to close Gap 1: δ|W|/|W| vs known factors
  VI  — Acceptance window width from Penrose scaling: does w = 2/(1+phi)?
  VII — Connection to f_frac and the four corroboration channels
  VIII — VERDICT

Run: python analysis/alpha/crys1_tool2_acceptance.py
Theory: alpha_theory.txt Part 0l — extend section 0l.9
"""

import math
import numpy as np

pi    = math.pi
sqrt5 = math.sqrt(5)
PHI   = (1 + sqrt5) / 2

# ── LOCKED CONSTANTS ─────────────────────────────────────────────────────────
alpha    = 7.2973525693e-3
eps_L5   = 3 / (8 * pi)
gj5      = 1 - math.cos(pi / 5)
Rs       = sqrt5 / (4 * pi)
Q        = 4 * pi**2 / PHI
dn       = 0.16857744391041
d2n      = 0.19763679211711
I_el     = 0.15614610339308
C_geo    = 10.33418281379304
delta_n  = 2.24745624e-6
delta_eps= delta_n / dn
f_frac   = delta_eps / eps_L5
n_EM     = 2.01868734358082
n_exact  = n_EM + delta_n
N_lock   = 2 * pi / (alpha * PHI)

SEP = '=' * 72


# ══════════════════════════════════════════════════════════════════════════════
#  PART I — Z^6 ICOSAHEDRAL EMBEDDING
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART I — Z^6 ICOSAHEDRAL EMBEDDING: PROJECTION MATRICES")
print(SEP)
print()
print("  The icosahedral quasicrystal is obtained by cut-and-project from Z^6.")
print("  6 basis vectors span the 6 five-fold axes of the icosahedron.")
print()
print("  For the (1,2) torus knot, we need a 2D cut from a higher-dimensional")
print("  lattice.  The relevant embedding is Z^4 (Hopf torus) or Z^6 (full I_h).")
print()

# ── 2D case first: Penrose tiling from Z^4 (the Ammann-Beenker / Z4 to R2) ──
# Actually for the (1,2) torus in 3D icosahedral: use Z^6 -> R^3
# 6 basis vectors (standard icosahedral quasicrystal):
# e_k = (cos(2*pi*k/5), sin(2*pi*k/5), zeta) for k=0..4, e_5 = (0,0,1)
# But different conventions exist. Use the Duneau-Katz convention:

# E_par: 3x6 matrix projecting Z^6 -> R^3 physical
# E_perp: 3x6 matrix projecting Z^6 -> R^3 perpendicular
zeta = 1 / sqrt5

# 5 five-fold axes (plus one vertical)
angles = [2 * pi * k / 5 for k in range(5)]
E_par  = np.zeros((3, 6))
E_perp = np.zeros((3, 6))

for k in range(5):
    a = angles[k]
    E_par[0, k]  =  math.cos(a)
    E_par[1, k]  =  math.sin(a)
    E_par[2, k]  =  zeta
    # Perpendicular: rotate by pi/5 and flip z
    E_perp[0, k] =  math.cos(a + pi/5)
    E_perp[1, k] =  math.sin(a + pi/5)
    E_perp[2, k] = -zeta

# 6th basis vector: z-axis (north pole)
E_par[2, 5]  = 1.0
E_perp[2, 5] = -1.0

# Normalization factor for orthonormality of E_par rows
# <e_i_par, e_j_par> = (5/2)*delta_ij ... rescale by sqrt(2/5)
norm_par  = math.sqrt(2 / 5)    # standard normalization
norm_perp = math.sqrt(2 / 5)

E_par_n  = E_par  * norm_par
E_perp_n = E_perp * norm_perp

# Check: E_par_n @ E_par_n.T should be proportional to identity
M_check = E_par_n @ E_par_n.T
print("  Projection check: E_par_n @ E_par_n.T (should be identity * const):")
for i in range(3):
    row = "    " + "  ".join(f"{M_check[i,j]:+.6f}" for j in range(3))
    print(row)
print()

# Perp check
M_perp = E_perp_n @ E_perp_n.T
print("  E_perp_n @ E_perp_n.T:")
for i in range(3):
    row = "    " + "  ".join(f"{M_perp[i,j]:+.6f}" for j in range(3))
    print(row)
print()

# Mixed: should be zero
M_mix = E_par_n @ E_perp_n.T
print("  E_par_n @ E_perp_n.T (should be ~0):")
for i in range(3):
    row = "    " + "  ".join(f"{M_mix[i,j]:+.6f}" for j in range(3))
    print(row)
print()

# Effective projection scaling
par_scale  = math.sqrt(np.trace(E_par_n @ E_par_n.T) / 3)
perp_scale = math.sqrt(np.trace(E_perp_n @ E_perp_n.T) / 3)
print(f"  Par  projection RMS scale: {par_scale:.8f}")
print(f"  Perp projection RMS scale: {perp_scale:.8f}")
print(f"  Ratio perp/par:            {perp_scale/par_scale:.8f}")
print(f"  Expected (= 1 for isometric): {1.0}")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART II — ACCEPTANCE WINDOW GEOMETRY
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART II — ACCEPTANCE WINDOW GEOMETRY FOR ICOSAHEDRAL QUASICRYSTAL")
print(SEP)
print()
print("  The acceptance window W ⊂ R^3_perp selects which Z^6 lattice points")
print("  project into physical space.  For the canonical icosahedral QC:")
print("    W = triacontahedron (30 faces, dual of icosidodecahedron)")
print("    OR: W = icosahedron (12 faces, depending on convention)")
print()

# Standard icosahedron vertices in perpendicular space
# Circumradius of regular icosahedron with edge a: R = a*sin(2*pi/5)
# Standard: unit edge a=1, R_ic = sqrt(1 + phi^2) / sqrt(2) ... or:
a_ic = 1.0   # unit edge length
R_ic = a_ic * math.sqrt(1 + PHI**2) / math.sqrt(2)  # circumradius formula? Let me recalculate
# Actually for regular icosahedron: circumradius = a * phi * sqrt(3/5)  ... various formulas
# Circumradius of icosahedron with edge a = a * sqrt(1 + phi^2) / sqrt(2)?
# More carefully: R = a/4 * sqrt(10 + 2*sqrt5)
R_ic_correct = a_ic / 4 * math.sqrt(10 + 2*sqrt5)
print(f"  Regular icosahedron, unit edge a = {a_ic}:")
print(f"    Circumradius R = a/4 * sqrt(10+2*sqrt5) = {R_ic_correct:.10f}")
# Check: R / a = phi * sin(pi/5) / cos(pi/10)... simpler:
# R = a * sin(2*pi/5) = a * sqrt(10+2*sqrt5)/4
R_check = a_ic * math.sin(2*pi/5)
print(f"    Alternate: a * sin(2*pi/5) = {R_check:.10f}  (should match)")
print()

# Volume of icosahedron with edge a:
V_ic = a_ic**3 * 5 * (3 + sqrt5) / 12
print(f"    Volume V_ic = 5(3+sqrt5)/12 * a^3 = {V_ic:.10f}")
print()

# The standard acceptance window for icosahedral QC from Z^6:
# In the Elser-Henley convention, W is the unit triacontahedron with
# volume V_W = (8/3)*phi^3 (in units where Z^6 unit cell has vol = 1)
# But more commonly: W = rhombic triacontahedron with edge a_W = 1/phi^2
a_W_std = 1 / PHI**2   # = 2 - phi = 1/phi^2
print("  Standard acceptance window (rhombic triacontahedron):")
print(f"    Characteristic edge a_W = 1/phi^2 = {a_W_std:.10f}")
print(f"    = 2-phi = phi^(-2)       = {2-PHI:.10f}  (check)")
print()

# Alternative: Elser's window a_W = 2/phi = 2*(phi-1) = 2/phi
a_W_elser = 2 / PHI
print(f"    Elser convention: a_W = 2/phi = {a_W_elser:.10f}")
print(f"    = 2*(phi-1)                   = {2*(PHI-1):.10f}")
print()

# Simple Penrose 2D analog:
# 2D acceptance window for Penrose tiling (Z^4 -> R^2):
# Regular pentagon with diagonal d = phi (edge = 1)
# Width: w = 2*sin(pi/5) * phi = phi * 2*sin(pi/5)
w_penrose = 2 * math.sin(pi/5) * PHI  # Not standard; standard is different
# Actually: for Penrose from Z^4, acceptance window in R^2_perp is a pentagon
# with edge length 1/phi (in lattice units), giving
# Area_W = (5 * phi^2) / (4 * tan(pi/5)) ... normalised
# The important ratio: span of window / span of lattice = 1/(1+phi) = 1/phi^2
w_std_2D = 2 / (1 + PHI)    # = 2*phi^(-2), standard Penrose window width
print(f"  2D Penrose analog acceptance window width:")
print(f"    w = 2/(1+phi) = 1/phi^2 * 2 = {w_std_2D:.10f}")
print(f"    = 2*phi^(-2)               = {2/PHI**2:.10f}")
print()

# The characteristic length scale of the window in the Z^6 embedding
# relative to the lattice spacing:
# For icosahedral QC the relevant scale is:
#   λ_W = (volume of window)^(1/3) / (Z^6 unit cell length in perp space)
# In standard normalization where Z^6 unit cell has volume 1:
lambda_W = a_W_std   # characteristic scale
print(f"  Window characteristic scale λ_W = {lambda_W:.10f}")
print(f"  Physical: λ_W sets the density of accepted lattice points.")
print(f"  Density ρ ~ (λ_W)^3 / (perp-space unit cell volume)")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART III — VERTEX DENSITY AND N_LOCK
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART III — VERTEX DENSITY AND CONNECTION TO N_LOCK")
print(SEP)
print()
print("  N_lock = 2*pi / (alpha * phi) = grain count the electron encounters")
print(f"         = {N_lock:.8f}")
print()
print("  In the quasicrystal picture:")
print("    Grain density along the torus knot ∝ (number of vertices per unit length)")
print("    = ρ_QC * d_bar (average nearest-neighbor distance)")
print()

# In 3D icosahedral QC the vertex density (number per unit volume) is:
# rho_QC = (V_W / V_perp_cell) * (1 / V_par_cell)
# Where:
#   V_W        = acceptance window volume (triacontahedron)
#   V_perp_cell = perpendicular unit cell volume
#   V_par_cell  = physical unit cell volume (diverges for QC, use per-unit-length)

# For our 1D problem (winding along torus knot), effective density is:
# rho_1D = (window cross-section area in perp 2D) / (Z^6 unit cell volume in perp 2D)

# In 2D Penrose:
# rho_2D = A_W_perp / A_cell_perp = (area of acceptance pentagon) / (Z^4 perp cell)
# = (5*phi^2)/(4*tan(pi/5)) / (some normalization)
# For Penrose with unit tile edge: vertex density = 2/(1+sqrt5) per unit area
# = 1/(phi+1) per unit area = phi^(-2)

rho_penrose_std = 1 / PHI**2
print(f"  2D Penrose vertex density (standard): {rho_penrose_std:.10f} per tile area")
print(f"  = phi^(-2) = {1/PHI**2:.10f}")
print()

# For 3D icosahedral QC the vertex density in Z^6 / Z^6 unit cell:
# rho_3D = (V_triacontahedron) / (V_par_cell * V_perp_cell_normalized)
# Numerically from simulations: rho_3D ≈ phi^(-3) / (1 + phi^(-2))
# = phi^(-3) * phi^2 / (phi^2 + 1) = phi^(-1) / (phi^2 + 1)
# For standard normalization: phi^(-3) ≈ 0.23607

rho_3D_approx = 1 / PHI**3
rho_3D_corr   = 1 / PHI**3 * (PHI**2 + 1) / PHI
print(f"  3D icosahedral QC density estimates:")
print(f"    phi^(-3) = {rho_3D_approx:.8f}")
print(f"    phi^(-3) * correction = {rho_3D_corr:.8f}")
print()

# N_lock is the number of grains along the torus knot of length L_knot:
# L_knot ~ 2*pi*R2 = (2*pi)^2 (since R2 = 2*pi)
# N_lock = rho_1D * L_knot
# So rho_1D = N_lock / L_knot
L_knot = (2 * pi)**2   # approximate arc length of (1,2) knot at eps~0
rho_1D_inferred = N_lock / L_knot
print(f"  Torus knot arc length (leading): L_knot ~ (2*pi)^2 = {L_knot:.8f}")
print(f"  Inferred 1D grain density:       rho_1D = {rho_1D_inferred:.8f}")
print(f"  Nearest simple form:             alpha*phi/(2*pi)^3 = {alpha*PHI/(2*pi)**3:.8f}")
print(f"  Check: 1/N_lock * 1/L_knot = {1/(N_lock*L_knot):.8f}")
print()

# If window shifts by δλ_W:
# rho ~ λ_W^d (d = 2 for 2D window cross-section)
# δrho/rho = d * δλ_W/λ_W
# δN_lock/N_lock = δrho/rho = d * δλ_W/λ_W
# δN_lock/N_lock = f_frac = 1.117e-4
# => δλ_W/λ_W = f_frac / d
d_window = 2   # effective dimension of window cross-section (2D window in 3D space)
dL_W_over_L_W = f_frac / d_window
print(f"  Window shift needed to produce δN_lock/N_lock = f_frac:")
print(f"    d (window dimension): {d_window}")
print(f"    δλ_W/λ_W = f_frac/d = {dL_W_over_L_W:.8e}")
print(f"    δλ_W     = λ_W * f_frac/d = {lambda_W * dL_W_over_L_W:.8e}")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART IV — PERPENDICULAR SPACE WINDING: PHASON CONTRIBUTION TO n
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART IV — PERPENDICULAR-SPACE WINDING: PHASON CONTRIBUTION")
print(SEP)
print()
print("  HYPOTHESIS: The torus knot path in 6D has both a physical-space")
print("  component (measured by n_EM) and a perpendicular-space component.")
print("  The phason winding contributes delta_n_phason = n_EM * |v_perp|/|v_par|.")
print()
print("  For a (1,2) torus knot: the Z^6 lattice vector is v = (1,0,0,0,0,2)")
print("  (in the icosahedral basis), with physical projection v_par and")
print("  perpendicular projection v_perp.")
print()

# For the (1,2) torus knot v = (m, n) with m=1 (poloidal) and n=2 (toroidal)
# In the 6D embedding this becomes a more complex vector.
# Let's compute the projection of a general icosahedral lattice vector onto
# par and perp spaces.

# Test with v = e_1 (first basis vector):
v1 = np.array([1, 0, 0, 0, 0, 0], dtype=float)
v1_par  = E_par_n @ v1
v1_perp = E_perp_n @ v1

print(f"  Projection of e_1 = (1,0,0,0,0,0):")
print(f"    |v_par|  = {np.linalg.norm(v1_par):.8f}")
print(f"    |v_perp| = {np.linalg.norm(v1_perp):.8f}")
print(f"    ratio    |v_perp|/|v_par| = {np.linalg.norm(v1_perp)/np.linalg.norm(v1_par):.8f}")
print()

# For the (1,2) winding vector in 6D:
# The (1,2) torus knot in the 5-fold layer: winding (1,2) maps to
# v = sum of 2 contributions: poloidal 1-turn + toroidal 2-turns
# In icosahedral basis, a natural vector with 1:2 ratio:
v12 = np.array([1, 2, 0, 0, 0, 0], dtype=float)
v12_par  = E_par_n @ v12
v12_perp = E_perp_n @ v12
r12 = np.linalg.norm(v12_perp) / np.linalg.norm(v12_par)

print(f"  Projection of v_(1,2) = (1,2,0,0,0,0):")
print(f"    |v_par|  = {np.linalg.norm(v12_par):.8f}")
print(f"    |v_perp| = {np.linalg.norm(v12_perp):.8f}")
print(f"    ratio    |v_perp|/|v_par| = {r12:.8f}")
print()

# Scan over all low-index (1,q) winding vectors:
print("  Scan over (1,q,0,0,0,0) vectors in Z^6, q = 1..6:")
print(f"  {'v':<20}  {'|par|':>8}  {'|perp|':>8}  {'perp/par':>10}")
for q in range(1, 7):
    v = np.zeros(6); v[0] = 1; v[1] = q
    vp = E_par_n @ v
    vq = E_perp_n @ v
    print(f"  (1,{q},0,0,0,0)          {np.linalg.norm(vp):>8.6f}  {np.linalg.norm(vq):>8.6f}  {np.linalg.norm(vq)/np.linalg.norm(vp):>10.8f}")
print()

# The phason contribution to delta_n:
# delta_n_phason = n_EM * (|v_perp|/|v_par|)^2  (from projection area ratio)
# OR: delta_n_phason = n_EM * |v_perp|/|v_par|  (from projection length ratio)
print("  Phason contribution to delta_n for (1,2) vector:")
r_12 = r12
delta_n_phason_sq  = n_EM * r_12**2
delta_n_phason_lin = n_EM * r_12
print(f"    n_EM * (perp/par)    = {delta_n_phason_lin:.8e}  vs delta_n = {delta_n:.8e}")
print(f"    n_EM * (perp/par)^2  = {delta_n_phason_sq:.8e}  vs delta_n = {delta_n:.8e}")
print(f"    (perp/par) ratio for linear match: {delta_n/n_EM:.8e}")
print(f"    (perp/par) needed for actual ratio: {r_12:.8e}")
print(f"    linear match error:  {(delta_n_phason_lin/delta_n - 1)*100:+.4f}%")
print()

# How many orders of magnitude is r_12 from what's needed?
needed = delta_n / n_EM
print(f"  SCALE CHECK:")
print(f"    Needed ratio:  {needed:.6e}")
print(f"    Actual r_12:   {r_12:.6e}")
print(f"    Ratio of ratios: {r_12/needed:.4f}  (how many times too large?)")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART V — REQUIRED WINDOW SHIFT AND KNOWN FACTORS
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART V — REQUIRED WINDOW SHIFT δλ_W AND KNOWN FACTORS")
print(SEP)
print()
print("  If the acceptance window shifts by δλ_W, the winding number shifts by:")
print("  δn = n_EM * (2 * δλ_W / λ_W)   (area-scaling in 2D window cross-section)")
print("  => δλ_W / λ_W = delta_n / (2 * n_EM) = f_frac / 2")
print()

need_ratio = delta_n / (2 * n_EM)
need_ratio_f = f_frac / 2
print(f"  δλ_W/λ_W needed = {need_ratio:.6e}")
print(f"  = f_frac/2      = {need_ratio_f:.6e}")
print(f"  Check: {(need_ratio/need_ratio_f):.8f}")
print()

print("  Candidate known factors that could produce this shift:")
candidates = [
    ("alpha^2/2",             alpha**2/2),
    ("alpha^2/eps_L5",        alpha**2/eps_L5),
    ("alpha^2*Q/2",           alpha**2*Q/2),
    ("gj5*alpha^2/eps_L5",    gj5*alpha**2/eps_L5),
    ("(2-PHI)*eps_L5",        (2-PHI)*eps_L5),
    ("1/N_lock",              1/N_lock),
    ("alpha*eps_L5",          alpha*eps_L5),
    ("delta_pi/pi",           3*(3+sqrt5)/5/pi - 1),
    ("excess5^2",             (2*eps_L5 - gj5)**2),
    ("eps_L5^2*alpha^2",      eps_L5**2*alpha**2),
    ("PHI^(-2)/Q",            1/PHI**2/Q),
    ("alpha/(2*Q)",           alpha/(2*Q)),
    ("alpha^2*d_dn",          alpha**2 * (d2n/dn)),
    ("Rs*alpha/2",            Rs*alpha/2),
    ("alpha^3/gj5",           alpha**3/gj5),
    ("w_std = 2/(1+PHI)",     2/(1+PHI)),
    ("w_std * eps_L5",        2/(1+PHI)*eps_L5),
    ("w_std * alpha",         2/(1+PHI)*alpha),
    ("w_std * f_frac",        2/(1+PHI)*f_frac),
    ("(1/PHI^2)/N_lock",      1/(PHI**2*N_lock)),
]

print(f"  {'Candidate':<36}  {'Value':>12}  {'/ (f_frac/2)':>12}  {'err from 1':>10}")
print("  " + "-" * 76)
for lbl, val in sorted(candidates, key=lambda x: abs(x[1]/need_ratio - 1)):
    ratio = val / need_ratio
    err = (ratio - 1) * 100
    flag = "  ***" if abs(err) < 0.5 else ("  **" if abs(err) < 5 else ("  *" if abs(err) < 15 else ""))
    print(f"  {lbl:<36}  {val:>12.6e}  {ratio:>12.6f}  {err:>+9.3f}%{flag}")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART VI — PENROSE SCALING AND ACCEPTANCE WINDOW WIDTH
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART VI — PENROSE SCALING: EXACT ACCEPTANCE WINDOW WIDTH")
print(SEP)
print()
print("  The 2D Penrose acceptance window (for Z^4 -> R^2 projection) is a")
print("  regular pentagon with CHARACTERISTIC WIDTH w = 2*sin(pi/5) / cos(pi/10).")
print()

# For Penrose tiling from Z^4:
# The perpendicular-space acceptance window is a regular pentagon
# with edge length a_pent = 1/PHI (in Z^4 lattice units)
a_pent_W = 1 / PHI
w_A = 2 * math.sin(pi/5) * a_pent_W    # inradius-based width
w_B = 2 * a_pent_W * math.sin(pi/5) / math.cos(pi/10)   # full width

print(f"  Pentagon edge a_W = 1/phi = {a_pent_W:.10f}")
print(f"  Pentagon width (2*r_in) = 2*sin(pi/5)*a_W = {w_A:.10f}")
print(f"  Pentagon full width     = {w_B:.10f}")
print()

# Standard width formula from literature:
# w = sqrt((5-sqrt5)/5) * phi^(-1) ... or simply w = 2/(1+phi)
w_std = 2 / (1 + PHI)
print(f"  Standard: w_std = 2/(1+phi) = {w_std:.10f}")
print(f"  phi^(-2) * 2    = {2/PHI**2:.10f}  (= w_std)")
print()

# The DEVIATION of the torus wave epsilon from the QC window scale:
print("  Connection to eps_L5:")
print(f"    eps_L5 = 3/(8*pi) = {eps_L5:.10f}")
print(f"    w_std              = {w_std:.10f}")
print(f"    eps_L5 / w_std     = {eps_L5/w_std:.10f}")
print(f"    eps_L5 * w_std     = {eps_L5*w_std:.10f}")
print()

# Higher-dimensional: the 3D icosahedral window has circumradius
# R_W = a_W * sin(2*pi/5) where a_W is the edge length
# For the canonical QC: a_W = 1/phi^2
R_W = a_W_std * math.sin(2*pi/5)
V_W_icos = 5*(3+sqrt5)/12 * a_W_std**3
print(f"  3D icosahedral acceptance window (edge = 1/phi^2 = {a_W_std:.6f}):")
print(f"    Circumradius R_W = {R_W:.10f}")
print(f"    Volume V_W       = {V_W_icos:.10f}")
print()

# How does a shift δa_W/a_W = f_frac affect V_W?
# V ~ a^3, so δV/V = 3*δa/a = 3*f_frac => δa/a = f_frac/3
daw_over_aw = f_frac / 3
print(f"  Window edge shift for δV_W/V_W = f_frac:")
print(f"    δa_W/a_W = f_frac/3 = {daw_over_aw:.8e}")
print(f"    δa_W     = a_W * f_frac/3 = {a_W_std*daw_over_aw:.8e}")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART VII — CONNECTION TO f_frac AND FOUR CHANNELS
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART VII — QC WINDOW SCALES AND THE FOUR CORROBORATION CHANNELS")
print(SEP)
print()

# Corroboration channel errors
err_A =  PHI / I_el / C_geo - 1
err_B =  (PHI**2 / sqrt5) / (d2n/dn) - 1
err_C_val = 0.232689e-2  # from crys1_tool2 computation not available; use stored
err_D =  f_frac

channels = [('A', err_A), ('B', err_B), ('D=f_frac', err_D)]

print("  QC scales vs channel residuals:")
print()
qc_scales = [
    ("w_std = 2/(1+phi)",           w_std),
    ("1/phi^2",                      1/PHI**2),
    ("w_std^2",                      w_std**2),
    ("w_std^3",                      w_std**3),
    ("w_std * eps_L5",               w_std * eps_L5),
    ("w_std^2 * eps_L5",             w_std**2 * eps_L5),
    ("R_W",                          R_W),
    ("R_W * alpha",                  R_W * alpha),
    ("R_W * f_frac",                 R_W * f_frac),
    ("V_W_icos",                     V_W_icos),
    ("V_W_icos * alpha",             V_W_icos * alpha),
    ("a_W_std * alpha^2 / eps_L5",   a_W_std * alpha**2 / eps_L5),
    ("a_W_std * alpha^2",            a_W_std * alpha**2),
    ("daw_over_aw",                   daw_over_aw),
    ("a_W_std * f_frac",             a_W_std * f_frac),
]

for qlbl, qval in qc_scales:
    for clabel, cval in channels:
        err = (qval / abs(cval) - 1) * 100
        if abs(err) < 10:
            print(f"  {qlbl:<36}  ({qval:.4e}) vs |err_{clabel}| ({abs(cval):.4e})  err={err:+.3f}%  ***")

print()

# Key dimensionless ratio: N_lock * f_frac
# = 2*pi/(alpha*phi) * f_frac = 2*pi*delta_eps/(alpha*phi*eps_L5)
N_f = N_lock * f_frac
print(f"  N_lock * f_frac = {N_f:.8f}")
print(f"  Compare to: phi/(2*pi*alpha) * f_frac... that is N_f itself = {N_f:.4f}")
print(f"  Compare to: 1/(4*pi) = {1/(4*pi):.8f}")
print(f"  Compare to: 1/(4*pi^2) = {1/(4*pi**2):.8f}")
print(f"  Compare to: phi/Q = {PHI/Q:.8f}")
print(f"  Compare to: 1/Q = {1/Q:.8f}")
print()

# Effective quasicrystal correction to n_EM from window boundary:
# When window edge is at R_W from origin in perp space,
# lattice points at distance d_perp ≈ R_W are MARGINAL (near-miss).
# Their fractional contribution to n_EM is ~ (d_perp/R_W) * (something).
# This is a SOFT correction, not a step function.
print("  Soft-boundary correction model:")
print("  Near-boundary lattice points within distance δ of window edge contribute")
print("  fractionally with weight ~ (1 - d_perp/R_W).")
print()
print(f"  The fraction of grains near the boundary is:")
frac_boundary = 3 * a_W_std / R_W   # surface/volume ~ 3*thickness/radius
print(f"    f_boundary = 3 * δ / R_W  (for shell of thickness δ = a_W_std)")
print(f"    = {frac_boundary:.8f}")
print(f"  Compare: f_frac = {f_frac:.8f}")
print(f"  Ratio:   {frac_boundary/f_frac:.4f}")
print()

# What δ would produce exactly f_frac?
delta_needed = f_frac * R_W / 3
print(f"  Shell thickness δ needed: {delta_needed:.8e}")
print(f"  vs a_W_std = {a_W_std:.8e}")
print(f"  δ / a_W_std = {delta_needed/a_W_std:.8e}")
print(f"  This ratio = f_frac * R_W / (3 * a_W_std^2) = {f_frac*R_W/(3*a_W_std**2):.8e}")
print()


# ══════════════════════════════════════════════════════════════════════════════
#  PART VIII — VERDICT
# ══════════════════════════════════════════════════════════════════════════════
print(SEP)
print("PART VIII — VERDICT: DOES CUT-AND-PROJECT OFFER A ROUTE TO GAP 1?")
print(SEP)
print()
print("  QUESTION 1: Is the phason winding (perp-space component) = delta_n?")
print(f"  n_EM * (perp/par) for (1,2,0,0,0,0) = {delta_n_phason_lin:.4e}")
print(f"  delta_n actual                        = {delta_n:.4e}")
print(f"  Ratio (how many times too large?):     {delta_n_phason_lin/delta_n:.2f}")
print()
if delta_n_phason_lin / delta_n > 10:
    print("  VERDICT 1: RULED OUT by scale mismatch.")
    print(f"  The phason winding contribution is {delta_n_phason_lin/delta_n:.0f}x too large.")
    print("  The geometric projection ratio is O(1), not O(f_frac).")
    print("  The perpendicular-space component is comparable to the physical-space")
    print("  component — not a small correction.")
else:
    print("  VERDICT 1: PLAUSIBLE — scale is consistent.")
print()
print("  QUESTION 2: Can window edge shift δλ_W/λ_W produce δn = delta_n?")
print(f"  Required: δλ_W/λ_W = f_frac/2 = {need_ratio:.4e}")
print("  Nearest candidate from Part V search: see above.")
print()
print("  QUESTION 3: Is the window size related to any fundamental scale?")
print(f"  a_W = 1/phi^2 = {a_W_std:.8f}")
print(f"  eps_L5 / a_W = {eps_L5/a_W_std:.8f}  (how many window-widths is eps_L5?)")
print(f"  gj5 / a_W    = {gj5/a_W_std:.8f}  (how many window-widths is gj5?)")
print()
print("  QUESTION 4: Does the QC vertex density give N_lock = 532?")
print("  The QC density from Z^6 is determined by V_W and the lattice determinant.")
print("  V_W ~ a_W^3 = phi^{-6} is a FIXED TOPOLOGICAL quantity.")
print("  N_lock = 2*pi/(alpha*phi) links grain count to ALPHA, not to V_W.")
print("  Alpha and the QC window are NOT directly related in this model.")
print()
print("  OVERALL VERDICT:")
print()
print("  The cut-and-project framework identifies the structure of the medium")
print("  (icosahedral quasicrystal) but does NOT directly give a mechanism for")
print("  Gap 1 unless the acceptance window size itself encodes alpha.")
print()
print("  The phason winding is O(n_EM), not O(f_frac * n_EM) — wrong scale.")
print("  The window-shift route requires a physical reason why delta_a_W/a_W = f_frac/d.")
print()
print("  NEW LEAD (not previously explored):")
print("  The acceptance window is FIXED at a_W = 1/phi^2 for the canonical QC.")
print("  But the ACTUAL grain medium may have a MODIFIED window due to the")
print("  electromagnetic wave strain: eps_L5 deforms the quasicrystal.")
print("  The phason strain K from eps_L5 shifts a_W by:")
print("    δa_W = K_eff * eps_L5 * a_W")
print(f"  If K_eff = f_frac / (2 * eps_L5) = {f_frac/(2*eps_L5):.6e},")
print(f"  that gives δa_W/a_W = f_frac/2 exactly.")
print(f"  K_eff = {f_frac/(2*eps_L5):.8e}")
print(f"  vs alpha^2/(2*phi) = {alpha**2/(2*PHI):.8e}  (err = {(alpha**2/(2*PHI)/(f_frac/(2*eps_L5))-1)*100:+.4f}%)")
print()
print("  KEY RESULT: K_eff * eps_L5 = f_frac/2 is equivalent to K_eff = alpha^2/(2*phi)")
print("  if alpha^2/(2*phi) = f_frac/(2*eps_L5).")
Keff = f_frac / (2 * eps_L5)
Keff_alpha = alpha**2 / (2 * PHI)
print(f"  CHECK: alpha^2/(2*phi) = {Keff_alpha:.8e}")
print(f"         f_frac/(2*eps_L5) = {Keff:.8e}")
print(f"         ratio              = {Keff_alpha/Keff:.8f}")
err_K = (Keff_alpha/Keff - 1) * 100
print(f"         error              = {err_K:+.4f}%")
print()
if abs(err_K) < 1.0:
    print("  *** NEAR HIT: alpha^2/(2*phi) ≈ f_frac/(2*eps_L5) to within 1% ***")
    print("  This means: if the phason strain K = alpha^2/(2*phi) from EM coupling,")
    print("  it produces exactly the required window shift to explain Gap 1.")
    print("  NEXT STEP: derive K_eff = alpha^2/(2*phi) from first principles.")
    print("  This is the EM-quasicrystal coupling — a NEW mechanism not yet explored.")
else:
    print(f"  MISS: {err_K:+.2f}% off — not a clean hit.")
    print("  Proceed to [crys1] Tool 4 (6D Brillouin zone / EW spectrum).")
print()
print("Script: analysis/alpha/crys1_tool2_acceptance.py")
print("Theory: alpha_theory.txt Part 0l — new section 0l.9")
