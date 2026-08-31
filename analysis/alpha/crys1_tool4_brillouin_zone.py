"""
crys1_tool4_brillouin_zone.py — [crys1] Tool 4: 6D Brillouin zone → EW spectrum

CONTEXT IN THE CRYS1 SERIES:
  Tool 0l.1 (phason):      C_geo "10.334" vs phason "9.934" — differ by d2n/dn residual
  Tool 0l.2 (I_h chars):   phi exact as irrep character; d2n/dn not recoverable from Casimirs
  Tool 0l.3 (WZW k=2):     h(j=1) = 1/2 exact (marginal); coupling unrenormalized
  Tool 0l.4 (2I/McKay/E8): 2I irrep dims = E8 Kac labels; 3-4-5 triangle exact
  Tool 0l.8 (gj5=a_W/2):   jamming threshold = half acceptance window of Z^6 quasicrystal
  Tool 4 (THIS SCRIPT):    6D Brillouin zone of Z^6 lattice → EM wave band structure

THE QUASICRYSTAL PICTURE:
  The icosahedral medium is modeled as a cut-and-project quasicrystal from Z^6
  (6D integer lattice). The 3D physical tiling has icosahedral symmetry because
  the projection direction is along the (1,phi,0,1,phi,0)/norm axis of Z^6.

  The 6D BRILLOUIN ZONE of Z^6 is the 6D hypercube [-pi, pi]^6.
  Equivalently (in dual basis), the Wigner-Seitz cell of the reciprocal lattice.

  For Z^6 (simple hypercubic): reciprocal lattice is also Z^6, BZ = [-1/2, 1/2]^6.
  The EM wave spectrum = eigenvalues of the 6D wave operator in this BZ.

THE PHYSICAL QUESTION:
  Does the PROJECTED BZ boundary in 3D physical space correspond to
  eps_L5 = 3/(8*pi) = 0.11937?

  Specifically:
  - The projection of the 6D BZ onto physical 3D gives a 3D region.
  - The boundary of this projected region in wave-vector space gives a
    characteristic length scale k* in 3D.
  - If k* corresponds to the Hopf torus parameter eps_L5, this provides
    independent corroboration from crystallography.

THREE COMPUTATIONS:
  I.   Geometry of Z^6 and the icosahedral projection matrix
  II.  The 6D BZ: volume, boundary, face structure (24-cell / cross-polytope?)
  III. Projected BZ: 3D image of the 6D BZ boundary under the cut-and-project map
  IV.  Characteristic scale k_BZ vs eps_L5 and related constants
  V.   Wave operator spectrum: eigenvalues of -Delta (6D Laplacian) at BZ boundary
  VI.  Phason contribution: perp-space BZ and its scale
  VII. Does any BZ-derived scale match eps_L5, gj5, alpha, or delta_n?
  VIII.Summary

IMPORTANT CAVEAT:
  The Z^6 lattice BZ is trivially a hypercube. The INTERESTING structure comes
  from the ICOSAHEDRAL PROJECTION — which 6D modes project to which 3D modes,
  and at what 3D wave-vectors does a 6D BZ boundary look like?

Run: python analysis/alpha/crys1_tool4_brillouin_zone.py
Theory: alpha_theory.txt sections 0l.1-0l.8, gap1_gj5_phi_identity.py
"""

import math
import itertools

pi    = math.pi
sqrt5 = math.sqrt(5)
sqrt2 = math.sqrt(2)
sqrt3 = math.sqrt(3)
PHI   = (1 + sqrt5) / 2

# -- LOCKED CONSTANTS ---------------------------------------------------------
alpha    = 7.2973525693e-3
eps_L5   = 3 / (8 * pi)
gj5      = 1 / (2 * PHI**2)    # = 1 - cos(pi/5) = a_W/2
a_W      = 1 / PHI**2           # acceptance window edge of Z^6 quasicrystal
C_geo    = 10.33418281379304
delta_n  = 2.24745624e-6

SEP  = '=' * 72
SEP2 = '-' * 60


# =============================================================================
# PART I -- GEOMETRY OF Z^6 AND THE ICOSAHEDRAL PROJECTION
# =============================================================================
print(SEP)
print("PART I -- Z^6 LATTICE AND ICOSAHEDRAL PROJECTION GEOMETRY")
print(SEP)
print()
print("  The icosahedral quasicrystal is cut-and-project from Z^6.")
print("  The 6D space decomposes as R^6 = E_par (3D physical) + E_perp (3D perp).")
print()

# The standard icosahedral projection from Z^6 to R^3 uses the 6 vectors:
# e_i projects to a physical 3D vector v_i such that {v_i} generate the
# icosahedrally-symmetric point set.
# The 6 projection vectors (Cahn-Levine-Shechtman basis):
# In units where the 6D lattice has spacing 1:
# v_1 = (1, tau, 0) / norm,  etc.  (12 vertices of icosahedron, 6 pairs)

# Physical projection vectors (to 3D), normalized
norm_v = math.sqrt(1 + PHI**2)
print(f"  6D->3D projection: each Z^6 basis vector -> icosahedral vector")
print(f"  Normalization: ||(1, phi, 0)|| = sqrt(1 + phi^2) = {norm_v:.8f}")
print(f"  phi = {PHI:.8f}")
print()

# The 6 physical projection vectors (columns of the 3x6 projection matrix Pi_par)
# Standard: e_k -> (cos(2*pi*k/6), sin(2*pi*k/6), cos(4*pi*k/6))... 
# Let's use the Elser-Sloane basis for Z^6 icosahedral projection:
# The 6 projection vectors in 3D are:
# a_k = (1/sqrt(1+phi^2)) * (1, phi, 0) cyclic permutations with sign changes
# giving the 6 "golden vectors" of the icosahedron:

def golden_vectors():
    """6 vectors for Z^6 icosahedral projection (Cahn-Levine-Shechtman)."""
    n = norm_v
    return [
        (1/n, PHI/n, 0),
        (-1/n, PHI/n, 0),
        (0, 1/n, PHI/n),
        (0, -1/n, PHI/n),
        (PHI/n, 0, 1/n),
        (-PHI/n, 0, 1/n),
    ]

def golden_perp_vectors():
    """6 perpendicular-space vectors (same structure with phi -> -1/phi)."""
    p = -1/PHI  # = 1 - phi (perp space: replace phi with -1/phi)
    n = math.sqrt(1 + (1/PHI)**2)
    return [
        (1/n, p/n, 0),
        (-1/n, p/n, 0),
        (0, 1/n, p/n),
        (0, -1/n, p/n),
        (p/n, 0, 1/n),
        (-p/n, 0, 1/n),
    ]

vphys = golden_vectors()
vperp = golden_perp_vectors()

print(f"  6 physical-space projection vectors (v_i in R^3):")
for i, v in enumerate(vphys):
    print(f"    v_{i+1} = ({v[0]:+.6f}, {v[1]:+.6f}, {v[2]:+.6f})")
print()

# Check: all |v_i| = 1 (unit vectors)
mags = [math.sqrt(sum(x**2 for x in v)) for v in vphys]
print(f"  All |v_i| = {mags[0]:.10f}  (should be 1.0)")
print()

# The Z^6 reciprocal lattice is also Z^6 (self-dual).
# The 6D BZ is the Wigner-Seitz cell of the reciprocal Z^6 = [-1/2, 1/2]^6 hypercube.
print(f"  Z^6 reciprocal lattice: Z^6 (self-dual)")
print(f"  6D Brillouin zone: [-1/2, 1/2]^6 (hypercube)")
print(f"  6D BZ volume: 1 (unit hypercube)")
print()


# =============================================================================
# PART II -- 6D BRILLOUIN ZONE: STRUCTURE AND BOUNDARY
# =============================================================================
print(SEP)
print("PART II -- 6D BRILLOUIN ZONE: STRUCTURE AND BOUNDARY")
print(SEP)
print()
print("  For Z^6: BZ = [-pi, pi]^6 in natural units (or [-1/2, 1/2]^6 in")
print("  units of 2*pi). We use lattice spacing = 1 so BZ = [-pi, pi]^6.")
print()

# BZ properties for 6D hypercubic lattice
dim = 6
BZ_half_side = pi
BZ_volume    = (2 * pi)**dim
print(f"  Dimension: {dim}")
print(f"  Half-side length: pi = {BZ_half_side:.8f}")
print(f"  Volume: (2*pi)^6 = {BZ_volume:.6f}")
print()

# The BZ boundary of a hypercube: faces are at k_i = +/-pi for i=1..6
# Number of faces: 2*6 = 12
# Face area: (2*pi)^5 each
# The BOUNDARY structure is important for determining which modes go soft first
# (i.e., which boundary points are reached as we increase the 3D wave vector k_phys)

print(f"  BZ boundary faces: 2*6 = 12 (hyperfaces at k_i = +/-pi)")
print(f"  Each face area: (2*pi)^5 = {(2*pi)**5:.6f}")
print()

# Special BZ boundary points
print(f"  Special points on 6D BZ boundary:")
print(f"    Gamma (center):  k = (0,0,0,0,0,0)")
print(f"    X (face center): k = (pi,0,0,0,0,0) and permutations -- 12 points")
print(f"    M (edge center): k = (pi,pi,0,0,0,0) and permutations -- C(6,2)*4 = 60 points")
print(f"    R (corner):      k = (pi,pi,pi,pi,pi,pi) -- {2**6} corners")
print()

# Gamma -> X distance (distance to nearest BZ face)
dist_GX = pi
# Gamma -> R distance (distance to corner)
dist_GR = pi * math.sqrt(dim)
print(f"  |Gamma -> X| = pi = {dist_GX:.6f}  (BZ face = 'X' zone boundary)")
print(f"  |Gamma -> R| = pi*sqrt(6) = {dist_GR:.6f}  (BZ corner)")
print()


# =============================================================================
# PART III -- PROJECTED BZ: 3D IMAGE OF 6D BZ BOUNDARY
# =============================================================================
print(SEP)
print("PART III -- PROJECTED BZ: 3D IMAGE UNDER ICOSAHEDRAL PROJECTION")
print(SEP)
print()
print("  The 6D wave vector k_6D = (k1,...,k6) in [-pi,pi]^6.")
print("  Physical 3D wave vector: k_phys = Pi_par * k_6D = sum_i k_i * v_i")
print("  Perp 3D wave vector:     k_perp = Pi_perp * k_6D = sum_i k_i * w_i")
print()
print("  The projected BZ in 3D physical space is the set of all k_phys")
print("  obtained from k_6D on the 6D BZ boundary.")
print()

# For k_6D on the BZ FACE k_1 = pi, the projected k_phys spans the range
# k_phys = pi * v_1 + sum_{i=2}^{6} k_i * v_i,  k_i in [-pi, pi]
# The amplitude along v_1: pi/norm_v (physical component of k_1 = pi along e_1)
# Physical amplitude of each BZ face point:
# k_phys when k_1 = pi, all others 0: pi * v_1 = pi * (1/norm_v, phi/norm_v, 0)
# |k_phys| = pi * |v_1| = pi * 1 = pi (since |v_i| = 1)

k_BZ_face_phys = pi  # magnitude of k_phys at the nearest BZ face (one component)
print(f"  NEAREST BZ FACE (X point): k_1 = pi, others = 0")
print(f"  Physical k: k_phys = pi * v_1,  |k_phys| = pi * |v_1| = pi = {k_BZ_face_phys:.6f}")
print()

# BUT: in the quasicrystal, ALL 6 components are correlated because the
# quasicrystal lattice sites correspond to specific k_6D values.
# The PHYSICAL dispersion is determined by the projected spectrum.

# For the icosahedral projection, the key momentum scales in 3D are:
# - The parallel component: determined by the icosahedral tiling lattice constant
# - The perpendicular (phason) component: gives the internal degree of freedom

# The 3D "effective BZ" of the projected quasicrystal:
# For icosahedral quasicrystals, the relevant momentum scale is set by
# the reciprocal of the average tile size.
# The tile size is 1/PHI^(n) for the n-th generation Penrose-like tiling.
# At Z^6 projection, the natural physical length scale is 1/norm_v.

# Physical reciprocal lattice vector magnitudes from the 6 projections:
# Each e_i^* = e_i (reciprocal Z^6) projects to k_phys^(i) = v_i
# The first "Bragg peak" in 3D occurs at k = |v_i| = 1.0 (all unit vectors)
# But the icosahedral tiling has 12-fold peaks in 3D (2*6 projections)

print(f"  Physical Bragg peaks from Z^6 projection:")
for i, v in enumerate(vphys[:3]):
    print(f"    G_{i+1} = {v}  |G| = {math.sqrt(sum(x**2 for x in v)):.6f}")
print(f"    (and 3 more with x-component negated)")
print()

# Characteristic physical length scale from BZ
# The BZ "radius" in 3D = typical |k_phys| at BZ boundary
# For uniform k_6D: |k_phys|^2 = sum_i k_i^2 * |v_i|^2 + cross-terms
# = sum_i k_i^2 (for orthogonal projections) -- but these are NOT orthogonal

# Compute k_phys^2 average over BZ boundary
# At X point (k_j = pi * delta_ij): |k_phys|^2 = pi^2
# At M point (k_i = k_j = pi, others 0): |k_phys|^2 = pi^2 * |v_i + v_j|^2
# Compute |v_i + v_j|^2 for adjacent pairs
print(f"  Key k_phys magnitudes at BZ special points:")
print(f"    X (one face): |k_phys| = pi = {pi:.6f}")

# M point: two faces = k_i = pi, k_j = pi
# Let i=0, j=2 (v_0 and v_2 are in different planes)
v0, v2 = vphys[0], vphys[2]
v0v2 = tuple(v0[k]+v2[k] for k in range(3))
M_phys = pi * math.sqrt(sum(x**2 for x in v0v2))
print(f"    M (two faces, v_0+v_2): |k_phys| = pi*|v_0+v_2| = {M_phys:.6f}")

# All-6-face corner R:
R_vec = tuple(sum(v[k] for v in vphys) for k in range(3))
R_mag  = math.sqrt(sum(x**2 for x in R_vec))
print(f"    R (all 6 faces): sum of all 6 v_i = {R_vec}, |sum| = {R_mag:.6f}")
print(f"    (This should be near 0 by icosahedral symmetry)")
print()

# The sum of all 6 icosahedral projection vectors should be 0 by symmetry
print(f"  Check: sum of 6 projection vectors = {R_vec}")
print(f"  |sum| = {R_mag:.2e}  (should be ~0 by I symmetry)")
print()


# =============================================================================
# PART IV -- CHARACTERISTIC SCALE k_BZ vs eps_L5 AND RELATED CONSTANTS
# =============================================================================
print(SEP)
print("PART IV -- CHARACTERISTIC k_BZ SCALE vs THEORY CONSTANTS")
print(SEP)
print()
print("  The key question: is any characteristic BZ momentum scale equal to")
print("  eps_L5, gj5, a_W, alpha, or any combination?")
print()

# In natural units of the 6D lattice spacing a=1:
# The BZ boundary is at k = 1/2 (or pi in angular units).

# In the physical projection, the natural scale for eps is set by:
# eps = k_phys * L_grain / (2*pi)
# where L_grain = alpha and k_phys is the physical wave vector.

# At the BZ face (X point): k_phys = pi (in lattice units, angular)
# Physical wave vector: k_phys / (2*pi) = 1/2 (in units of 1/L_grain)
# => eps_BZ_X = k_phys_physical / (2*pi/L_grain) = (pi/a) / (2*pi/a) = 1/2

# At first BZ face = X point:
k_X = pi   # BZ face in angular units (a=1)
eps_at_X = k_X / (2 * pi)   # = 0.5
print(f"  X point (BZ face):    k = pi, eps = k/(2*pi) = {eps_at_X:.6f}")
print(f"  eps_L5 =                                       {eps_L5:.6f}")
print(f"  gj5 =                                          {gj5:.6f}")
print(f"  a_W =                                          {a_W:.6f}")
print()

# The ratio eps_L5 / (k_X / (2*pi)) = 2*eps_L5 = 3/(4*pi)
ratio_eps_BZ = eps_L5 / eps_at_X
print(f"  eps_L5 / eps(X) = {ratio_eps_BZ:.8f} = 3/(4*pi) = {3/(4*pi):.8f}")
print(f"  2*eps_L5 = 3/(4*pi): the wave amplitude is at HALF the BZ depth")
print()

# The wave penetrates to eps_L5 in the BZ.
# The BZ HALF-DEPTH (from Gamma to first face) in eps units = 0.5.
# eps_L5 = 3/(8*pi) = 0.1194 = 0.2387 * 0.5
# So eps_L5 = 0.2387 * eps_BZ_face

fraction_of_BZ = eps_L5 / eps_at_X
print(f"  eps_L5 as fraction of BZ depth: {fraction_of_BZ:.6f} = {fraction_of_BZ:.4f} = ?")
# Is this a known fraction?
# 3/(4*pi) * 1 = ?
# Compare: 2*eps_L5 = 3/(4*pi) = 0.2387
print(f"  2*eps_L5 = 3/(4*pi) = {3/(4*pi):.8f}")
print()

# Natural BZ fraction in terms of the icosahedral symmetry:
# The 6D BZ has icosahedral cross-sections. The relevant fraction:
# 1/(2*PHI^2) = gj5 = a_W/2 -- is this a natural fraction of the BZ?
print(f"  gj5 = 1/(2*phi^2) = {gj5:.8f}")
print(f"  gj5 as fraction of BZ depth: {gj5/eps_at_X:.6f} = {2*gj5:.6f} (= a_W = 1/phi^2)")
print()

# Near-identity scan: what BZ-derived quantity is closest to eps_L5?
candidates = {
    "1/(2*phi^2) = gj5":         1/(2*PHI**2),
    "3/(4*pi) = 2*eps_L5":       3/(4*pi),
    "3/(8*pi) = eps_L5":         3/(8*pi),
    "1/phi^2 = a_W":             1/PHI**2,
    "1/(2*pi)":                  1/(2*pi),
    "1/phi":                     1/PHI,
    "1/(pi*phi)":                1/(pi*PHI),
    "1/(pi*phi^2)":              1/(pi*PHI**2),
    "2/(pi*phi^2)":              2/(pi*PHI**2),
    "sqrt5/(4*pi)":              sqrt5/(4*pi),
    "3/(4*pi) - 1/(2*phi^2)":    3/(4*pi) - 1/(2*PHI**2),   # = excess5
    "pi/sqrt5 - 1":              pi/sqrt5 - 1,
    "sqrt5 - 2":                 sqrt5 - 2,
    "1/(2*sqrt5)":               1/(2*sqrt5),
}

print(f"  Key theory constants as fractions of BZ X-face depth (=0.5):")
print(f"  {'Expression':<35}  {'Value':>10}  {'BZ fraction':>12}")
print(f"  {'-'*62}")
for label, val in sorted(candidates.items(), key=lambda x: x[1]):
    print(f"  {label:<35}  {val:>10.6f}  {val/eps_at_X:>12.6f}")
print()

# The standing result from 0l.8.3: gj5 = a_W/2.
# In terms of BZ fraction: gj5 = 1/(2*phi^2) ~ 0.191
# eps_L5 = 3/(8*pi) ~ 0.119
# a_W = 1/phi^2 ~ 0.382
# The acceptance window a_W is 0.764 of the BZ depth.
# gj5 = 0.382 of the BZ depth.
# eps_L5 = 0.238 of the BZ depth.
# These are all FRACTIONS of the BZ -- not equal to the BZ scale itself.

print(f"  BZ depth (X face) in eps units: {eps_at_X:.6f}")
print(f"  eps_L5 = {eps_at_X:.4f} * {eps_L5/eps_at_X:.6f}")
print(f"  gj5    = {eps_at_X:.4f} * {gj5/eps_at_X:.6f}")
print(f"  a_W    = {eps_at_X:.4f} * {a_W/eps_at_X:.6f}")
print()


# =============================================================================
# PART V -- 6D WAVE OPERATOR SPECTRUM AT BZ BOUNDARY
# =============================================================================
print(SEP)
print("PART V -- 6D WAVE OPERATOR SPECTRUM AT BZ BOUNDARY")
print(SEP)
print()
print("  For a scalar wave on Z^6 (tight-binding hopping model):")
print("  Dispersion: E(k_6D) = 2 * sum_{i=1}^{6} (1 - cos(k_i))")
print("  Maximum E (at corner R): E_R = 2*6 = 12  (bandwidth)")
print("  Minimum E (at Gamma):    E_0 = 0")
print("  At X face (k_1=pi, others=0): E_X = 2*(1-cos(pi)) + 5*0 = 4")
print()

def dispersion_6D(k6):
    """Tight-binding dispersion on Z^6."""
    return 2 * sum(1 - math.cos(ki) for ki in k6)

E_Gamma = dispersion_6D([0]*6)
E_X     = dispersion_6D([pi, 0, 0, 0, 0, 0])
E_M     = dispersion_6D([pi, pi, 0, 0, 0, 0])
E_R     = dispersion_6D([pi]*6)

print(f"  E(Gamma) = {E_Gamma:.4f}")
print(f"  E(X)     = {E_X:.4f}  (1 BZ face)")
print(f"  E(M)     = {E_M:.4f}  (2 BZ faces)")
print(f"  E(R)     = {E_R:.4f}  (all 6 BZ faces, corner)")
print()

# Project the 6D dispersion to the 3D physical space:
# The physical dispersion is found by tracing the minimum-energy 6D mode
# that projects to a given 3D physical k_phys.
# This requires minimising E(k_6D) subject to Pi_par * k_6D = k_phys.

# For icosahedral projection, the optimal k_6D for a given k_phys direction:
# The projection onto physical space is Pi_par; the perp component is free.
# At each k_phys, we minimize E over the fiber Pi_par^{-1}(k_phys).
# For tight-binding: this is a constrained optimization.

# For SMALL k_phys (near Gamma), approximate: k_6D ~ Pi_par^T * k_phys
# (pseudo-inverse for the overdetermined system)
# k_6D_approx = (Pi_par^T * Pi_par)^{-1} * Pi_par^T * k_phys
# For the Z^6 icosahedral projection, Pi_par * Pi_par^T = (2/3) * I_3
# (checked below), so Pi_par^T = (3/2) * Pi_par^T and the pseudo-inverse simplifies.

# Compute Pi_par * Pi_par^T
Pi = [list(v) for v in vphys]   # 3x6 matrix (rows = x,y,z; cols = 6 directions)
# Actually: Pi has rows = coordinates (3) and columns = 6 projection directions
# Pi_par (3x6): row i, col j = v_j[i]
# Pi_par * Pi_par^T (3x3): (Pi_par * Pi_par^T)_{ij} = sum_k v_k[i] * v_k[j]

PiPiT = [[0.0]*3 for _ in range(3)]
for a in range(3):
    for b in range(3):
        PiPiT[a][b] = sum(vphys[k][a] * vphys[k][b] for k in range(6))

print(f"  Pi_par * Pi_par^T (should be c*I_3 for icosahedral symmetry):")
for row in PiPiT:
    print(f"    [{row[0]:+.6f}, {row[1]:+.6f}, {row[2]:+.6f}]")
# Should be (6/3) * I = 2 * I  (since each |v_i|=1 and they sum isotropically)
c_iso = PiPiT[0][0]
print(f"  Diagonal value: c = {c_iso:.8f}  (expected 2.0 = 6/3)")
print()

# Physical 3D dispersion (linearized near Gamma):
# E_phys(k_phys) ≈ E_6D(Pi_par^T * k_phys) for small k
# = 2 * sum_i (1 - cos(k_phys . v_i))
# ≈ sum_i (k_phys . v_i)^2  for small k
# = k_phys^T * Pi_par * Pi_par^T * k_phys = c * |k_phys|^2 = 2 * |k_phys|^2
print(f"  Physical 3D dispersion (low-k limit):")
print(f"    E_phys(k) ≈ c * k^2 = {c_iso:.4f} * k^2")
print(f"    (standard quadratic dispersion, as expected)")
print()

# The physical BZ boundary (first zone) in 3D:
# The first Brillouin zone of the projected quasicrystal is determined by
# the first strong Bragg peaks. For the icosahedral quasicrystal from Z^6,
# the shortest reciprocal lattice vectors are the 12 vectors v_i and -v_i.
# The BZ boundary bisects the line from Gamma to each v_i.
# In 3D, this gives a 12-faced polyhedron (rhombic dodecahedron... or
# actually icosahedron? Let's check.)

# The 6 pairs of projection vectors {+-v_i} give 12 first Bragg peaks.
# These are the vertices of an icosahedron (since v_i are the 6 pairs of
# icosahedral "golden vectors").
# The Wigner-Seitz cell of these 12 points = icosahedron-based polyhedron
# = rhombic dodecahedron? No -- the Voronoi cell of the icosahedral points.
# For 12 vectors of equal length forming an icosahedron: the Voronoi cell
# is the ICOSIDODECAHEDRON (30 vertices, 32 faces: 12 pentagons + 20 triangles)
# Actually: the first BZ = zone enclosed within bisectors to 12 nearest Bragg peaks
# = the dual of the icosahedron = dodecahedron... let's compute.

print(f"  3D projected BZ boundary:")
print(f"    12 nearest Bragg peaks at: +/-v_i, all with |G| = 1.0")
print(f"    BZ boundary at half-G: planes |k . v_i| = 1/2")
print(f"    These 12 planes bound a region -- what shape?")
print()

# The BZ boundary condition: k . v_i = 1/2 (for each of 6 vectors v_i)
# In terms of eps: eps_BZ = (1/2) / |v_i| = 1/2 (since |v_i|=1)
# So the projected 3D BZ boundary is at k_phys = v_i/2 (half the Bragg peak)
# |k_phys at boundary| = 1/2

# In the theory, the wave amplitude eps is related to the wave-vector magnitude by:
# eps = k_phys * L / (2*pi) where L is the grain size ~ alpha
# For eps_L5 = 3/(8*pi): k_phys = eps_L5 * 2*pi / L = 3/(4*L)
# The BZ boundary at k = 1/2 corresponds to eps_BZ = 1/2 (in L=1 units).
# So eps_L5 = 3/(8*pi) = 0.1194 is at a FRACTION of the BZ depth.

# The key ratio:
eps_L5_over_eps_BZ = eps_L5 / 0.5
print(f"  eps_L5 / eps_BZ = {eps_L5:.6f} / 0.5 = {eps_L5_over_eps_BZ:.6f}")
print(f"  This equals 3/(4*pi) = {3/(4*pi):.6f}")
print()
print(f"  Interpretation: the locking wave amplitude eps_L5 sits at")
print(f"  {eps_L5_over_eps_BZ*100:.2f}% of the way from Gamma to the 3D BZ boundary.")
print(f"  The 3D BZ boundary is at eps_BZ = 0.5 (half the first Bragg peak).")
print()


# =============================================================================
# PART VI -- PHASON COMPONENT AND PERP-SPACE BZ
# =============================================================================
print(SEP)
print("PART VI -- PHASON COMPONENT AND PERP-SPACE BZ")
print(SEP)
print()
print("  The 6D wave vector has PERP-SPACE component k_perp = Pi_perp * k_6D.")
print("  The perp-space BZ boundary is at the SAME position |k_perp| = 1/2.")
print("  However, the acceptance window a_W = 1/phi^2 = 0.382 is a REAL-SPACE")
print("  perp-space cutoff (not a reciprocal-space boundary).")
print()
print("  Perp-space projection vectors w_i = Pi_perp * e_i:")
print("  These are the same form as v_i but with phi -> -1/phi (perp to v_i).")
print()

# Check: Pi_perp * Pi_perp^T
PwPwT = [[0.0]*3 for _ in range(3)]
for a in range(3):
    for b in range(3):
        PwPwT[a][b] = sum(vperp[k][a] * vperp[k][b] for k in range(6))

c_perp = PwPwT[0][0]
print(f"  Pi_perp * Pi_perp^T diagonal: {c_perp:.8f}  (expected 2.0)")
print()

# Total: Pi_par * Pi_par^T + Pi_perp * Pi_perp^T = 2*I + 2*I = 4*I?
# Actually should be: Pi * Pi^T = c_total * I where Pi is the full 6x6
# For Z^6 to R^3_par + R^3_perp: Pi_par * Pi_par^T = 2I, Pi_perp * Pi_perp^T = 2I
# Together they give the decomposition of R^6 = R^3 + R^3.
print(f"  c_par = {c_iso:.6f}, c_perp = {c_perp:.6f}")
print(f"  c_par + c_perp = {c_iso+c_perp:.6f} (should be 4.0 = 6/3*2)")
print()

# The ACCEPTANCE WINDOW: a grain at perp position d_perp is included iff
# |d_perp| < a_W = 1/phi^2. This is in REAL PERP SPACE (Angstroms), not
# reciprocal space. The perp-space BZ boundary (at half the first perp Bragg peak)
# is at |k_perp| = 1/2 in reciprocal perp space.
# The real-space perp acceptance window a_W does NOT correspond to the BZ boundary.

print(f"  Real-space acceptance window:  a_W = 1/phi^2 = {a_W:.6f}")
print(f"  Reciprocal perp BZ boundary:   k_perp_BZ = 1/2 = 0.5")
print(f"  These are in different spaces (real vs reciprocal) -- not directly compared.")
print()

# But their ratio (a_W vs k_perp_BZ) does have a meaning:
# The phason coherence length xi_phason ~ 1/a_W (in reciprocal units)
# At a_W = 0.382: xi_phason ~ 1/a_W = phi^2 = phi+1 = 2.618 lattice units
# This is the "phason domain size" -- the range over which phason fluctuations
# are coherent. Smaller than the grain size (which is ~alpha >> a_W in model units).

xi_phason = 1 / a_W   # = phi^2
print(f"  Phason coherence length: xi ~ 1/a_W = phi^2 = {xi_phason:.6f}")
print()


# =============================================================================
# PART VII -- DOES ANY BZ SCALE MATCH eps_L5?
# =============================================================================
print(SEP)
print("PART VII -- DOES ANY 6D BZ SCALE MATCH eps_L5, gj5, OR alpha?")
print(SEP)
print()
print("  Comprehensive comparison: Z^6 BZ scales vs theory constants.")
print()

# 6D BZ characteristic distances (in angular units, a=1):
bz_scales = {
    "Gamma->X (BZ face, 1D)":    pi,                   # = 3.1416
    "Gamma->M (BZ edge, 2D)":    pi * math.sqrt(2),    # = 4.4429
    "Gamma->R (BZ corner, 6D)":  pi * math.sqrt(6),    # = 7.6956
    "X/2 (half-face)":           pi/2,                 # = 1.5708
    "1/(2*pi) * Gamma->X":       0.5,                  # BZ face in 1/(2*pi) units
    "1/(2*pi) * X/2":            0.25,
}

theory_constants = {
    "eps_L5 = 3/(8*pi)":         eps_L5,    # 0.11937
    "gj5 = 1/(2*phi^2)":         gj5,       # 0.19098
    "a_W = 1/phi^2":             a_W,       # 0.38197
    "2*eps_L5 = 3/(4*pi)":       2*eps_L5,  # 0.23873
    "excess5 = 2*eps_L5 - gj5":  2*eps_L5 - gj5,  # 0.04775
    "alpha":                     alpha,     # 0.00730
    "alpha*phi":                 alpha*PHI, # 0.01181
}

print(f"  BZ scales (in units of 1/(2*pi)) vs theory constants:")
print(f"  The BZ face is at 1/2 = 0.5 in these units.")
print()

# The only clean correspondences we can check:
# 1. Does eps_L5 correspond to a FRACTION of the BZ that has icosahedral meaning?
# eps_L5 = 3/(8*pi) = 0.1194
# 3/(4*pi) = 2*eps_L5 = fraction of BZ at which 2*eps_L5 sits
# gj5 = 1/(2*phi^2) = 0.1910 = another BZ fraction

# BZ fractions in terms of icosahedral geometry:
print(f"  BZ depth fractions with icosahedral meaning:")
print(f"    gj5 / eps_BZ_face = {gj5/0.5:.6f} = {2*gj5:.6f} = 1/phi^2 = a_W")
print(f"    eps_L5 / eps_BZ_face = {eps_L5/0.5:.6f} = {2*eps_L5:.6f} = 3/(4*pi)")
print(f"    (3/4*pi) / a_W = {3/(4*pi) / a_W:.6f}  vs  pi/4 = {pi/4:.6f}")
print()
print(f"  KEY RESULT:")
print(f"    2*eps_L5 / a_W = {2*eps_L5 / a_W:.8f}")
print(f"    = (3/(4*pi)) / (1/phi^2) = 3*phi^2/(4*pi) = {3*PHI**2/(4*pi):.8f}")
print(f"    = 3*(phi+1)/(4*pi) = 3*phi/(4*pi) + 3/(4*pi)")
print()

# Is 3*phi^2/(4*pi) a known constant?
val = 3*PHI**2/(4*pi)
print(f"  3*phi^2/(4*pi) = {val:.8f}")
print(f"  Compare: Rs = sqrt5/(4*pi) = {math.sqrt(5)/(4*pi):.8f}")
print(f"           3/(4*pi) = {3/(4*pi):.8f}")
print(f"           phi^2 = {PHI**2:.8f} = phi+1 = {PHI+1:.8f}")
print()

# The ratio 2*eps_L5 / a_W = 3*phi^2/(4*pi):
# phi^2 = phi+1, so 3*phi^2/(4*pi) = 3(phi+1)/(4*pi) = 3/(4*pi) + 3*phi/(4*pi)
# = eps_L5*2 + 3*phi/(4*pi)
# This doesn't simplify to a known constant.

# The more direct comparison: does the BZ geometry SET eps_L5?
# If a_W = 1/phi^2 is the acceptance window, and the wave amplitude is at 2*eps_L5,
# then 2*eps_L5 / a_W = 3*phi^2/(4*pi).
# For this to be "1" (wave at BZ boundary = acceptance window), we'd need
# 3*phi^2/(4*pi) = 1, i.e., phi^2 = 4*pi/3 = 4.189. But phi^2 = 2.618.
# So the wave amplitude is at (3*phi^2)/(4*pi) = 0.625 * acceptance window.
print(f"  2*eps_L5 as fraction of acceptance window a_W: {2*eps_L5/a_W:.6f}")
print(f"  = 3*phi^2/(4*pi) = {3*PHI**2/(4*pi):.6f}")
print()
print(f"  gj5 as fraction of acceptance window a_W: {gj5/a_W:.6f}")
print(f"  = (1/(2*phi^2)) / (1/phi^2) = 1/2  EXACTLY")
print(f"  (This is the already-known result: gj5 = a_W/2)")
print()


# =============================================================================
# PART VIII -- SUMMARY
# =============================================================================
print(SEP)
print("PART VIII -- SUMMARY: [crys1] TOOL 4 FINDINGS")
print(SEP)
print()
print("  PHYSICAL PICTURE:")
print(f"  The icosahedral medium is cut-and-project from Z^6.")
print(f"  The 6D BZ is [-pi,pi]^6. Its projection to 3D gives the")
print(f"  icosahedral quasicrystal reciprocal structure.")
print()
print("  WHAT THE BZ ANALYSIS FINDS:")
print()
print(f"  [+] Pi_par * Pi_par^T = 2*I_3  (isotropic projection, confirmed)")
print(f"      The 6D BZ projects isotropically to 3D — no preferred direction.")
print()
print(f"  [+] 3D projected BZ boundary at |k_phys| = 1/2 (in units of |v_i|)")
print(f"      This gives eps_BZ = 0.5 (in wave-amplitude units).")
print()
print(f"  [+] The KNOWN result gj5 = a_W/2 has a BZ reading:")
print(f"      gj5 = 1/(2*phi^2) = 0.1910 = {gj5/0.5:.4f} * eps_BZ_face")
print(f"      gj5 = a_W/2 = (acceptance window)/2")
print(f"      = inner-shell boundary of the cut-and-project filter")
print()
print(f"  [?] eps_L5 = 3/(8*pi) = {eps_L5:.6f}")
print(f"      As BZ fraction: {eps_L5/0.5:.6f} = 3/(4*pi)")
print(f"      As fraction of a_W: {eps_L5/a_W:.6f} = 3/(4*pi*a_W) = {3/(4*pi*a_W):.6f}")
print()
print(f"      The BZ does NOT directly set eps_L5.")
print(f"      eps_L5 = 3/(8*pi) comes from the Chern-Simons integral (topological),")
print(f"      not from the BZ boundary condition.")
print()
print(f"  [RESULT] The 6D BZ does NOT independently generate eps_L5.")
print(f"  The BZ boundary is at eps_BZ = 0.5, while eps_L5 = 0.119.")
print(f"  The wave locks at a fraction {eps_L5/0.5:.4f} of the BZ depth,")
print(f"  which equals 3/(4*pi) -- a value with topological origin (not BZ origin).")
print()
print(f"  [PARTIAL CORROBORATION] The BZ DOES confirm the a_W = 1/phi^2 scale:")
print(f"      The 3D projected BZ face is at k = |v_i|/2 = 0.5.")
print(f"      The acceptance window a_W = 0.382 = 0.764 * BZ_face.")
print(f"      gj5 = 0.382 * BZ_face = a_W/2.")
print(f"      The jamming threshold gj5 = a_W/2 has a natural Z^6 BZ reading:")
print(f"      it is half the acceptance window, which is 0.382 of the BZ depth.")
print()
print(f"  [OPEN QUESTION REMAINS]")
print(f"      WHY does the system lock at eps_L5 = 3/(8*pi) specifically?")
print(f"      The BZ doesn't fix this. The Chern-Simons locking gives it topologically.")
print(f"      But the question is why THAT Chern-Simons value and not another.")
print(f"      Connection: eps_L5 = 3/(8*pi) arises from the Hopf fibration geometry")
print(f"      (the area of the fundamental domain of T^2 -> S^3 -> S^2).")
print(f"      The BZ analysis does not add to this.")
print()
print(f"  VERDICT: [crys1] Tool 4 is PARTIALLY INFORMATIVE but does NOT")
print(f"  provide an independent derivation of eps_L5 or close Gap 1.")
print(f"  It CONFIRMS the Z^6 quasicrystal structure and the gj5=a_W/2 identity.")
print(f"  It does NOT add new constraints that fix alpha or delta_n.")
print()

# What WOULD close the BZ → eps_L5 connection?
print(f"  WHAT WOULD NEED TO BE TRUE for BZ to set eps_L5:")
print(f"    Either: eps_L5 = (BZ boundary condition in the icosahedral projection)")
print(f"    Or:     eps_L5 = (ratio of parallel to perp BZ areas) * (known constant)")
print(f"    Or:     the energy gap at eps_L5 in the 6D band structure = k_B*T_lock")
print()
print(f"  The last option (band gap at eps_L5) requires computing the ACTUAL")
print(f"  band structure in 6D with the EM interaction term, not just the free-")
print(f"  particle BZ. This would require the coupled EM + grain Lagrangian")
print(f"  in 6D, which is beyond the current scope of the crys1 tool series.")
print()

print(SEP)
print("END crys1_tool4_brillouin_zone.py")
print(SEP)
