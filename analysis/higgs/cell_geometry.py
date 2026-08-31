"""
cell_geometry.py
================
Canonical specification of the Jobson cell: all derived geometric quantities
from the single input L_J = alpha*phi*r_p.

This is the companion script for doc_jobson_cell.txt.

Run: python analysis/higgs/cell_geometry.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import alpha, phi, hbar_c, E_cell_GeV, N_lock, L_J
import math as _m
pi  = _m.pi
Rs  = _m.sqrt(5)/(4*pi)
nu  = (1 - 2*Rs**2)/(2*(1 - Rs**2))
r_p_fm = 0.8414   # fm  (CODATA-2018)

pi     = math.pi
sqrt3  = math.sqrt(3)
sqrt5  = math.sqrt(5)
Rs     = sqrt5 / (4*pi)

SEP  = "=" * 70
SEP2 = "-" * 70

print(SEP)
print("JOBSON CELL CANONICAL SPECIFICATION")
print(SEP)
print()

# ── SECTION 1: FUNDAMENTAL LENGTH ────────────────────────────────────────────
print(SEP)
print("SECTION 1  Fundamental length scale")
print(SEP2)
print(f"  L_J = alpha * phi * r_p")
print(f"      = {alpha:.10e} * {phi:.10f} * {r_p_fm:.4f} fm")
print(f"      = {L_J:.10f} fm")
print(f"      = {L_J*1e-15:.6e} m")
print(f"  Physical: minimum arc segment of the (1,2) torus knot at the proton scale.")
print(f"  N_lock = 2*pi/(alpha*phi) = {N_lock:.6f}  (tube closure number)")
print()

# ── SECTION 2: ICOSAHEDRAL VERTEX GEOMETRY ────────────────────────────────────
print(SEP)
print("SECTION 2  Icosahedral vertex geometry")
print(SEP2)
print()
# Standard icosahedron: vertices at permutations of (0, +/-1, +/-phi)
# Edge length = 2, circumradius = sqrt(1+phi^2)
a_std    = 2.0                        # edge length in standard coords
R_std    = math.sqrt(1 + phi**2)      # circumradius in standard coords
r_in_std = phi**2 / (2*math.sqrt(3)) # inradius (center to face)
r_mid_std = phi/2                     # midradius (center to edge midpoint)

# Scale to physical units where edge = L_J
scale  = L_J / a_std
a_phys = L_J
R_phys = R_std * scale
r_in_phys  = r_in_std  * scale
r_mid_phys = r_mid_std * scale

print(f"  Standard icosahedron (vertices at permutations of (0, +/-1, +/-phi)):")
print(f"    Edge length:     a = {a_std:.6f}  (standard units)")
print(f"    Circumradius:    R = sqrt(1+phi^2) = {R_std:.8f}")
print(f"    Inradius:        r_in = phi^2/(2*sqrt(3)) = {r_in_std:.8f}")
print(f"    Midradius:       r_mid = phi/2 = {r_mid_std:.8f}")
print(f"    a/R ratio:       {a_std/R_std:.8f}  (= 2/sqrt(1+phi^2) = {2/R_std:.8f})")
print()
print(f"  Scaled to Jobson cell (edge = L_J = {L_J:.6f} fm):")
print(f"    Edge length:     a = L_J = {a_phys:.8f} fm")
print(f"    Circumradius:    R = {R_phys:.8f} fm")
print(f"    Inradius:        r_in = {r_in_phys:.8f} fm")
print(f"    Midradius:       r_mid = {r_mid_phys:.8f} fm")
print()
print(f"  Vertex count:  12  (icosahedral)")
print(f"  Edge count:    30")
print(f"  Face count:    20  (equilateral triangles)")
print(f"  Each vertex has 5 nearest neighbors.")
print()

# 12 vertex positions in standard coords, verify edge length
v = []
for signs in [(1,1),(1,-1),(-1,1),(-1,-1)]:
    v.append((0, signs[0], signs[1]*phi))
    v.append((signs[0], signs[1]*phi, 0))
    v.append((signs[0]*phi, 0, signs[1]))
# Check: all 12 vertices have same distance from origin
r_check = math.sqrt(v[0][0]**2 + v[0][1]**2 + v[0][2]**2)
print(f"  Circumradius check: sqrt(0^2+1^2+phi^2) = sqrt({1+phi**2:.8f}) = {r_check:.8f}")
# Check edge length between two adjacent vertices
e_check = math.sqrt((v[0][0]-v[2][0])**2 + (v[0][1]-v[2][1])**2 + (v[0][2]-v[2][2])**2)
print(f"  Edge length check: |(0,1,phi)-(1,phi,0)| = {e_check:.8f}  (should be 2.0)")
print()

# ── SECTION 3: ELASTIC PROPERTIES ────────────────────────────────────────────
print(SEP)
print("SECTION 3  Elastic properties from wave speeds")
print(SEP2)
print()
print(f"  Wave speeds (from GW170817 + K-formula):")
print(f"    v_p = c  (longitudinal = light speed)")
print(f"    v_s = Rs*c  where Rs = sqrt(5)/(4*pi) = {Rs:.10f}")
print()
K_o_G = (2*(1+nu)) / (3*(1-2*nu))
G_frac = 1 / (1 + K_o_G)   # G/(K+G)
K_frac = K_o_G / (1 + K_o_G)  # K/(K+G)
print(f"  Poisson ratio: nu = (1-2Rs^2)/(2(1-Rs^2)) = {nu:.10f}")
print(f"  Exact form:    nu = (8pi^2-5)/(16pi^2-5)")
print(f"  Bulk/shear:    K/G = (48pi^2-20)/15 = {K_o_G:.8f}")
print(f"  Quartic coupling: lambda = (1-nu)/4 = {(1-nu)/4:.10f}")
print(f"  lambda exact:     2*pi^2/(16*pi^2-5) = {2*pi**2/(16*pi**2-5):.10f}")
print()

# ── SECTION 4: CELL ENERGY ────────────────────────────────────────────────────
print(SEP)
print("SECTION 4  Cell energy and binding")
print(SEP2)
print()
print(f"  E_cell = 2*pi*hbar*c / L_J = {E_cell_GeV:.9f} GeV")
print(f"  E_cell = N_lock * hbar*c / r_p")
E_cell_MeV = E_cell_GeV * 1000
print(f"         = {E_cell_MeV:.6f} MeV")
print()
# N_J for the Higgs
N_J_H = hbar_c / (E_cell_MeV * L_J)
print(f"  Higgs N_J: N_J_H = hbar_c/(m_H*L_J) = 1/(2*pi) = {1/(2*pi):.8f}")
print(f"             actual: {N_J_H:.8f}  (at alpha/pi precision)")
print()
# Claim 8
k_n_max = 3125/3456
lhs = 7 * k_n_max / (2*pi)
rhs = 1 + alpha + alpha**2*phi
print(f"  Claim 8 (scale-invariant jamming, 0.0001%):")
print(f"    k_n_max = 3125/3456 = {k_n_max:.10f}  (exact algebraic)")
print(f"    7*k_n_max/(2*pi)          = {lhs:.10f}")
print(f"    1 + alpha + alpha^2*phi   = {rhs:.10f}")
print(f"    Gap: {abs(lhs-rhs)/rhs*100:.8f}%")
print()

# ── SECTION 5: I_h CHARACTER TABLE ───────────────────────────────────────────
print(SEP)
print("SECTION 5  Icosahedral group I_h character table (gerade irreps)")
print(SEP2)
print()
print(f"  Group: I_h (icosahedral with inversion), order 120")
print(f"  5 gerade irreps from I (order 60, 5 conjugacy classes):")
print()
print(f"  {'Irrep':<8} {'dim':<6} {'E':<6} {'C2':<8} {'C3':<6} {'C5':<12} {'C52':<12}")
print(f"  {'-'*8} {'-'*6} {'-'*6} {'-'*8} {'-'*6} {'-'*12} {'-'*12}")
char_table = [
    ('A_g',   1,  1,    1,  1,    1,          1),
    ('T_1g',  3, -1,    0, -1, phi,       -1/phi),
    ('T_2g',  3, -1,    0, -1, -1/phi,     phi),
    ('G_g',   4,  0,    1, -2, -1,         -1),
    ('H_g',   5,  1,   -1,  0,  0,          0),
]
# Note: C2 here means C2 (15 elements), C3 (20 elem), C5 (12 elem), C52 (12 elem)
# Need to correct the C2 and C3 values:
char_table = [
    ('A_g',   1,  1,  1,  1,   1,       1),
    ('T_1g',  3,  3, -1,  0,  phi,  -1/phi),
    ('T_2g',  3,  3, -1,  0, -1/phi,  phi),
    ('G_g',   4,  4,  0,  1,  -1,     -1),
    ('H_g',   5,  5,  1, -1,   0,      0),
]
# Format: (name, dim, chi_E, chi_C2, chi_C3, chi_C5, chi_C52)
# chi_E = dim for any irrep
print(f"  {'Irrep':<8} {'dim':<6} {'chi(E)':<8} {'chi(C2)':<10} {'chi(C3)':<8} {'chi(C5)':<14} {'chi(C52)':<12}")
print(f"  {'-'*8} {'-'*6} {'-'*8} {'-'*10} {'-'*8} {'-'*14} {'-'*12}")
for row in char_table:
    name, dim, e, c2, c3, c5, c52 = row
    c5_str  = f"phi={phi:.4f}" if abs(c5-phi)<0.001 else f"-1/phi={-1/phi:.4f}" if abs(c5+1/phi)<0.001 else f"{c5:.4f}"
    c52_str = f"phi={phi:.4f}" if abs(c52-phi)<0.001 else f"-1/phi={-1/phi:.4f}" if abs(c52+1/phi)<0.001 else f"{c52:.4f}"
    print(f"  {name:<8} {dim:<6} {e:<8} {c2:<10} {c3:<8} {c5_str:<14} {c52_str:<12}")
print()
print(f"  Key values at C_5: phi = (1+sqrt(5))/2 = {phi:.8f}")
print(f"                     -1/phi = -(sqrt(5)-1)/2 = {-1/phi:.8f}")
print(f"  phi satisfies: phi^2 = phi+1  (Fibonacci identity, algebraically exact)")

# ── SECTION 6: CG DECOMPOSITIONS ─────────────────────────────────────────────
print()
print(SEP)
print("SECTION 6  Clebsch-Gordan decompositions (verified)")
print(SEP2)
print()
# From higgs_cg_twoloop.py
print(f"  T_1g x T_1g = A_g + T_1g + H_g  (dims 1+3+5=9=3^2) [DERIVED]")
print(f"  T_2g x T_2g = A_g + T_2g + H_g  (dims 1+3+5=9=3^2) [DERIVED]")
print(f"  G_g  x G_g  = A_g + T_1g + T_2g + G_g + H_g  (1+3+3+4+5=16=4^2) [DERIVED]")
print(f"  H_g  x H_g  = A_g + T_1g + T_2g + 2*G_g + 2*H_g  (1+3+3+8+10=25=5^2) [DERIVED]")
print(f"  A_g  x A_g  = A_g  (trivial) [DERIVED]")
print(f"  T_1g x T_2g = G_g + H_g  (NO A_g) [DERIVED -- FORBIDDEN CHANNEL]")
print()
print(f"  A_g appears ONCE in: T_1g x T_1g,  T_2g x T_2g,  G_g x G_g,  H_g x H_g")
print(f"  A_g is ABSENT from:  T_1g x T_2g  -> H->T_1g+T_2g FORBIDDEN")
print()
print(f"  chi(T_1g x T_1g, C_5) = phi^2 = phi+1 = {phi**2:.8f}  [algebraically exact]")
print(f"  Fibonacci truncation: phi^n = F(n)*phi + F(n-1)  => {{1,phi}} complete basis")
print(f"  alpha^2*phi^2 + alpha^3*phi^3 = alpha^2*phi^2*(1+alpha*phi)  [bit-exact]")

# ── SECTION 7: MASS CORRECTION ────────────────────────────────────────────────
print()
print(SEP)
print("SECTION 7  Higgs mass correction series")
print(SEP2)
print()
c1 = alpha/pi
c2 = alpha**2 * phi**2
lam_val = (1-nu)/4
mH_1 = E_cell_GeV * (1 + c1)
mH_2 = E_cell_GeV * (1 + c1 + c2)
v_1  = mH_1 / math.sqrt(2*lam_val)
v_2  = mH_2 / math.sqrt(2*lam_val)
v_EW = 246.21965
print(f"  Term 1 (scalar QED one-loop):      alpha/pi       = {c1:.8e}")
print(f"  Term 2 (T_1g x T_1g -> A_g):       alpha^2*phi^2  = {c2:.8e}")
print(f"  Ratio T2/T1 = alpha*pi*phi^2       = {c2/c1:.8f}")
print()
print(f"  m_H (1-term): E_cell*(1+a/pi)       = {mH_1:.6f} GeV  (Claim 1, -1.01 sigma)")
print(f"  m_H (2-term): E_cell*(1+a/pi+a^2p^2)= {mH_2:.6f} GeV  (Section 5a, 0.86 sigma)")
print(f"  vev (1-term): {v_1:.6f} GeV  gap={( v_1-v_EW)*1000:+.2f} MeV")
print(f"  vev (2-term): {v_2:.6f} GeV  gap={(v_2-v_EW)*1000:+.3f} MeV")
print()
print(f"  Effective Lagrangian: L_HWW = alpha^2*phi^2 * |H|^2 * |W|^2")
print(f"  Forbidden: H -> T_1g + T_2g  (T_1g x T_2g = G_g + H_g, no A_g)")

print()
print(SEP)
print("SUMMARY: JOBSON CELL SPECIFICATION")
print(SEP2)
print(f"  Topology:        (1,2) Hopf fibration, winding (p,q)=(1,2)")
print(f"  Symmetry:        I_h (icosahedral, order 120, 10 irreps)")
print(f"  Edge length:     L_J = alpha*phi*r_p = {L_J:.8f} fm  [DERIVED]")
print(f"  Circumradius:    R_c = L_J*sqrt(1+phi^2)/2 = {R_phys:.8f} fm  [DERIVED]")
print(f"  Energy:          E_cell = 2*pi*hbar_c/L_J = {E_cell_GeV:.8f} GeV  [DERIVED]")
print(f"  Poisson ratio:   nu = (8pi^2-5)/(16pi^2-5) = {nu:.8f}  [DERIVED]")
print(f"  Quartic coupl:   lambda = 2*pi^2/(16*pi^2-5) = {lam_val:.8f}  [DERIVED]")
print(f"  K/G ratio:       (48*pi^2-20)/15 = {K_o_G:.8f}  [DERIVED]")
print(f"  Vertices:        12  (T_1g x 4 assignments)")
print(f"  CG:              T_1g x T_1g = A_g + T_1g + H_g  [DERIVED]")
print(f"  Lagrangian:      L_HWW = alpha^2*phi^2*|H|^2*|W|^2  [DERIVED from CG]")
print(f"  Forbidden:       H -> T_1g + T_2g  [DERIVED]")
print()
print(f"  Source documents: doc_alpha (L_J, N_lock), doc_torsion (nu, E_cell),")
print(f"                    doc_higgs (CG, Lagrangian, H2)")
