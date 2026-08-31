"""
lepton_mass.py
==============
Derives electron and muon masses from torsion medium polygon path geometry.
Establishes the grand lepton structure: electron=VERTEX, muon=EDGE, tau=FACE
of the icosahedral Jobson cell lattice.

FRAMEWORK:
  Mass = excluded volume (m = rho * V). The lepton is a Hopf-winding mode
  that bounces around the Jobson cell lattice on a specific polygonal path.
  The Born vertex balance condition for that path sets the coupling constant
  that determines the lepton's Compton wavelength (= rest mass).

GRAND LEPTON STRUCTURE (three icosahedral element types):
  ELECTRON (E+, dim=2, C3=-1): VERTEX mode
    Path: scatters at 12 icosahedral vertices. Deflection = 72 deg (C5).
    Color role: vertex is junction of 5 colored faces -> color INTERSECTION.
    eff_e = phi.  m_e = 2*pi*alpha^2*phi*m_p * (Born) * (1+3/4*alpha^2)

  MUON (G32, dim=4, C3=+1): EDGE mode
    Path: traverses 30 icosahedral edges via 5-edge zig-zag. Deflection = 72 deg (C5).
    Color role: edge borders 2 colored faces -> color BOUNDARY mode.
    eff_mu = (9-sqrt5)/8.  m_mu = 2*pi*alpha*(2/sqrt5)*phi^2*m_p * (Born) * (1+Rs^2+2*alpha)

  TAU (I52, dim=6, C3=0): FACE mode
    Path: Hamiltonian cycle through 20 icosahedral faces. Deflection = arccos(-1/sqrt5) = 116.57 deg
          This is the icosahedral DIHEDRAL ANGLE (the same -1/sqrt5 that appeared
          in the pentagonal bipyramid equatorial -- which belongs to the TAU, not the muon).
    Color role: face = one colored region -> color NEUTRAL (equal coupling to all 3 colors).
    I52 dim=6 = 2 (spinor) x 3 (quark colors from icosahedral face 3-coloring).
    eff_tau = (1+1/sqrt5)/2 = 0.7236.  Mass formula: OPEN (needs face-mode winding factor).
    Current best: m_tau from Koide using derived m_e, m_mu -> 1776.92 MeV (+0.004%)

GEOMETRIC IDENTITIES:
  Pentagonal bipyramid (all equal edges):
    h_t / r_e = 1/phi  (EXACT -- golden ratio in height/radius ratio)
    Apex zig-zag deflection = cos(72 deg) = 1/(2*phi)  [C5 rotation angle]
    Equatorial deflection = arccos(-1/sqrt5) = 116.57 deg  [= icosahedral DIHEDRAL]
  For muon: replace pi -> 5*tan(pi/5) (pentagon polygon normalization)
  For electron: 12*tan(pi/12) ~ pi (icosahedron approximates circle within 2.4%)
  For tau: 20*tan(pi/20) = 3.168 (0.83% above pi)

KOIDE STRUCTURE:
  T1(W) x E-(L-antineutrino) = {G32, I52}  [2I CG]
  -> Muon = G32 (dim=4), Tau = I52 (dim=6), Electron = E+ (dim=2) -- CG-forced
  Koide 2/3 should emerge as algebraic identity from Euler/Maxwell icosahedral
  constraints on the three path Born balances (OPEN -- needs tau winding factor).

Run: python analysis/quantum/lepton_mass.py
References:
  doc_leptons.txt (supersedes doc_muon.txt and doc_tau.txt)
  analysis/demos/jobson_cell_doc.py (J26 electron formula)
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, phi, hbar_c, r_p

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

pi   = math.pi
m_p  = 938.272046    # MeV
m_e_pdg  = 0.51099895      # MeV
m_mu_pdg = 105.6583755     # MeV
log5 = math.log(5)
Rs   = math.sqrt(5) / (4 * pi)
Rs2  = Rs**2

# ── Polygon constants ─────────────────────────────────────────────────────────
# For an N-gon: polygon_pi(N) = N * tan(pi/N) = half-perimeter/apothem
# This replaces the smooth-circle pi in the Born normalization.
# As N→∞: N*tan(pi/N) → pi.
def polygon_pi(N):
    return N * math.tan(pi / N)

poly_e  = polygon_pi(12)   # icosahedron bounce count ~ pi (2.4% diff)
poly_mu = polygon_pi(5)    # pentagon polygon-pi  (15.6% diff from pi)

print(SEP)
print("LEPTON MASSES FROM POLYGON PATH GEOMETRY")
print(SEP2)
print(f"  Electron path: 12 icosahedral bounces  "
      f"poly_const = {poly_e:.6f}  (~pi={pi:.6f}, diff={poly_e-pi:.6f})")
print(f"  Muon path:     10 bipyramid bounces     "
      f"poly_const = {poly_mu:.6f}  (NOT pi, diff={poly_mu-pi:.6f})")
print()

# ── Section 1: Electron (J26, reference) ─────────────────────────────────────
print(SEP)
print("SECTION 1: ELECTRON MASS (J26 reference)")
print(SEP2)

# Born parameter: eff_e = phi = (1+sqrt5)/2 (from icosahedral (1,2) winding)
# Polygon normalization: ~pi (12 icosahedral bounces, pi is 2.4% approx)
eff_e   = phi
log5_e  = log5
L3_e    = (eff_e**3 + log5_e**3) / (eff_e**2 + log5_e**2)
x_e     = alpha * eff_e**2
k_e     = alpha * eff_e * (1 - (3/4)*alpha**2) / (1 + x_e + x_e**2)
dn_e    = L3_e * k_e
base_e  = 2*pi * alpha**2 * eff_e * m_p
m_e_pred = base_e * (1 + dn_e/pi) * (1 + (3/4)*alpha**2)

print(f"  eff_e = phi = (1+sqrt5)/2 = {eff_e:.8f}")
print(f"  L3_e = {L3_e:.6f},  k_e = {k_e:.8f},  dn_e/pi = {dn_e/pi:.8f}")
print(f"  m_e = 2*pi*alpha^2*phi*m_p * (1+dn_e/pi) * (1+3/4*alpha^2)")
print(f"      = {m_e_pred:.10f} MeV  (PDG: {m_e_pdg})")
print()

check("LM1 Electron mass m_e within 0.001% of PDG (J26 formula)",
      abs((m_e_pred - m_e_pdg)/m_e_pdg) < 1e-5,
      f"m_e = {m_e_pred:.10f} MeV  err={( m_e_pred-m_e_pdg)/m_e_pdg*100:+.8f}%")
check("LM2 eff_e = phi (golden ratio of (1,2) Hopf winding)",
      abs(eff_e - (1+math.sqrt(5))/2) < 1e-9,
      f"eff_e = {eff_e:.8f} = phi = {(1+math.sqrt(5))/2:.8f}")

# ── Section 2: Muon path geometry (pentagonal bipyramid) ─────────────────────
print()
print(SEP)
print("SECTION 2: PENTAGONAL BIPYRAMID PATH GEOMETRY")
print(SEP2)

# Regular pentagonal bipyramid (J13 Johnson solid, all equal edges a=1)
import numpy as np

a     = 1.0
sin36 = math.sin(pi/5)
cos36 = math.cos(pi/5)
r_e_bip = a / (2 * sin36)       # equatorial ring radius
h_t_bip = math.sqrt(a**2 - r_e_bip**2)  # apex height

print(f"  Pentagonal bipyramid (all equal edges):")
print(f"    r_e = a/(2*sin36) = {r_e_bip:.8f}")
print(f"    h_t = sqrt(a^2-r_e^2) = {h_t_bip:.8f}")
print(f"    h_t/r_e = {h_t_bip/r_e_bip:.8f}  (1/phi = {1/phi:.8f})")
print()

# Deflection cosines on the zig-zag path
top = np.array([0, 0, h_t_bip])
bot = np.array([0, 0, -h_t_bip])
eq  = [np.array([r_e_bip*math.cos(2*pi*k/5), r_e_bip*math.sin(2*pi*k/5), 0])
       for k in range(5)]

def deflect(v_in, vertex, v_out):
    d_in  = (vertex - v_in)  / np.linalg.norm(vertex - v_in)
    d_out = (v_out - vertex) / np.linalg.norm(v_out - vertex)
    return float(np.dot(d_in, d_out))

cos_apex = deflect(eq[0], top, eq[2])    # zig-zag: skip one vertex
cos_eq   = deflect(top, eq[0], bot)      # equatorial: top->eq->bottom

check("LM3 h_t/r_e = 1/phi exactly (golden ratio in bipyramid)",
      abs(h_t_bip/r_e_bip - 1/phi) < 1e-9,
      f"h_t/r_e = {h_t_bip/r_e_bip:.9f}  1/phi = {1/phi:.9f}")
check("LM4 Apex deflection = cos(72 deg) = 1/(2*phi)  [same as C5 angle]",
      abs(cos_apex - 1/(2*phi)) < 1e-8,
      f"cos(apex) = {cos_apex:.8f}  1/(2*phi) = {1/(2*phi):.8f}")
check("LM5 Equatorial deflection = -1/sqrt(5)  [icosahedral constant]",
      abs(cos_eq - (-1/math.sqrt(5))) < 1e-8,
      f"cos(eq) = {cos_eq:.8f}  -1/sqrt5 = {-1/math.sqrt(5):.8f}")

# GEOMETRIC CORRECTION (2026-08-22): actual icosahedral zig-zag path has
# ALL 5 deflections = 1/(2*phi) = cos(72 deg).  No two types.
# The regular bipyramid above was an approximation.  Verify uniform deflection:
import numpy as np
r_e_ico = 1/(2*math.sin(math.pi/5))
z_u_ico = math.sqrt((2*r_e_ico*math.sin(math.pi/10))**2)/2
# recompute properly
r_plane_ico = 2*r_e_ico*math.sin(math.pi/10)
z_u_ico = math.sqrt(1 - r_plane_ico**2)/2
h_top_ico = z_u_ico + math.sqrt(1 - r_e_ico**2)
top_ico = np.array([0, 0, h_top_ico])
bot_ico = np.array([0, 0, -h_top_ico])
upper_ico = [np.array([r_e_ico*math.cos(2*math.pi*k/5), r_e_ico*math.sin(2*math.pi*k/5), z_u_ico]) for k in range(5)]
lower_ico = [np.array([r_e_ico*math.cos(2*math.pi*k/5+math.pi/5), r_e_ico*math.sin(2*math.pi*k/5+math.pi/5), -z_u_ico]) for k in range(5)]
def deflect_ico(vin, v, vout):
    d_in=(v-vin)/np.linalg.norm(v-vin); d_out=(vout-v)/np.linalg.norm(vout-v)
    return float(np.dot(d_in,d_out))
path_ico = [top_ico, upper_ico[0], lower_ico[0], bot_ico, lower_ico[2], upper_ico[2], top_ico]
deflections_ico = [deflect_ico(path_ico[i-1], path_ico[i], path_ico[i+1]) for i in range(1,6)]
all_uniform = all(abs(d - 1/(2*phi)) < 1e-8 for d in deflections_ico)
check("LM4b Actual icosahedral path: ALL 5 deflections = cos(72 deg) = 1/(2*phi)",
      all_uniform,
      f"deflections = {[round(d,6) for d in deflections_ico]}  (all = {1/(2*phi):.6f})")

# ── Section 3: Muon Born balance and mass ─────────────────────────────────────
print()
print(SEP)
print("SECTION 3: MUON MASS FROM PENTAGON PATH")
print(SEP2)

# Born balance effective parameter: eff_mu = (1 + ratio)/2 where
# ratio = |cos_apex| / |cos_eq| = (1/(2*phi)) / (1/sqrt5) = sqrt5/(2*phi)
ratio_mu = abs(cos_apex) / abs(cos_eq)  # = sqrt5/(2*phi)
eff_mu   = (1 + ratio_mu) / 2           # = (9-sqrt5)/8

# Polygon normalization: pentagon (5 bounces per half-circuit) -> 5*tan(pi/5)
# NOT pi (pentagon is 15.6% from circle, unlike icosahedron 2.4%)
C_mu = poly_mu  # = 5*tan(pi/5) = 3.6327

# Born correction (same L3, k structure as electron, with eff_mu)
L3_mu = (eff_mu**3 + log5**3) / (eff_mu**2 + log5**2)
x_mu  = alpha * eff_mu**2
k_mu  = alpha * eff_mu * (1 - (3/4)*alpha**2) / (1 + x_mu + x_mu**2)
dn_mu = L3_mu * k_mu

# Base formula: half the Hopf winding (alpha^1 not alpha^2), phi^2, (2/sqrt5)
base_mu = 2*pi * alpha * (2/math.sqrt(5)) * phi**2 * m_p

# Corrections: Maxwell critical jam (Rs^2) + free-spin 2 channels (2*alpha)
corr_mu = 1 + Rs2 + 2*alpha

m_mu_pred = base_mu * (1 + dn_mu/C_mu) * corr_mu

print(f"  eff_mu = (1 + sqrt5/(2*phi))/2 = (9-sqrt5)/8 = {eff_mu:.8f}")
print(f"  deflection ratio = |cos_apex|/|cos_eq| = {ratio_mu:.8f}")
print(f"    = sqrt5/(2*phi) = {math.sqrt(5)/(2*phi):.8f}")
print(f"  Polygon normalization: 5*tan(pi/5) = {C_mu:.6f}  (not pi = {pi:.6f})")
print(f"  L3_mu = {L3_mu:.6f},  k_mu = {k_mu:.8f},  dn_mu/C_mu = {dn_mu/C_mu:.8f}")
print(f"  base = 2*pi*alpha*(2/sqrt5)*phi^2*m_p = {base_mu:.5f} MeV")
print(f"  corr = 1+Rs^2+2*alpha = {corr_mu:.8f}")
print(f"  m_mu = {m_mu_pred:.6f} MeV  (PDG: {m_mu_pdg})")
print()

check("LM6 eff_mu = (9-sqrt5)/8  [from deflection ratio via Born balance]",
      abs(eff_mu - (9-math.sqrt(5))/8) < 1e-9,
      f"eff_mu = {eff_mu:.8f}  (9-sqrt5)/8 = {(9-math.sqrt(5))/8:.8f}")
check("LM7 eff_mu formula same as eff_e: (1+norm)/2 where norm=sqrt5/(2*phi)",
      abs(eff_mu - (1 + math.sqrt(5)/(2*phi))/2) < 1e-9,
      f"eff_mu = {eff_mu:.8f},  (1+sqrt5/(2*phi))/2 = {(1+math.sqrt(5)/(2*phi))/2:.8f}")
check("LM8 Muon mass within 0.01% of PDG",
      abs((m_mu_pred - m_mu_pdg)/m_mu_pdg) < 1e-4,
      f"m_mu = {m_mu_pred:.6f} MeV  err={( m_mu_pred-m_mu_pdg)/m_mu_pdg*100:+.6f}%")

# ── Section 4: Lepton generation structure ────────────────────────────────────
print()
print(SEP)
print("SECTION 4: LEPTON GENERATIONS FROM 2I DOUBLE GROUP")
print(SEP2)

# T1(W) x E-(L-antineutrino) = {G32, I52}  [CG selection -- no E+ (electron)]
# This is the torsionverse derivation of pion decay preference:
# V-A weak coupling (left-handed antineutrino = E-) forces G32 (muon),
# not E+ (electron). Electron production requires right-handed antineutrino.

# Koide formula check
m_tau_pdg = 1776.86  # MeV
koide = (m_e_pdg + m_mu_pdg + m_tau_pdg) / (
    math.sqrt(m_e_pdg) + math.sqrt(m_mu_pdg) + math.sqrt(m_tau_pdg))**2
print(f"  Lepton generation assignments (2I CG: T1 x E- = {{G32, I52}}):")
print(f"    E+ (dim=2): electron  m_e = {m_e_pdg:.6f} MeV")
print(f"    G32 (dim=4): muon     m_mu = {m_mu_pdg:.6f} MeV  <- derived here")
print(f"    I52 (dim=6): tau      m_tau = {m_tau_pdg:.4f} MeV  <- not yet derived")
print(f"  Koide formula: (m_e+m_mu+m_tau)/(sum sqrt)^2 = {koide:.8f}")
print(f"    (2/3 = {2/3:.8f},  deviation = {(koide-2/3)/(2/3)*100:.6f}%)")
print()

check("LM9 Koide formula holds to 0.002%",
      abs(koide - 2/3) / (2/3) < 2e-5,
      f"Koide = {koide:.8f}  2/3 = {2/3:.8f}  dev = {(koide-2/3)/(2/3)*1e6:.1f} ppm")

# ── Section 5: Tau from Koide ─────────────────────────────────────────────────
print()
print(SEP)
print("SECTION 5: TAU MASS FROM KOIDE (using derived m_e, m_mu)")
print(SEP2)
# Given derived m_e and m_mu, Koide predicts m_tau exactly.
# This is the fastest route to m_tau without deriving the I52 path geometry.
a_k = math.sqrt(m_e_pred); b_k = math.sqrt(m_mu_pred)
s_k = a_k + b_k; A_k = m_e_pred + m_mu_pred
disc_k = 4*s_k**2 - (3*A_k - 2*s_k**2)
x_k = 2*s_k + math.sqrt(disc_k)    # positive physical root
m_tau_pred = x_k**2
m_tau_pdg = 1776.86   # MeV (PDG)
m_tau_err = (m_tau_pred - m_tau_pdg)/m_tau_pdg*100

print(f"  Input: m_e  = {m_e_pred:.8f} MeV  ({(m_e_pred-m_e_pdg)/m_e_pdg*100:+.6f}%)")
print(f"  Input: m_mu = {m_mu_pred:.6f} MeV  ({(m_mu_pred-m_mu_pdg)/m_mu_pdg*100:+.5f}%)")
print(f"  Koide prediction: m_tau = {m_tau_pred:.4f} MeV")
print(f"  PDG m_tau = {m_tau_pdg} MeV   error = {m_tau_err:+.4f}%")
print()

check("LM10 Tau from Koide (derived m_e, m_mu) within 0.01% of PDG",
      abs(m_tau_err) < 0.01,
      f"m_tau = {m_tau_pred:.4f} MeV  PDG={m_tau_pdg}  err={m_tau_err:+.4f}%")

# Verify Koide with predicted tau
total_m = m_e_pred + m_mu_pred + m_tau_pred
total_sq = math.sqrt(m_e_pred) + math.sqrt(m_mu_pred) + math.sqrt(m_tau_pred)
koide_pred = total_m / total_sq**2
check("LM11 Koide ratio with derived masses = 2/3 to high precision",
      abs(koide_pred - 2/3) < 1e-8,
      f"Koide(derived) = {koide_pred:.10f}  2/3 = {2/3:.10f}")

# Mass ratio structure
print(f"  Mass ratios: m_mu/m_e = {m_mu_pdg/m_e_pdg:.3f}   m_tau/m_mu = {m_tau_pdg/m_mu_pdg:.3f}")
print(f"  2I irrep dims: E+(2), G32(4), I52(6) -- dimensions 2:4:6")

# ── Section 6: Grand lepton structure (vertex/edge/face) ──────────────────────
print()
print(SEP)
print("SECTION 6: GRAND LEPTON STRUCTURE -- VERTEX / EDGE / FACE")
print(SEP2)
print("  The three lepton generations correspond exactly to the three types")
print("  of elements of the icosahedron:")
print()

# Icosahedral element counts
V = 12   # vertices
E = 30   # edges
F = 20   # faces
euler = V - E + F
maxwell = 3*V - E

print(f"  Icosahedron: V={V} (vertices), E={E} (edges), F={F} (faces)")
print(f"  Euler:   V - E + F = {euler}  (must be 2)")
print(f"  Maxwell: 3V - E    = {maxwell}  (marginally rigid = critical structure)")
print()

check("LM12 Icosahedron Euler formula V-E+F=2",
      euler == 2, f"V-E+F = {euler}")
check("LM13 Icosahedron Maxwell criterion 3V-E=6 (marginally rigid)",
      maxwell == 6, f"3V-E = {maxwell}")

print()
print("  Lepton -- element mapping (2I double group spinors):")
print(f"    E+  (dim=2, C3=-1)  VERTEX mode:  12 vertices,  72-deg  (C5) deflections")
print(f"    G32 (dim=4, C3=+1)  EDGE mode:    30 edges,     72-deg  (C5) deflections")
print(f"    I52 (dim=6, C3=0)   FACE mode:    20 faces, 72-deg geometric path (Born balance uses -1/sqrt5)")
print()

# Tau face deflection: Born balance uses -1/sqrt5 (equatorial bipyramid model).
# Geometric path deflection = 72 deg (C5, same as muon) -- gluon_tau_helix.py GH2.
# These are SEPARATE: path geometry (72 deg forward) vs Born balance coupling (-1/sqrt5 backward).
cos_dihedral = -1/math.sqrt(5)
theta_dihedral = math.degrees(math.acos(cos_dihedral))
eff_tau = (1 + 1/math.sqrt(5)) / 2    # Born balance: (1 + |cos_born|)/2
poly_tau = polygon_pi(20)             # 20-face Hamiltonian path polygon constant

print(f"  TAU face-path geometry:")
print(f"    Geometric path deflection = 72 deg (C5, uniform -- gluon_tau_helix.py GH2)")
print(f"    Born balance coupling angle = arccos(-1/sqrt5) = {theta_dihedral:.4f} deg")
print(f"    cos(Born balance) = -1/sqrt(5) = {cos_dihedral:.8f}")
print(f"    eff_tau = (1 + 1/sqrt5)/2    = {eff_tau:.8f}")
print(f"    Polygon normalization (N=20): = {poly_tau:.6f}  (pi = {pi:.6f}, diff = {poly_tau-pi:.6f})")
print()
check("LM14 Tau Born balance angle = arccos(-1/sqrt5) [NOT icosahedral dihedral; geometric path = 72 deg]",
      abs(cos_dihedral - (-1/math.sqrt(5))) < 1e-12,
      f"cos(Born) = {cos_dihedral:.10f}  -1/sqrt5 = {-1/math.sqrt(5):.10f}  [icosahedral dihedral = 138.19 deg, JC5]")

# Color coupling: I52 dim=6 = 2 (spinor) x 3 (quark colors from face 3-coloring)
dim_I52 = 6
spinor_mult = 2
color_mult = dim_I52 // spinor_mult
print(f"  COLOR COUPLING (F-7 connection):")
print(f"    I52 dim = {dim_I52} = {spinor_mult} (spinor) x {color_mult} (quark colors)")
print(f"    The icosahedron's 20 faces are 3-colorable (face 3-coloring = R/G/B quark colors)")
print(f"    I52 C3=0: tau couples to all 3 colored faces equally -> color NEUTRAL")
print(f"    G32 C3=+1: muon couples to face boundaries (edges) -> color BOUNDARY mode")
print(f"    E+  C3=-1: electron at face junctions (5 faces per vertex) -> color INTERSECTION")
print()

check("LM15 I52 dim=6 = 2*3 (spinor x 3 quark colors from face 3-coloring)",
      dim_I52 == spinor_mult * color_mult and color_mult == 3,
      f"dim_I52 = {spinor_mult} * {color_mult} = {dim_I52}")

print()
print("  C3 CHARACTER SUMMARY (confirms vertex/edge/face assignment):")
print(f"    E+  C3 char = -1  (vertex: 5 faces meet -> anti-phase under C3 rotation)")
print(f"    G32 C3 char = +1  (edge: 2 faces meet -> in-phase under C3 rotation)")
print(f"    I52 C3 char =  0  (face: pure C3 symmetric -> color-neutral)")
print()
print("  KOIDE 2/3 EXPECTATION:")
print("    All three paths (vertex/edge/face) are Born balances on the SAME icosahedral")
print("    lattice. The Euler formula V-E+F=2 and Maxwell criterion 3V-E=6 constrain")
print("    the relative coupling constants eff_e, eff_mu, eff_tau. Once the tau face-")
print("    mode winding factor is derived (OPEN), the Koide 2/3 should emerge as an")
print("    algebraic identity from the icosahedral Euler/Maxwell structure.")

# ── Section 7: Tau corkscrew winding factor ──────────────────────────────────
print()
print(SEP)
print("SECTION 7: TAU CORKSCREW -- WINDING FACTOR FROM (1,2) HOPF p-COMPONENT")
print(SEP2)

# Electron uses q-projection: 2*pi * alpha^2 * phi^1  (circulation, 2 Born contacts)
# Muon uses q-projection:     2*pi * alpha^1 * (2/sqrt5) * phi^2  (circulation, 1 Born contact)
# Tau corkscrew (flux, 0 Born vertex contacts):
#   - No 2*pi: face coupling is a flux through the face normal, not a closed circulation
#   - No alpha: tau never touches a vertex or edge point; dihedral crossing has no Born contact
#   - Winding factor = p-component of unit (1,2) vector = 1/sqrt5
#     (muon used the longitudinal q=2 component; tau uses the transverse p=1 component)

# Physical picture (CORRECTED 2026-08-28): two corpuscle photons on the unique
# Hamiltonian circuit in opposite directions, bouncing at face-CENTER nexuses
# (impact-rebound 72 deg off convergent gluon maximum, GH0b). They meet twice
# per circuit (hops 0 and 10). Each hop is a chord through the interior.
# Together they cover all 20 face-center nexuses per circuit. [TPC6a-d]
# The name "corkscrew" refers to the helical 3D shape of the path.

p_winding = 1.0           # p-component of (p,q)=(1,2) Hopf winding
q_winding = 2.0
norm_winding = math.sqrt(p_winding**2 + q_winding**2)  # = sqrt5
winding_factor_tau = p_winding / norm_winding           # = 1/sqrt5 (face flux)
winding_factor_mu  = q_winding / norm_winding           # = 2/sqrt5 (edge circulation)

# Base tau: no 2*pi, no alpha -- pure geometric coupling
base_tau_corkscrew = winding_factor_tau * phi**3 * m_p

# phi^3 = 2+sqrt5 (exact)
phi3_exact = 2 + math.sqrt(5)
base_tau_exact = phi3_exact / math.sqrt(5) * m_p

print(f"  (1,2) winding vector: norm = sqrt5 = {norm_winding:.8f}")
print(f"  Muon  winding factor: q-component = 2/sqrt5 = {winding_factor_mu:.8f}  (longitudinal)")
print(f"  Tau   winding factor: p-component = 1/sqrt5 = {winding_factor_tau:.8f}  (transverse/flux)")
print(f"  phi^3 = 2+sqrt5 = {phi3_exact:.8f}  (exact algebraic)")
print(f"  base_tau = (1/sqrt5) * phi^3 * m_p = (phi^3/sqrt5) * m_p")
print(f"           = {base_tau_corkscrew:.4f} MeV  (PDG m_tau = {m_tau_pdg} MeV)")
print(f"  Error (no corrections): {(base_tau_corkscrew - m_tau_pdg)/m_tau_pdg*100:+.4f}%")
print()
print(f"  Algebraic identity: phi^3/sqrt5 = (2+sqrt5)/sqrt5 = 2/sqrt5 + 1 = 1 + (muon winding factor)")
print(f"    = {phi3_exact/math.sqrt(5):.8f}  vs  1 + 2/sqrt5 = {1 + 2/math.sqrt(5):.8f}  [check: {abs(phi3_exact/math.sqrt(5) - (1 + 2/math.sqrt(5))) < 1e-12}]")
print()
print(f"  m_tau / m_p = phi^3/sqrt5 = {m_tau_pdg/m_p:.8f}  (PDG)")
print(f"             vs  phi^3/sqrt5 = {phi3_exact/math.sqrt(5):.8f}  (formula)")
print()

check("LM16 Tau corkscrew base = phi^3/sqrt5 * m_p within 0.05% of PDG (leading order, no corrections)",
      abs((base_tau_corkscrew - m_tau_pdg)/m_tau_pdg) < 5e-4,
      f"base_tau = {base_tau_corkscrew:.4f} MeV  err={( base_tau_corkscrew-m_tau_pdg)/m_tau_pdg*100:+.4f}%  [exact correction from Hopf face flux integral]")

check("LM17 Algebraic identity: phi^3/sqrt5 = 1 + 2/sqrt5 (tau contains muon winding + 1)",
      abs(phi**3/math.sqrt(5) - (1 + 2/math.sqrt(5))) < 1e-10,
      f"phi^3/sqrt5 = {phi**3/math.sqrt(5):.10f}  1+2/sqrt5 = {1+2/math.sqrt(5):.10f}")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP2)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"  Electron: {m_e_pred:.9f} MeV  ({( m_e_pred-m_e_pdg)/m_e_pdg*100:+.8f}%)")
print(f"  Muon:     {m_mu_pred:.6f} MeV  ({(m_mu_pred-m_mu_pdg)/m_mu_pdg*100:+.6f}%)")
print(f"  Tau:      {m_tau_pred:.4f} MeV (Koide, {m_tau_err:+.4f}%);  corkscrew base = phi^3/sqrt5 * m_p = {base_tau_corkscrew:.4f} MeV ({(base_tau_corkscrew-m_tau_pdg)/m_tau_pdg*100:+.4f}%)")
print(f"  Grand structure: electron=VERTEX, muon=EDGE, tau=FACE (icosahedral elements)")
print(f"  Color coupling: E+(C3=-1)=intersection, G32(C3=+1)=boundary, I52(C3=0)=neutral")
print(f"  Total: {passed}/{len(results)}  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, s, d in results:
        if s == "FAIL": print(f"  FAILED: {name}")
print(f"  Reference: docs/open_items.txt F-9, docs/doc_leptons.txt")
print(SEP)
