"""
higgs_bond_geometry.py
======================
Tests the hypothesis: is the Higgs boson the bond/excitation between
adjacent Jobson cells?

Key geometric question: what is the distance between adjacent
icosahedral cell centers, and how does the Higgs Compton wavelength
compare to the inter-cell geometry?

Three types of adjacency for icosahedral cells:
  (A) Edge-sharing: two cells share one edge (2 vertices)
  (B) Face-sharing: two cells share one face (3 vertices, 3 edges)
  (C) Vertex-sharing: two cells share one vertex

Run: python analysis/higgs/higgs_bond_geometry.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2

SEP  = "=" * 65
SEP2 = "-" * 65

print(SEP)
print("JOBSON CELL BOND GEOMETRY")
print("Testing Higgs-as-inter-cell-bond hypothesis")
print(SEP)
print()

# ── Icosahedral geometry (edge = 2, vertices at (0,±1,±phi) etc.) ─────────────
# Standard icosahedron: 12 vertices at (0, ±1, ±phi) and cyclic permutations.
# Edge = 2 (distance between adjacent vertices), scales to edge = L_J by * L_J/2.

edge_model = 2.0   # model edge length (vertices at (0, ±1, ±phi) and permutations)

# Circumradius: distance from center to vertex
# Vertex (0, 1, phi): distance = sqrt(1 + phi^2) = sqrt(phi+2)
circumradius_model = math.sqrt(1 + phi**2)
circumradius = circumradius_model * L_J / edge_model

# Face inradius: distance from center to face center
# Face containing (0,1,phi), (1,phi,0), (phi,0,1)
# Face center = ((0+1+phi)/3, (1+phi+0)/3, (phi+0+1)/3) = ((1+phi)/3, (1+phi)/3, (1+phi)/3)
# = (phi^2/3, phi^2/3, phi^2/3)  since 1+phi = phi^2
face_center_model = phi**2 / 3
face_inradius_model = math.sqrt(3) * face_center_model  # distance from origin
face_inradius = face_inradius_model * L_J / edge_model

# Edge midpoint inradius: distance from center to edge midpoint
# Edge between (0,1,phi) and (0,-1,phi): midpoint = (0,0,phi)
# Distance from origin = phi
edge_mid_inradius_model = phi
edge_mid_inradius = edge_mid_inradius_model * L_J / edge_model

print("ICOSAHEDRAL CELL GEOMETRY (edge = L_J)")
print(SEP2)
print(f"  L_J              = {L_J:.6f} fm = {L_J*1000:.4f} am")
print(f"  Circumradius R   = sqrt(phi^2+1)/2 * L_J = {circumradius:.6f} fm")
print(f"                   = {circumradius/L_J:.6f} * L_J")
print(f"                   = sqrt(phi+2)/2 * L_J  [phi+2 = {phi+2:.6f}]")
print(f"  Face inradius    = phi^2/sqrt(3) * L_J/2 = {face_inradius:.6f} fm")
print(f"                   = {face_inradius/L_J:.6f} * L_J")
print(f"  Edge-mid inradius = phi/2 * L_J = {edge_mid_inradius:.6f} fm")
print(f"                   = {edge_mid_inradius/L_J:.6f} * L_J")
print()

# ── Center-to-center distances for each adjacency type ───────────────────────
print("CENTER-TO-CENTER DISTANCES (cells touching)")
print(SEP2)
print()

# (A) Edge-sharing: centers on opposite sides of shared edge midpoint
# c-to-c = 2 * edge_mid_inradius = phi * L_J  (EXACT)
d_edge = 2 * edge_mid_inradius
print(f"(A) Edge-sharing cells:")
print(f"    center-to-center = 2 * phi/2 * L_J = phi * L_J")
print(f"                     = {d_edge/L_J:.8f} * L_J")
print(f"                     = {d_edge:.6f} fm = {d_edge*1000:.4f} am")
print(f"    phi = {phi:.8f}")
print(f"    Exact: d = phi * L_J  [no approximation -- pure icosahedral geometry]")
print()

# (B) Face-sharing: centers on opposite sides of shared face
# c-to-c = 2 * face_inradius = phi^2/sqrt(3) * L_J
d_face = 2 * face_inradius
print(f"(B) Face-sharing cells:")
print(f"    center-to-center = 2 * phi^2/sqrt(3)/2 * L_J = phi^2/sqrt(3) * L_J")
print(f"                     = {d_face/L_J:.8f} * L_J")
print(f"                     = {d_face:.6f} fm")
print()

# (C) Vertex-sharing: centers separated by 2*circumradius (touching at vertex)
d_vertex = 2 * circumradius
print(f"(C) Vertex-sharing cells:")
print(f"    center-to-center = 2 * sqrt(phi+2)/2 * L_J = sqrt(phi+2) * L_J")
print(f"                     = {d_vertex/L_J:.8f} * L_J")
print(f"                     = {d_vertex:.6f} fm")
print()

# ── Higgs Compton wavelength vs cell geometry ─────────────────────────────────
hbar_c_fm = hbar_c  # MeV*fm
m_H_MeV   = m_H_pdg22 * 1000  # MeV
lambda_H  = hbar_c_fm / m_H_MeV  # fm = Compton wavelength of Higgs

print("HIGGS COMPTON WAVELENGTH vs CELL GEOMETRY")
print(SEP2)
print()
print(f"  lambda_H = hbar*c / m_H = {lambda_H:.8f} fm = {lambda_H*1000:.4f} am")
print(f"           = {lambda_H/L_J:.8f} * L_J")
print(f"           = {lambda_H/L_J:.8f}  [should be ~ 1/(2*pi) if m_H = E_cell]")
print(f"  1/(2*pi) = {1/(2*pi):.8f}")
print(f"  Deviation from 1/(2*pi): {abs(lambda_H/L_J - 1/(2*pi))/(1/(2*pi))*100:.4f}%")
print()
print(f"  lambda_H / d_edge  = lambda_H / (phi * L_J) = {lambda_H/d_edge:.8f}")
print(f"                     = 1/(2*pi*phi)           = {1/(2*pi*phi):.8f}")
print(f"  The Higgs fills 1/(2*pi*phi) ~ {1/(2*pi*phi)*100:.2f}% of the edge-sharing bond")
print()

# ── The inter-cell gap hypothesis ─────────────────────────────────────────────
print("INTER-CELL GAP HYPOTHESIS")
print(SEP2)
print()
print("  Hypothesis: the equilibrium gap between adjacent cells = lambda_H")
print("  Physical reasoning: the inter-cell mode zero-point energy = m_H*c^2 = E_cell")
print("  The gap is set by the Higgs zero-point oscillation -- not zero, not arbitrary.")
print()
print(f"  Predicted gap = lambda_H = {lambda_H:.8f} fm = {lambda_H*1000:.4f} am")
print(f"               = L_J / (2*pi) = {L_J/(2*pi):.8f} fm")
print(f"  Gap / L_J = {lambda_H/L_J:.6f} = 1/(2*pi)")
print()
print("  The inter-cell gap is L_J/(2*pi) -- one radian of the cell circumference.")
print("  This is NOT zero. The medium is NOT space-filling in the rigid sense.")
print("  The cells are packed with a thin elastic film of thickness L_J/(2*pi).")
print()

# ── Can we corroborate with W and Z masses? ───────────────────────────────────
print("W AND Z MASSES vs CELL GEOMETRY")
print(SEP2)
print()
m_W = 80400  # MeV
m_Z = 91200  # MeV

lambda_W = hbar_c_fm / m_W
lambda_Z = hbar_c_fm / m_Z

print(f"  lambda_W = {lambda_W/L_J:.6f} * L_J  (W boson Compton wavelength)")
print(f"  lambda_Z = {lambda_Z/L_J:.6f} * L_J  (Z boson Compton wavelength)")
print(f"  lambda_H = {lambda_H/L_J:.6f} * L_J")
print()
print("  If Higgs = scalar bond mode (breathing), W/Z = transverse bond modes:")
print(f"  lambda_W / lambda_H = m_H / m_W = {m_H_MeV/m_W:.6f}")
print(f"  phi/2               =             {phi/2:.6f}  (ratio {abs(m_H_MeV/m_W - phi/2)/(phi/2)*100:.2f}% off)")
print(f"  pi/2                =             {pi/2:.6f}  (ratio {abs(m_H_MeV/m_W - pi/2)/(pi/2)*100:.2f}% off)")
print(f"  sqrt(phi^2+1)/phi   =             {math.sqrt(phi**2+1)/phi:.6f}")
print()
print(f"  lambda_Z / lambda_H = m_H / m_Z = {m_H_MeV/m_Z:.6f}")
print(f"  phi/2               =             {phi/2:.6f}  ({abs(m_H_MeV/m_Z - phi/2)/(phi/2)*100:.2f}% off)")
print(f"  Closest I_h number: {min([phi/2, pi/2, math.sqrt(phi**2+1)/phi, phi-0.3], key=lambda x: abs(m_H_MeV/m_Z-x)):.4f}")
print()

# ── Higgs as breathing mode vs bond mode ──────────────────────────────────────
print("BREATHING MODE vs BOND MODE DISTINCTION")
print(SEP2)
print()
print("  BREATHING MODE picture (intra-cell):")
print("    The Higgs = A_g (l=0) excitation of the cell itself.")
print("    All 12 vertices move radially in/out.")
print("    Energy = E_cell, spin = 0 (l=0 mode), CP-even (gerade).")
print("    This is consistent with ALL observed Higgs quantum numbers.")
print()
print("  BOND MODE picture (inter-cell):")
print("    The Higgs = excitation of the elastic film between edge-sharing cells.")
print("    Film thickness = L_J/(2*pi), bond length = phi * L_J.")
print("    The film is compressed when cells are pushed together.")
print("    This picture explains WHY the Higgs is sub-cell:")
print("      -- the bond is between cells, so the Higgs lives OUTSIDE any single cell")
print("      -- its Compton wavelength = gap thickness = L_J/(2*pi) < L_J")
print()
print("  KEY DISTINCTION: these two pictures may be DUAL.")
print("    The breathing mode OF the cell IS the compression of adjacent bonds.")
print("    When a cell breathes (vertices move out), all 12 edges to neighbors compress.")
print("    The number of edges per vertex in the lattice...")
icosahedron_edges = 30
icosahedron_vertices = 12
edges_per_vertex = icosahedron_edges * 2 / icosahedron_vertices  # each edge shared by 2 vertices
print(f"    Icosahedron: {icosahedron_edges} edges, {icosahedron_vertices} vertices")
print(f"    Edges per vertex: {edges_per_vertex:.1f}")
print(f"    When cell breathes, {edges_per_vertex:.0f} bonds are simultaneously compressed.")
print()

# ── The phi * L_J result corroboration ────────────────────────────────────────
print("CORROBORATION: phi * L_J as the natural INTER-CELL scale")
print(SEP2)
print()
d_bond = phi * L_J  # edge-sharing center-to-center
E_bond = hbar_c_fm / d_bond * 1e-3  # GeV  (energy associated with bond length)
print(f"  Bond length (edge-sharing) = phi * L_J = {d_bond:.8f} fm = {d_bond*1000:.4f} am")
print(f"  Energy scale of bond: E_bond = hbar*c / (phi*L_J) = {E_bond:.4f} GeV")
print(f"  E_bond / E_cell = {E_bond/E_cell_GeV:.6f} = 1/phi = {1/phi:.6f}")
print()
print(f"  The bond energy scale is E_cell / phi -- a sub-cell particle from")
print(f"  the bond perspective would have N_J_bond = R_Compton / (phi*L_J)")
print()
m_bond_particle = hbar_c_fm / d_bond  # MeV
print(f"  A particle 'filling' one bond: mass = hbar*c/(phi*L_J) = {m_bond_particle:.1f} MeV")
print(f"  Known particles near {m_bond_particle:.0f} MeV:")
print(f"    tau lepton: 1777 MeV  (ratio {1777/m_bond_particle:.3f})")
print(f"    charm quark: 1280 MeV  (ratio {1280/m_bond_particle:.3f})")
print(f"    proton: 938 MeV  (ratio {938/m_bond_particle:.3f})")
print()

# ── Summary table ─────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()
print(f"  Cell edge L_J                = {L_J:.6f} fm = {L_J*1000:.4f} am")
print(f"  Higgs Compton wavelength     = L_J/(2*pi) = {lambda_H:.6f} fm = {lambda_H*1000:.4f} am  (N_J = 0.159)")
print(f"  Edge-sharing bond length     = phi*L_J    = {d_bond:.6f} fm = {d_bond*1000:.4f} am  (EXACT)")
print(f"  Face-sharing bond length     = phi^2/sqrt(3)*L_J = {d_face:.6f} fm")
print()
print("  KEY RESULT: edge-sharing center-to-center = phi * L_J  (exact, icosahedral)")
print("  The Higgs Compton wavelength is 1/(2*pi*phi) of this bond length.")
print()
print("  PICTURE: The torsion medium has icosahedral cells packed with equilibrium")
print("  inter-cell gaps of thickness lambda_H = L_J/(2*pi). The Higgs boson is the")
print("  quantum of the inter-cell breathing mode. The Higgs FIELD is the collective")
print("  displacement field of those inter-cell gaps -- distinct from the torsion")
print("  medium (the cells themselves), but embedded within it.")
print()
print("  HIGGS FIELD vs TORSION MEDIUM:")
print("    Torsion medium = the icosahedral grains (cells)")
print("    Higgs field    = the inter-cell elastic displacement field (the 'film')")
print("    Higgs boson    = quantum of that film's excitation at energy E_cell")
print()
print("  The Higgs FIELD and the torsion medium FIELD are complementary:")
print("  one is what the cells ARE, the other is how cells INTERACT.")
print(SEP)
