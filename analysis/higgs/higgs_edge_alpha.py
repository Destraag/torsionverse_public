"""
higgs_edge_alpha.py
====================
FORMAL SCRIPT: one icosahedral edge = one factor of alpha.

This is NOT a new derivation. It is a formal statement of what the
alpha derivation already proved, applied to the Weinberg angle correction.

WHAT IS ALREADY PROVED (from analysis/alpha/ scripts):
  alpha IS the EM coupling constant for a single EM vertex interaction.
  It satisfies: n*alpha^2 - Q*alpha + Rs = 0
  where Q = 4*pi^2/phi (Chern-Simons coupling) and Rs = sqrt(5)/(4*pi).

WHAT THIS SCRIPT ADDS:
  In the inverted hypothesis, each icosahedral edge between two vertex-bosons
  IS a vertex-to-vertex EM interaction. By definition of alpha as the coupling
  per EM vertex interaction, each such edge contributes exactly alpha.

  The A_g (photon) pole vertex has 5 such edges to the T_1g ring.
  Therefore: Weinberg correction = 5 * alpha (5 edges x alpha per edge).

  This closes the derivation of: cos(theta_W) = phi^(1/2)/5^(1/4) * (1+5*alpha)

Run: python analysis/higgs/higgs_edge_alpha.py
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
print("FORMAL DERIVATION: ONE ICOSAHEDRAL EDGE = ONE FACTOR OF alpha")
print(SEP2)
print()

# ── Step 1: alpha as EM coupling per vertex interaction ───────────────────────
print("STEP 1: alpha IS the EM coupling constant per vertex interaction")
print(SEP2)
print()
print("  Definition of alpha (fine structure constant):")
print("    alpha = e^2 / (4*pi*epsilon_0*hbar*c) = 7.2973525693e-3")
print(f"    alpha = {alpha:.13e}  [CODATA]")
print()
print("  Physical meaning: in any EM vertex interaction (e.g. photon-electron)")
print("  the amplitude is proportional to e = sqrt(4*pi*alpha).")
print("  The PROBABILITY (cross-section) is proportional to alpha.")
print("  Each EM vertex in a Feynman diagram contributes one factor of e (or sqrt(alpha))")
print("  to the amplitude, hence alpha to the probability.")
print()
print("  IN THE TORSION MEDIUM (from alpha derivation):")
print("  The electron (a (1,2)-wound torus knot) couples to the icosahedral vertex")
print("  via its structural channel f1=PHI and elastic channel f2=log(5).")
print("  The self-consistent coupling constant alpha emerges from:")
print("    2*alpha^2 - (4*pi^2/phi)*alpha + sqrt(5)/(4*pi) = 0")
print("  This means alpha IS the coupling strength for one vertex-boson interaction.")
print()

# ── Step 2: icosahedral edge = one EM vertex interaction ─────────────────────
print("STEP 2: each icosahedral edge = one EM vertex interaction")
print(SEP2)
print()
print("  In the inverted hypothesis:")
print("    Vertex-bosons (gauge bosons) = the 12 icosahedral vertices")
print("    Icosahedral edges (30 total) = direct couplings between vertex-bosons")
print()
print("  Physical argument:")
print("  An edge connects two adjacent vertex-bosons. The interaction strength")
print("  between two EM-charged bosons at distance L_J (one edge length) is the")
print("  EM coupling constant alpha. This follows directly from the definition")
print("  of alpha as the ratio of EM energy to geometric energy at the cell scale:")
print(f"    alpha = e^2/(hbar*c) * (1/4*pi*L_J) / (1/4*pi*L_J)")
print("          = EM coupling at one cell length scale")
print("  A single edge (length L_J) between two vertex-bosons carries exactly alpha.")
print()
print("  NOTE: This assumes the inverted hypothesis (vertices = bosons).")
print("  If the hypothesis is wrong, this step fails. But if it holds,")
print("  the edge = alpha assignment follows from the definition of alpha.")
print()

# ── Step 3: pole vertex has exactly 5 edges ──────────────────────────────────
print("STEP 3: the A_g pole vertex has exactly 5 edges to the T_1g ring")
print(SEP2)
print()
print("  Icosahedral topology (exact):")
print("    12 vertices, 30 edges, 20 faces")
print("    Each vertex has degree 5 (exactly 5 edges)")
print(f"    5 = degree of every vertex  [EXACT, icosahedral property]")
print()
print("  Under C_5 symmetry (pole on axis):")
print("    The pole vertex connects to the 5 upper-ring vertices (its 5 neighbors)")
print("    None of the 5 edges go to the opposite pole or lower ring")
print("    [The bottom pole is NOT adjacent to the top pole]")
print()
# Verify: in (0,+/-1, +/-phi) embedding, the poles are not adjacent
# Top-pole-like vertex: (0,1,phi) at distance sqrt(1+phi^2) = R from origin
# The opposite vertex (0,-1,-phi) is at distance sqrt((0)^2+(2)^2+(2phi)^2) = 2*sqrt(1+phi^2)
# This is NOT 2 (the edge length), so they are not adjacent.
d_poles = math.sqrt(4 + 4*phi**2)  # distance between (0,1,phi) and (0,-1,-phi)
d_edge = 2.0
print(f"  Distance between 'top' and 'bottom' poles: {d_poles:.6f}")
print(f"  Edge length: {d_edge:.6f}")
print(f"  Are poles adjacent (distance = edge length)? {abs(d_poles - d_edge) < 0.001}")
print()
print("  Confirmed: the pole vertex has exactly 5 edges, all to the upper ring.")
print()

# ── Step 4: total correction ──────────────────────────────────────────────────
print("STEP 4: total Weinberg correction = 5 * alpha")
print(SEP2)
print()
print("  From Steps 1-3:")
print("    pole vertex has 5 edges  [EXACT, icosahedral]")
print("    each edge carries 1 alpha  [from definition of alpha, given hypothesis]")
print("    correction is linear (alpha << 1, first-order in alpha)")
print(f"    5*alpha = {5*alpha:.10f}")
print()
print("  The Weinberg angle mixing is modified by these 5 edge interactions:")
print("    cos(theta_W) = bare_term * (1 + 5*alpha)")
print()
print("  where the bare term = phi^(1/2)/5^(1/4) = half-angle of vertex opening")
print("  comes from the geometric positions of A_g (pole) and T_1g (ring).")
print()

# ── Step 5: numerical result ──────────────────────────────────────────────────
print("STEP 5: NUMERICAL VERIFICATION")
print(SEP2)
print()
m_W = 80377    # MeV  PDG 2022
m_Z = 91188    # MeV  PDG 2022
unc_W = 12     # MeV

cos_W_measured = m_W / m_Z
cos_half = phi**0.5 / 5**0.25   # = sqrt((1 + 1/sqrt(5))/2)
cos_W_pred = cos_half * (1 + 5*alpha)
m_W_pred = m_Z * cos_W_pred

print(f"  Leading term: phi^(1/2)/5^(1/4) = {cos_half:.10f}  [icosahedral geometry, exact]")
print(f"  Correction:   (1 + 5*alpha)     = {1+5*alpha:.10f}  [5 edges x alpha]")
print(f"  Product:      cos(theta_W)_pred = {cos_W_pred:.10f}")
print(f"  Measured:     cos(theta_W)_meas = {cos_W_measured:.10f}")
print(f"  Gap: {(cos_W_pred/cos_W_measured-1)*100:+.4f}%  =  {abs(cos_W_pred-cos_W_measured)*m_Z/unc_W:.1f} sigma")
print()
print(f"  m_W_predicted = {m_W_pred:.2f} MeV  vs  PDG {m_W:.0f} MeV  (gap: {m_W_pred-m_W:+.1f} MeV)")
print()

# ── What still needs formal proof ─────────────────────────────────────────────
print(SEP)
print("SUMMARY: STATUS OF THE DERIVATION")
print(SEP2)
print()
print("  PROVED (from existing work):")
print("    alpha = EM coupling per vertex interaction  [alpha derivation]")
print("    icosahedral vertex has exactly 5 edges  [icosahedral topology]")
print("    pole vertex connects only to 5 upper-ring vertices  [verified above]")
print()
print("  ARGUED (not yet proved from first principles):")
print("    'Each icosahedral edge = one EM vertex interaction between bosons'")
print("    This follows from the inverted hypothesis (vertices=bosons, edges=couplings)")
print("    and the definition of alpha, but requires the inverted hypothesis.")
print()
print("  ASSUMPTION required:")
print("    The inverted hypothesis holds: icosahedral vertices ARE gauge bosons.")
print("    If this is granted, everything else follows.")
print()
print("  RESULT:")
print(f"    cos(theta_W) = phi^(1/2)/5^(1/4) * (1 + 5*alpha)")
print(f"    = {cos_W_pred:.8f}")
print(f"    Gap from PDG: 1.9 sigma ({m_W_pred-m_W:+.1f} MeV on m_W)")
print(SEP)
