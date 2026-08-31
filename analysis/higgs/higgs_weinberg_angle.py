"""
higgs_weinberg_angle.py
=======================
Tests the Weinberg angle formula derived from icosahedral vertex geometry.

SETUP (from inverted hypothesis: bosons = vertices, cell = emergent Higgs):
  - Bosons interact ONLY with other vertices (via edges) and with particle waves
  - The Weinberg angle = mixing between A_g vertex (pole) and T_1g vertex (equatorial)
  - This mixing is determined by the ANGLE between vertex positions

ICOSAHEDRAL GEOMETRY:
  The angle between two adjacent vertices (subtended at cell center):
    cos(theta_vertex) = 1/sqrt(5)  [exact]
    theta_vertex = arccos(1/sqrt(5)) = 63.43 deg

  The Weinberg angle ~ half this angle:
    cos(theta_W) ~ sqrt((1 + cos(theta_vertex))/2)  [half-angle formula]
               = sqrt((1 + 1/sqrt(5))/2)
               = 0.8507  (3.5% off measurement)

  CORRECTION: each vertex has 5 edges (5-fold symmetry).
  The same 5-fold vertex count drove f1=PHI in the alpha derivation.
  Correction: (1 + 5*alpha) where 5 = edges per vertex.

  CANDIDATE FORMULA:
    cos(theta_W) = sqrt((1 + 1/sqrt(5))/2) * (1 + 5*alpha)

Run: python analysis/higgs/higgs_weinberg_angle.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi   = math.pi
phi  = (1 + math.sqrt(5)) / 2
Rs   = math.sqrt(5) / (4*pi)

SEP  = "=" * 65
SEP2 = "-" * 65

m_W    = 80377   # MeV  PDG 2022
m_Z    = 91188   # MeV  PDG 2022
unc_W  = 12      # MeV  PDG uncertainty

cos_W_measured = m_W / m_Z
theta_W_measured = math.acos(cos_W_measured)
sin2_W_measured  = math.sin(theta_W_measured)**2

print(SEP)
print("WEINBERG ANGLE FROM ICOSAHEDRAL VERTEX GEOMETRY")
print(SEP2)
print()
print(f"  Measured: m_W/m_Z = {m_W}/{m_Z} = {cos_W_measured:.8f} = cos(theta_W)")
print(f"  theta_W  = {theta_W_measured*180/pi:.4f} deg")
print(f"  sin^2(theta_W) = {sin2_W_measured:.6f}")
print()

# ── Icosahedral vertex geometry ───────────────────────────────────────────────
print("ICOSAHEDRAL VERTEX GEOMETRY")
print(SEP2)
print()
print("  Adjacent vertex opening angle (subtended at cell center):")
cos_vertex = 1/math.sqrt(5)
theta_vertex = math.acos(cos_vertex)
print(f"    cos(theta_vertex) = 1/sqrt(5) = {cos_vertex:.8f}  [EXACT]")
print(f"    theta_vertex      = {theta_vertex*180/pi:.4f} deg")
print()
print("  This angle is related to the C_5 axis structure:")
print(f"    theta_vertex = arccos(1/sqrt(5)) = arccos(phi/(phi+2)^(1/2) * (phi+2)^(1/2)/phi^(1/2)^2)...")
print(f"    More directly: 1/sqrt(5) = 1/||(1,2)||  [winding vector norm]  [EXACT]")
print()
print("  Half-angle formula: cos(theta_vertex/2) = sqrt((1+cos(theta_vertex))/2)")
cos_half_vertex = math.sqrt((1 + cos_vertex)/2)
theta_half = math.acos(cos_half_vertex)
print(f"    cos(theta_vertex/2) = sqrt((1 + 1/sqrt(5))/2)")
print(f"                       = {cos_half_vertex:.8f}")
print(f"    theta_vertex/2     = {theta_half*180/pi:.4f} deg")
print(f"    vs theta_W measured  = {theta_W_measured*180/pi:.4f} deg")
print(f"    Bare half-angle gap: {(cos_half_vertex/cos_W_measured - 1)*100:+.4f}%  (3.5% off)")
print()

# ── Alpha correction with 5-fold vertex count ─────────────────────────────────
print(SEP)
print("ALPHA CORRECTION: (1 + 5*alpha)")
print(SEP2)
print()
print("  Each icosahedral vertex has 5 edges (5-fold local symmetry).")
print("  This is the SAME 5-fold structure that gives f1=PHI in the alpha derivation.")
print("  In alpha: the 5 edges provide the coupling factor, giving PHI at leading order.")
print("  Analogously: the Weinberg mixing gets a correction of 5*alpha (5 edges * coupling).")
print()
print(f"  5 * alpha = {5*alpha:.8f} = {5*alpha*100:.4f}%")
correction = 1 + 5*alpha
print(f"  (1 + 5*alpha) = {correction:.8f}")
print()

cos_W_pred = cos_half_vertex * correction
theta_W_pred = math.acos(cos_W_pred)
sin2_W_pred  = math.sin(theta_W_pred)**2

print(f"  CANDIDATE FORMULA:")
print(f"    cos(theta_W) = sqrt((1 + 1/sqrt(5))/2) * (1 + 5*alpha)")
print(f"                 = {cos_half_vertex:.8f} * {correction:.8f}")
print(f"                 = {cos_W_pred:.8f}")
print(f"    theta_W_pred = {theta_W_pred*180/pi:.4f} deg")
print(f"    sin^2(theta_W_pred) = {sin2_W_pred:.6f}")
print()
print(f"  vs measured:  cos(theta_W) = {cos_W_measured:.8f}")
gap_pct = (cos_W_pred/cos_W_measured - 1)*100
sigma_from_W = abs(cos_W_pred - cos_W_measured) * m_Z / unc_W
print(f"  Gap: {gap_pct:+.4f}%  = {sigma_from_W:.2f} sigma from PDG m_W uncertainty")
print()
print(f"  Improvement over bare half-angle: {3.5:.1f}% -> {abs(gap_pct):.3f}%")
print(f"  Factor improvement: {3.5/abs(gap_pct):.0f}x")
print()

# ── Physical interpretation ───────────────────────────────────────────────────
print(SEP)
print("PHYSICAL INTERPRETATION")
print(SEP2)
print()
print("  WHY the half-angle of the vertex opening angle?")
print("  The photon (A_g, pole vertex) and W_3 (T_1g, equatorial vertex) mix.")
print("  The OPENING ANGLE between them (theta_vertex = 63.43 deg) determines")
print("  the mixing. The HALF-ANGLE enters because Weinberg mixing is a")
print("  coherent superposition (like |theta_W> = cos|B> + sin|W_3>) where")
print("  the amplitude is related to sqrt((1+cos(full_angle))/2).")
print()
print("  WHY (1 + 5*alpha)?")
print("  Each vertex has 5 neighboring vertices (edges). When a particle wave")
print("  (like the electron) interacts with the vertex, all 5 edge contacts")
print("  contribute, each at strength alpha. This shifts the effective")
print("  mixing angle from the bare geometric value by 5*alpha.")
print("  Identical structure to the vertex stiffness correction in alpha:")
print("    alpha derivation: n_exact = n_topo + L3(PHI,log5)*delta_k")
print("                              = integer + (5-fold vertex correction)")
print("    Weinberg mixing:  cos(theta_W) = bare * (1 + 5*alpha)")
print()
print("  The 5-fold vertex count is the bridge between alpha and theta_W.")
print()

# ── Algebraic form ────────────────────────────────────────────────────────────
print(SEP)
print("ALGEBRAIC SIMPLIFICATION")
print(SEP2)
print()
print("  sqrt((1 + 1/sqrt(5))/2) in terms of phi:")
print("    1/sqrt(5) = sqrt(5)/5 = (phi^2-phi)/(phi^2-phi+... ) ...")
print("    Actually: 1/sqrt(5) = phi/(phi+2) * ??? -- let me compute directly")
val = (1 + 1/math.sqrt(5))/2
print(f"    (1+1/sqrt(5))/2 = {val:.10f}")
print(f"    = (sqrt(5)+1)/(2*sqrt(5))")
print(f"    = phi*2/(2*sqrt(5))  [since phi = (1+sqrt(5))/2, so 2*phi = 1+sqrt(5)]")
print(f"    = phi/sqrt(5)")
val2 = phi/math.sqrt(5)
print(f"    phi/sqrt(5) = {val2:.10f}  vs (1+1/sqrt(5))/2 = {val:.10f}")
print(f"    Equal? {abs(val - val2) < 1e-10}")
print()
print("  THEREFORE:")
print("    sqrt((1+1/sqrt(5))/2) = sqrt(phi/sqrt(5))")
print("                          = phi^(1/2) / 5^(1/4)")
print(f"    = phi^(1/2) / 5^(1/4) = {phi**0.5 / 5**0.25:.10f}")
print(f"    vs direct            = {cos_half_vertex:.10f}")
print(f"    Equal? {abs(phi**0.5 / 5**0.25 - cos_half_vertex) < 1e-10}")
print()
print("  FINAL FORM of candidate formula:")
print("    cos(theta_W) = (phi/sqrt(5))^(1/2) * (1 + 5*alpha)")
print("                 = phi^(1/2) / 5^(1/4) * (1 + 5*alpha)")
print()
print(f"  Numerical check: {phi**0.5 / 5**0.25 * correction:.8f}")
print(f"  vs measured:     {cos_W_measured:.8f}  (gap: {(phi**0.5/5**0.25*correction/cos_W_measured-1)*100:+.4f}%)")
print()

# ── Comparison with other candidates ──────────────────────────────────────────
print(SEP)
print("COMPARISON WITH OTHER CANDIDATES")
print(SEP2)
print()
print(f"  {'Formula':<45} {'cos(theta_W)':>12}  {'gap%':>8}")
print(SEP2)
candidates = [
    ("measured m_W/m_Z",                    cos_W_measured),
    ("sqrt((1+1/sqrt(5))/2) bare",          cos_half_vertex),
    ("sqrt((1+1/sqrt(5))/2)*(1+5*alpha)",   cos_half_vertex*(1+5*alpha)),
    ("phi/sqrt(phi+2)  [vertex polar]",     phi/math.sqrt(phi+2)),
    ("phi/sqrt(phi+2)*(1+5*alpha)",         phi/math.sqrt(phi+2)*(1+5*alpha)),
    ("arctan(1/phi) -> cos",                math.cos(math.atan(1/phi))),
    ("dim(A_g)/sqrt(dim(A_g)^2+dim(T1g)^2)",1/math.sqrt(1+9)),
    ("GUT: sin^2=3/8 -> cos",               math.sqrt(1-3/8)),
]
for name, val in candidates:
    if val == cos_W_measured:
        pct = 0.0
    else:
        pct = (val/cos_W_measured - 1)*100
    print(f"  {name:<45} {val:>12.6f}  {pct:>+8.4f}%")
print()

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
print()
print("  cos(theta_W) = phi^(1/2) / 5^(1/4) * (1 + 5*alpha)")
print(f"              = {phi**0.5/5**0.25:.8f} * {correction:.8f}")
print(f"              = {cos_W_pred:.8f}")
print(f"  vs measured   {cos_W_measured:.8f}")
print(f"  Gap: {gap_pct:+.4f}% = {sigma_from_W:.1f} sigma")
print()
print("  FORMULA STRUCTURE:")
print("  Leading term: phi^(1/2) / 5^(1/4) = icosahedral vertex half-angle cosine")
print("    - phi = golden ratio from (1,2) winding [from alpha derivation]")
print("    - 5   = 5-fold vertex count of icosahedron")
print("  Correction:   (1 + 5*alpha) = 5-edge vertex coupling")
print("    - 5   = edges per vertex")
print("    - alpha = EM coupling [from alpha derivation]")
print()
print("  STATUS: promising (0.03% off), same structural origin as alpha correction.")
print("  Needs: a DERIVED reason why the correction is exactly 5*alpha (not fitted).")
print("  If exact: predicts m_W = m_Z * phi^(1/2)/5^(1/4) * (1+5*alpha)")
m_W_pred = m_Z * cos_W_pred
print(f"  m_W_predicted = {m_W_pred:.3f} MeV  vs PDG {m_W} MeV")
print(f"  Gap: {m_W_pred - m_W:+.3f} MeV = {abs(m_W_pred-m_W)/unc_W:.1f} sigma")
print(SEP)
