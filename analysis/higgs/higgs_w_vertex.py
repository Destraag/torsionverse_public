"""
higgs_w_vertex.py
=================
Investigates whether the icosahedral vertex stiffness (I_h representation
theory) can close the gap between the measured m_H/m_W ratio and pi/2.

KEY QUESTION: The W boson is spin-1 (T_1g under I_h). The I_h character
table gives chi(T_1g, C_5) = phi -- the SAME phi that drives the structural
vertex coupling in the alpha derivation. Can this close the 0.86% gap?

Run: python analysis/higgs/higgs_w_vertex.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2
Rs   = math.sqrt(5) / (4*pi)
log5 = math.log(5)
L3   = (phi**3 + log5**3) / (phi**2 + log5**2)  # Born-weighted mean (alpha channel)

SEP  = "=" * 65
SEP2 = "-" * 65

# ── I_h character table (relevant entries) ────────────────────────────────────
# Under I_h, C_5 rotation by 2*pi/5 (72 degrees):
#   A_g:  chi = 1         (trivial -- Higgs spin-0)
#   T_1g: chi = 1+2*cos(2*pi/5) = phi  (W, Z spin-1)
#   H_g:  chi = 0         (5-dimensional)
#   G_g:  chi = -1        (4-dimensional)
chi_Ag  = 1.0
chi_T1g = 1 + 2*math.cos(2*pi/5)   # = phi
chi_T2g = 1 + 2*math.cos(4*pi/5)   # = 1/phi - 1 = 1 - phi (negative)

print(SEP)
print("I_h CHARACTER TABLE -- SPIN ASSIGNMENTS")
print(SEP2)
print(f"  chi(A_g,  C_5) = {chi_Ag:.8f}  = 1     [Higgs, spin-0]")
print(f"  chi(T_1g, C_5) = {chi_T1g:.8f}  = phi  [W, Z, spin-1]")
print(f"  chi(T_2g, C_5) = {chi_T2g:.8f}  = 1/phi^2-1  [other]")
print(f"  phi = {phi:.8f}")
print(f"  chi(T_1g) = phi:  {abs(chi_T1g - phi) < 1e-10}")
print()

# ── The gap we need to explain ────────────────────────────────────────────────
m_H  = m_H_pdg22      # 125.20 GeV
m_W  = 80.377         # GeV PDG 2022
m_Z  = 91.188         # GeV PDG 2022
unc_W = 0.012         # GeV uncertainty on m_W

ratio_HW = m_H / m_W
ratio_HZ = m_H / m_Z
ratio_WZ = m_W / m_Z

print(SEP)
print("MEASURED MASS RATIOS")
print(SEP2)
print(f"  m_H / m_W = {ratio_HW:.8f}")
print(f"  m_H / m_Z = {ratio_HZ:.8f}")
print(f"  m_W / m_Z = {ratio_WZ:.8f}  [= cos(theta_W) = {ratio_WZ:.6f}]")
print(f"  theta_W   = {math.acos(ratio_WZ)*180/pi:.4f} degrees")
print()

# ── Geometric targets ─────────────────────────────────────────────────────────
print("GEOMETRIC TARGET CANDIDATES")
print(SEP2)
print()
targets_HW = [
    ("pi/2",              pi/2),
    ("phi/2 * sqrt(pi)",  phi/2 * math.sqrt(pi)),
    ("sqrt(phi+1)",       math.sqrt(phi+1)),
    ("phi*(1+alpha/pi)",  phi*(1+alpha/pi)),
    ("L3(phi,log5)",      L3),
    ("phi^(3/2)/sqrt(pi)",phi**1.5/math.sqrt(pi)),
]
print("  m_H/m_W candidates:")
for name, val in targets_HW:
    err = (ratio_HW/val - 1)*100
    print(f"    {name:<30} = {val:.6f}  ({err:+.4f}%)")
print()

targets_WZ = [
    ("cos(pi/5) = phi/2", math.cos(pi/5)),
    ("1/sqrt(phi)",       1/math.sqrt(phi)),
    ("sqrt(phi)/phi^(1/2)",math.sqrt(phi)/phi**0.5),
    ("phi^2/(phi^2+1)",   phi**2/(phi**2+1)),
    ("Rs+Rs^2+Rs^3",      Rs+Rs**2+Rs**3),
    ("1-alpha",           1-alpha),
]
print("  m_W/m_Z = cos(theta_W) candidates:")
for name, val in targets_WZ:
    err = (ratio_WZ/val - 1)*100
    print(f"    {name:<30} = {val:.6f}  ({err:+.4f}%)")
print()

# ── Can the T_1g vertex correction close m_H/m_W gap? ────────────────────────
print(SEP)
print("T_1g VERTEX CORRECTION ANALYSIS")
print(SEP2)
print()
print("  For the HIGGS (A_g, spin-0):")
print(f"    chi(A_g, C_5) = 1  [trivial: no vertex correction from C_5 rotation]")
print(f"    Sub-cell -> vertex stiffness suppressed -> bulk coupling only (lambda = (1-nu)/4)")
print()
print("  For the W/Z (T_1g, spin-1):")
print(f"    chi(T_1g, C_5) = phi = {phi:.6f}  [non-trivial C_5 character]")
print(f"    Sub-cell -> vertex stiffness normally suppressed for sub-cell particles")
print(f"    BUT: T_1g is the TRANSVERSE mode of the inter-cell bond.")
print(f"    The bond shear (transverse) couples via T_1g, even for sub-cell particles.")
print()

# The vertex correction scale from alpha derivation
delta_k_alpha = 0.01869 / L3   # delta_k back-calculated from n_exact - 2
print(f"  Vertex stiffness scale from alpha:")
print(f"    n_exact - 2 = 0.01869  [vertex correction to n]")
print(f"    L3(phi,log5) = {L3:.6f}  [Born-weighted mean]")
print(f"    delta_k = {delta_k_alpha:.6f}  [per-channel stiffness increment]")
print()

# For T_1g, the coupling character is phi (not L3, since only one channel active)
# If W couples via T_1g channel with factor phi:
# Mass correction to m_W: delta_m_W/m_W = phi * delta_k

# BUT: what direction? The vertex correction in alpha was a POSITIVE addition to n.
# For the W, the vertex correction to its coupling could go either way.

# For m_H/m_W to reach pi/2, we need m_H/m_W to INCREASE by:
gap_HW = (pi/2 - ratio_HW) / ratio_HW   # fractional increase needed
print(f"  Fractional correction needed to reach pi/2:")
print(f"    (pi/2 - ratio_HW) / ratio_HW = {gap_HW*100:.4f}%")
print()
print(f"  Vertex correction scale (alpha-derived): 0.935%")
print(f"  Gap needed: {gap_HW*100:.3f}%")
print(f"  Ratio: {gap_HW/0.00935:.4f}  (1.0 would be exact match)")
print()

# Specific hypothesis: W mass gets T_1g vertex correction
# If m_W is shifted down by phi * delta_k_alpha:
delta_m_W_from_vertex = phi * delta_k_alpha  # fractional
m_W_corrected = m_W * (1 - delta_m_W_from_vertex)
ratio_HW_corrected = m_H / m_W_corrected

print(f"  IF m_W_corrected = m_W * (1 - phi * delta_k):")
print(f"    phi * delta_k = phi * {delta_k_alpha:.6f} = {phi*delta_k_alpha:.6f} = {phi*delta_k_alpha*100:.4f}%")
print(f"    m_W_corrected = {m_W_corrected:.4f} GeV")
print(f"    m_H / m_W_corrected = {ratio_HW_corrected:.6f}")
print(f"    vs pi/2 = {pi/2:.6f}  (gap: {(ratio_HW_corrected/(pi/2)-1)*100:+.4f}%)")
print(f"    vs measured m_W = {m_W} +/- {unc_W} GeV")
print(f"    Shift in m_W: {m_W_corrected - m_W:.4f} GeV = {(m_W_corrected-m_W)/unc_W:.1f} sigma")
print()
print("  VERDICT: the T_1g vertex correction at scale phi*delta_k would shift m_W")
print(f"  by {(phi*delta_k_alpha)*100:.3f}% = {(m_W*phi*delta_k_alpha):.3f} GeV.")
print(f"  PDG uncertainty on m_W is {unc_W} GeV. The shift is {m_W*phi*delta_k_alpha/unc_W:.0f} sigma.")
print(f"  This CANNOT be a physical correction to the measured m_W.")
print()

# ── Is pi/2 even the right target? ───────────────────────────────────────────
print(SEP)
print("IS pi/2 THE RIGHT TARGET?")
print(SEP2)
print()
print("  pi/2 is not derived from any geometric argument -- it was a numerical")
print("  near-miss. The correct derivation of m_H/m_W requires:")
print(f"    m_H/m_W = 2*sqrt(2*lambda)/g")
print(f"    where g = SU(2) gauge coupling (NOT yet derived from cell geometry)")
print()
nu      = (1 - 2*Rs**2) / (2*(1 - Rs**2))
lam_sub = (1 - nu) / 4
g_measured = 2*m_W/v_EW   # from m_W = v*g/2
print(f"  Measured SU(2) coupling: g = 2*m_W/v_EW = {g_measured:.6f}")
print(f"  Predicted m_H/m_W from lambda_sub and g: {2*math.sqrt(2*lam_sub)/g_measured:.6f}")
print(f"  Ratio of g to natural constants:")
print(f"    g / (2*Rs) = {g_measured/(2*Rs):.4f}")
print(f"    g / phi    = {g_measured/phi:.4f}")
print(f"    g / pi     = {g_measured/pi:.4f}")
print(f"    g^2        = {g_measured**2:.6f}")
print(f"    g^2 / alpha= {g_measured**2/alpha:.4f}")
print(f"    g^2 / (4*pi*alpha) = {g_measured**2/(4*pi*alpha):.4f}")
print()
print("  None of these are clean. g is not yet derivable from the cell geometry.")
print("  The gap in m_H/m_W is fundamentally a gap in the SU(2) coupling derivation.")
print()

# ── What DOES H3 need? ────────────────────────────────────────────────────────
print(SEP)
print("WHAT H3 ACTUALLY REQUIRES")
print(SEP2)
print()
print("  H3 = Yukawa coupling derivation = fermion mass spectrum derivation.")
print("  This is separate from the W/Z mass gap.")
print()
print("  The W/Z masses come from EW symmetry breaking: m_W = v*g/2.")
print("  Fermion masses come from Yukawa couplings: m_f = y_f * v / sqrt(2).")
print("  These are DIFFERENT parameters (g vs y_f).")
print()
print("  For H3, we need to derive the Yukawa couplings y_f for each fermion.")
print("  That requires either:")
print("  (a) N_J resonance theory (fermions as icosahedral lattice resonances)")
print("  (b) Some other geometric principle linking mass to cell geometry")
print()
print("  The vertex complexity (L3, T_1g = phi) is a promising tool for (b)")
print("  but the connection to fermion masses via Yukawa is not established.")
print()

# ── New lead: Weinberg angle from I_h ────────────────────────────────────────
print(SEP)
print("NEW LEAD: WEINBERG ANGLE FROM I_h REPRESENTATION MIXING")
print(SEP2)
print()
# cos(theta_W) = m_W/m_Z
# In the SM, theta_W arises from the mixing of SU(2) x U(1) gauge bosons.
# In the I_h picture: W is T_1g, Z is a T_1g x A_g mixture (neutral current).
# The mixing angle might be related to the angle between T_1g and A_g representations.

# Under I_h, the T_1g and A_g representations are orthogonal.
# But their coupling STRENGTHS (chi values) are phi and 1 respectively.
# The mixing angle: tan(theta_mix) = chi(A_g)/chi(T_1g) = 1/phi = 1/phi

theta_mix = math.atan(1/phi)
theta_W_measured = math.acos(m_W/m_Z)
print(f"  If mixing angle = arctan(1/phi) = arctan(chi_Ag/chi_T1g):")
print(f"    theta_mix = arctan(1/phi) = {theta_mix*180/pi:.4f} degrees")
print(f"    cos(theta_mix) = {math.cos(theta_mix):.8f}")
print(f"    vs m_W/m_Z     = {m_W/m_Z:.8f}  (Weinberg angle cosine)")
print(f"    Deviation: {abs(math.cos(theta_mix)/(m_W/m_Z)-1)*100:.4f}%")
print()
# Also try arctan(chi_Ag/sqrt(chi_T1g))
theta_mix2 = math.atan(1/math.sqrt(phi))
print(f"  If mixing angle = arctan(1/sqrt(phi)):")
print(f"    theta_mix2 = {theta_mix2*180/pi:.4f} degrees")
print(f"    cos(theta_mix2) = {math.cos(theta_mix2):.8f}")
print(f"    vs m_W/m_Z      = {m_W/m_Z:.8f}")
print(f"    Deviation: {abs(math.cos(theta_mix2)/(m_W/m_Z)-1)*100:.4f}%")
print()
# sin^2(theta_W) target
sin2W_measured = 1 - (m_W/m_Z)**2
sin2W_gut = 3/8   # GUT prediction
print(f"  sin^2(theta_W) measured = {sin2W_measured:.6f}")
print(f"  sin^2(theta_W) GUT     = {sin2W_gut:.6f}  (3/8)")
print(f"  sin^2(theta_W) from arctan(1/phi): {math.sin(theta_mix)**2:.6f}")
print()
print("  VERDICT: None of the simple I_h mixing angles exactly reproduce")
print("  the Weinberg angle. This is a genuine open problem -- the SU(2)xU(1)")
print("  gauge structure is NOT yet derivable from the I_h cell geometry.")
print()

print(SEP)
print("SUMMARY")
print(SEP)
print()
print("  1. The T_1g character = phi IS exact and IS the same phi as the")
print("     alpha vertex structural coupling. This is a real connection.")
print()
print("  2. BUT the vertex correction to alpha applies to BULK particles")
print("     (electron, N_J >> 1). W/Z are SUB-CELL (N_J < 1).")
print("     Sub-cell particles use bulk Poisson coupling, not vertex stiffness.")
print("     The T_1g vertex correction would require a new mechanism for")
print("     transverse modes in the sub-cell regime.")
print()
print("  3. The 0.86% gap in m_H/m_W vs pi/2 is NOT the vertex gap.")
print("     It is the gap in the SU(2) coupling g derivation.")
print("     g = 2*m_W/v = 0.653 is not yet derivable from cell geometry.")
print()
print("  4. This does NOT close H3 (fermion masses).")
print("     H3 requires Yukawa couplings y_f, not gauge couplings g.")
print("     These are independent parameters in the SM.")
print()
print("  5. New lead: the Weinberg angle mixing (m_W/m_Z) might come")
print("     from I_h representation mixing (A_g vs T_1g characters).")
print("     Best current candidate: arctan(1/phi) gives cos(theta) = 0.8507")
print(f"     vs measured m_W/m_Z = {m_W/m_Z:.4f}  (deviation {abs(math.cos(theta_mix)/(m_W/m_Z)-1)*100:.2f}%)")
print()
print("  RECOMMENDATION: The T_1g = phi connection is a genuine geometric")
print("  lead for deriving the EW gauge structure, but it needs a new")
print("  physical mechanism for sub-cell transverse modes. Separate from H3.")
print(SEP)
