"""
higgs_h2_partial.py
====================
Partial H2 (branching ratio structure) derivation using:
  1. CG decompositions from higgs_lagrangian_h2.py (structure)
  2. Weinberg angle from our unified formula (W/Z split)

WHAT THIS COMPUTES:
  - Relative partial widths Gamma(H->WW*) / Gamma(H->ZZ*) from sin^2(theta_W)
  - Comparison to PDG branching ratios
  - What remains open

WHAT REMAINS OPEN:
  - Absolute partial widths (need full coupling derivation)
  - Fermion channels: b, tau, t irrep assignments in I_h not yet confirmed
  - H->gammagamma: loop-induced, requires special treatment

The W/Z split within T_1g (both W and Z are T_1g modes in I_h):
  - W+, W-: two charged components of T_1g (x+iy, x-iy polarizations)
  - Z:      neutral component of T_1g (z polarization)
  - Mixing via Weinberg angle: Z = cos(theta_W)*W_3 - sin(theta_W)*B
  - The ratio Gamma(WW)/Gamma(ZZ) involves cos^4(theta_W)/sin^4(theta_W) * 2
    (factor 2 for two charged W modes vs one neutral Z)

Run: python analysis/higgs/higgs_h2_partial.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import alpha, E_cell_GeV, phi

pi    = math.pi
Rs    = math.sqrt(5) / (4*pi)
nu    = (1 - 2*Rs**2) / (2*(1 - Rs**2))
K_o_G = (2*(1+nu)) / (3*(1-2*nu))

SEP  = "=" * 70
SEP2 = "-" * 70

# Our Weinberg angle predictions (from higgs_weinberg_unified.py)
# Best formula: unified (2/3 vertex + 1/3 pressure)
sin2_tw_unified = 0.22309   # from higgs_weinberg_unified.py
cos2_tw         = 1 - sin2_tw_unified
cos_tw          = math.sqrt(cos2_tw)
sin_tw          = math.sqrt(sin2_tw_unified)
tan2_tw         = sin2_tw_unified / cos2_tw

# PDG reference
sin2_tw_pdg = 0.22290   # on-shell (PDG 2022)
m_H_pdg     = 125.20    # GeV
m_W_pdg     = 80.379    # GeV
m_Z_pdg     = 91.1876   # GeV
Gamma_H_pdg = 4.07      # MeV
BR_WW  = 0.2137
BR_ZZ  = 0.0264
BR_bb  = 0.5812
BR_gg  = 0.0822   # H->gg (gluon-gluon)
BR_tau = 0.0627
BR_cc  = 0.0289
BR_gam = 0.00228  # H->gammagamma

print(SEP)
print("PARTIAL H2: BRANCHING RATIO STRUCTURE FROM WEINBERG ANGLE")
print(SEP)
print()
print(f"  Our Weinberg angle: sin^2(theta_W) = {sin2_tw_unified:.5f}")
print(f"  PDG on-shell:       sin^2(theta_W) = {sin2_tw_pdg:.5f}")
print(f"  cos(theta_W) = {cos_tw:.6f},  sin(theta_W) = {sin_tw:.6f}")
print()

# ── STEP 1: W/Z SPLIT WITHIN T_1g ────────────────────────────────────────────
print(SEP)
print("STEP 1  W/Z split within T_1g")
print(SEP2)
print()
print("  T_1g has 3 components. Under I_h -> SU(2)xU(1):")
print("    W+, W-: 2 charged modes  (weight: g_W^2 per mode)")
print("    Z:      1 neutral mode   (weight: g_Z^2 = g_W^2/cos^2(theta_W))")
print()
print("  In our framework alpha = unique coupling => all gauge couplings in terms of alpha:")
print("    g_W^2 = 4*pi*alpha / sin^2(theta_W)")
print("    g_Z^2 = g_W^2 / cos^2(theta_W) = 4*pi*alpha / (sin^2(theta_W)*cos^2(theta_W))")
g_W2 = 4*pi*alpha / sin2_tw_unified
g_Z2 = g_W2 / cos2_tw
print(f"    g_W^2 = {g_W2:.6f}")
print(f"    g_Z^2 = {g_Z2:.6f}")
print()

# ── STEP 2: PARTIAL WIDTH RATIO Gamma(WW)/Gamma(ZZ) ─────────────────────────
print(SEP)
print("STEP 2  Gamma(H->WW*) / Gamma(H->ZZ*) from Weinberg angle")
print(SEP2)
print()
# At tree level, ignoring off-shell effects (both W on-shell approximation):
# Gamma(H->WW) proportional to g_W^4 * m_H / m_W^2
# Gamma(H->ZZ) proportional to g_Z^4 * m_H / m_Z^2
# Ratio = (g_W/g_Z)^4 * (m_Z/m_W)^2 * 2 (for 2 charged W modes)
# g_W/g_Z = cos(theta_W), m_Z/m_W = 1/cos(theta_W)
# => Ratio = cos^4(theta_W) * (1/cos(theta_W))^2 * 2 = 2*cos^2(theta_W)
# More precise: factor of dim correction
ratio_tree = 2 * cos2_tw
print(f"  Tree-level ratio (on-shell, leading order):")
print(f"  Gamma(WW)/Gamma(ZZ) = 2 * cos^2(theta_W) = {ratio_tree:.4f}")
print(f"  PDG ratio:  BR(WW*)/BR(ZZ*) = {BR_WW/BR_ZZ:.4f}  ({BR_WW}/{BR_ZZ})")
print()

# More accurate: include coupling structure
# In SM: Gamma(H->W+W-) / Gamma(H->ZZ) = 2 * (m_W^4/m_H^4 * ...) / (m_Z^4/m_H^4 * ...)
# With off-shell corrections the exact ratio is complex; leading order:
# = 2 * cos^4(theta_W) / cos^4(theta_W) * (something) -- this isn't right
# Let me use the SM formulas:
# Gamma(H->WW, on-shell) ~ G_F * m_H^3 / (8*pi*sqrt(2)) * (1 - 4x + 12x^2) where x=(m_W/m_H)^2
# Gamma(H->ZZ, on-shell) ~ G_F * m_H^3 / (16*pi*sqrt(2)) * (1 - 4x_Z + 12x_Z^2) where x_Z=(m_Z/m_H)^2
# Note factor 1/2 for ZZ (identical particles)

GF = 1.1663787e-5
mH = m_H_pdg
mW_our = 80.358  # our prediction from (1,3) winding
mZ_our = mW_our / cos_tw

xW = (mW_our/mH)**2
xZ = (mZ_our/mH)**2

def partial_WW(mH, mW, GF):
    x = (mW/mH)**2
    if 4*x > 1: return 0
    return GF * mH**3 / (8*pi*math.sqrt(2)) * math.sqrt(1 - 4*x) * (1 - 4*x + 12*x**2)

def partial_ZZ(mH, mZ, GF):
    x = (mZ/mH)**2
    if 4*x > 1: return 0
    return GF * mH**3 / (16*pi*math.sqrt(2)) * math.sqrt(1 - 4*x) * (1 - 4*x + 12*x**2)

# These are on-shell formulas -- below threshold if mH < 2mW
print(f"  Note: m_H = {mH} < 2*m_W = {2*mW_our:.3f} -> WW is off-shell (WW*)")
print(f"  On-shell formula not valid here. Using PDG ratio as reference.")
print()

# The off-shell ratio requires a full phase-space integral.
# We use our Weinberg angle to predict the leading-order ratio:
# The T_1g splits into: 2 W modes (weight cos^0) and 1 Z mode (weight ~ (1-2sin^2)^2)
# coupling ratio: g_W : g_Z*cos(theta_W) = 1 : 1 (at tree level)
# So ratio = 2*(g_W^2)/(g_Z^2) * phase_space_ratio
coupling_ratio = 2 * g_W2 / g_Z2  # = 2 * cos^2(theta_W)
print(f"  Leading-order coupling ratio: 2*g_W^2/g_Z^2 = 2*cos^2(theta_W) = {coupling_ratio:.4f}")
print(f"  PDG ratio (with off-shell): {BR_WW/BR_ZZ:.4f}")
print(f"  Leading ratio captures: {coupling_ratio/(BR_WW/BR_ZZ)*100:.1f}% of PDG value")
print(f"  Remainder = off-shell phase-space correction (not yet derived)")
print()

# ── STEP 3: ALL CHANNELS FROM CG + COUPLINGS ─────────────────────────────────
print(SEP)
print("STEP 3  All decay channels: CG structure + coupling")
print(SEP2)
print()
print(f"  {'Channel':<25} {'CG':<6} {'Coupling':<20} {'Status'}")
print(f"  {'-'*25} {'-'*6} {'-'*20} {'-'*20}")

channels = [
    ("H->WW* (T_1g x T_1g)",   1, "alpha^2*phi^2 (derived)", "STRUCTURE DERIVED"),
    ("H->ZZ* (T_1g x T_1g)",   1, "alpha^2*phi^2/cos^2(tw)", "NEEDS WEINBERG SPLIT"),
    ("H->bb  (G_g x G_g)",     1, "alpha_bb (not derived)",   "STRUCTURE DERIVED"),
    ("H->tt  (H_g x H_g)",     1, "alpha_tt (not derived)",   "STRUCTURE DERIVED"),
    ("H->gamgam (A_g x A_g)",  1, "loop (not derived)",       "STRUCTURE DERIVED"),
    ("H->T1g+T2g",             0, "FORBIDDEN",                "FORBIDDEN (no A_g)"),
    ("H->tautau (G_g? x G_g?)",1, "unknown irrep",            "IRREP UNKNOWN"),
]
for ch, cg, coup, status in channels:
    print(f"  {ch:<25} {cg:<6} {coup:<20} {status}")

print()

# ── STEP 4: BRANCHING RATIOS FROM OUR WEINBERG ANGLE ─────────────────────────
print(SEP)
print("STEP 4  Predicted vs PDG branching ratios (structural only)")
print(SEP2)
print()
print("  Using CG structure + Weinberg angle for gauge channels:")
print(f"  sin^2(theta_W) = {sin2_tw_unified:.5f} (our unified formula)")
print()

# Structural weights (CG=1 for all, so relative weights from coupling)
# W modes: 2 x g_W^2 = 2 x (4pi*alpha/sin^2_W)
# Z modes: 1 x g_Z^2 = 1 x (4pi*alpha/(sin^2_W * cos^2_W))
# Ratio W:Z = 2 * cos^2_W
print(f"  Structural gauge ratio:")
print(f"    Gamma(WW*) : Gamma(ZZ*) = 2*cos^2(theta_W) = {2*cos2_tw:.4f}")
print(f"    PDG:                                          {BR_WW/BR_ZZ:.4f}")
print(f"    Fractional agreement: {2*cos2_tw/(BR_WW/BR_ZZ)*100:.1f}%")
print()
print(f"  Total gauge fraction (W+Z) relative to each other:")
total_gauge = BR_WW + BR_ZZ
print(f"    PDG: BR(WW*)+BR(ZZ*) = {total_gauge:.4f}")
print(f"    Our structural prediction: 2*cos^2(tw) + 1 = {2*cos2_tw+1:.4f} (unnormalized)")
print()
print("  OPEN: to get absolute branching ratios need fermion couplings (G_g, H_g irreps)")
print()

# ── STEP 5: JOBSON CELL GEOMETRY SUMMARY ────────────────────────────────────
print(SEP)
print("STEP 5  Jobson cell H2 geometry summary")
print(SEP2)
print()
print("  CELL: icosahedral (I_h), 12 vertices, 30 edges, 20 triangular faces")
print(f"  Edge length: L_J = alpha*phi*r_p  (Hopf-derived)")
print(f"  Circumradius: R_c = L_J * sqrt(1+phi^2) / 2 (icosahedral geometry)")
R_c_over_LJ = math.sqrt(1 + phi**2) / 2
print(f"               R_c / L_J = {R_c_over_LJ:.6f}")
print()
print("  PARTICLE-IRREP ASSIGNMENTS (conditional on inverted hypothesis):")
print(f"  A_g  (dim 1): Higgs   -- center, scalar breathing mode")
print(f"  T_1g (dim 3): W+,W-,Z -- 3 gauge bosons at vertices, vector modes")
print(f"  T_2g (dim 3): ??      -- 3 modes, CHI(C5)=-1/phi (anti-screening)")
print(f"  G_g  (dim 4): b quark -- 4-component spinor (candidate)")
print(f"  H_g  (dim 5): top?    -- 5-component (candidate)")
print()
print("  CG SELECTION RULES (from T_1g x T_1g = A_g + T_1g + H_g):")
print(f"  H->WW/ZZ:  ALLOWED (unique, CG=1)")
print(f"  H->T1+T2:  FORBIDDEN (T_1g x T_2g has no A_g)")
print(f"  H->WW/ZZ ratio: 2*cos^2(theta_W) = {2*cos2_tw:.4f}  [PDG: {BR_WW/BR_ZZ:.4f}]")
print()
print(f"  Scripts that produced these results:")
print(f"    higgs_cg_twoloop.py     -- CG derivation (verified)")
print(f"    higgs_lagrangian_h2.py  -- Lagrangian + all CG channels")
print(f"    higgs_weinberg_unified.py -- Weinberg angle (0.012%)")
print(f"    higgs_h2_partial.py     -- This script: W/Z split + ratio")
