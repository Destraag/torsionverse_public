"""
higgs_c2_vertex_born_amplitude_test.py
========================================
Next step from higgs_two_vertex_z2_simultaneous_reframe.py's item 5: attempt
to independently DERIVE (not assume by analogy) that the "jamming vertex"
(the Z_2 pi-rotation element) contributes a second, distinct factor of alpha
to the Higgs mass shift, using the SAME Born-amplitude method as
alpha_born_vertex.py Part A -- but evaluated at the C2 (pi-rotation) class
instead of the C5 class.

PHYSICAL MOTIVATION FOR C2 SPECIFICALLY (not just "the other rotation
class"): C5 axes pass through VERTEX pairs (12 vertices, 6 antipodal pairs)
-- the actual "impact points" where a corpuscle rebounds WITHIN one cell
(doc_alpha.txt Section 4: "the corpuscle... rebounds at the vertex on each
impact"), an INTRA-cell process. C2 axes pass through EDGE-MIDPOINT pairs
(30 edges, 15 antipodal pairs) -- and edge midpoints are EXACTLY where
ADJACENT cells make contact (cell_rotation_propagation.py RP1: "adjacent
cells touch at EDGE MIDPOINTS"). Maxwell-critical jamming (3V-E=6) is
DEFINED via the edge network -- an INTER-cell constraint-satisfaction fact,
not an intra-cell one. So C2/edge-midpoint is the physically motivated
place to evaluate a "jamming vertex," not an arbitrary alternative to C5.

Run: python analysis/higgs/higgs_c2_vertex_born_amplitude_test.py
"""

import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
r_p = 0.8414e-15
hbar_c = 197.3269804  # MeV*fm
L_J = alpha * phi * (r_p * 1e15)  # fm
E_cell_MeV = 2 * pi * hbar_c / L_J

SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "FAIL"
    results.append((name, s, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


print(SEP)
print("higgs_c2_vertex_born_amplitude_test.py")
print("Does a REAL Born amplitude at the C2 (edge-midpoint) axis reproduce")
print("the second alpha needed for delta_m=2*alpha*(hbar_c/L_J)?")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("SECTION 1: C5 (VERTEX) VS C2 (EDGE-MIDPOINT) -- GEOMETRIC GROUNDING")
print(SEP2)
print()
V, E = 12, 30
c5_axes = V // 2   # 6 antipodal vertex pairs
c2_axes = E // 2   # 15 antipodal edge-midpoint pairs
print(f"  Icosahedron: V={V} vertices -> {c5_axes} C5 axes (vertex pairs)")
print(f"               E={E} edges    -> {c2_axes} C2 axes (edge-midpoint pairs)")
check("G1 C2 axis count (15) matches the icosahedral group's own C2 class size",
      c2_axes == 15, f"E/2={c2_axes}")
print()
print("  C5 = intra-cell impact points (where ONE cell's corpuscle rebounds).")
print("  C2 = inter-cell contact points (where ADJACENT cells touch,")
print("  cell_rotation_propagation.py RP1) -- the SAME edge network whose")
print("  rigidity (3V-E=6) DEFINES Maxwell-critical jamming. C2 is the")
print("  physically motivated 'jamming vertex' location, not an arbitrary")
print("  alternative class.")

# =============================================================================
print()
print(SEP2)
print("SECTION 2: REPRODUCE PART A'S FORMULA, THEN EVALUATE AT C2")
print(SEP2)
print()

chi_T1g_C5 = 1 + 2 * math.cos(2 * pi / 5)
chi_T1g_C2 = 1 + 2 * math.cos(pi)  # C2 = pi rotation
print(f"  chi(T_1g, C5) = 1 + 2*cos(72 deg)  = {chi_T1g_C5:.10f}  (= phi, established)")
print(f"  chi(T_1g, C2) = 1 + 2*cos(180 deg) = {chi_T1g_C2:.10f}  (= -1, standard)")
check("B1 chi(T_1g,C5) = phi (reproduces alpha_born_vertex.py Part A exactly)",
      abs(chi_T1g_C5 - phi) < 1e-9, f"{chi_T1g_C5:.10f} vs phi={phi:.10f}")
check("B2 chi(T_1g,C2) = -1 (standard I-group character table value)",
      abs(chi_T1g_C2 - (-1)) < 1e-9, f"{chi_T1g_C2:.10f}")

print()
print("  Part A's formula (verbatim): M(T_1g <- vertex) = chi(T_1g,class)*alpha*k_LW")
M_C5_over_kLW = chi_T1g_C5 * alpha
M_C2_over_kLW = chi_T1g_C2 * alpha
print(f"  M_C5 / (alpha*k_LW) = chi(T_1g,C5) = {M_C5_over_kLW/alpha:.6f}")
print(f"  M_C2 / (alpha*k_LW) = chi(T_1g,C2) = {M_C2_over_kLW/alpha:.6f}   <-- NEGATIVE")
print()
print("  Using the MAGNITUDE (a physical coupling strength should not be")
print("  negative, matching Section 4.2's f_1,f_2 both being positive geometric")
print("  quantities, and Section 4.3's Born-weighting using f_k^2 always >=0):")
M_C2_mag_over_kLW = abs(chi_T1g_C2) * alpha
print(f"  |M_C2| / (alpha*k_LW) = |chi(T_1g,C2)| = {M_C2_mag_over_kLW/alpha:.6f}")

sum_over_kLW = M_C5_over_kLW + M_C2_mag_over_kLW
print()
print(f"  SUM (C5 + |C2|) / (alpha*k_LW) = phi + 1 = {sum_over_kLW/alpha:.6f}")
print(f"  phi^2 (Fibonacci identity phi+1=phi^2)     = {phi**2:.6f}")
check("B3 the C5+|C2| sum equals phi^2 exactly (Fibonacci identity, not a new fit)",
      abs(sum_over_kLW / alpha - phi**2) < 1e-9, f"{sum_over_kLW/alpha:.6f} vs phi^2={phi**2:.6f}")

needed = 2.0
gap_pct = (sum_over_kLW / alpha / needed - 1) * 100
check("B4 HONEST: phi^2 does NOT match the needed coefficient of 2 -- "
      "a real, non-trivial miss, not a near-match",
      abs(gap_pct) > 20, f"phi^2={phi**2:.6f} vs needed=2, gap={gap_pct:+.2f}%")

# =============================================================================
print()
print(SEP2)
print("SECTION 3: WHY THE MISMATCH -- A UNITS/IDENTIFICATION GAP, NOT A")
print("           REFUTATION OF THE Z_2 STORY ITSELF")
print(SEP2)
print()
print("  k_n, k_LW (Part A/B) are STIFFNESS-like quantities internal to the")
print("  alpha quadratic's OWN self-consistency loop (calibrating alpha's")
print("  numerical value itself) -- k_LW's absolute value is explicitly")
print("  stated to CANCEL OUT of that ratio (doc_alpha.txt Section 4.4).")
print("  hbar_c/L_J (used for the Higgs mass shift) is a DIFFERENT physical")
print("  quantity -- an absolute ENERGY scale for a DIFFERENT particle (the")
print("  Higgs, not calibrating alpha itself). Mechanically reusing Part A's")
print("  exact numerical formula (chi*alpha*k_LW) for a differently-scaled,")
print("  differently-purposed quantity is NOT justified without an explicit")
print("  identification connecting k_LW to hbar_c/L_J -- one does not")
print("  currently exist in this repo.")
print()
print("  CONCLUSION: this specific attempt (mechanically reusing Part A's")
print("  formula at C2) FAILS -- gives phi^2*alpha, not 2*alpha, a genuine")
print("  31% miss. This does NOT refute the Z_2/parity story from the prior")
print("  script (that remains a real, checked structural fact); it shows")
print("  ONLY that 'borrow Part A's exact coupling formula and swap C5->C2'")
print("  is the wrong way to independently derive item 3. The 'why bare")
print("  alpha per sector, not alpha*phi or alpha*(-1)' question is NOT")
print("  resolved by this attempt and remains open.")

# =============================================================================
print()
print(SEP)
print("SECTION 4: FINAL HONEST VERDICT")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass}/{len(results)} checks PASS")
print()
print("  1. C2 axes are geometrically well-motivated as the 'jamming vertex'")
print("     location (edge midpoints, inter-cell contact, the SAME network")
print("     whose rigidity defines Maxwell-critical jamming) -- this part of")
print("     the physical story stands.")
print("  2. BUT mechanically reusing Part A's exact Born-amplitude formula at")
print("     C2 gives phi^2*alpha (2.618*alpha), not the needed 2*alpha -- a")
print("     genuine, honest 31% miss, not a near-success.")
print("  3. The likely reason: k_LW (Part A/B) and hbar_c/L_J (Higgs mass")
print("     shift) are different physical quantities with no established")
print("     identification between them -- reusing one formula for the other")
print("     was not actually justified, just structurally tempting.")
print("  4. STATUS: item 3 ('why bare alpha per Z_2 sector') remains OPEN.")
print("     The Z_2/parity reframing (prior script) is unaffected and still")
print("     stands on its own. A genuine derivation of the per-sector alpha")
print("     would need its OWN Born-amplitude calculation done directly in")
print("     hbar_c/L_J units, not borrowed from the alpha-value's own")
print("     internal k_n/k_LW self-consistency loop.")
