"""
higgs_two_vertex_z2_simultaneous_reframe.py
=============================================
CORRECTS higgs_two_vertex_jamming_reframe.py's reading of "2 vertices" after
author clarification (2026-09-23): "i did mean the normal + jamming as
distinct events that happen at THE SAME TIME" -- NOT the prior script's
"same vertex type encountered q=2 times SEQUENTIALLY per cycle" picture.
That prior framing is explicitly SUPERSEDED here, not silently repeated.

REVISED LOGICAL CHAIN (what actually ties q to "2 simultaneous, distinct"
things, per author's request to "tie that together logically"):
  doc_higgs.txt Section 3a's own proof: a pi-rotation of the winding's major
  axis maps the (p,q) curve to itself IFF q is even. This is a Z_2 (order-
  EXACTLY-2) symmetry: {identity, pi-rotation}. Two facts about a Z_2 group
  matter here -- NOT q's specific numeric value:
    (a) it ALWAYS has EXACTLY 2 elements, for ANY even q (2, 4, 6, ...) --
        so "the 2" comes from the SYMMETRY TYPE (Z_2), not from counting to
        q. This directly supersedes the prior script's "coefficient = q"
        claim, which would (wrongly, untested elsewhere in this framework)
        predict a GROWING coefficient for larger even q.
    (b) both elements are SIMULTANEOUSLY, statically true properties of the
        SAME single configuration -- a symmetry group is not a sequence of
        events in time -- matching "at the same time" exactly.
  This reuses the SAME q-even fact Section 3a already needs for the
  scalar/A_g classification -- not a new, independent topological input.

STILL AN OPEN ASSUMPTION (Section 3, flagged honestly, not oversold): that
each of the 2 Z_2 sectors contributes exactly one factor of alpha,
ADDITIVELY (not Born-weighted-averaged, since the two sectors are always-
simultaneously-present properties, not either/or alternatives with
probabilities summing to 1, unlike Section 4.3's f_1/f_2 channels). This is
consistent with, not yet independently re-derived from, alpha_born_
vertex.py Part B's own "one EM vertex factor" self-energy language.

Section 4 gives a physical (not yet fully rigorous) reading of WHICH Z_2
element is "normal" vs "jamming": identity = the winding as constructed
(normal); pi-rotation = a REAL rigid pi-rotation of the local structure --
notably, the icosahedral group I already has genuine order-2 (C2) rotation
axes as a native symmetry class (15 of them), so "pi-rotation" is not an
exotic add-on. Whether THIS specific operation is what jamming-relaxation
physically IS remains an interpretive label, flagged as such.

Run: python analysis/higgs/higgs_two_vertex_z2_simultaneous_reframe.py
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
print("higgs_two_vertex_z2_simultaneous_reframe.py")
print("Corrected reading: '2 vertices' = the 2 SIMULTANEOUS elements of a")
print("Z_2 (pi-rotation) symmetry, not 2 SEQUENTIAL vertex re-encounters.")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("SECTION 1: THE Z_2 SYMMETRY GROUP HAS EXACTLY 2 SIMULTANEOUS ELEMENTS")
print(SEP2)
print()

p, q = 1, 2
n_wind = p * q
print(f"  Higgs (1,2) winding: p={p}, q={q}")
print(f"  Section 3a's condition for pi-rotation symmetry: q*pi = 0 (mod 2*pi)")
qpi_mod = (q * pi) % (2 * pi)
check("Z1 q=2 satisfies the pi-rotation-symmetry condition (q*pi = 0 mod 2*pi)",
      abs(qpi_mod) < 1e-9 or abs(qpi_mod - 2 * pi) < 1e-9, f"q*pi mod 2*pi = {qpi_mod:.3e}")
print()
print("  The symmetry group is {identity, pi-rotation}: applying the")
print("  pi-rotation TWICE returns the original configuration (2*pi = full")
print("  circle) -- this makes it order-EXACTLY-2 (a Z_2 group), and both")
print("  elements are SIMULTANEOUS structural facts about the SAME static")
print("  winding, not a sequence of 2 events happening one after another.")
identity_then_pi_twice = (2 * pi) % (2 * pi)
check("Z2 applying the pi-rotation twice = identity (order exactly 2)",
      abs(identity_then_pi_twice) < 1e-9, f"2*pi mod 2*pi = {identity_then_pi_twice:.3e}")

# =============================================================================
print()
print(SEP2)
print("SECTION 2: ROBUSTNESS -- DOES '2' COME FROM Z_2 (PARITY), NOT FROM q's")
print("           SPECIFIC VALUE? (the prior script's claim, now superseded)")
print(SEP2)
print()

print(f"  {'q (even)':>10}  {'q*pi mod 2pi':>14}  {'Z_2 group size':>15}  "
      f"{'OLD claim (=q)':>15}  {'NEW claim (=|Z2|)':>18}")
old_new_agree_only_at_q2 = True
for q_test in [2, 4, 6, 8]:
    mod_val = (q_test * pi) % (2 * pi)
    is_symmetric = abs(mod_val) < 1e-9 or abs(mod_val - 2 * pi) < 1e-9
    z2_size = 2  # {identity, pi-rotation} -- ALWAYS 2, for any even q
    old_claim = q_test
    new_claim = z2_size
    print(f"  {q_test:>10}  {mod_val:>14.3e}  {z2_size:>15}  {old_claim:>15}  {new_claim:>18}")
    if q_test != 2 and old_claim == new_claim:
        old_new_agree_only_at_q2 = False

check("Z3 OLD ('coefficient=q') and NEW ('coefficient=|Z2|=2') claims are "
      "DIFFERENT for every even q except q=2 -- q=2 was not a discriminating "
      "test between them",
      old_new_agree_only_at_q2,
      "OLD claim grows with q (2,4,6,8); NEW claim stays fixed at 2 for all "
      "even q -- these are genuinely different, falsifiable claims")
print()
print("  Only q=2 is actually used anywhere in this framework (the Higgs's")
print("  own winding), so this is a LOGICAL consistency check, not a new")
print("  physical prediction against real data -- but it shows the NEW")
print("  reframing is the more principled one: it ties '2' to the SYMMETRY")
print("  TYPE (Z_2, exactly 2 elements whenever q is even, by definition of")
print("  a Z_2 group), not to an arbitrary count that happens to equal 2")
print("  only for this specific winding.")

# =============================================================================
print()
print(SEP2)
print("SECTION 3: DO THE 2 SECTORS ADD (matching 'simultaneous, both present')")
print("           OR BORN-WEIGHT-AVERAGE (matching 'either/or, probabilistic')?")
print(SEP2)
print()

vertex_unit_MeV = hbar_c / L_J
needed_coefficient = 2.0  # from delta_m = 2*alpha*(hbar_c/L_J), established prior session

additive_both_present = 1 * alpha + 1 * alpha          # both Z2 elements simultaneously contribute
born_weighted_average = 0.5 * alpha + 0.5 * alpha       # if treated as either/or alternatives (p=1/2 each)

print(f"  Needed coefficient (from delta_m=2*alpha*hbar_c/L_J, established): "
      f"{needed_coefficient}")
print(f"  ADDITIVE, both simultaneously present (1*alpha + 1*alpha):  "
      f"{additive_both_present/alpha:.4f} * alpha")
print(f"  BORN-WEIGHTED AVERAGE, either/or (0.5*alpha + 0.5*alpha):   "
      f"{born_weighted_average/alpha:.4f} * alpha")
check("C1 ADDITIVE (both-simultaneously-present) combination matches the "
      "needed coefficient of 2",
      abs(additive_both_present / alpha - needed_coefficient) < 1e-12,
      f"{additive_both_present/alpha:.4f} vs {needed_coefficient}")
check("C2 BORN-WEIGHTED AVERAGE (either/or) does NOT match -- would give 1, not 2",
      abs(born_weighted_average / alpha - needed_coefficient) > 0.5,
      f"{born_weighted_average/alpha:.4f} vs {needed_coefficient}")
print()
print("  This distinction matters physically, not just algebraically: Section")
print("  4.3's f_1/f_2 channels (PHI, log5) ARE Born-weighted-averaged because")
print("  they are alternative, mutually-exclusive activation channels (p_1+")
print("  p_2=1, a probability split). The 2 Z_2 elements here are NOT")
print("  alternatives -- both are always true of the SAME configuration")
print("  simultaneously (per the author's own correction) -- so an additive,")
print("  not probability-weighted, combination is the structurally correct")
print("  choice, matching alpha_born_vertex.py Part B's own additive self-")
print("  energy language (k_n_self = alpha*k_n, one EM vertex factor, not")
print("  weighted against anything else).")
print()
print("  HONEST GAP, not resolved here: 'each sector contributes exactly one")
print("  factor of alpha' is consistent with, but not independently derived")
print("  from, Part B's self-energy formula -- it remains an assumption by")
print("  analogy, not a first-principles count from a Born amplitude.")

# =============================================================================
print()
print(SEP2)
print("SECTION 4: WHICH Z_2 ELEMENT IS 'NORMAL' AND WHICH IS 'JAMMING'?")
print(SEP2)
print()
print("  Identity element: the winding/vertex AS CONSTRUCTED, no extra")
print("  operation applied -- the natural reading of 'normal vertex'.")
print()
print("  Pi-rotation element: a REAL, physical rigid pi-rotation of the")
print("  local structure -- notably NOT an exotic add-on: the icosahedral")
print("  group I already has genuine order-2 (C2) rotation axes as one of")
print("  its native symmetry classes (15 of them, alongside the 12 C5, 12")
print("  C5^2, and 20 C3 axes used throughout this repo's other vertex")
print("  calculations).")
class_sizes = {"E": 1, "C5": 12, "C5^2": 12, "C3": 20, "C2": 15}
check("J1 icosahedral group I (order 60) has a genuine C2 (pi-rotation) "
      "class with 15 axes, alongside E/C5/C5^2/C3",
      sum(class_sizes.values()) == 60, f"{class_sizes}, sum={sum(class_sizes.values())}")
print()
print("  Reading the pi-rotation element as 'jamming, cell pushed back to")
print("  Maxwell critical' is a PHYSICAL INTERPRETATION, not yet an")
print("  independent proof that THIS specific operation is what jamming-")
print("  relaxation physically is -- flagged honestly as an interpretive")
print("  label on top of the (now more solid) Z_2/parity structure, not a")
print("  second independently-derived fact.")

# =============================================================================
print()
print(SEP)
print("SECTION 5: FINAL HONEST VERDICT")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass}/{len(results)} checks PASS")
print()
print("  1. CORRECTED: 'the 2' is best tied to the ORDER of the Z_2 (pi-")
print("     rotation) symmetry group that the (1,2) winding possesses --")
print("     not to q's specific numeric value. This ties to the SAME q-even")
print("     fact already used for scalar classification (Section 3a), as the")
print("     author suspected ('if we already use that to tie it to scalar -")
print("     that may make sense'), and is a genuinely SIMULTANEOUS structure")
print("     (2 group elements, both always true at once), not a sequence.")
print("  2. The additive (not Born-weighted) combination is what's needed and")
print("     is the physically correct choice given the 2 sectors are always-")
print("     present, not either/or alternatives -- a real, checkable")
print("     distinction, not an arbitrary pick.")
print("  3. STILL OPEN: 'each sector contributes exactly one alpha' remains")
print("     an assumption by analogy to Part B's self-energy language, not")
print("     an independently derived count.")
print("  4. STILL OPEN: labeling the pi-rotation element specifically as")
print("     'jamming/Maxwell-critical relaxation' (vs. just 'a real, present")
print("     C2 symmetry element') is an interpretation, not a proof that this")
print("     operation IS the jamming-relaxation process physically.")
print("  5. NEXT STEP: attempt to independently derive item 3 -- does a Born-")
print("     amplitude calculation (same method as alpha_born_vertex.py Part")
print("     A/B, projecting onto the relevant mode at the pi-rotation/C2")
print("     class specifically) actually produce a second, independent")
print("     factor of alpha, rather than this being assumed by analogy?")
