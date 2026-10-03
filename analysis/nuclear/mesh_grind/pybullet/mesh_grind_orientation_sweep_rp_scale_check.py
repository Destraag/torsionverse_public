#!/usr/bin/env python3
"""
mesh_grind_orientation_sweep_rp_scale_check.py

Answers the "NOT yet connected" gap flagged after the orientation-sweep N=20
run (/memories/repo/mesh-grind-chirality-derivation.md, 2026-10-01 entries):
that run's "mult" is in units of the toy 12-cell driver's OWN tiny touching
distance (TOUCH_SEP_FM), not r_p -- so "signal vanishes by mult=7" could not
yet be read as a statement about real physical reach relative to r_p (the
distance scale the bound-neutron exact-distance investigation,
bound_neutron_exact_distance_solve.py, independently found is the physically
relevant range for a real bound pair, and the scale the author explicitly
flagged as the target reach: "the pull should extend meaningfully all the
way out to the EDGE OF ZONE 3 (r_p scale)").

This script does ONLY the pure-algebra conversion -- no new simulation, no
new physics, just expressing the ALREADY-RUN sweep's own geometry constants
(copied verbatim from mesh_grind_pybullet_two_driver_orientation_sweep_test.py)
in units of r_p, to see whether the observed signal-vanishing point (between
mult=3 and mult=7) is a small or large fraction of r_p.

HONEST SCOPE LIMIT (stated up front, not a late caveat): DRIVER_FOOTPRINT_FM
is the 12-cell TOY driver's own footprint, not a realistically-sized Zone 2
structure (~1.8e4-4e4 cells per zone2_cell_count_estimate.py). TOUCH_SEP_FM
scales with DRIVER_FOOTPRINT_FM, so a realistically-sized driver would push
TOUCH_SEP_FM itself out further in absolute fm -- this script answers "how
far did the ALREADY-RUN toy test actually reach, in r_p units," not "how far
would a realistically-sized driver's force reach." That bigger question is a
separate, much more expensive rebuild, not attempted here.

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_orientation_sweep_rp_scale_check.py
"""
import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP = "=" * 72
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, bool(cond)))
    print(f"  [{s}] {name}")
    if detail:
        print(f"         {detail}")

pi = math.pi

# ── SAME constants, copied verbatim from the orientation sweep script ───────
alpha    = 7.2973525693e-3
hbar_c   = 197.3269804
m_p      = 938.272
lambda_p = hbar_c / m_p
r_p_fm   = 4 * lambda_p
phi      = (1 + math.sqrt(5)) / 2
L_J      = alpha * phi * r_p_fm
R_c      = L_J * math.sqrt(1 + phi**2) / 2
r_in     = L_J * phi**2 / (2 * math.sqrt(3))

DRIVER_FOOTPRINT_FM = 2 * r_in + R_c + r_in
TOUCH_SEP_FM = 2 * DRIVER_FOOTPRINT_FM

print(SEP)
print("TOY-DRIVER SCALE vs r_p -- WHAT DOES 'mult' ACTUALLY MEAN IN r_p UNITS?")
print(SEP)
print(f"  r_p_fm              = {r_p_fm:.6f} fm")
print(f"  L_J                 = {L_J:.6f} fm")
print(f"  DRIVER_FOOTPRINT_FM = {DRIVER_FOOTPRINT_FM:.6f} fm")
print(f"  TOUCH_SEP_FM        = {TOUCH_SEP_FM:.6f} fm  (mult=1 separation)")
print()

footprint_over_rp = DRIVER_FOOTPRINT_FM / r_p_fm
touch_over_rp = TOUCH_SEP_FM / r_p_fm
print(f"  DRIVER_FOOTPRINT_FM / r_p = {footprint_over_rp:.4f}")
print(f"  TOUCH_SEP_FM / r_p        = {touch_over_rp:.4f}")
check("S1 the toy driver's own footprint is much smaller than r_p (expected: "
      "a 12-cell toy proxy, not a realistic ~1.8e4-cell Zone 2 driver)",
      footprint_over_rp < 0.5, f"footprint/r_p = {footprint_over_rp:.4f}")

print()
print("  Separation actually tested, in BOTH fm and r_p units:")
for mult in (2.0, 3.0, 7.0):
    sep_fm = mult * TOUCH_SEP_FM
    sep_over_rp = sep_fm / r_p_fm
    print(f"    mult={mult:4.1f}:  separation = {sep_fm:.5f} fm  = {sep_over_rp:.4f} * r_p")

mult_at_one_rp = r_p_fm / TOUCH_SEP_FM
print()
print(f"  mult value corresponding to EXACTLY 1*r_p separation: {mult_at_one_rp:.3f}")
check("S2 mult=7 (where the hetero touch-rate signal vanished, 100% -> 0%) "
      "corresponds to a separation that is LESS than 1*r_p -- i.e. the signal "
      "was already gone before reaching the Zone-3-edge reach target, not far "
      "beyond it",
      7.0 * TOUCH_SEP_FM < r_p_fm,
      f"mult=7 separation = {7.0*TOUCH_SEP_FM/r_p_fm:.4f}*r_p; "
      f"mult=1*r_p would be mult={mult_at_one_rp:.3f}")
check("S3 even mult=3 (where hetero touch-rate was CLEANEST, 100%) is itself "
      "still well inside r_p, not already at/beyond it",
      3.0 * TOUCH_SEP_FM < r_p_fm,
      f"mult=3 separation = {3.0*TOUCH_SEP_FM/r_p_fm:.4f}*r_p")

print()
print(SEP)
print("PART 2: DOES SWITCHING n=2 -> n=3 (THE OTHER ESTABLISHED ZONE-3 ")
print("CANDIDATE) HELP, OR IS THE CURRENT FORCE ALREADY THE MORE GENEROUS ONE?")
print(SEP)
print("  The current approach force IS already an inverse-SQUARE law")
print("  (force ~ (TOUCH_SEP_FM/d)^2) -- i.e. already 'n=2', just anchored at")
print("  the toy driver's own tiny touching distance rather than independently")
print("  re-derived from r_p-scale physics. Comparing candidate exponents")
print("  n=1,2,3, all normalized to the SAME magnitude at mult=1 (the touch")
print("  distance), shows how much force would remain at the r_p-equivalent")
print("  separation (mult=13.158) under each:")
print()
for n in (1, 2, 3):
    frac_at_rp = 1.0 / (mult_at_one_rp ** n)
    frac_at_mult7 = 1.0 / (7.0 ** n)
    print(f"    n={n}:  force(r_p)/force(touch) = {frac_at_rp:.6f}   "
          f"force(mult=7)/force(touch) = {frac_at_mult7:.6f}")
check("S4 n=2 (what's actually implemented) retains MORE force at r_p than "
      "n=3 (the other established candidate) would -- switching to n=3 would "
      "make the reach problem WORSE, not better",
      (1.0 / (mult_at_one_rp ** 2)) > (1.0 / (mult_at_one_rp ** 3)),
      f"n=2: {1.0/(mult_at_one_rp**2):.6f}  n=3: {1.0/(mult_at_one_rp**3):.6f}")
check("S5 n=1 (slower than either established Zone-3 candidate) would retain "
      "roughly an order of magnitude MORE force at r_p than the current n=2",
      (1.0 / (mult_at_one_rp ** 1)) / (1.0 / (mult_at_one_rp ** 2)) > 5,
      f"ratio n=1/n=2 at r_p = {(1.0/mult_at_one_rp)/(1.0/mult_at_one_rp**2):.2f}x")

print()
print(SEP)
passed = sum(1 for _, ok in results if ok)
failed = len(results) - passed
print(f"RESULT: {len(results)} checks  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
print()
print("  READING: the entire tested range (mult=2 to mult=7) sits WELL INSIDE")
print("  r_p, not straddling or exceeding it. The signal-vanishing point found")
print("  (between mult=3 and mult=7) is therefore NOT yet informative about")
print("  whether the approach force has real reach out to the r_p/Zone-3-edge")
print("  scale that matters physically -- the toy driver's force was already")
print("  dying out at a separation several times SMALLER than r_p itself, using")
print(f"  mult={mult_at_one_rp:.2f} as the reference for '1 full r_p' in this")
print("  toy driver's own units. This sharpens (does not resolve) the author's")
print("  falloff-law/reach critique: it is not just that the signal fades by")
print("  some large multiple of the driver's own touching distance -- it fades")
print("  well BEFORE that touching distance has even grown to match the real")
print("  physical scale (r_p) the mechanism needs to reach. A realistically-")
print("  sized driver (not this 12-cell toy) would have a larger own")
print("  DRIVER_FOOTPRINT_FM/TOUCH_SEP_FM to begin with, which could push the")
print("  absolute reach further out in fm -- that is the separate, NOT yet")
print("  attempted, much more expensive rebuild this script does not address.")
print()
print("  PART 2 READING: switching the approach force's exponent from n=2 to")
print("  n=3 (the OTHER established Zone-3 candidate, bound_neutron_g_factor_")
print("  distance_dependence.py Section 8) would NOT help -- n=3 falls off")
print("  FASTER with distance, retaining even LESS force at r_p than the")
print("  current n=2 already does. Of the two candidates already established")
print("  elsewhere in this framework, the current implementation already uses")
print("  the MORE generous one for long-range reach. A genuinely slower law")
print("  (e.g. n=1) would retain substantially more force at r_p, but n=1 has")
print("  no existing physical motivation anywhere in this framework -- it")
print("  would be a new, ungrounded assumption, not a reuse of an established")
print("  result, unlike n=2/n=3. This narrows the 'wrong functional form'")
print("  angle: it is NOT a simple swap between two already-derived")
print("  candidates; either the MAGNITUDE needs retuning (a real option, not")
print("  yet tried), or the falloff-law critique points toward a mechanism")
print("  with NO existing derivation in this repo (speculative), or one of")
print("  the OTHER two open angles (cascade/escape validity; whether pressure-")
print("  well attraction is the right mechanism at all) is the more productive")
print("  place to look next.")
print(SEP)
