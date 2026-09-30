"""
sphere_contact_patch_geometry_v2.py
=====================================
Correction to sphere_contact_patch_geometry.py (author caught a real zone-
identification error, 2026-09-11): that script used R=lambda_p (the Zone
1/2 jamming BOUNDARY) as "the sphere of influence." But the canonical
four-zone structure (doc_nucleus.txt Section 1.1, PS1-PS4) defines:
  ZONE 1 -- sub-cell core            (r < lambda_p)
  ZONE 2 -- Maxwell jamming BOUNDARY  (r ~ lambda_p, a thin transition)
  ZONE 3 -- CO-ROTATING shell, cells SPIN HERE  (lambda_p < r < r_p)
Zone 3 -- the actually-spinning, chirality-carrying region -- is bounded
by r_p (4*lambda_p), not lambda_p. R=lambda_p was the wrong radius for
"where do two nuclei's spinning regions first make contact" -- that
question needs R=r_p (Zone 3's OUTER edge), giving first contact at
separation 2*r_p, four times farther out than r_grind.

TWO DISTINCT CHECKPOINTS ALONG THE SAME APPROACH (not two competing
answers -- two different physical events at two different separations):
  d = 2*r_p     : Zone 3 (co-rotating, spinning cells) shells first touch.
                  Chirality-dependent mesh/grind effects can only start
                  here at the earliest -- nothing before this separation.
  d = 2*lambda_p = r_grind : Zone 2 (frozen jamming boundary) shells touch.
                  The established "hard core" number -- a LATER, more
                  extreme stage of the same approach, well inside the
                  Zone-3-overlap regime.
This script computes the contact-patch radius a(d) and cell count N(d)
using R=r_p (correct, Zone 3) across the FULL range from first Zone-3
contact down through r_grind, replacing the incomplete R=lambda_p-only
picture.

Run: python analysis/nuclear/sphere_contact_patch_geometry_v2.py
Reference: /memories/repo/mesh-grind-chirality-derivation.md,
  doc_nucleus.txt Section 1.1 (canonical four-zone structure)
"""

import sys, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

alpha = 7.2973525693e-3
phi = (1 + math.sqrt(5)) / 2
hbar_c = 197.3269804
m_p = 938.272
lambda_p_fm = hbar_c / m_p           # Zone 1/2 boundary
r_p_fm = 4 * lambda_p_fm              # Zone 3 OUTER edge (PS4)
r_grind_fm = 2 * lambda_p_fm          # established "hard core" number (Zone 2 boundaries touch)
L_J_fm = alpha * phi * r_p_fm
r_mid_fm = phi * L_J_fm / 2

SEP = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

def contact_patch_radius(d, R):
    if d >= 2 * R:
        return 0.0
    return math.sqrt(max(R**2 - (d / 2) ** 2, 0.0))

# ── Section 1: the two distinct checkpoints, both computed ───────────────────
print(SEP)
print("SECTION 1: TWO CHECKPOINTS ON ONE APPROACH -- ZONE 3 FIRST TOUCH, THEN ZONE 2")
print(SEP2)

d_zone3_touch = 2 * r_p_fm
d_zone2_touch = 2 * lambda_p_fm   # = r_grind

print(f"  lambda_p (Zone 1/2 boundary)          = {lambda_p_fm:.4f} fm")
print(f"  r_p (Zone 3 OUTER edge, = 4*lambda_p) = {r_p_fm:.4f} fm")
print()
print(f"  d = 2*r_p      = {d_zone3_touch:.4f} fm  -- Zone 3 (spinning cells) FIRST touch")
print(f"  d = 2*lambda_p = {d_zone2_touch:.4f} fm  -- Zone 2 (frozen boundary) touch = r_grind")
print(f"  Ratio: Zone-3-first-touch separation is {d_zone3_touch/d_zone2_touch:.1f}x farther out than r_grind")

check("V2-1: Zone 3 first-touch separation (2*r_p) is exactly 4x r_grind (since r_p=4*lambda_p)",
      abs(d_zone3_touch / d_zone2_touch - 4.0) < 1e-9,
      f"ratio = {d_zone3_touch/d_zone2_touch:.9f}")

# ── Section 2: contact patch radius and cell count, using the CORRECT R=r_p ──
print()
print(SEP)
print("SECTION 2: CONTACT PATCH a(d), N(d) USING R=r_p (Zone 3, the correct sphere)")
print(SEP2)
print(f"  {'d (fm)':>10}  {'regime':>28}  {'a(d) (fm)':>12}  {'N(d) ~ (a/r_mid)^2':>20}")

# Sample from first Zone-3 touch, down through r_grind, to full overlap
d_samples = [d_zone3_touch, 1.5, 1.2, 1.0, d_zone2_touch, 0.3, 0.15, 0.0]
Nd = []
for d in d_samples:
    a_d = contact_patch_radius(d, r_p_fm)
    N_d = (a_d / r_mid_fm) ** 2
    Nd.append((d, N_d))
    regime = "Zone 3 first touch" if abs(d - d_zone3_touch) < 1e-9 else (
             "= r_grind (Zone 2 touches)" if abs(d - d_zone2_touch) < 1e-9 else
             ("inside r_grind" if d < d_zone2_touch else "Zone-3-overlap, pre-r_grind"))
    print(f"  {d:>10.4f}  {regime:>28}  {a_d:>12.6f}  {N_d:>20.2f}")

check("V2-2: contact patch is exactly 0 at Zone 3 first touch (2*r_p), matching Stage 1's "
      "original single-point-contact picture, now at the CORRECT (larger) separation",
      Nd[0][1] == 0.0,
      f"N(d=2*r_p) = {Nd[0][1]}")

check("V2-3: at r_grind itself (the OLD script's 'first touch'), the Zone-3-based patch is "
      "ALREADY substantially grown -- r_grind is deep inside the Zone-3-overlap regime, not the start of it",
      Nd[4][1] > 100,
      f"N(d=r_grind={d_zone2_touch:.4f} fm) using R=r_p = {Nd[4][1]:.1f} cells "
      f"(vs the old script's N=0 there, which used the wrong radius)")

# ── Section 3: what changes vs the withdrawn v1 script ───────────────────────
print()
print(SEP)
print("SECTION 3: WHAT THIS CHANGES")
print(SEP2)
print(f"""
  v1 (WRONG radius, R=lambda_p): treated r_grind itself as the FIRST-TOUCH
  point (a=0 there), with the patch only growing for separations BELOW
  r_grind. This made r_grind look like the ONSET of contact.

  v2 (CORRECTED, R=r_p): r_grind is NOT the onset -- it is a checkpoint
  deep INSIDE an already-substantial Zone-3 contact patch that started
  forming much earlier, at d=2*r_p={d_zone3_touch:.4f} fm (4x farther out).
  By the time the approach reaches r_grind, roughly {Nd[4][1]:.0f} cells are
  already estimated to be in the Zone-3 contact patch (using the same
  area/area proxy as v1), not zero.

  This matters for the "cells are few, then scaling up as the spheres are
  forced closer" picture (author's own framing): the SCALING starts at
  Zone 3 first touch (2*r_p), not at r_grind -- r_grind is close to the
  END of that scaling range explored here, not the beginning.
""")

# ── Summary ────────────────────────────────────────────────────────────────────
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print(SEP)
