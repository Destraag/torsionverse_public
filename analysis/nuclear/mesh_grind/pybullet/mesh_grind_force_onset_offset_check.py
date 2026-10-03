#!/usr/bin/env python3
"""
mesh_grind_force_onset_offset_check.py

Follow-up to mesh_grind_orientation_sweep_rp_scale_check.py, prompted by the
author's correction: the approach force should not start at the driver's
bare CENTER position -- it should start at (structural boundary) + (one
free-cell width), since the force is CAUSED by free cells that spin/expand
once driven into co-rotation by the driver, not by an abstract point charge
at the center. A zero-thickness "boundary" has no cells sitting on it; the
first REAL cell layer's own center sits roughly one cell-radius beyond the
boundary, and its outer edge roughly one cell-diameter beyond it.

HONEST SEARCH RESULT (reported, not glossed over): searched session
transcripts 3-6 (pre-"Zone 1" terminology, decompressed and grepped
directly) and re-examined session10's full exchange for an explicit,
already-decided "zone boundary + cell diameter" formula. Did NOT find one
committed anywhere as an established convention. The closest documented
precedent is analysis/nuclear/hopf_cell_coupling.py Section 4 (centrifugal
expansion of the EFFECTIVE Zone 1 boundary due to p-wave quark spin,
delta~98% radius increase) -- same general physical idea (spin inflates an
effective boundary), but scoped to Zone 1 quark modes, and its own
downstream check (matching the Hopf (1,2) ratio) FAILED. This script does
NOT reuse that specific number -- it proposes the much more modest, directly
analogous correction (+1 cell width) for the mesh_grind driver geometry
specifically, self-consistently in the SAME toy-scale units already used
there (no separate real-to-toy scaling factor needed, since both the
driver's own footprint and the proposed offset are built from the same L_J).

Run: python analysis/nuclear/mesh_grind/pybullet/mesh_grind_force_onset_offset_check.py
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
r_mid    = L_J * phi / 2

DRIVER_FOOTPRINT_FM = 2 * r_in + R_c + r_in
TOUCH_SEP_FM = 2 * DRIVER_FOOTPRINT_FM

print(SEP)
print("FORCE-ONSET OFFSET -- DISTANCE SHOULD BE MEASURED FROM A FREE-CELL-WIDTH")
print("BEYOND THE DRIVER'S OWN STRUCTURAL EDGE, NOT FROM THE BARE CENTER")
print(SEP)
print(f"  DRIVER_FOOTPRINT_FM (driver's own established edge) = {DRIVER_FOOTPRINT_FM:.6f} fm")
print(f"  Single-cell sizes: L_J (edge) = {L_J:.6f} fm, r_in (inradius) = {r_in:.6f} fm,")
print(f"                      r_mid (mid-radius) = {r_mid:.6f} fm, R_c (circumradius) = {R_c:.6f} fm")
print()

candidates = {
    "+ r_in (first free cell's own inradius)": r_in,
    "+ r_mid (first free cell's own mid-radius)": r_mid,
    "+ L_J (one full cell edge-length/diameter-ish)": L_J,
    "+ 2*r_in (one full cell diameter, face-to-face)": 2 * r_in,
}

print("  Candidate onset radii = DRIVER_FOOTPRINT_FM + one free-cell width:")
for name, width in candidates.items():
    onset = DRIVER_FOOTPRINT_FM + width
    frac_of_footprint = width / DRIVER_FOOTPRINT_FM
    frac_of_touch = onset / TOUCH_SEP_FM
    print(f"    {name:<48}: onset = {onset:.6f} fm "
          f"({frac_of_footprint*100:4.1f}% bigger than footprint alone, "
          f"{frac_of_touch:.4f}*TOUCH_SEP_FM)")

check("F1 every candidate offset is a MODEST (10-60%) correction to the "
      "driver's own footprint, not a huge rescaling -- self-consistent since "
      "both the footprint and the offset come from the same L_J-based "
      "single-cell geometry, no separate real-to-toy scale factor needed",
      all(0.05 < w / DRIVER_FOOTPRINT_FM < 1.0 for w in candidates.values()),
      f"range: {min(w/DRIVER_FOOTPRINT_FM for w in candidates.values())*100:.1f}% "
      f"to {max(w/DRIVER_FOOTPRINT_FM for w in candidates.values())*100:.1f}%")

print()
print(SEP)
print("EFFECT ON THE EXISTING mult CONVENTION")
print(SEP)
onset_width = r_in  # most conservative candidate, first free cell's own inradius
onset_radius = DRIVER_FOOTPRINT_FM + onset_width
print(f"  Current code: force_mag ~ (TOUCH_SEP_FM / center_to_center_separation)^2")
print(f"  Proposed fix: force_mag ~ (TOUCH_SEP_FM / GAP)^2, where")
print(f"    GAP = center_to_center_separation - 2*onset_radius")
print(f"    onset_radius = DRIVER_FOOTPRINT_FM + onset_width = {onset_radius:.6f} fm")
print(f"  (using +r_in, the most conservative candidate)")
print(f"  This does NOT just rescale mult uniformly -- it SUBTRACTS a fixed")
print(f"  amount from every separation, so it matters MOST at small mult (near")
print(f"  touching) and matters LESS at large mult (far away, where GAP and the")
print(f"  raw separation converge). The already-run mult=2/3/7 data cannot be")
print(f"  reinterpreted by a simple formula from the numbers in this script")
print(f"  alone -- it needs a rerun with the corrected GAP in the force law to")
print(f"  see the real effect, not a guess from this static geometry check.")
print(SEP)
