"""
mercury_precession.py
=====================
Tests the hypothesis that orbits = pressure-minimum paths in the torsion medium,
via the framework's Einstein-Cartan identification (stated in doc_orbit_pressure.txt,
Lense-Thirring / Einstein-Cartan note).

CHAIN:
  torsion medium pressure well = spacetime curvature  [framework = Einstein-Cartan]
  GR geodesic = pressure-isobar path in torsion medium  [by identification]
  Mercury perihelion precession = GR prediction  [post-Newtonian isobar rotation]
  If GR formula matches observation => pressure-isobar orbit hypothesis supported

This is INDIRECT evidence (uses GR as intermediate step). Direct evidence would
require deriving the post-Newtonian correction from the torsion medium equations
of motion. That derivation is OPEN.

CONTEXT: Three tiers of evidence for "orbit = pressure-minimum path":
  TIER 1 (PROVEN): Bohr atom -- electron settles to EM pressure minimum.
                   G3 in gravity_doc.py proves to 0.000060%.
  TIER 2 (THIS):   Mercury precession -- GR geodesic = medium pressure isobar.
                   Indirect; requires framework = Einstein-Cartan.
  TIER 3 (OPEN):   Direct derivation from torsion medium equations of motion.

Usage: python analysis/gravity/mercury_precession.py
Reference: docs/doc_orbit_pressure.txt, HYPOTHESIS block in Section 4
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# ── constants ─────────────────────────────────────────────────────────────────
# Use GM_sun directly: known ~10,000x more precisely than G or M_sun separately
# (planetary ephemerides are built on GM, not G and M independently)
GM_sun    = 1.32712440018e20   # m^3/s^2  (IAU 2012, uncertainty ~1e9 m^3/s^2)
c_SI      = 2.99792458e8       # m/s      (exact by SI)

# Mercury orbital elements (IAU / DE430 ephemeris)
a_merc    = 5.7909050e10       # m        semi-major axis
e_merc    = 0.20563069         # (dimensionless) eccentricity
T_merc_d  = 87.969257          # days     sidereal orbital period
days_per_century = 36524.25    # days

# Measured GR precession (after removing all Newtonian perturbations)
# Le Verrier (1859) residual, confirmed by modern radar ranging
precession_measured_lo = 42.56   # arcsec/century  (3-sigma lower)
precession_measured    = 43.1    # arcsec/century  (best value)
precession_measured_hi = 43.64   # arcsec/century  (3-sigma upper)

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

# ── Section 1: GR precession formula ─────────────────────────────────────────
print(SEP)
print("MERCURY PERIHELION PRECESSION — PRESSURE-ISOBAR ORBIT TEST")
print(SEP2)
print()
print("  Framework chain:")
print("    torsion medium = Einstein-Cartan [doc_orbit_pressure.txt, L-T note]")
print("    => GR geodesic = pressure-isobar path in torsion medium")
print("    => GR precession formula applies directly")
print()

# GR precession per orbit: Δφ = 6π·GM / (a·c²·(1-e²))  [radians]
one_minus_e2   = 1.0 - e_merc**2
delta_phi_rad  = 6 * math.pi * GM_sun / (a_merc * c_SI**2 * one_minus_e2)
delta_phi_arcs = delta_phi_rad * (180.0 / math.pi) * 3600.0   # arcsec per orbit

orbits_per_century = days_per_century / T_merc_d
precession_pred    = delta_phi_arcs * orbits_per_century        # arcsec/century

print(f"  1 - e^2          = {one_minus_e2:.8f}")
print(f"  Δφ per orbit     = {delta_phi_rad:.6e} rad  =  {delta_phi_arcs:.6f} arcsec/orbit")
print(f"  Orbits/century   = {orbits_per_century:.4f}")
print(f"  Precession pred  = {precession_pred:.3f} arcsec/century")
print(f"  Precession meas  = {precession_measured:.1f}  ±  0.5 arcsec/century")
err_pct = (precession_pred - precession_measured) / precession_measured * 100
print(f"  Error            = {err_pct:+.2f}%")
print()

check("M1 Mercury precession = GR formula [framework=Einstein-Cartan => medium isobar]",
      precession_measured_lo <= precession_pred <= precession_measured_hi,
      f"predicted={precession_pred:.3f}  measured={precession_measured:.1f}±0.5  err={err_pct:+.2f}%")

# ── Section 2: Connection to orbit hypothesis ─────────────────────────────────
print()
print(SEP)
print("WHAT THIS PROVES (AND DOES NOT PROVE)")
print(SEP2)
print()
print("  PROVEN by M1:")
print("    The orbit precession rate is consistent with GR geodesic motion.")
print("    Since framework = Einstein-Cartan, this is consistent with the body")
print("    tracking the torsion medium pressure isobar (post-Newtonian isobars")
print("    are slightly non-spherical, causing the isobar to rotate = precession).")
print()
print("  NOT PROVEN by M1:")
print("    The formula is GR (borrowed), not derived from torsion medium EOMs.")
print("    Direct derivation: show the torsion medium pressure field around M_sun")
print("    produces the Schwarzschild metric to post-Newtonian order.")
print("    That derivation is OPEN.")
print()

# ── Section 3: Comparison with other evidence tiers ──────────────────────────
print(SEP)
print("EVIDENCE TIERS FOR 'ORBIT = PRESSURE-MINIMUM PATH'")
print(SEP2)
print()
print("  TIER 1 -- Bohr atom [PROVEN, G3 in gravity_doc.py]:")
print("    Electron settles to EM pressure minimum = Bohr radius.")
print("    Reproduced to 0.000060% from derived m_e. No free parameters.")
print("    Direct proof for EM case; extends to gravity by scale invariance.")
print()
print("  TIER 2 -- Mercury precession [M1, THIS SCRIPT]:")
print(f"    GR geodesic = pressure-isobar path (framework = Einstein-Cartan).")
print(f"    {precession_pred:.3f} vs {precession_measured:.1f} arcsec/century ({err_pct:+.2f}%).")
print("    Indirect; relies on GR as intermediate step. Consistent, not derived.")
print()
print("  TIER 3 -- Flat rotation curves [G11 + SPARC, doc_torsion]:")
print("    Galactic bodies at different radii all move at same v (orbit flow).")
print("    Direct evidence for coherent medium entrainment = isobar tracking.")
print("    Confirmed for 153 SPARC galaxies.")
print()
print("  TIER 4 -- Direct derivation [OPEN]:")
print("    Show Keplerian orbit = pressure-minimum variational path in torsion EOMs.")
print("    Requires post-Newtonian expansion of torsion medium equations of motion.")
print()

# ── Summary ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY")
print(SEP)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"  Total checks: {len(results)}   PASS: {passed}   FAIL: {failed}")
print()
if failed == 0:
    print("  M1 PASS: Mercury precession consistent with pressure-isobar hypothesis.")
    print("  Evidence tier 2 of 4 confirmed. Direct derivation (tier 4) remains open.")
    print()
    print("  CANDIDATE FOR gravity_doc.py once direct derivation is complete.")
    print("  DO NOT add to gravity_doc.py until torsion-medium EOMs derivation exists.")
