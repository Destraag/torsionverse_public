"""
higgs_unilateral_contact_rectification_check.py
==================================================
Multi-model subagent review of item 6 (2026-09-23) produced a CONVERGENT
proposal from two independent models (GPT-6 Astra, Gemini 3.8 Flash): the
medium's granular contacts are UNILATERAL (Signorini condition -- transmit
compression, never tension), so phase-averaging a breathing mode's contact
engagement over a full 2*pi cycle gives exactly max(0,cos(theta)) averaged
= 1/pi (the standard half-wave-rectified-cosine "DC value", the same fact
used to compute a rectifier circuit's average output voltage). Both models
proposed this as the source of the mass-correction's 1/pi, with the
already-established bare-alpha vertex correction (alpha_born_vertex.py Part
B: k_n_self = alpha*k_n, explicitly "no 1/pi... a CONTACT correction") as
the un-rectified starting point.

Both models self-rated LOW-MODERATE confidence: this establishes the
KINEMATIC/FORCE gate = 1/pi, but does NOT by itself prove the ENERGY
observable (a mass shift IS an energy, not a force) is weighted the same
way. This script checks that specific gap directly and honestly.

THREE THINGS CHECKED:
  1. The elementary half-wave-rectified-cosine average IS exactly 1/pi
     (checked both analytically and by numerical quadrature -- trivial,
     but should not be taken on faith either).
  2. The ENERGY-based version of the same mechanism: for a compression-
     only (Signorini) contact undergoing a KINEMATICALLY-IMPOSED cos(theta)
     oscillation (matching this framework's own established "geometry, not
     dynamics, sets the phase" convention -- cell_rotation_propagation.py
     RP3's "72-deg rotation forced by C5 geometry, no free parameter"),
     stored elastic energy U~x^2 is quadratic, not linear -- does its
     full-cycle-averaged ratio (unilateral/bilateral) ALSO come out to
     1/pi, or does the quadratic (always non-negative) nature of energy
     change the answer?
  3. Whether EITHER ratio (force-based 1/pi, or whatever the energy-based
     check gives) composes naturally with the ALREADY-established k_n_self
     = alpha*k_n bare-alpha correction to reproduce alpha/pi specifically,
     without needing any further unexplained adjustment.

Run: python analysis/higgs/higgs_unilateral_contact_rectification_check.py
"""

import math
import numpy as np

pi = math.pi
SEP = "=" * 78
SEP2 = "-" * 78
results = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


print(SEP)
print("UNILATERAL CONTACT RECTIFICATION CHECK (2026-09-23)")
print("GPT-6 Astra / Gemini 3.8 Flash subagent proposal for the alpha/pi 1/pi")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("CHECK 1: elementary half-wave-rectified-cosine average = 1/pi")
print(SEP2)
analytic_avg = 1.0 / pi

# numerical cross-check via quadrature (fine grid, mean over evenly-spaced samples)
theta = np.linspace(0, 2 * pi, 2_000_001)
rectified = np.maximum(0.0, np.cos(theta))
numeric_avg = np.mean(rectified)
print(f"  Analytic:  (1/(2*pi)) * integral_0^2pi max(0,cos(theta)) dtheta = 1/pi = {analytic_avg:.8f}")
print(f"  Numeric (2e6-point quadrature):                                        {numeric_avg:.8f}")
check("UC1: analytic 1/pi matches numerical quadrature to 1e-6",
      abs(analytic_avg - numeric_avg) < 1e-6,
      f"diff = {abs(analytic_avg - numeric_avg):.2e}")

# =============================================================================
print()
print(SEP2)
print("CHECK 2: does the ENERGY (quadratic, U~x^2) version give 1/pi too,")
print("         or does quadratic (always-non-negative) energy behave")
print("         differently from the linear (signed) force/amplitude case?")
print(SEP2)
print()
print("  Setup: theta=omega*t is a KINEMATICALLY-IMPOSED phase (matching this")
print("  framework's own 'geometry forces the phase, not dynamics' convention).")
print("  Bilateral (ordinary Hookean) contact: U_bi(theta) = U0*cos^2(theta) for")
print("  ALL theta in [0,2*pi) -- always engaged, both compression and tension.")
print("  Unilateral (Signorini, no tension): U_uni(theta) = U0*cos^2(theta) for")
print("  theta in the compressive half [-pi/2,+pi/2], ZERO otherwise (the")
print("  contact is simply not engaged during the 'tension' half -- no energy")
print("  stored, not negative energy).")
print()

# Full-cycle average of bilateral energy: <cos^2> over 2pi = 1/2 (standard)
avg_bi_analytic = 0.5
theta_full = np.linspace(0, 2 * pi, 2_000_001)
avg_bi_numeric = np.mean(np.cos(theta_full) ** 2)
print(f"  <U_bi>/U0 analytic = 1/2 = {avg_bi_analytic}")
print(f"  <U_bi>/U0 numeric  = {avg_bi_numeric:.8f}")
check("UC2a: bilateral full-cycle energy average = 1/2 (standard <cos^2>=1/2)",
      abs(avg_bi_numeric - 0.5) < 1e-6, f"numeric={avg_bi_numeric:.8f}")

# Full-cycle average of unilateral energy: only integrate over the active half,
# still normalize by the FULL 2*pi period (the contact is simply silent, not
# absent from the reference cycle).
theta_active = np.linspace(-pi / 2, pi / 2, 1_000_001)
integral_active = np.mean(np.cos(theta_active) ** 2) * pi
avg_uni_numeric = integral_active / (2 * pi)
avg_uni_analytic = (pi / 2) / (2 * pi)  # integral_{-pi/2}^{pi/2} cos^2 = pi/2 ; /(2*pi) = 1/4
print(f"\n  integral_[-pi/2,pi/2] cos^2(theta) dtheta = pi/2 (standard)")
print(f"  <U_uni>/U0 analytic = (pi/2)/(2*pi) = 1/4 = {avg_uni_analytic:.8f}")
print(f"  <U_uni>/U0 numeric  = {avg_uni_numeric:.8f}")
check("UC2b: unilateral (Signorini) full-cycle energy average = 1/4, NOT the force-case's 1/pi",
      abs(avg_uni_numeric - 0.25) < 1e-6, f"numeric={avg_uni_numeric:.8f}")

ratio_energy = avg_uni_numeric / avg_bi_numeric
print(f"\n  RATIO <U_uni>/<U_bi> = {ratio_energy:.8f}")
print(f"  Compare: force-based duty-cycle ratio (GPT-6/Gemini's claim) = 1/pi = {1/pi:.8f}")
check("UC2c: the ENERGY ratio (1/2) is NOT 1/pi -- the elegant force-based duty-cycle result does NOT carry over to the energy observable a mass correction actually needs",
      abs(ratio_energy - 0.5) < 1e-6 and abs(ratio_energy - 1 / pi) > 0.05,
      f"energy ratio = {ratio_energy:.6f} (=1/2, exactly, not 1/pi={1/pi:.6f})")

# =============================================================================
print()
print(SEP2)
print("CHECK 3: does EITHER ratio compose naturally with the ALREADY-")
print("         established bare-alpha correction to give alpha/pi cleanly?")
print(SEP2)
alpha = 7.2973525693e-3
alpha_pi_target = alpha / pi
print(f"  Target: delta_m/E_cell = alpha/pi = {alpha_pi_target:.8e}")
print()
print("  Route A (force-based duty cycle, GPT-6/Gemini's literal claim):")
print(f"    bare alpha (alpha_born_vertex.py Part B, no 1/pi) * (1/pi duty cycle)")
route_a = alpha * (1 / pi)
print(f"    = alpha/pi = {route_a:.8e}  [MATCHES target EXACTLY, by construction --")
print(f"      this is definitionally alpha/pi if the duty-cycle factor IS 1/pi,")
print(f"      which CHECK 1 confirmed for the FORCE quantity]")
check("UC3a: bare-alpha x (force duty cycle 1/pi) reproduces alpha/pi exactly (tautological given CHECK 1, not new evidence by itself)",
      abs(route_a - alpha_pi_target) < 1e-12, f"route_a={route_a:.10e}")

print()
print("  Route B (energy-based duty cycle, the physically-appropriate")
print("  observable type for a MASS/ENERGY shift):")
route_b = alpha * ratio_energy
print(f"    bare alpha * (1/2 energy duty cycle) = {route_b:.8e}")
sigma_equiv_b = abs(route_b - alpha_pi_target) / alpha_pi_target
print(f"    vs target alpha/pi = {alpha_pi_target:.8e}  ({100*sigma_equiv_b:.1f}% off)")
check("UC3b: bare-alpha x (energy duty cycle 1/2) does NOT reproduce alpha/pi (off by pi/2 ~ 57%, not a rounding-level miss)",
      sigma_equiv_b > 0.3, f"{100*sigma_equiv_b:.1f}% relative miss")

# =============================================================================
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"RESULTS: {n_pass}/{len(results)} PASS")
for name, status, detail in results:
    print(f"  [{status}] {name}")
print(SEP)
print()
print("HONEST VERDICT: the GPT-6 Astra / Gemini 3.8 Flash 'unilateral contact")
print("rectification' mechanism is mathematically real for a FORCE-like (linear,")
print("signed) quantity -- max(0,cos(theta)) really does average to exactly 1/pi,")
print("checked both analytically and numerically (CHECK 1). But a mass/energy")
print("shift is fundamentally an ENERGY (quadratic, always non-negative)")
print("observable, and CHECK 2 shows the SAME unilateral/bilateral construction")
print("applied to energy gives a full-cycle ratio of exactly 1/2 (a completely")
print("mundane number, no pi), NOT 1/pi. Route A (CHECK 3a) reproducing alpha/pi")
print("is therefore closer to definitional/tautological (it assumes the very")
print("thing being tested applies to the mass correction) than a genuine")
print("independent derivation; Route B (the more defensible energy-based")
print("version) misses by ~57%, not a rounding-level miss.")
print("CONCLUSION: NOT CONFIRMED as stated. The 1/pi duty-cycle idea is a real,")
print("elegant, and previously-unconsidered mechanism for FORCE/impulse-type")
print("quantities, but does not, as tested, transfer to the ENERGY observable")
print("alpha/pi actually needs. This narrows (does not fully rule out -- a more")
print("careful QFT-level treatment of how a duty-cycle-gated force feeds into a")
print("self-energy diagram might still recover 1/pi through a different route)")
print("but should NOT be reported as a working native derivation without that")
print("further step, which has not been attempted here.")
