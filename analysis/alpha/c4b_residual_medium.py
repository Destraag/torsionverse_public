"""
c4b_residual_medium.py — Does the local medium loading state explain the C4b residual?

Context:
  C4b quadratic n_exact = 2.01869; residual = n_exact - 2 = 0.01869.
  analysis/c4b_residual_scale.py found two geometric near-hits:
    Rs/pi^2       = 0.01803  (-3.5%)
    alpha*phi^2   = 0.01910  (+2.2%)

New hypothesis from medium_chains.txt saturation table:
  The Sun's saturation fraction sat_frac_Sun = (v_rot/v_esc) / R_s = 0.018.
  This is the fraction of the medium's saturation capacity loaded by the Sun's
  rotation in our local region. It is the framework's measure of the medium's
  "back-pressure" against the Sun's rotational loading — the coupling efficiency
  loss, or what the user called "kickback in the transmission."

  If the local medium state (set by the Sun) imprints on the electron's
  torus knot winding: n_effective = 2 + sat_frac_Sun, then the C4b residual
  comes from the medium's local loading level, not pure geometry.

This script:
  PART A — Compute sat_frac_Sun precisely from solar constants.
  PART B — Compare to residual: is this the "kickback"?
  PART C — Check other solar system bodies: is the Sun unique?
            (If Mercury or Venus also appear, the match is less specific.)
  PART D — Two-scale summary: PSR B1828-11 (macroscopic) + C4b (quantum)
            both exhibit n=2. What the pulsar data says vs. what is needed.
  PART E — What a forward derivation of n = 2 + sat_frac would require.

Run: python analysis/c4b_residual_medium.py
"""

import math

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
sqrt5 = math.sqrt(5)
Rs    = sqrt5 / (4 * pi)
alpha = 7.2973525693e-3
G     = 6.67430e-11          # m^3 kg^-1 s^-2

n_exact  = (4 * pi**2 / phi * alpha - Rs) / alpha**2
residual = n_exact - 2

SEP  = "=" * 65
SEP2 = "-" * 65


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("C4B RESIDUAL — LOCAL MEDIUM LOADING HYPOTHESIS")
print("Is the residual n_exact-2 = 0.01869 the Sun's saturation fraction?")
print(SEP)
print()
print(f"  Residual = n_exact - 2 = {residual:.8f}")
print()

# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART A — SUN'S SATURATION FRACTION (precise computation)")
print(SEP)
print()

# Solar constants (IAU 2015 / CODATA)
R_sun  = 6.957e8          # m  (solar radius)
M_sun  = 1.98892e30       # kg (solar mass)
P_sun  = 25.4 * 86400     # s  (sidereal rotation period at equator, days->s)

v_rot_sun  = 2 * pi * R_sun / P_sun
v_esc_sun  = math.sqrt(2 * G * M_sun / R_sun)
ratio_sun  = v_rot_sun / v_esc_sun
sat_frac_sun = ratio_sun / Rs

print(f"  Solar constants (IAU 2015):")
print(f"    R_sun  = {R_sun:.4e} m")
print(f"    M_sun  = {M_sun:.4e} kg")
print(f"    P_sun  = {P_sun:.4e} s  ({P_sun/86400:.1f} days, equatorial sidereal)")
print()
print(f"  Derived:")
print(f"    v_rot  = 2*pi*R_sun / P_sun = {v_rot_sun:.4f} m/s")
print(f"    v_esc  = sqrt(2*G*M/R)     = {v_esc_sun:.4f} m/s  ({v_esc_sun/1000:.2f} km/s)")
print(f"    v_rot/v_esc = {ratio_sun:.6f}")
print(f"    sat_frac = (v_rot/v_esc) / R_s = {ratio_sun:.6f} / {Rs:.6f}")
print(f"             = {sat_frac_sun:.8f}")
print()
print(f"  Comparison:")
print(f"    sat_frac_Sun = {sat_frac_sun:.8f}")
print(f"    residual     = {residual:.8f}")
print(f"    difference   = {sat_frac_sun - residual:+.8f}")
print(f"    % off        = {(sat_frac_sun - residual)/residual*100:+.4f}%")
print()

# Sensitivity to solar period (differential rotation: equatorial 25.4d vs polar 35d)
P_polar = 35.0 * 86400
sat_frac_polar = (2*pi*R_sun/P_polar) / v_esc_sun / Rs
print(f"  Solar differential rotation range:")
print(f"    Equatorial (25.4d): sat_frac = {sat_frac_sun:.6f}")
print(f"    Polar      (35.0d): sat_frac = {sat_frac_polar:.6f}")
print(f"    The match to residual is at the equatorial surface value.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART B — IS THIS THE 'KICKBACK IN THE TRANSMISSION'?")
print(SEP)
print()
print("  In the torsion medium framework, the saturation fraction is the")
print("  coupling efficiency: what fraction of the body's rotational loading")
print("  the medium actually absorbs vs. resists. For sub-threshold bodies")
print("  (sat_frac < 1), the medium pushes back — most of the rotational")
print("  loading doesn't fully couple through. The Sun's sat_frac = 0.018")
print("  means only 1.8% of its rotation loads the medium to saturation.")
print("  The remaining 98.2% is the medium's back-pressure — the 'kickback'.")
print()
print(f"  sat_frac_Sun (the kickback fraction)  = {sat_frac_sun:.6f}")
print(f"  residual n_exact - 2 (from C4b)       = {residual:.6f}")
print(f"  % match                                = {(1 - abs(sat_frac_sun-residual)/residual)*100:.2f}%")
print()

# Compare to the best hits from c4b_residual_scale.py
print(f"  Comparison to previous candidates (from c4b_residual_scale.py):")
print(f"    Rs/pi^2       = {Rs/pi**2:.8f}  ({(Rs/pi**2-residual)/residual*100:+.2f}%)")
print(f"    alpha*phi^2   = {alpha*phi**2:.8f}  ({(alpha*phi**2-residual)/residual*100:+.2f}%)")
print(f"    sat_frac_Sun  = {sat_frac_sun:.8f}  ({(sat_frac_sun-residual)/residual*100:+.2f}%)")
print()
sat_pct = abs(sat_frac_sun - residual)/residual*100
print(f"  sat_frac_Sun is {sat_pct:.2f}% off. For context:")
if sat_pct < abs((Rs/pi**2-residual)/residual*100):
    print(f"  This is BETTER than Rs/pi^2 ({abs((Rs/pi**2-residual)/residual*100):.2f}% off).")
elif sat_pct < abs((alpha*phi**2-residual)/residual*100):
    print(f"  This is BETTER than alpha*phi^2 ({abs((alpha*phi**2-residual)/residual*100):.2f}% off).")
else:
    print(f"  This is NOT better than the best previous candidate.")
print()
print(f"  IMPORTANT CAVEAT — what sat_frac_Sun depends on:")
print(f"    - Solar period P_sun: known to ~0.1% (differential rotation is ~40%)")
print(f"    - Solar mass M_sun: known to 1 ppm")
print(f"    - Solar radius R_sun: known to 0.01%")
print(f"    - Rs = sqrt(5)/(4*pi): exact in the framework")
print(f"    The dominant uncertainty is which rotation period to use.")
print(f"    Using the equatorial period gives the closest match.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART C — IS THE SUN SPECIAL? OTHER SOLAR SYSTEM BODIES")
print(SEP)
print()
print("  If the residual = sat_frac_local, it should be THE Sun's value,")
print("  not coincidentally matching Mercury or Jupiter. Check all bodies.")
print()

# (name, M_kg, R_m, P_days)
bodies = [
    ("Sun",          1.98892e30, 6.957e8,   25.4),
    ("Mercury",      3.3011e23,  2.4397e6,  58.65),
    ("Venus",        4.8675e24,  6.0518e6,  243.0),
    ("Earth",        5.9722e24,  6.3781e6,  1.0),
    ("Mars",         6.4171e23,  3.3895e6,  1.026),
    ("Jupiter",      1.8982e27,  7.1492e7,  0.414),
    ("Saturn",       5.6834e26,  6.0268e7,  0.444),
    ("Uranus",       8.6810e25,  2.5362e7,  0.718),
    ("Neptune",      1.0241e26,  2.4622e7,  0.671),
]

print(f"  {'Body':<12} {'v_rot/v_esc':>12}  {'sat_frac':>10}  {'diff from residual':>20}  {'% off':>8}")
print(f"  {'-'*12} {'-'*12}  {'-'*10}  {'-'*20}  {'-'*8}")

for name, M, R, P_days in bodies:
    P_s   = P_days * 86400
    v_r   = 2 * pi * R / P_s
    v_e   = math.sqrt(2 * G * M / R)
    ratio = v_r / v_e
    sf    = ratio / Rs
    diff  = sf - residual
    pct   = diff / residual * 100
    flag  = "  <--" if abs(pct) < 5 else ""
    print(f"  {name:<12} {ratio:>12.6f}  {sf:>10.6f}  {diff:>+20.6f}  {pct:>+8.2f}%{flag}")

print()
print(f"  Only the Sun's sat_frac falls within ~3% of the residual.")
print(f"  This specificity supports the Sun being the relevant body, not")
print(f"  an arbitrary coincidence across the solar system table.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART D — TWO-SCALE SUMMARY: WHERE n=2 APPEARS IN NATURE")
print(SEP)
print()
print("  SCALE 1 — Quantum (electron Hopf torus):")
print(f"    C4b quadratic: n*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0")
print(f"    Best fit integer: n=2, giving alpha error = -0.000560%")
print(f"    n_exact = {n_exact:.6f}  (integer 2 is {abs(n_exact-2)/n_exact*100:.3f}% below)")
print(f"    Physical picture: (1,2) torus knot crossing ring on Hopf torus")
print()
print("  SCALE 2 — Astrophysical (freely precessing pulsar):")
print(f"    PSR B1828-11 (Stairs, Lyne & Shemar 2000, Nature 406, 484):")
print(f"    Dominant modulation periods: P1 = ~511 days, P2 = ~256 days")
print(f"    Ratio P1/P2 = {511/256:.6f}  (within 0.07 sigma of n=2)")
print(f"    Physical picture: beam traces (1,2) torus knot in angle-space")
print(f"    per precession cycle — same topology as electron crossing ring")
print()
print("  AGREEMENT:")
print(f"    Both scales observe n=2 as the dominant integer, to the precision")
print(f"    currently measurable. The pulsar ratio uncertainty (~2.8%) is")
print(f"    3x too coarse to resolve n_exact vs n=2. The quantum value")
print(f"    (n_exact = {n_exact:.4f}) is the sharper constraint.")
print()
print("  RESIDUAL (the 0.019 above n=2):")
print(f"    Quantum:    n_exact - 2 = {residual:.6f}")
print(f"    Pulsar:     P1/P2 - 2   = {511/256-2:.6f}  (low precision end)")
print(f"    Sun sat_frac (kickback): {sat_frac_sun:.6f}")
print()
print(f"  The residual is not yet explained. The three nearest candidates:")
print(f"    Sun sat_frac  = {sat_frac_sun:.6f}  ({(sat_frac_sun-residual)/residual*100:+.2f}%)  [local medium loading]")
print(f"    alpha*phi^2   = {alpha*phi**2:.6f}  ({(alpha*phi**2-residual)/residual*100:+.2f}%)  [geometric]")
print(f"    Rs/pi^2       = {Rs/pi**2:.6f}  ({(Rs/pi**2-residual)/residual*100:+.2f}%)  [geometric]")
print(f"  None is clean enough to be definitive.")
print()


# ─────────────────────────────────────────────────────────────────────────────
print(SEP)
print("PART E — WHAT A DERIVATION OF n = 2 + sat_frac WOULD REQUIRE")
print(SEP)
print()
print("  If n_effective = 2 + sat_frac_local (where sat_frac_local is the")
print("  medium's loading state from the dominant nearby body), then:")
print()
print("  E1 — Location dependence of alpha:")
print("    n_effective varies with sat_frac_local, which varies between")
print("    solar system bodies. A measurement of alpha on Jupiter (sat_frac=1.19)")
print("    would give a DIFFERENT n_effective than on Earth (sat_frac=0.234).")
print("    This is a testable prediction if the theory is correct:")
print("      n_eff(Jupiter) = 2 + 1.187 = 3.187  (very different)")
print("      n_eff(Earth)   = 2 + 0.234 = 2.234  (measurably different?)")
print("      n_eff(Sun loc) = 2 + 0.018 = 2.018")
print()
print("    BUT: We measure alpha on Earth with precision 1.5e-10 (relative).")
print("    Has anyone measured alpha in different gravity environments?")
print("    Atomic clock comparisons at altitude test gravitational redshift,")
print("    not alpha directly. This would require an independent alpha")
print("    measurement (g-2 experiment or atom recoil) in space.")
print()
print("  E2 — Why the Sun, not Earth?")
print("    If alpha measured on Earth uses Earth's sat_frac = 0.234, then")
print("    n_eff = 2.234 and the C4b error would be much larger.")
print("    But CODATA alpha is Earth-based and matches C4b with n=2.")
print("    This implies either:")
print("      (a) The relevant loading is the Sun's (dominant body in system),")
print("          not Earth's. The electron's topology samples the medium on")
print("          the scale of the solar system's dominant torsion loading.")
print("      (b) sat_frac does NOT enter alpha this way and the Sun match")
print("          is a coincidence at the 3% level.")
print("    Option (a) would predict: alpha measured in isolation (far from")
print("    any star) would have n_eff = 2 exactly, giving C4b = CODATA.")
print()
print("  E3 — Forward derivation:")
print("    Showing n_effective = 2 + sat_frac would require:")
print("    - The Hopf linking integral in Part D2 to depend on the local")
print("      medium loading state, not just the torus geometry alone.")
print("    - This is structurally plausible: if the medium's stiffness")
print("      enters the linking integral as a boundary condition, the")
print("      local sat_frac modifies the effective winding by that amount.")
print()
print(SEP)
print("SUMMARY")
print(SEP)
print()
print("  Q: Is the C4b residual 0.01869 the 'kickback in the transmission'?")
print(f"  A: It could be. sat_frac_Sun = {sat_frac_sun:.5f} is {abs(sat_frac_sun-residual)/residual*100:.2f}%")
print(f"     off from the residual — the best physical candidate found.")
print(f"     Conceptually consistent: the medium's local loading (set by")
print(f"     the Sun's rotation) would imprint on the electron's torus knot")
print(f"     winding by the saturation fraction.")
print()
print(f"  Q: Is this completely disjoint from the 'kickback' concept?")
print(f"  A: No. The saturation fraction IS the kickback — it measures how")
print(f"     much of the rotational coupling 'doesn't get through'. The Sun's")
print(f"     value (1.8% coupling efficiency) matches the winding residual to 3%.")
print(f"     This is not clean enough to be definitive but too close to ignore.")
print()
print(f"  n=2 appears at two scales (quantum C4b + pulsar B1828-11).")
print(f"  The residual 0.019 has three candidates within 4%:")
print(f"    sat_frac_Sun ({(sat_frac_sun-residual)/residual*100:+.2f}%), alpha*phi^2 (+2.22%), Rs/pi^2 (-3.53%).")
print(f"  None is confirmed. All three warrant a forward derivation test.")
print(SEP)
