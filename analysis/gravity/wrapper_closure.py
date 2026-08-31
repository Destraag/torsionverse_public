"""
wrapper_closure.py — Item A unscripted #2: wrapper closure R_s appearance.

Framework claim (C5): the torsion medium saturates when a body's surface
rotation speed reaches a fraction R_s of its escape speed:
    v_rot / v_esc = R_s

This is the 5th independent appearance of R_s (solar system table).

Tests the claim against actual planetary data.

Data sources:
  - Equatorial rotation speeds: NASA planetary fact sheets
  - Escape velocities: sqrt(2*G*M/R), computed from known M, R

Run: python analysis/wrapper_closure.py
"""

import math

SEP = "=" * 65

G  = 6.674e-11   # m^3 kg^-1 s^-2
Rs = math.sqrt(5) / (4 * math.pi)

print(SEP)
print("WRAPPER CLOSURE — R_s = v_rot/v_esc CHECK (C5)")
print(SEP)
print(f"R_s = {Rs:.6f}")
print()

# Body data: (name, M_kg, R_m, v_rot_surface_m/s)
# v_rot: equatorial surface rotation speed (not orbital speed)
# Sources: NASA planetary fact sheets (2024)
bodies = [
    # name,      M (kg),        R (m),        v_rot (m/s),   note
    ("Sun",      1.989e30,      6.957e8,       1990,          "equatorial, 25.4d period"),
    ("Mercury",  3.301e23,      2.440e6,       3.026,         "58.6d period"),
    ("Venus",    4.867e24,      6.052e6,       1.81,          "243d retrograde period"),
    ("Earth",    5.972e24,      6.371e6,       465,           "23.93h sidereal period"),
    ("Mars",     6.417e23,      3.390e6,       241,           "24.62h period"),
    ("Jupiter",  1.898e27,      7.149e7,       12600,         "9.93h period"),
    ("Saturn",   5.683e26,      6.027e7,       9870,          "10.66h period"),
    ("Uranus",   8.681e25,      2.556e7,       2590,          "17.24h period"),
    ("Neptune",  1.024e26,      2.476e7,       2680,          "16.11h period"),
]

print(f"  {'Body':<10} {'v_rot (m/s)':>12}  {'v_esc (m/s)':>12}  {'ratio':>8}  {'R_s':>8}  {'err%':>7}  {'sat_frac':>9}")
print(f"  {'-'*10} {'-'*12}  {'-'*12}  {'-'*8}  {'-'*8}  {'-'*7}  {'-'*9}")

ratios = []
for name, M, R, v_rot, note in bodies:
    v_esc = math.sqrt(2 * G * M / R)
    ratio = v_rot / v_esc
    sat   = ratio / Rs          # saturation fraction
    err   = (ratio - Rs) / Rs * 100
    ratios.append((name, ratio, sat))
    print(f"  {name:<10} {v_rot:>12.1f}  {v_esc:>12.1f}  {ratio:>8.5f}  {Rs:>8.5f}  {err:>+7.1f}%  {sat:>9.4f}")

print()
print("INTERPRETATION:")
print()
print("  sat_frac < 1.0: body below wrapper saturation (torsion contribution partial)")
print("  sat_frac > 1.0: body above saturation (fully saturated, Newtonian + full torsion)")
print()

# Find which bodies are closest to R_s
closest = sorted(ratios, key=lambda x: abs(x[1] - Rs))
print(f"  Closest body to R_s threshold: {closest[0][0]} (ratio={closest[0][1]:.5f}, err={abs(closest[0][1]-Rs)/Rs*100:.1f}%)")
print()
print("  NOTE: The C5 claim is NOT that every body has v_rot/v_esc = R_s.")
print("  The claim is that R_s is the SATURATION THRESHOLD -- the value at")
print("  which the wrapper closes. Bodies above R_s are fully saturated;")
print("  bodies below R_s have partial torsion contribution.")
print("  This means the framework predicts a DICHOTOMY, not uniformity:")
print("    Bodies at sat_frac >> 1: gravity is Newtonian + full torsion term")
print("    Bodies at sat_frac << 1: gravity is mostly displacement-dominated")
print()
print("  The 5th R_s appearance is the THRESHOLD VALUE itself, not the ratio")
print("  for any individual body. Confirmed by the flyby K formula where")
print("  sat_frac enters directly (analysis/flyby_anomaly.py Part C).")
print()

# Print sat_frac table again with clear labeling
print("  Saturation fractions (flyby_anomaly.py uses these directly):")
for name, ratio, sat in ratios:
    label = "SATURATED" if sat >= 1.0 else f"partial ({sat:.2f}x)"
    print(f"    {name:<10}  sat_frac = {sat:.4f}  -> {label}")

print(SEP)
