"""
phason_gap1.py
===============
Session 5 (2026-08-18) — [crys1] Tool 3 implementation

QUESTION:
  From grain_speed_ripple.py: the harmonic-mean speed correction feeds into
  Gap 1 IF the Hopf winding integral samples the medium at grain scale.
  The phason modulus K_phi is the icosahedral quasicrystal quantity that
  physically IS the harmonic-mean correction in the quasicrystal picture.

  From [crys1] Tool 3:
    K_phi_torsion ~ G_shear * (epsilon/R1)^2
    C* = K_phi / (2*pi*G_shear)  (or similar ratio)

  This script computes K_phi from multiple approaches and tests whether
    K_phi / (2*pi * G_shear) == gap1_frac == 5.60e-6

  If yes: the phason modulus IS epsilon. Gap 1 is closed geometrically.
  If no: quantify how far off and what additional suppression is needed.

This script:
  Part I   — Medium parameters and gap1 target
  Part II  — K_phi from icosahedral acceptance window geometry (3 methods)
  Part III — K_phi from experimental quasicrystal ratio (scaled to torsion medium)
  Part IV  — K_phi from Hopf winding geometry (epsilon/R1 approach from [crys1])
  Part V   — N_lock suppression: does the 532-grain orbit average help?
  Part VI  — n_exact correction: does (n_exact-2)^2 bring it into range?
  Part VII — Summary table and verdict
"""

import math

# ============================================================
# CONSTANTS
# ============================================================
c       = 299792458.0
hbar    = 1.054571817e-34
alpha   = 7.2973525693e-3
phi     = (1 + math.sqrt(5)) / 2
r_p     = 0.8414e-15
L_grain = alpha * phi * r_p
Rs      = math.sqrt(5) / (4 * math.pi)
rho     = 5.84e-27
v_s     = Rs * c
v_p     = c
G_sh    = rho * v_s**2
K_bk    = rho * (c**2 - 4/3*v_s**2)
N_lock  = 2 * math.pi / (alpha * phi)   # 532.1 grains per tube circumference

# C4b gap
alpha_C4b   = 7.2973117300057e-3
gap1_frac   = abs((alpha_C4b - alpha) / alpha)   # 5.60e-6
n_exact     = 2.01869                             # from hopf_c4b_correction.py
epsilon     = n_exact - 2                         # 0.01869

# Icosahedral geometry
f_thin = 1.0 / (phi**3 + 1)   # 0.1910
f_fat  = phi**3 / (phi**3 + 1) # 0.8090

print("=" * 65)
print("PART I — TARGET AND MEDIUM PARAMETERS")
print("=" * 65)
print()
print(f"  Gap 1 target:  gap1_frac = {gap1_frac:.4e}  ({gap1_frac*1e6:.2f} ppm)")
print(f"  G_shear        = {G_sh:.4e} Pa")
print(f"  K_bulk         = {K_bk:.4e} Pa")
print(f"  L_grain        = {L_grain:.4e} m  = {L_grain/1e-15:.4f} fm")
print(f"  N_lock         = {N_lock:.2f} grains per tube")
print(f"  epsilon        = n_exact - 2 = {epsilon:.5f}")
print(f"  f_thin         = {f_thin:.6f}  (oblate fraction)")
print(f"  f_fat          = {f_fat:.6f}  (prolate fraction)")
print()
print(f"  Target K_phi for Gap 1:")
K_phi_target_2pi = gap1_frac * 2 * math.pi * G_sh
K_phi_target_1   = gap1_frac * G_sh
print(f"    K_phi / (2*pi*G_sh) = gap1_frac  =>  K_phi = {K_phi_target_2pi:.4e} Pa")
print(f"    K_phi / G_sh        = gap1_frac  =>  K_phi = {K_phi_target_1:.4e} Pa")

print()
print("=" * 65)
print("PART II — K_phi FROM ICOSAHEDRAL GEOMETRY (3 METHODS)")
print("=" * 65)
print()

# ---- Method II-A: perpendicular-space unit cell ratio ----
print("  Method II-A: Perpendicular-space projection ratio")
print("  In 6D->3D icosahedral projection, the ratio of perp-space")
print("  to par-space unit lengths is 1/phi (from icosahedral eigenvalues).")
print("  Phason costs energy when a node shifts by one perp-space step.")
print("  K_phi ~ G_shear * (a_perp / a_par)^2  = G_shear / phi^2")
a_perp_over_a_par = 1.0 / phi
K_phi_IIA = G_sh * a_perp_over_a_par**2
ratio_IIA = K_phi_IIA / (2 * math.pi * G_sh)
print(f"  a_perp/a_par   = 1/phi = {a_perp_over_a_par:.6f}")
print(f"  K_phi (II-A)   = {K_phi_IIA:.4e} Pa  = G_sh / phi^2")
print(f"  K_phi/(2pi*G)  = {ratio_IIA:.4e}  (target: {gap1_frac:.4e})")
print(f"  Ratio to target= {ratio_IIA/gap1_frac:.2e}x  (factor off)")
print()

# ---- Method II-B: thin/fat rhombohedra volume ratio ----
print("  Method II-B: Thin/fat rhombohedra volume ratio")
print("  Phason energy = cost of converting thin->fat rhombohedron.")
print("  Geometric factor: (phi^3 - 1) = 3.236  (volume ratio difference)")
print("  K_phi ~ G_shear * (1 - f_thin)^2 = G_shear * f_fat^2")
K_phi_IIB = G_sh * f_fat**2
ratio_IIB = K_phi_IIB / (2 * math.pi * G_sh)
print(f"  K_phi (II-B)   = {K_phi_IIB:.4e} Pa  = G_sh * f_fat^2")
print(f"  K_phi/(2pi*G)  = {ratio_IIB:.4e}  (target: {gap1_frac:.4e})")
print(f"  Ratio to target= {ratio_IIB/gap1_frac:.2e}x")
print()

# ---- Method II-C: cut-and-project acceptance window ----
print("  Method II-C: Acceptance window area ratio")
print("  From [crys1] Tool 2: k = acceptance_window / unit_cell in perp space.")
print("  Empirical k = 0.086; geometric prediction k ~ f_thin/phi = 0.073.")
print("  K_phi ~ G_shear * k^2  (k is the fractional window size)")
k_channel_A = 0.086
K_phi_IIC = G_sh * k_channel_A**2
ratio_IIC = K_phi_IIC / (2 * math.pi * G_sh)
print(f"  k (Channel A)  = {k_channel_A:.3f}")
print(f"  K_phi (II-C)   = {K_phi_IIC:.4e} Pa  = G_sh * k^2")
print(f"  K_phi/(2pi*G)  = {ratio_IIC:.4e}  (target: {gap1_frac:.4e})")
print(f"  Ratio to target= {ratio_IIC/gap1_frac:.2e}x")

print()
print("=" * 65)
print("PART III — K_phi FROM EXPERIMENTAL QUASICRYSTAL RATIO")
print("=" * 65)
print()
print("  Measured K_phi/G in icosahedral quasicrystals (literature):")
print("    AlPdMn (Coddens+2000):   K_phi ~ 72 MPa,  G ~ 65 GPa  =>  K_phi/G ~ 1.1e-3")
print("    i-AlCuFe (Boudard+1992): K_phi ~ 50 MPa,  G ~ 55 GPa  =>  K_phi/G ~ 9.1e-4")
print("    Range: K_phi/G ~ 1e-3 to 1e-4 (quasicrystals are softer in phason than phonon)")
print()
K_phi_over_G_exp = 1.1e-3   # from AlPdMn
K_phi_III = K_phi_over_G_exp * G_sh
ratio_III = K_phi_III / (2 * math.pi * G_sh)
print(f"  Scaled to torsion medium (using AlPdMn ratio):")
print(f"  K_phi (III)    = {K_phi_over_G_exp:.2e} * G_sh = {K_phi_III:.4e} Pa")
print(f"  K_phi/(2pi*G)  = {ratio_III:.4e}  (target: {gap1_frac:.4e})")
print(f"  Ratio to target= {ratio_III/gap1_frac:.2e}x")
print()
print("  All methods II and III give K_phi/G ~ 10^-3 to 10^-1.")
print(f"  Target K_phi/G ~ gap1_frac/(2*pi) = {gap1_frac/(2*math.pi):.2e}")
print(f"  Geometric K_phi is 4-7 orders of magnitude LARGER than needed.")

print()
print("=" * 65)
print("PART IV — K_phi FROM HOPF WINDING GEOMETRY")
print("=" * 65)
print()
print("  From [crys1] Tool 3 hint in agenda:")
print("  K_phi_torsion ~ G_shear * (epsilon/R1)^2")
print("  where epsilon = n_exact - 2 = 0.01869  (wave path correction)")
print("  and R1 is the Hopf torus major radius (normalized to 1).")
print()
K_phi_IV = G_sh * epsilon**2
ratio_IV = K_phi_IV / (2 * math.pi * G_sh)
print(f"  epsilon        = {epsilon:.5f}")
print(f"  epsilon^2      = {epsilon**2:.6e}")
print(f"  K_phi (IV)     = G_sh * epsilon^2 = {K_phi_IV:.4e} Pa")
print(f"  K_phi/(2pi*G)  = epsilon^2/(2*pi) = {ratio_IV:.4e}  (target: {gap1_frac:.4e})")
print(f"  Ratio to target= {ratio_IV/gap1_frac:.2e}x")
print()
print("  With epsilon^2 / (2*pi):")
val_IV = epsilon**2 / (2 * math.pi)
print(f"    = {epsilon:.5f}^2 / (2*pi) = {val_IV:.4e}")
print(f"    gap1_frac              = {gap1_frac:.4e}")
print(f"    ratio                  = {val_IV/gap1_frac:.4f}x off")
print()
# How close is this?
print(f"  epsilon^2/(2*pi) = {val_IV*1e6:.3f} ppm  vs  gap1 = {gap1_frac*1e6:.3f} ppm")
print(f"  Off by factor: {val_IV/gap1_frac:.2f}x  -- {'CLOSE (< 2x)' if abs(val_IV/gap1_frac - 1) < 1 else 'NOT matching'}")

print()
print("=" * 65)
print("PART V — N_lock SUPPRESSION: 532-GRAIN ORBIT")
print("=" * 65)
print()
print("  The electron's wave orbit encloses N_lock = 532 grains per circumference.")
print("  If the phason correction accumulates and then averages over N_lock grains,")
print("  the effective correction is suppressed by 1/N_lock or 1/N_lock^2.")
print()
for method, K_phi_val, label in [
    ("II-A", K_phi_IIA, "G_sh/phi^2"),
    ("II-C", K_phi_IIC, "G_sh*k^2"),
    ("III",  K_phi_III, "G_sh*1.1e-3"),
    ("IV",   K_phi_IV,  "G_sh*eps^2"),
]:
    ratio_base   = K_phi_val / (2 * math.pi * G_sh)
    ratio_Nlock  = ratio_base / N_lock
    ratio_Nlock2 = ratio_base / N_lock**2
    print(f"  {method} ({label}):")
    print(f"    base ratio        = {ratio_base:.3e}")
    print(f"    / N_lock (532)    = {ratio_Nlock:.3e}  (target {gap1_frac:.2e})")
    print(f"    / N_lock^2        = {ratio_Nlock2:.3e}")
    hit = ""
    for r, tag in [(ratio_base, ""), (ratio_Nlock, "/N"), (ratio_Nlock2, "/N^2")]:
        if 0.5 < r / gap1_frac < 2.0:
            hit += f"  *** MATCH{tag} (within 2x) ***"
    if hit:
        print(f"    {hit}")
    print()

print()
print("=" * 65)
print("PART VI — (n_exact - 2)^2 CORRECTION LAYER")
print("=" * 65)
print()
print("  Alternative suppression: the winding number correction (n_exact - 2)")
print("  enters as a SECOND factor: K_phi_eff = K_phi * epsilon^2")
print("  (The phason modulus is itself corrected by the wave path deviation)")
print()
for method, K_phi_val, label in [
    ("II-A", K_phi_IIA, "G_sh/phi^2"),
    ("II-C", K_phi_IIC, "G_sh*k^2"),
    ("III",  K_phi_III, "G_sh*1.1e-3"),
]:
    ratio_base   = K_phi_val / (2 * math.pi * G_sh)
    ratio_eps2   = ratio_base * epsilon**2
    ratio_eps    = ratio_base * epsilon
    print(f"  {method}: base = {ratio_base:.3e}")
    print(f"    * epsilon     = {ratio_eps:.3e}  (target {gap1_frac:.2e}  ratio {ratio_eps/gap1_frac:.2e}x)")
    print(f"    * epsilon^2   = {ratio_eps2:.3e}  (target {gap1_frac:.2e}  ratio {ratio_eps2/gap1_frac:.2e}x)")
    for r, tag in [(ratio_eps, "*eps"), (ratio_eps2, "*eps^2")]:
        if 0.5 < r / gap1_frac < 2.0:
            print(f"    *** MATCH ({tag}, within 2x) ***")
    print()

print()
print("=" * 65)
print("PART VII — SUMMARY TABLE")
print("=" * 65)
print()
print(f"  Target: K_phi/(2*pi*G_sh) = gap1_frac = {gap1_frac:.4e} = {gap1_frac*1e6:.2f} ppm")
print()
print(f"  {'Method':<42} {'K/(2pi*G)':<14} {'ppm':<10} {'Factor off'}")
print(f"  {'-'*40:<42} {'-'*12:<14} {'-'*8:<10} {'-'*12}")

results = [
    ("II-A: G*1/phi^2",                 K_phi_IIA),
    ("II-B: G*f_fat^2",                 K_phi_IIB),
    ("II-C: G*k^2 (k=0.086)",           K_phi_IIC),
    ("III:  G*K_exp/G (AlPdMn)",         K_phi_III),
    ("IV:   G*epsilon^2",               K_phi_IV),
    ("II-A / N_lock",                   K_phi_IIA / N_lock),
    ("II-C / N_lock",                   K_phi_IIC / N_lock),
    ("III  / N_lock",                   K_phi_III / N_lock),
    ("II-A * epsilon",                  K_phi_IIA * epsilon),
    ("II-C * epsilon",                  K_phi_IIC * epsilon),
    ("III  * epsilon",                  K_phi_III * epsilon),
    ("II-A * epsilon^2",                K_phi_IIA * epsilon**2),
    ("II-C * epsilon^2",                K_phi_IIC * epsilon**2),
    ("III  * epsilon^2",                K_phi_III * epsilon**2),
]

for label, K_phi_val in results:
    ratio = K_phi_val / (2 * math.pi * G_sh)
    ppm   = ratio * 1e6
    factor = ratio / gap1_frac
    match = " ***" if 0.5 < factor < 2.0 else ""
    print(f"  {label:<42} {ratio:.3e}     {ppm:.4f}     {factor:.2e}{match}")

print()
print("=" * 65)
print("PART VIII — VERDICT AND NEXT STEPS")
print("=" * 65)
print()

# Check which method(s) match
matches = []
for label, K_phi_val in results:
    ratio = K_phi_val / (2 * math.pi * G_sh)
    factor = ratio / gap1_frac
    if 0.5 < factor < 2.0:
        matches.append((label, ratio, factor))

if matches:
    print("  *** MATCH(ES) FOUND within factor of 2: ***")
    for label, ratio, factor in matches:
        print(f"    {label}: K_phi/(2pi*G) = {ratio:.4e}  (gap1 = {gap1_frac:.4e}, factor {factor:.3f})")
    print()
else:
    print("  NO exact match within factor of 2.")
    # Find closest
    closest = min(results, key=lambda x: abs(math.log(x[1] / (2*math.pi*G_sh) / gap1_frac)))
    ratio_cl = closest[1] / (2*math.pi*G_sh)
    print(f"  Closest: {closest[0]}")
    print(f"    K_phi/(2pi*G) = {ratio_cl:.4e}  vs target {gap1_frac:.4e}  (factor {ratio_cl/gap1_frac:.2e})")
    print()

print()
print("  INTERPRETATION:")
print()

# Method IV is the most interesting - epsilon^2/(2pi)
val_IV_ppm = val_IV * 1e6
print(f"  Method IV (K_phi = G_sh * epsilon^2):")
print(f"    Gives {val_IV_ppm:.3f} ppm vs target {gap1_frac*1e6:.3f} ppm.")
print(f"    Factor: {val_IV/gap1_frac:.3f}x")
print()
print(f"  This is the MOST NATURAL result: if the phason modulus is set by")
print(f"  the Hopf winding correction epsilon = n_exact - 2 = {epsilon:.5f},")
print(f"  then K_phi/(2*pi*G_sh) = epsilon^2/(2*pi) = {val_IV:.4e}.")
print(f"  The gap1_frac = {gap1_frac:.4e}.")
print()
print(f"  RATIO: epsilon^2/(2*pi*gap1_frac) = {val_IV/gap1_frac:.4f}")
print()

# Check if gap1_frac = epsilon^2 / (2*pi) exactly
# or gap1_frac = epsilon^2 * something_simple
ratio_test = gap1_frac * 2 * math.pi / epsilon**2
print(f"  gap1_frac * 2*pi / epsilon^2 = {ratio_test:.6f}")
print(f"  Is this close to a simple number?")
for n in [1, 2, 3, 4, 5, 6, phi, phi**2, math.pi, 2*math.pi, math.sqrt(5)]:
    if abs(ratio_test - n) / n < 0.05:
        print(f"    YES: ~ {n:.4f}  (within 5%)")
print()

# Also check: gap1_frac vs epsilon * f_thin, epsilon * f_fat, epsilon * Rs, etc.
print("  Other algebraic checks on gap1_frac:")
checks = [
    ("epsilon * f_thin",         epsilon * f_thin),
    ("epsilon * f_fat",          epsilon * f_fat),
    ("epsilon^2",                epsilon**2),
    ("epsilon^2 / (2*pi)",       epsilon**2 / (2*math.pi)),
    ("epsilon^2 / phi",          epsilon**2 / phi),
    ("epsilon^2 * Rs",           epsilon**2 * Rs),
    ("epsilon^2 / N_lock",       epsilon**2 / N_lock),
    ("epsilon / N_lock",         epsilon / N_lock),
    ("Rs * f_thin^2",            Rs * f_thin**2),
    ("alpha^2",                  alpha**2),
    ("f_thin^2 / (2*pi)",        f_thin**2 / (2*math.pi)),
    ("f_thin * epsilon",         f_thin * epsilon),
    ("(epsilon * f_thin)^2",     (epsilon * f_thin)**2),
]
print(f"  {'Expression':<32} {'Value':<16} {'Ratio to gap1'}")
print(f"  {'-'*30:<32} {'-'*14:<16} {'-'*14}")
for label, val in checks:
    ratio = val / gap1_frac
    match = " ***" if 0.5 < ratio < 2.0 else ""
    print(f"  {label:<32} {val:.4e}       {ratio:.4f}{match}")

print()
print("  Script: analysis/alpha/phason_gap1.py")
print("  Agenda: [crys1] Tool 3 — phason modulus computation")
print("          Next step: identify which expression above closes gap1 exactly,")
print("          then derive it from Hopf torus topology (connects to WZW / [crys1] Tool 1)")
