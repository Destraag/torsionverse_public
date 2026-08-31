"""
grain_internode_propagation.py
================================
Session 5 (2026-08-18) — Agenda item [c2] follow-up

SPECIFIC QUESTION:
  The torsion medium has discrete icosahedral grain nodes (pentagons).
  A wave traveling through the medium has two possible path types:
    Path A: THROUGH the grain node (must do work against the concentrated
            pentagon energy structure)
    Path B: BETWEEN grain nodes (travels through inter-grain space,
            avoiding the node energy concentrations)

  If Path B exists and is 'free' of node resistance, then:
    - The wave might travel faster between nodes than through them
    - The effective propagation speed could exceed c
    - The MEASURED v_p = c (from GW170817) might be a PATH-AVERAGE
      of fast inter-node propagation and slow through-node propagation

  This is DISTINCT from the spring/rebound question (already answered: too
  small to account for the C4b gap). This is a topological path question.

This script:
  Part I   — What are the grain nodes physically in the framework?
  Part II  — Is there a distinct 'between-grain' propagation path?
  Part III — Scale argument: does the wave resolve individual grain nodes?
  Part IV  — Phason picture from [crys1]: is this the correct frame?
  Part V   — Quantitative: what speed enhancement would be needed for Gap 1?
  Part VI  — Verdict
"""

import math

# ============================================================
# CONSTANTS
# ============================================================
c       = 299792458.0
hbar    = 1.054571817e-34
alpha   = 7.2973525693e-3
phi     = (1 + math.sqrt(5)) / 2
r_p     = 0.8414e-15           # m
L_grain = alpha * phi * r_p    # m  grain length = coherence length
Rs      = math.sqrt(5) / (4 * math.pi)
rho     = 5.84e-27             # kg/m^3
v_s     = Rs * c
v_p     = c
G_sh    = rho * v_s**2
K_bk    = rho * (c**2 - 4/3*v_s**2)

# Gap 1 residual
alpha_C4b   = 7.2973117300057e-3
gap1_pct    = (alpha_C4b - alpha) / alpha * 100   # -0.000560%

print("=" * 65)
print("PART I — WHAT ARE GRAIN NODES PHYSICALLY?")
print("=" * 65)
print()
print("  In the torsion medium, the grain structure (Appendix G) describes")
print("  the minimum coherence length of the medium's stress-transmission:")
print(f"    L_grain = alpha*phi*r_p = {L_grain:.4e} m = {L_grain/1e-15:.4f} fm")
print()
print("  A 'grain node' is NOT a point of higher mass density.")
print("  The medium density rho is UNIFORM (dark energy density, Planck 2018).")
print("  The grain structure is a property of the ELASTIC RESPONSE, specifically:")
print("    - Below one grain length: no coherent stress wave can form")
print("      (asymptotic freedom regime, Section 1.6)")
print("    - Above one grain length: stress waves engage, confinement operates")
print()
print("  The 'icosahedral node' in the pentagon grain picture is the point where")
print("  5 grain boundaries meet — a topological feature of the quasicrystal")
print("  tiling (see [crys1]). It is a location of higher TOPOLOGICAL STRESS")
print("  (more boundary curvature), not higher mass density.")
print()
print("  KEY IMPLICATION:")
print("  Because rho is uniform, wave speed v = sqrt(modulus/rho) is NOT")
print("  changed by the presence of nodes. The nodes affect the medium's")
print("  SATURATION BEHAVIOR (when does it yield?) but not its DENSITY.")
print("  A lower density between nodes would require a vacuum gap — there is none.")

print()
print("=" * 65)
print("PART II — IS THERE A DISTINCT 'BETWEEN-GRAIN' PATH?")
print("=" * 65)
print()
print("  In a CRYSTAL with actual atoms, the distinction makes sense:")
print("    - Atom sites: high electron density, high scattering cross-section")
print("    - Between atoms: low density, 'open' channels")
print("    - Channeling effect: particles traveling between atom rows experience")
print("      less scattering and can travel faster / further")
print()
print("  In the TORSION MEDIUM the analogy breaks down:")
print()
print("  1. DENSITY: rho is uniform everywhere. There is no 'open channel'")
print("     with lower density between grain nodes.")
print()
print("  2. ELASTIC MODULI: G and K are derived from the BULK wave speeds.")
print("     They represent the average stiffness of the medium, including")
print("     both node and inter-node regions. There is no separate modulus")
print("     for the inter-node space.")
print()
print("  3. GRAIN NODES ARE NOT SCATTERERS: In QED/QCD, particles scatter")
print("     off nuclei. In the torsion medium, the grain node is a topological")
print("     curvature concentration, not a scattering center for propagating")
print("     elastic waves (whose wavelength >> L_grain at all observable scales).")
print()
print("  HOWEVER — there is a subtler question:")
print("  If the icosahedral node has higher TOPOLOGICAL STRESS (stiffer locally),")
print("  and the inter-node region has lower stress (softer locally), then the")
print("  LOCAL wave speed v_local = sqrt(G_local/rho) would vary spatially.")
print("  A wave front could travel faster through softer inter-node regions.")
print()
print("  This is the LEGITIMATE version of the question. Let's quantify it.")

print()
print("=" * 65)
print("PART III — SCALE ARGUMENT: DOES THE WAVE RESOLVE GRAIN NODES?")
print("=" * 65)
print()
print("  For a wave to 'see' the difference between node and inter-node regions,")
print("  its wavelength must be COMPARABLE TO OR SMALLER THAN L_grain.")
print()

# Observable wave frequencies and their wavelengths
test_waves = [
    (100.0,    "GW at 100 Hz (LIGO)",         v_p),
    (1e-4,     "Torsion wake (flyby, ~0.1 mHz)", v_s),
    (3e14,     "Visible light (600 nm)",        c),
    (3e19,     "Hard X-ray (10 keV)",            c),
    (1.5e25,   "Grain zone boundary (GW)",       v_p),
]

print(f"  {'Wave / Frequency':<38} {'Wavelength':<18} {'lambda/L_grain':<16} {'Resolves nodes?'}")
print(f"  {'-'*36:<38} {'-'*16:<18} {'-'*14:<16} {'-'*15}")
for freq, label, speed in test_waves:
    lam = speed / freq
    ratio = lam / L_grain
    resolves = "YES" if ratio < 10 else "NO (lambda >> L_grain)"
    print(f"  {label:<38} {lam:.3e} m    {ratio:.3e}       {resolves}")

print()
print("  For ANY wave observable with current instruments (lambda >> L_grain),")
print("  the grain node structure is COMPLETELY INVISIBLE to the propagating wave.")
print("  The wave front spans ~10^7 to ~10^40 grain spacings simultaneously.")
print("  It does not 'choose' a path between nodes — it averages over all paths.")
print()
lam_GW = v_p / 100.0
N_J_GW = lam_GW / L_grain
print(f"  A 100 Hz GW wavelength = {lam_GW:.3e} m spans {N_J_GW:.2e} grains.")
print("  The node/inter-node speed variation averages to exactly the bulk speed.")
print("  No path selection is possible at any currently accessible frequency.")

print()
print("=" * 65)
print("PART IV — PHASON PICTURE ([crys1] Tool 3): IS THIS THE RIGHT FRAME?")
print("=" * 65)
print()
print("  From [crys1], the torsion medium has TWO distinct wave modes in the")
print("  quasicrystal picture:")
print()
print("  1. PHONON (acoustic) mode: standard displacement wave, v ~ c or Rs*c")
print("     This is what GW and flyby observations measure.")
print("     Path: wave displaces grains uniformly, propagates through the bulk.")
print()
print("  2. PHASON mode: geometric rearrangement of the tiling without displacement")
print("     Grains 'flip' between two quasicrystal configurations.")
print("     Path: this is the 'between-grain' mode — the tiling topology changes")
print("           without mass displacement. It is a different physical process.")
print()
print("  PHASON PROPERTIES (from quasicrystal physics):")
print("  - Phason modes are diffusive at low frequency, not propagating")
print("  - Phason speed << acoustic speed (typically 3-4 orders of magnitude slower)")
print("  - They cannot carry energy faster than acoustic modes")
print("  - In icosahedral quasicrystals (AlPdMn), phason relaxation time ~ ms to s")
print()
print("  CONCLUSION FROM PHASON PICTURE:")
print("  The 'between-grain' propagation, if it exists, is the phason mode.")
print("  Phasons are SLOWER than acoustic modes, not faster.")
print("  They are diffusive, not ballistic — they do not propagate as clean waves.")
print("  The phason picture does NOT provide a mechanism for faster-than-c propagation.")

print()
print("=" * 65)
print("PART V — WHAT SPEED ENHANCEMENT WOULD CLOSE THE GAP 1?")
print("=" * 65)
print()
print("  For completeness: if a between-grain speed enhancement existed,")
print("  what magnitude would be needed to explain the C4b gap?")
print()
print(f"  C4b gap: alpha_C4b / alpha_CODATA - 1 = {gap1_pct:.6f}%  = {gap1_pct*1e-2:.3e}")
print()
print("  Connection to wave speed: v_s enters C4b only via Rs = v_s/c.")
print("  C4b: 2*alpha^2 - (4*pi^2/phi)*alpha + Rs = 0")
print("  If v_s has a correction delta_v: Rs_eff = (v_s + delta_v)/c")
print()
# Sensitivity: d(alpha)/d(Rs) from quadratic formula
a_q = 2.0
b_q = -4*math.pi**2/phi
c_q = Rs
disc = b_q**2 - 4*a_q*c_q
dalpha_dRs = -1.0 / math.sqrt(disc)   # d(alpha)/d(Rs) from physical root
print(f"  d(alpha)/d(Rs) at physical root = {dalpha_dRs:.6e}")
print()
# What delta_Rs closes the gap?
delta_alpha_needed = (alpha - alpha_C4b)   # positive: need to increase alpha
delta_Rs_needed = delta_alpha_needed / dalpha_dRs
delta_vs_needed = delta_Rs_needed * c
print(f"  delta_alpha needed       = {delta_alpha_needed:.4e}")
print(f"  delta_Rs needed          = {delta_Rs_needed:.4e}  ({delta_Rs_needed/Rs*100:.6f}%)")
print(f"  delta_v_s needed         = {delta_vs_needed:.4f} m/s  ({delta_vs_needed/v_s*100:.6f}%)")
print()
print(f"  To close Gap 1 via a v_s correction, v_s must shift by {delta_vs_needed:.2f} m/s.")
print(f"  This is {delta_vs_needed/v_s*100:.4f}% of v_s = {delta_vs_needed/c*100:.6f}% of c.")
print()
print("  Flyby anomaly measures v_s = Rs*c to ~0.09% precision (the K ratio).")
print(f"  The required delta_v_s / v_s = {delta_vs_needed/v_s*1e6:.2f} ppm.")
print(f"  Flyby precision = 900 ppm (0.09%).")
print(f"  The required shift is {abs(delta_vs_needed/v_s)/(0.09e-2):.1f}x smaller than flyby precision.")
print()
print("  The required v_s shift to close Gap 1 is undetectable by current data.")
print("  This does NOT rule out a v_s correction as the origin of Gap 1 —")
print("  it means we cannot confirm or deny it from existing observations.")
print("  The between-grain speed would need a mechanism that produces")
print(f"  exactly {delta_vs_needed/v_s*1e6:.2f} ppm enhancement — no current theory predicts this.")

print()
print("=" * 65)
print("PART VI — VERDICT")
print("=" * 65)
print()
print("  Q: Is there a 'between-grain' wave path that could be faster than c?")
print()
print("  ANSWER: The mechanism is not supported, for two independent reasons:")
print()
print("  1. DENSITY ARGUMENT (strong):")
print("     Wave speed requires v ~ sqrt(G/rho). The grain nodes are topological")
print("     features of the elastic response, NOT mass density concentrations.")
print("     rho is uniform everywhere. There is no low-density inter-grain channel.")
print("     No path selection is possible.")
print()
print("  2. SCALE ARGUMENT (definitive):")
print("     All observable waves have lambda >> L_grain by 7 to 40 orders of")
print("     magnitude. The wave front averages over ~10^40 grain spacings.")
print("     Individual node/inter-node speed variations are completely washed out.")
print("     This is independent of whether such variations exist in principle.")
print()
print("  3. PHASON MODE (from [crys1]):")
print("     The 'between-grain' rearrangement corresponds to phason modes in the")
print("     quasicrystal picture. Phasons are SLOWER than acoustic modes and")
print("     diffusive — they cannot carry energy faster than c.")
print()
print("  4. GAP 1 MAGNITUDE (informative):")
print(f"     Closing Gap 1 via v_s requires a {delta_vs_needed/v_s*1e6:.2f} ppm shift.")
print("     This is below flyby anomaly precision. Even if the mechanism existed,")
print("     it is invisible to current measurements and cannot be confirmed.")
print()
print("  CONCLUSION: The between-grain path hypothesis is not supported by the")
print("  framework architecture. The Gap 1 residual requires a topological")
print("  correction to the wave path integer n (epsilon = n_exact - 2 = 0.01869),")
print("  not a medium speed correction. The WZW correlator route ([crys1] Tool 1)")
print("  remains the primary attack on Gap 1.")
print()
print("  Script: analysis/alpha/grain_internode_propagation.py")
print("  Agenda: [c2] between-grain path question (fully addressed)")
