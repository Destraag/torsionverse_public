"""
maxwell_corpuscle_relay_explore.py

EXPLORATORY ONLY -- no doc claims depend on this yet. Tests three ideas
raised in discussion, all checkable without touching v_p, v_s, K/G, or any
other already-load-bearing gravitational-wave-derived quantity:

  (1) v_p and v_s were NEVER direct EM/photon measurements. Per doc_torsion.txt
      Section 3.1: v_p=c is measured via GW170817 (a GRAVITATIONAL wave arrival
      time), v_s=Rs*c via the flyby anomaly K-formula (an ORBITAL MECHANICS /
      gravitational effect). Neither is a photon experiment. So there is no
      textual requirement that "photon speed = c" must be explained via the
      SAME elastic-continuum branch (v_p) as the GW measurement -- it could be
      a SEPARATE, independently-true fact with a different mechanism, that
      simply happens to equal the same number c.

  (2) CORPUSCLE-INTRINSIC-c HYPOTHESIS: every corpuscle already established in
      this framework (gluon, tau photon-pair, W/Z before SSB) moves AT c
      intrinsically (not derived from sqrt(K/rho) or sqrt(G/rho) at all -- c
      is simply the corpuscle's own built-in speed, same as it is for gluons
      and tau photons elsewhere in the same document). A photon relaying
      vertex-to-vertex through a medium made of c-moving corpuscles inherits
      speed c for a completely different, non-conflicting reason than the
      elastic P-wave formula. This does not require resolving the P/S
      classification at all -- it sidesteps it.

  (3) EXACT DECOMPOSITION CHECK: chi(T_1g, C5) = phi is algebraically EXACTLY
      1 (longitudinal, unchanged along the C5 axis) + 2*cos(72 deg) (from the
      2 transverse components' standard 2D rotation trace). This is a genuine
      structural fact about the 3x3 rotation matrix, not an assumption. Tests
      whether using ONLY the transverse piece (2*cos72, dropping the "1")
      in the alpha vertex-coupling formula still gives something close to
      the real alpha, or whether it's wildly wrong -- diagnostic for whether
      alpha's OWN derivation structurally requires the full (longitudinal+
      transverse) character, or would tolerate a transverse-only interpretation.

Run: python analysis/quantum/maxwell_corpuscle_relay_explore.py
"""
import math

SEP = "=" * 70
phi = (1 + math.sqrt(5)) / 2
alpha_CODATA = 7.2973525693e-3

print(SEP)
print("PART 1: WERE v_p / v_s EVER DIRECT EM MEASUREMENTS?")
print(SEP)
print("  Per doc_torsion.txt Section 3.1 (verbatim, already in repo):")
print("    v_p = c        [GW170817: gravitational wave arrival, direct observation]")
print("    v_s = Rs * c   [K-formula + flyby anomaly, direct observation]")
print()
print("  Neither is a photon/EM experiment. GW170817 = gravitational wave timing.")
print("  Flyby anomaly = spacecraft trajectory (orbital mechanics/gravitational).")
print("  CONCLUSION: no existing measurement forces 'photon must be v_p branch'.")
print("  This is a genuine textual fact, not a reinterpretation.")
print()

print(SEP)
print("PART 2: EXACT DECOMPOSITION OF chi(T_1g, C5) = phi")
print(SEP)
theta = math.radians(72)
longitudinal_piece = 1.0                  # unchanged along the C5 axis
transverse_piece = 2 * math.cos(theta)    # standard 2D rotation trace, C5 plane
chi_reconstructed = longitudinal_piece + transverse_piece
print(f"  Rotation matrix about the C5 axis, in (along-axis, perp1, perp2) basis:")
print(f"    R = diag(1, [[cos72,-sin72],[sin72,cos72]])")
print(f"  trace(R) = 1 (longitudinal) + 2*cos(72deg) (2x transverse)")
print(f"    1                    = {longitudinal_piece:.6f}")
print(f"    2*cos(72deg)         = {transverse_piece:.6f}")
print(f"    sum                  = {chi_reconstructed:.6f}")
print(f"    phi (established)    = {phi:.6f}")
print(f"    match: {abs(chi_reconstructed - phi) < 1e-12}")
print()
print("  CONFIRMED: chi(T_1g,C5)=phi is EXACTLY '1 (longitudinal) + 2cos72 (transverse)'.")
print("  This is an algebraic fact about the matrix trace, not an assumption --")
print("  phi structurally MIXES a longitudinal and a transverse contribution.")
print()

print(SEP)
print("PART 3: DOES ALPHA NEED THE FULL phi, OR JUST THE TRANSVERSE PIECE?")
print(SEP)
# Reproduce the established alpha formula (doc_alpha Sec 4.5 / jobson_cell_doc J17/J24)
# using phi (full character) vs transverse_piece (2*cos72 alone) in place of the
# Born coupling constant, holding everything else identical.
def k_n_over_k_eff(f):
    """f = the coupling constant used in place of phi in the Born balance."""
    x = alpha_CODATA * f**2
    return alpha_CODATA * f * (1 - 0.75*alpha_CODATA**2) / (1 + x + x**2)

log5 = math.log(5)
def alpha_from_f(f):
    L3 = (f**3 + log5**3) / (f**2 + log5**2)
    kk = k_n_over_k_eff(f)
    # n_exact = 2 + L3*k_n/k_eff  (doc_alpha n=2 base + vertex correction)
    n_exact = 2 + L3 * kk
    return n_exact  # just report n_exact shift; comparing structurally, not refitting alpha itself

n_full = alpha_from_f(phi)
n_trans = alpha_from_f(transverse_piece)
print(f"  Using FULL phi = {phi:.6f}:              n_exact = {n_full:.8f}")
print(f"  Using TRANSVERSE-ONLY 2cos72 = {transverse_piece:.6f}:  n_exact = {n_trans:.8f}")
print(f"  Empirical target n_exact ~ 2.01869 (established, doc_alpha.txt)")
print(f"  Full-phi gap:        {(n_full-2.01869)/2.01869*100:+.4f}%")
print(f"  Transverse-only gap: {(n_trans-2.01869)/2.01869*100:+.4f}%")
print()
print("  If transverse-only is wildly worse: alpha's OWN vertex coupling needs")
print("  the FULL (longitudinal+transverse) character -- meaning alpha's mechanism")
print("  is structurally about something with all 3 T_1g components active (consistent")
print("  with it being the ELECTRON's own vertex process, not a statement about a free")
print("  massless photon's propagation specifically).")
print("  If transverse-only is comparably good: the full-phi requirement is not sharp,")
print("  and a transverse-only (2-polarization) interpretation remains open.")
