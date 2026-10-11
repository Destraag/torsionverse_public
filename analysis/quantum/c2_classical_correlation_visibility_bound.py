"""
C2 hypothesis (i): rigorous visibility bound for a classical (pseudo-thermal)
correlation substitute in Kim et al.'s OWN apparatus topology, superseding
the earlier "~12%, transferred by near-field-imaging analogy" estimate.

See notes/open_items/kim_eraser_explanation_proposal.txt, C2 section,
remaining-steps item 4.

WHY THIS SUPERSEDES THE OLD ~12% ANALOGY: Valencia et al. 2005 (quant-ph/
0408001, read in full) measured ~12% for a DIFFERENT optical configuration
-- near-field GHOST IMAGING, where a bucket detector sums incoherently over
ALL illuminated object features, diluting any one feature's point-to-point
correlation by ~1/N (N = number of features; their own Eq. 4: G ~ N +
|T(..)|^2). Kim et al.'s own idler-arm detectors D1/D2 do NOT work this way:
they coherently recombine EXACTLY the two paths (slit A, slit B) via a
dedicated final beamsplitter BSc (Kim's own eq 4: Psi(t0,t1) = A(t0,t1A) +
A(t0,t1B)), regardless of total slit count. This is a structurally
different "reference arm," so Valencia's N-dependent bucket formula does
not directly apply -- Kim's own topology is modeled here on its own terms,
directly from Kim et al.'s own published equations (confirmed via
ar5iv.labs.arxiv.org/html/quant-ph/9903047), not imported from elsewhere.

MODEL: each slit region (A, B) is treated as an independent chaotic
("thermal") amplitude, resampled every coincidence window -- the natural
classical-correlation substitute for SPDC's own "pair generated from
EITHER region A OR region B" exclusivity (mirrors the standard pseudo-
thermal technique: laser + rotating ground-glass diffuser, Martienssen &
Spiller, Am. J. Phys. 32, 919 (1964)). Both arms' measured fields are
COHERENT sums over these shared-per-shot random amplitudes (signal arm:
far-field position-dependent phase per Kim's eq 10 geometry; idler arm D1:
Kim's own eq 4 "+" recombination; D3: Kim's own eq 5 single-path tap).

MATHEMATICAL BASIS: for zero-mean, circular, jointly-Gaussian complex
amplitudes (standard thermal/pseudo-thermal statistics), the complex-
Gaussian moment theorem (Wick/Isserlis, specialized to this case) gives:
    <|E_s|^2 |E_i|^2> = <|E_s|^2><|E_i|^2> + |<E_s* E_i>|^2
(the Siegert relation). Verified numerically below (TC1-TC2) before being
used to derive Kim's own D1/D3 visibility (TC3-TC6) -- a direct Monte
Carlo coincidence-counting simulation, not just a plugged-in formula.
"""
import numpy as np

rng = np.random.default_rng(20261010)


def thermal_modes(n_samples, n_modes, rng):
    """n_modes independent, zero-mean, unit-variance complex circular
    Gaussian amplitudes per sample -- standard pseudo-thermal statistics."""
    re = rng.normal(size=(n_samples, n_modes))
    im = rng.normal(size=(n_samples, n_modes))
    return (re + 1j * im) / np.sqrt(2)


checks = []

# ---------------------------------------------------------------------
# TC1: Siegert relation, single mode -- bedrock thermal-statistics check.
# g2(0) = <I^2>/<I>^2 = 2 is the defining signature of chaotic/thermal
# light (vs. g2(0)=1 for coherent/laser light) -- confirms the random
# number generation genuinely reproduces thermal statistics.
# ---------------------------------------------------------------------
a = thermal_modes(5_000_000, 1, rng)[:, 0]
I = np.abs(a) ** 2
g2_0 = np.mean(I ** 2) / np.mean(I) ** 2
tc1 = abs(g2_0 - 2.0) < 0.02
checks.append(("TC1", "Siegert g2(0)=2 (single thermal mode)", g2_0, 2.0, tc1))

# ---------------------------------------------------------------------
# TC2: general complex-Gaussian moment factorization,
#   <|Es|^2 |Ei|^2> = <|Es|^2><|Ei|^2> + |<Es* Ei>|^2
# checked directly against Monte Carlo for an arbitrary 2-mode combination
# (not yet tied to Kim's specific geometry -- a generic sanity check).
# ---------------------------------------------------------------------
EA, EB = thermal_modes(5_000_000, 2, rng).T
Es = EA + 0.6 * EB
Ei = 0.3 * EA + EB
lhs = np.mean(np.abs(Es) ** 2 * np.abs(Ei) ** 2)
rhs = np.mean(np.abs(Es) ** 2) * np.mean(np.abs(Ei) ** 2) + np.abs(np.mean(np.conj(Es) * Ei)) ** 2
tc2 = abs(lhs - rhs) / rhs < 0.03
checks.append(("TC2", "complex-Gaussian moment factorization", lhs, rhs, tc2))

# ---------------------------------------------------------------------
# Kim et al.'s own apparatus geometry (read directly from the paper):
# a = 0.3mm (slit width), d = 0.7mm (slit separation), lambda=702.2nm.
# Visibility is independent of the overall fringe pitch (set by d/lambda/f)
# -- only the PHASE DIFFERENCE between the two paths matters -- so we scan
# the dimensionless phase directly rather than needing f explicitly.
# ---------------------------------------------------------------------
d_over_lambda = 0.7e-3 / 702.2e-9  # Kim et al.'s own d/lambda, for reference
phase = np.linspace(0, 4 * np.pi, 400)  # scans several fringe periods
N_SIM = 500_000


def coincidence_rate(idler_weights, n_samples, rng, phase):
    """
    idler_weights: complex coefficients for how the idler-arm detector in
    question combines region A, region B (e.g. [1,1] for D1's "+"
    recombination per Kim's eq 4; [1,0] for D3's which-path tap per eq 5).
    Returns G2(phase) via direct Monte Carlo coincidence counting (the
    actual 4th-moment product, averaged over samples) -- not a shortcut
    formula plug-in.
    """
    EA, EB = thermal_modes(n_samples, 2, rng).T
    Ei = idler_weights[0] * EA + idler_weights[1] * EB
    G2 = np.empty_like(phase)
    for k, ph in enumerate(phase):
        Es = EA * np.exp(-1j * ph / 2) + EB * np.exp(1j * ph / 2)
        G2[k] = np.mean(np.abs(Es) ** 2 * np.abs(Ei) ** 2)
    return G2


# ---------------------------------------------------------------------
# TC3: D1 channel (erasure, "+" recombination) -- Kim's own eq 4 topology.
# Analytic prediction V=1/3 (see tracking-file writeup for the closed-form
# algebra); confirmed here by direct Monte Carlo coincidence counting.
# ---------------------------------------------------------------------
G2_D1 = coincidence_rate([1, 1], N_SIM, rng, phase)
V_D1 = (G2_D1.max() - G2_D1.min()) / (G2_D1.max() + G2_D1.min())
tc3 = abs(V_D1 - 1 / 3) < 0.02
checks.append(("TC3", "D1 (erasure) visibility = 1/3", V_D1, 1 / 3, tc3))

# ---------------------------------------------------------------------
# TC4: D3 channel (which-path, single-path tap) -- Kim's own eq 5
# topology. Predicts FLAT (V=0) -- no fringes, matching both standard QM
# and TV's existing prediction for the which-path channels.
# ---------------------------------------------------------------------
G2_D3 = coincidence_rate([1, 0], N_SIM, rng, phase)
V_D3 = (G2_D3.max() - G2_D3.min()) / (G2_D3.max() + G2_D3.min())
tc4 = V_D3 < 0.02
checks.append(("TC4", "D3 (which-path) visibility = 0 (flat)", V_D3, 0.0, tc4))


# ---------------------------------------------------------------------
# TC5: the 1/3 ceiling is topology-dependent, not slit-count-dependent --
# confirmed by re-running the SAME "+...+" coherent-recombination topology
# with M=3 independent regions instead of 2. If V stays fixed at 1/3
# regardless of M, this confirms the ceiling comes from "the idler arm
# coherently recombines all paths at one final beamsplitter," not from
# "how many total slits exist" -- the key structural difference from
# Valencia's N-dependent (bucket-detector) topology.
# ---------------------------------------------------------------------
def coincidence_rate_M(M, n_samples, rng, phase):
    E = thermal_modes(n_samples, M, rng)
    Ei = E.sum(axis=1)  # coherent "+...+" recombination of all M regions
    k_idx = np.arange(M) - (M - 1) / 2  # symmetric spacing
    G2 = np.empty_like(phase)
    for k, ph in enumerate(phase):
        Es = (E * np.exp(-1j * ph * k_idx)).sum(axis=1)
        G2[k] = np.mean(np.abs(Es) ** 2 * np.abs(Ei) ** 2)
    return G2


G2_M3 = coincidence_rate_M(3, N_SIM, rng, phase)
V_M3 = (G2_M3.max() - G2_M3.min()) / (G2_M3.max() + G2_M3.min())
tc5 = abs(V_M3 - 1 / 3) < 0.03
checks.append(("TC5", "M=3 regions, same topology, still V=1/3", V_M3, 1 / 3, tc5))

# ---------------------------------------------------------------------
# TC6: genuine SPDC/entangled contrast -- a SINGLE quantum state's
# coherent superposition (not independent per-shot random draws) gives
# full (up to the sinc envelope) visibility, confirming WHY entanglement
# does better than classical correlation in this SAME apparatus topology.
# Modeled as a single fixed complex amplitude pair (c_A, c_B) -- the same
# amplitude pair every "shot" (the entangled state's own fixed, non-random
# phase relationship), rather than independently resampled per shot.
# ---------------------------------------------------------------------
c_A, c_B = 1.0 + 0j, 1.0 + 0j  # fixed, not random -- one quantum state
Ei_spdc = c_A + c_B
G2_spdc = np.array([
    np.abs(c_A * np.exp(-1j * ph / 2) + c_B * np.exp(1j * ph / 2)) ** 2 * np.abs(Ei_spdc) ** 2
    for ph in phase
])
V_spdc = (G2_spdc.max() - G2_spdc.min()) / (G2_spdc.max() + G2_spdc.min())
tc6 = V_spdc > 0.99
checks.append(("TC6", "SPDC/entangled case: full visibility (>99%)", V_spdc, 1.0, tc6))

print(f"{'ID':4} {'Check':46} {'measured':>10} {'expected':>10} {'result':>6}")
for cid, desc, meas, exp, ok in checks:
    print(f"{cid:4} {desc:46} {meas:10.4f} {exp:10.4f} {'PASS' if ok else 'FAIL':>6}")

n_pass = sum(1 for *_, ok in checks if ok)
print(f"\n{n_pass}/{len(checks)} PASS")
print(f"\nKim et al.'s own d/lambda = {d_over_lambda:.1f} (for reference only --")
print("visibility is independent of this ratio; it sets the fringe pitch in x,")
print("not the contrast derived above).")
print("\nRESULT: classical-correlation substitute (C2 hypothesis i), modeled on")
print("Kim et al.'s OWN D1/D3 apparatus topology (not transferred from a")
print("different experiment's geometry), gives a RIGOROUS ceiling of")
print("V = 1/3 ~= 33.3% for the erasure channels (D1/D2), independent of how")
print("many slits exist -- supersedes the earlier ~12%-by-analogy estimate.")
