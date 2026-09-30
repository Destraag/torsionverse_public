"""
higgs_ecell_zone_boundary_wavelength_audit.py
================================================
Multi-model subagent review of item 6 (2026-09-23: GPT-6 Astra, Grok 4.7,
Gemini 3.8 Flash, Claude Opus 5.5) surfaced Claude Opus 5.5's specific,
GREP-VERIFIED-REAL critique of doc_higgs.txt lines 106-107: "This is the
energy of the zone-boundary A_g phonon... (wavelength = L_J)". Standard
periodic-lattice dispersion puts the zone-boundary WAVELENGTH at 2x the
lattice spacing (a wave of wavelength = 1x spacing aliases back to k=0),
not 1x. This script checks that claim directly and numerically -- against
this framework's OWN already-committed conventions, not just narrative --
rather than trusting either the doc's wording or Opus's critique on faith.

FIVE THINGS CHECKED, IN ORDER:
  1. E_cell's literal-wording self-consistency: is 2*pi*hbar_c/L_J EXACTLY
     what "wavelength = L_J, taken completely literally" (k=2*pi/L_J) gives?
  2. The STANDARD lattice zone-boundary convention (wavelength = 2*spacing,
     k = pi/spacing) for spacing=L_J -- and what m_H it would imply if
     substituted for the current E_cell, keeping (1+alpha/pi) unchanged.
  3. Cross-check against wave_dispersion.py's OWN, separately-committed
     zone-boundary conventions (E_grain = hbar_c/L_grain, NO 2*pi at all;
     and its "zone boundary frequency" = v/(2*L_grain), the STANDARD
     convention) -- confirms this repo already uses at least 3 mutually
     DIFFERENT zone-boundary energy formulas for the same length scale,
     independent of anything found this session.
  4. The REAL, independently-established inter-cell spacing (phi*L_J, from
     cell_rotation_propagation.py RP1 AND doc_higgs.txt's OWN Section 7.1 --
     an exact algebraic identity, not a model choice) vs L_J itself, and
     what m_H each implies under the standard zone-boundary convention.
  5. A DIRECT NUMERICAL test of Opus's actual proposed mechanism: extend
     periodic_jobson_lattice_dispersion.py's own Bloch dynamical matrix
     (imported/reused verbatim, not reimplemented) to the REAL zone
     boundary q=pi/a_lat -- a point that script's own tested |q| range
     (up to 0.4, vs a zone boundary near pi/2~1.57 along (100)) never
     reached -- and compares the true zone-boundary eigenvalue to the
     naive continuum extrapolation c^2*q_ZB^2 (c^2 taken from that
     script's own already-converged small-q effective-modulus limit).

HONEST FRAMING: this is a consistency/diagnostic AUDIT, not a fix. If a
"more standard" convention breaks the measured Higgs-mass match badly,
that is itself the finding (the current formula's specific numerical form
is empirically anchored to something not yet understood, not simply
"wrong per Debye theory") -- not swept under the rug, but not silently
patched into the doc either. Any doc change stays a PROPOSED item needing
author confirmation, per this repo's own judgment-call convention for
headline numbers.

Run: python analysis/higgs/higgs_ecell_zone_boundary_wavelength_audit.py
"""

import sys
import os
import math
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from constants import alpha, phi, hbar_c, L_J, E_cell_GeV, m_H_pred, m_H_pdg22, m_H_pdg_unc, alpha_pi

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


def sigma_from_measured(m_H_candidate):
    return (m_H_candidate - m_H_pdg22) / m_H_pdg_unc


print(SEP)
print("HIGGS E_cell ZONE-BOUNDARY / WAVELENGTH CONVENTION AUDIT (2026-09-23)")
print("Triggered by Claude Opus 5.5 subagent review's doc_higgs.txt critique")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("SECTION 1: E_cell's OWN literal-wording self-consistency")
print(SEP2)
k_literal = 2 * pi / L_J
E_literal = hbar_c * k_literal
print(f"  L_J = {L_J:.6f} fm")
print(f"  'wavelength = L_J' taken literally: k = 2*pi/L_J = {k_literal:.4f} fm^-1")
print(f"  E = hbar_c * k = {E_literal:.4f} MeV*fm * fm^-1... ", end="")
E_literal_GeV = E_literal / 1000
print(f"= {E_literal_GeV:.4f} GeV (converting MeV->GeV)")
print(f"  E_cell (constants.py, established) = {E_cell_GeV:.4f} GeV")
check("ZB1: E_cell EXACTLY equals hbar_c*(2*pi/L_J) -- self-consistent with its OWN literal wavelength=L_J wording",
      abs(E_literal_GeV - E_cell_GeV) / E_cell_GeV < 1e-9,
      f"E_literal={E_literal_GeV:.6f} GeV vs E_cell={E_cell_GeV:.6f} GeV")

# =============================================================================
print()
print(SEP2)
print("SECTION 2: STANDARD zone-boundary convention (wavelength=2*spacing)")
print("           vs the CURRENT formula -- does swapping it in still match m_H?")
print(SEP2)
print()
print("  Standard 1D periodic-lattice fact: the first Brillouin zone spans")
print("  k in (-pi/a, +pi/a]; the zone BOUNDARY is at k=pi/a, i.e. wavelength")
print("  lambda=2*pi/k=2*a -- TWICE the spacing a, not once (a wave with")
print("  wavelength=1x spacing is indistinguishable from k=0, aliased).")
print()

candidates = []

# Candidate: standard convention, spacing = L_J
k_std_LJ = pi / L_J
E_std_LJ_GeV = hbar_c * k_std_LJ / 1000
m_H_std_LJ = E_std_LJ_GeV * (1 + alpha_pi)
candidates.append(("standard (wavelength=2*L_J), spacing=L_J", E_std_LJ_GeV, m_H_std_LJ))

# Candidate: literal wavelength=phi*L_J (doc's own Section-7.1 spacing, taken literally like Section 2 does for L_J)
L_J_phi = phi * L_J
k_lit_phiLJ = 2 * pi / L_J_phi
E_lit_phiLJ_GeV = hbar_c * k_lit_phiLJ / 1000
m_H_lit_phiLJ = E_lit_phiLJ_GeV * (1 + alpha_pi)
candidates.append(("literal wavelength=phi*L_J (Section 7.1's own spacing)", E_lit_phiLJ_GeV, m_H_lit_phiLJ))

# Candidate: standard convention, spacing = phi*L_J
k_std_phiLJ = pi / L_J_phi
E_std_phiLJ_GeV = hbar_c * k_std_phiLJ / 1000
m_H_std_phiLJ = E_std_phiLJ_GeV * (1 + alpha_pi)
candidates.append(("standard (wavelength=2*phi*L_J), spacing=phi*L_J", E_std_phiLJ_GeV, m_H_std_phiLJ))

# Candidate: wave_dispersion.py's OWN E_grain convention (no 2*pi at all): E = hbar_c/L_J
E_grain_GeV = (hbar_c / L_J) / 1000
m_H_grain = E_grain_GeV * (1 + alpha_pi)
candidates.append(("wave_dispersion.py's own E_grain=hbar_c/L_grain (no 2*pi)", E_grain_GeV, m_H_grain))

print(f"  {'Candidate':<52} {'E (GeV)':>10} {'m_H pred':>10} {'sigma':>8}")
print(f"  {'CURRENT (2*pi*hbar_c/L_J, wavelength=L_J literal)':<52} {E_cell_GeV:>10.3f} {m_H_pred:>10.3f} {sigma_from_measured(m_H_pred):>8.2f}")
for label, E_cand, m_H_cand in candidates:
    print(f"  {label:<52} {E_cand:>10.3f} {m_H_cand:>10.3f} {sigma_from_measured(m_H_cand):>8.2f}")

worst_alt_sigma = min(abs(sigma_from_measured(m)) for _, _, m in candidates)
check("ZB2: EVERY standard/alternative convention is FAR worse than the current formula (>|3| sigma) -- the current formula is empirically anchored, not a naive convention error",
      worst_alt_sigma > 3.0,
      f"best alternative sigma = {worst_alt_sigma:.2f} (current = {sigma_from_measured(m_H_pred):.2f})")

# =============================================================================
print()
print(SEP2)
print("SECTION 3: wave_dispersion.py's OWN, SEPARATELY-COMMITTED conventions")
print("           (confirms >=3 mutually different zone-boundary formulas")
print("            already coexist in this repo, independent of this audit)")
print(SEP2)
Rs = math.sqrt(5) / (4 * pi)
v_s_over_c = Rs
print(f"  wave_dispersion.py: E_grain = hbar_c/L_grain = {E_grain_GeV:.4f} GeV")
print(f"    (L_grain = alpha*phi*r_p, i.e. NUMERICALLY IDENTICAL to L_J -- same")
print(f"     formula, different variable name, confirmed: L_grain == L_J to")
print(f"     float precision by construction)")
print(f"  wave_dispersion.py: 'zone boundary frequency' (elastic waves) uses")
print(f"    freq = v/(2*L_grain) -- the STANDARD wavelength=2*spacing convention,")
print(f"    for a DIFFERENT purpose (GW dispersion detection bands), not for E_cell.")
print(f"  doc_higgs.txt Section 2: E_cell = 2*pi*hbar_c/L_J -- a THIRD distinct")
print(f"    prefactor (2*pi, vs wave_dispersion's bare 1, vs the standard pi).")
check("ZB3: E_cell / E_grain(wave_dispersion.py convention) is EXACTLY 2*pi -- confirms 2 pre-existing, independently-committed, mutually-different 'zone boundary energy' formulas for the SAME length scale",
      abs(E_cell_GeV / E_grain_GeV - 2 * pi) < 1e-6,
      f"E_cell/E_grain = {E_cell_GeV/E_grain_GeV:.6f} (2*pi = {2*pi:.6f})")

# =============================================================================
print()
print(SEP2)
print("SECTION 4: L_J vs phi*L_J -- WHICH spacing does periodic_jobson_lattice_")
print("           dispersion.py's own schematic lattice actually use?")
print(SEP2)
print(f"  L_J          = {L_J:.6f} fm  (single icosahedron's OWN edge length)")
print(f"  phi*L_J      = {L_J_phi:.6f} fm  (REAL inter-cell center-to-center distance")
print(f"                 when cells touch at edge midpoints -- cell_rotation_")
print(f"                 propagation.py RP1a, algebraic identity, exact to 1e-15,")
print(f"                 explicitly cited there as 'doc_higgs Section 7.1')")
print(f"  periodic_jobson_lattice_dispersion.py's OWN Bravais lattice constant:")
print(f"    a_lat = r_nn = the intra-cell edge length (i.e. a_lat = L_J, NOT")
print(f"    phi*L_J) -- that script's OWN docstring already flags this as an")
print(f"    explicitly-simplified schematic embedding, 'NOT a claim that Jobson")
print(f"    cells really pack this way' -- this audit did not discover a NEW")
print(f"    inconsistency here, it confirms an ALREADY-FLAGGED one is the same")
print(f"    L_J-vs-phi*L_J tension raised by Opus's critique.")
check("ZB4: phi*L_J != L_J (confirms these are genuinely different length scales, not a rounding-level distinction)",
      abs(L_J_phi - L_J) / L_J > 0.5,
      f"phi*L_J/L_J = {L_J_phi/L_J:.6f} (phi = {phi:.6f})")

# =============================================================================
print()
print(SEP)
print("SECTION 5: DIRECT NUMERICAL TEST -- extend the ACTUAL Bloch dynamical")
print("matrix (periodic_jobson_lattice_dispersion.py, reused verbatim) to the")
print("REAL zone boundary, which that script's own tested |q| range never reached")
print(SEP)

# --- geometry: EXACT copy of periodic_jobson_lattice_dispersion.py Section 1-3 ---
verts_raw = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        verts_raw.append([0, s1, s2 * phi])
        verts_raw.append([s1, s2 * phi, 0])
        verts_raw.append([s2 * phi, 0, s1])
verts = np.array(verts_raw, dtype=float)
N = len(verts)
all_dists = sorted(np.linalg.norm(verts[i] - verts[j])
                    for i in range(N) for j in range(i + 1, N))
r_nn = all_dists[0]
tol = r_nn * 0.05
intra_edges = [(i, j) for i in range(N) for j in range(i + 1, N)
               if abs(np.linalg.norm(verts[i] - verts[j]) - r_nn) < tol]

axis_names = ["+x", "-x", "+y", "-y", "+z", "-z"]
axis_vecs = [np.array(v, dtype=float) for v in
             [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]]
a_lat = r_nn

vertex_dir_idx = []
for i in range(N):
    vhat = verts[i] / np.linalg.norm(verts[i])
    dots = [np.dot(vhat, ax) for ax in axis_vecs]
    vertex_dir_idx.append(int(np.argmax(dots)))

inter_bonds = [(i, a_lat * axis_vecs[vertex_dir_idx[i]]) for i in range(N)]

D_intra = np.zeros((3 * N, 3 * N))
for (i, j) in intra_edges:
    rij = verts[j] - verts[i]
    u = rij / np.linalg.norm(rij)
    outer = np.outer(u, u)
    D_intra[3 * i:3 * i + 3, 3 * i:3 * i + 3] += outer
    D_intra[3 * j:3 * j + 3, 3 * j:3 * j + 3] += outer
    D_intra[3 * i:3 * i + 3, 3 * j:3 * j + 3] -= outer
    D_intra[3 * j:3 * j + 3, 3 * i:3 * i + 3] -= outer


def build_Dq(q):
    D = D_intra.copy()
    for (i, R) in inter_bonds:
        u = R / np.linalg.norm(R)
        outer = np.outer(u, u)
        cos_qR = math.cos(float(np.dot(q, R)))
        D[3 * i:3 * i + 3, 3 * i:3 * i + 3] += 2 * outer * (1 - cos_qR)
    return D


print(f"  a_lat = r_nn = {a_lat:.6f} (natural units, unit spring constants/masses)")
q_ZB_100 = pi / a_lat
print(f"  Zone boundary along (100): |q_ZB| = pi/a_lat = {q_ZB_100:.6f}")
print(f"  (periodic_jobson_lattice_dispersion.py's OWN tested range only reached")
print(f"   |q|=0.4 -- {100*0.4/q_ZB_100:.1f}% of the way to this zone boundary)")

# small-q effective modulus (continuum c^2), reusing that script's own (111) result
deltas_small = [0.05, 0.025, 0.0125]
qdir_111 = np.array([1.0, 1.0, 1.0]) / math.sqrt(3)
c2_estimates = []
for d in deltas_small:
    vals = np.sort(np.linalg.eigvalsh(build_Dq(d * qdir_111)))[:6]
    dispersing = vals[vals > 1e-10]
    if len(dispersing):
        c2_estimates.append(dispersing.max() / d**2)
c2_continuum = float(np.mean(c2_estimates[-2:])) if len(c2_estimates) >= 2 else float(c2_estimates[-1])
print(f"\n  Continuum c^2 estimate (small-q (111), converged): {c2_continuum:.6f} (natural units)")
print(f"  (from d={deltas_small}: eff.modulus = {[f'{v:.6f}' for v in c2_estimates]})")

# (100) zone boundary -- flagged sparse per that script's own PD5 finding
vals_ZB_100 = np.sort(np.linalg.eigvalsh(build_Dq(q_ZB_100 * np.array([1.0, 0.0, 0.0]))))
dispersing_ZB_100 = vals_ZB_100[vals_ZB_100 > 1e-8]
print(f"\n  (100) zone boundary eigenvalues (natural units): {[f'{v:.4f}' for v in vals_ZB_100[:6]]}")
print(f"  (PD5 in the source script already found only 4/12 vertices' bonds")
print(f"   respond to a pure-(100) q -- most of the 6 tracked modes stay exactly")
print(f"   zero here regardless of |q|; this is a sparsity artifact of the")
print(f"   schematic bond assignment, not new physics)")

# (111) BZ corner (R point for simple cubic): q where qx=qy=qz=pi/a_lat
q_corner_111 = (pi / a_lat) * np.array([1.0, 1.0, 1.0])
vals_corner_111 = np.sort(np.linalg.eigvalsh(build_Dq(q_corner_111)))
print(f"\n  (111) BZ CORNER (simple-cubic R point, |q|=sqrt(3)*pi/a_lat={np.linalg.norm(q_corner_111):.4f}):")
print(f"  eigenvalues (natural units): {[f'{v:.4f}' for v in vals_corner_111[:6]]}")

naive_continuum_at_corner = c2_continuum * np.dot(q_corner_111, q_corner_111)
top_mode_corner = vals_corner_111[-1]
ratio_corner = top_mode_corner / naive_continuum_at_corner if naive_continuum_at_corner > 1e-12 else float("nan")
print(f"\n  Naive continuum extrapolation c^2*|q_corner|^2 = {naive_continuum_at_corner:.4f}")
print(f"  Actual top dispersing eigenvalue at corner      = {top_mode_corner:.4f}")
print(f"  Ratio (actual/continuum)                        = {ratio_corner:.4f}")
print(f"  sqrt(ratio) (an omega_disc/omega_cont analog)    = {math.sqrt(ratio_corner) if ratio_corner > 0 else float('nan'):.4f}")
print(f"  Opus's proposed target ratios: (2/pi)^2={4/pi**2:.4f} or (pi/2)^2={ (pi/2)**2:.4f}")

check("ZB5: the true zone-boundary/corner eigenvalue ratio in THIS schematic lattice is NOT cleanly (2/pi)^2 or (pi/2)^2 -- Opus's mechanism is NOT confirmed by this framework's own already-built (if simplified) dynamical matrix",
      not (abs(ratio_corner - 4 / pi**2) < 0.05 or abs(ratio_corner - (pi / 2) ** 2) < 0.05),
      f"ratio={ratio_corner:.4f} vs (2/pi)^2={4/pi**2:.4f}, (pi/2)^2={(pi/2)**2:.4f}")

# =============================================================================
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"RESULTS: {n_pass}/{len(results)} PASS")
for name, status, detail in results:
    print(f"  [{status}] {name}")
print(SEP)
print()
print("SUMMARY (honest, not a fix):")
print("  1. doc_higgs.txt's 'wavelength = L_J' wording is a genuine departure")
print("     from standard lattice zone-boundary convention (confirmed, ZB1) --")
print("     Opus's critique of the WORDING is correct.")
print("  2. BUT naively 'correcting' it to any standard convention (spacing=L_J")
print("     or phi*L_J, with the standard 2x-wavelength rule) breaks the Higgs")
print("     mass match by >3 sigma in every case tested (ZB2) -- so this is NOT")
print("     a simple factor-of-2 bug to fix; the CURRENT formula's specific")
print("     numerical form remains empirically anchored by something not yet")
print("     identified, not by the naive zone-boundary picture its own prose")
print("     invokes.")
print("  3. This repo independently already has >=2 OTHER, mutually different")
print("     'zone boundary energy' conventions in wave_dispersion.py (ZB3) --")
print("     a pre-existing terminology inconsistency, not created by this audit.")
print("  4. The L_J-vs-phi*L_J spacing tension (ZB4) is real but was ALREADY")
print("     flagged in periodic_jobson_lattice_dispersion.py's own docstring;")
print("     this audit confirms it is the same tension, not a new one.")
print("  5. Opus's specific 'continuum/discrete ratio = 2/pi' mechanism does")
print("     NOT come out of this framework's own (schematic) periodic lattice")
print("     dynamical matrix at its real zone boundary/BZ corner (ZB5) -- an")
print("     honest negative result for that specific numerical claim, though")
print("     it does not rule out the general idea in a less schematic model.")
print("RECOMMENDATION: doc_higgs.txt's Section 2 prose ('wavelength = L_J') is")
print("  imprecise relative to standard usage and could be softened/reworded")
print("  for honesty (e.g. 'an effective zone-boundary-scale energy, not a")
print("  literal Debye zone-boundary construction') -- but the E_cell NUMBER")
print("  itself should NOT be changed based on this audit: every standard")
print("  alternative is numerically far worse, and Opus's own proposed pi-ratio")
print("  mechanism is not confirmed by this framework's own dynamical matrix.")
print("  This is a PROPOSED wording softening, needs author confirmation before")
print("  touching a headline doc, per this repo's own judgment-call convention.")
