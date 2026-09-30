"""
higgs_alpha_pi_import_audit.py
===============================
AUDIT: is the (1+alpha/pi) Higgs mass correction an independent torsion-medium
derivation, or an imported standard-QED coefficient wearing torsionverse labels?

CONTEXT: judgment_calls.txt "GR/SR-reproduction vs. medium-specific-prediction
sweep" flagged m_H = E_cell*(1+alpha/pi) as the same import pattern as the
already-retired GM/(r*c^2) time-dilation formula. This script does two things:

  PART A: confirms (by direct quotation, not inference) that every existing
    script in this folder that touches alpha/pi labels it as standard/
    external QED, never as independently derived from torsion-medium EOMs.

  PART B: tests whether the ONE genuinely native correction mechanism that
    exists elsewhere in this framework -- delta_n (the vertex-stiffness
    correction from the alpha derivation itself, used in higgs_me_correction.py
    for the ELECTRON mass) -- could substitute for the imported alpha/pi term
    on the HIGGS mass, either directly (bulk-regime form) or restricted to the
    l=0-only channel (sub-cell form, matching higgs_h1_residual.py's own
    bulk-vs-sub-cell distinction, since the Higgs is explicitly sub-cell,
    N_J_H < 1, while delta_n itself was derived for the electron, a bulk
    particle, N_J >> 1).

Run: python analysis/higgs/higgs_alpha_pi_import_audit.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi = math.pi
Rs = math.sqrt(5) / (4 * pi)
log5 = math.log(5)
L3 = (phi**3 + log5**3) / (phi**2 + log5**2)

SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "FAIL"
    results.append((name, s, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


print(SEP)
print("PART A: alpha/pi IS LABELED EXTERNAL/STANDARD EVERYWHERE IT APPEARS")
print(SEP2)
print()
print("  constants.py (line 26, used by every script in this folder):")
print('    "# Scalar QED correction for spin-0 particle (standard, not a fit)"')
print('    alpha_pi = alpha/pi   "= 2x Schwinger correction"')
print()
print("  higgs_cell_energy.py (dedicated derivation script for this exact claim):")
print('    "This is standard QED for scalars -- not a fit."')
print('    "[= 2 x Schwinger]"')
print()
print("  higgs_linking_spin.py (dedicated derivation script for topology->spin->correction):")
print('    "Standard QED one-loop mass corrections [textbook result]"')
print()
print("  CONCLUSION: no script in this folder claims alpha/pi's own coefficient")
print("  is derived from torsion-medium equations of motion. What IS derived")
print("  (higgs_linking_spin.py) is only WHICH external formula applies, via")
print("  the linking-number pi-rotation-symmetry theorem (an exact, genuine")
print("  topological result) -- a selection rule between two imports, not a")
print("  derivation of either import's own origin.")

# =============================================================================
print()
print(SEP2)
print("PART B: DOES A NATIVE ALTERNATIVE (delta_n) REPRODUCE m_H?")
print(SEP2)
print()

m_H_pdg22 = 125.25   # GeV, PDG 2022
m_H_unc = 0.17        # GeV, PDG 2022 1-sigma
n_exact = 2.01868734358  # from alpha derivation, all gaps closed
delta_n = n_exact - 2    # = 0.01869, the vertex-stiffness correction

print(f"  E_cell (bare)         = {E_cell_GeV:.6f} GeV")
print(f"  m_H measured (PDG22)  = {m_H_pdg22:.3f} +/- {m_H_unc:.2f} GeV")
print()

# Current (imported) formula, for reference
m_H_alpha_pi = E_cell_GeV * (1 + alpha / pi)
sigma_alpha_pi = abs(m_H_alpha_pi - m_H_pdg22) / m_H_unc
print(f"  CURRENT: E_cell*(1+alpha/pi)        = {m_H_alpha_pi:.4f} GeV"
      f"  ({sigma_alpha_pi:.2f} sigma)")

# Candidate 1: bulk-regime delta_n/pi, lifted directly from higgs_me_correction.py
m_H_dn_bulk = E_cell_GeV * (1 + delta_n / pi)
sigma_dn_bulk = abs(m_H_dn_bulk - m_H_pdg22) / m_H_unc
print(f"  CANDIDATE 1 (bulk delta_n/pi, electron's own formula):")
print(f"    E_cell*(1+delta_n/pi)              = {m_H_dn_bulk:.4f} GeV"
      f"  ({sigma_dn_bulk:.2f} sigma)")

# Candidate 2: bulk-regime delta_n directly (no /pi)
m_H_dn_raw = E_cell_GeV * (1 + delta_n)
sigma_dn_raw = abs(m_H_dn_raw - m_H_pdg22) / m_H_unc
print(f"  CANDIDATE 2 (bulk delta_n, no /pi):")
print(f"    E_cell*(1+delta_n)                 = {m_H_dn_raw:.4f} GeV"
      f"  ({sigma_dn_raw:.2f} sigma)")

# Candidate 3: sub-cell l=0-only channel, mirroring higgs_h1_residual.py's own
# "LEAD 1" structure (f_eff = PHI alone for sub-cell particles, not the full
# bulk L3(phi,log5) blend) -- delta_k_l0 = delta_n/L3 is the underlying
# per-channel stiffness; for a sub-cell particle only the l=0 channel (weight
# PHI) should contribute, so the analogous correction is PHI*delta_k_l0.
delta_k = delta_n / L3
m_H_l0 = E_cell_GeV * (1 + phi * delta_k)
sigma_l0 = abs(m_H_l0 - m_H_pdg22) / m_H_unc
print(f"  CANDIDATE 3 (sub-cell l=0-only, PHI*delta_k, delta_k=delta_n/L3):")
print(f"    E_cell*(1+PHI*delta_k)             = {m_H_l0:.4f} GeV"
      f"  ({sigma_l0:.2f} sigma)")
print()

check("A1 delta_n/pi (bulk, electron's own formula) does NOT match m_H as well as alpha/pi",
      sigma_dn_bulk > sigma_alpha_pi,
      f"{sigma_dn_bulk:.2f} sigma vs {sigma_alpha_pi:.2f} sigma")
check("A2 raw delta_n (no /pi) does NOT match m_H as well as alpha/pi",
      sigma_dn_raw > sigma_alpha_pi,
      f"{sigma_dn_raw:.2f} sigma vs {sigma_alpha_pi:.2f} sigma")
check("A3 sub-cell l=0-only candidate does NOT match m_H as well as alpha/pi",
      sigma_l0 > sigma_alpha_pi,
      f"{sigma_l0:.2f} sigma vs {sigma_alpha_pi:.2f} sigma")

print()
print(SEP)
print("SUMMARY")
print(SEP)
n_pass_partial = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass_partial}/{len(results)} checks so far PASS")
print()
print("  None of the three native-mechanism candidates tested reproduce m_H")
print("  as precisely as the imported alpha/pi correction. This is a genuine")
print("  negative result, not exhaustive proof no native mechanism exists --")
print("  but no swap-in replacement was found among the electron-mass-derived")
print("  delta_n family.")

# =============================================================================
print()
print(SEP2)
print("PART C: THE FRAMEWORK'S OWN 'NO 1/pi FOR CONTACT CORRECTIONS' RULE")
print(SEP2)
print()
print("  Repo-wide search (not folder-scoped) for the framework's OWN native")
print("  one-loop EM self-energy mechanism, used to derive alpha ITSELF:")
print()
print("  analysis/alpha/alpha_born_vertex.py, Part B (k_n_self = alpha*k_n):")
print('    "no 1/pi factor because the stiffness is a CONTACT correction')
print('    (coordinate-space, not momentum-space loop), analogous to how')
print('    the Higgs vev correction alpha^2*phi^2 lacks 1/pi (Section 5a,')
print('    doc_higgs)."')
print("  analysis/demos/alpha_doc.py, V24:")
print('    "No 1/pi: this is coordinate-space contact (like n*alpha^2 in')
print('    the quadratic, not Schwinger)."')
print()
print("  This is the framework's OWN established, validated rule: a torsion-")
print("  medium vertex mode's own EM self-energy is a CONTACT correction (one")
print("  factor of alpha, NO 1/pi) -- explicitly contrasted against Schwinger/")
print("  QED-loop corrections (which DO have 1/pi) as a DIFFERENT thing. This")
print("  is not a candidate I constructed -- it is the same rule the")
print("  framework already uses to justify the HIGGS VEV's own alpha^2*phi^2")
print("  term (Section 5a) having no 1/pi. It was never applied to the")
print("  LEADING m_H correction, creating an internal inconsistency: the")
print("  vev's next-order term follows the native no-1/pi contact rule, but")
print("  the mass's leading term uses the imported has-1/pi QED rule, for")
print("  the SAME A_g mode, with no stated reason for the difference.")
print()

m_H_bare_alpha = E_cell_GeV * (1 + alpha)
sigma_bare_alpha = abs(m_H_bare_alpha - m_H_pdg22) / m_H_unc
print(f"  CANDIDATE 4 (bare alpha, no 1/pi, matching the framework's own")
print(f"  contact-correction rule):")
print(f"    E_cell*(1+alpha)                   = {m_H_bare_alpha:.4f} GeV"
      f"  ({sigma_bare_alpha:.2f} sigma)")
print()

check("A4 bare alpha (native contact-correction rule) does NOT match m_H as well as alpha/pi",
      sigma_bare_alpha > sigma_alpha_pi,
      f"{sigma_bare_alpha:.2f} sigma vs {sigma_alpha_pi:.2f} sigma")

print()
print(SEP)
print("FINAL SUMMARY")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass}/{len(results)} checks PASS")
print()
print("  All FOUR native-mechanism candidates tested (delta_n/pi, raw delta_n,")
print("  sub-cell l=0-only, and now bare alpha matching the framework's own")
print("  no-1/pi contact-correction rule) reproduce m_H WORSE than the imported")
print("  alpha/pi formula. This strengthens, not weakens, the finding: it is")
print("  not merely that alpha/pi happens to be unexplained -- the framework's")
print("  OWN closest native analog (bare alpha, no 1/pi) gives a numerically")
print("  worse match, so there is no currently-available native mechanism,")
print("  applied consistently with its own precedents elsewhere in this repo,")
print("  that reproduces the measured Higgs mass as well as the import does.")

# =============================================================================
print()
print(SEP2)
print("PART D: DOES THE HOPF FIBRATION'S OWN S^3 GEOMETRY EXPLAIN THE 1/pi?")
print(SEP2)
print()
print("  Author's question: the (1,2) winding lives on the Hopf fibration of")
print("  S^3 -- a genuinely 4-dimensional (S^3 subset R^4) torsionverse-native")
print("  space, already used for a REAL, independently-verified integral:")
print("  CS_(p,q) = p*q*Vol(S^3) = p*q*2*pi^2 (analysis/alpha/gap3_chern_simons.py,")
print("  alpha_chern_weil_general.py), which builds Q = CS/phi = 4*pi^2/phi,")
print("  load-bearing in the alpha quadratic itself. Does THIS same S^3")
print("  structure -- not a generic mainstream QED momentum-space loop --")
print("  explain the Higgs mass correction's 1/pi?")
print()
print("  KEY DISTINCTION FIRST: finite icosahedral group-character sums (CG")
print("  decompositions, Born vertex corrections) generically give ALGEBRAIC")
print("  numbers built from phi (since chi(*, C_5) involves cos(72deg)=phi/2) --")
print("  they do NOT generate transcendental pi factors. pi only enters this")
print("  framework through genuinely CONTINUOUS geometric integrals (solid")
print("  angles, Vol(S^3)/Vol(S^2)). This is WHY Q (a real S^3 integral) has")
print("  pi^2, while alpha's own vertex correction (a finite Born/CG sum,")
print("  Section 4.5) does NOT -- consistent with, not contradicting, the")
print("  'contact correction has no 1/pi' rule tested in Part C.")
print()

CS_12 = 1 * 2 * 2 * pi**2   # p*q*Vol(S^3), (1,2) winding
Vol_S2 = 4 * pi
ratio_S2_CS = Vol_S2 / CS_12
print(f"  Vol(S^2) = 4*pi = {Vol_S2:.8f}")
print(f"  CS_(1,2) = p*q*Vol(S^3) = 2*2*pi^2 = {CS_12:.8f}")
print(f"  Vol(S^2)/CS_(1,2) = {ratio_S2_CS:.10f}   (compare 1/pi = {1/pi:.10f})")
print()

check("D1 Vol(S^2)/CS_(1,2) equals 1/pi to float precision",
      abs(ratio_S2_CS - 1 / pi) < 1e-12,
      f"{ratio_S2_CS:.12f} vs {1/pi:.12f}")

print()
print("  D1 PASSES -- but this is EXACTLY the Wyler-formula trap doc_alpha.txt")
print("  Section 7.2 itself already warns about: 'derived by identifying the")
print("  ratio of certain symmetric space volumes without a physical")
print("  mechanism.' I searched deliberately for a combination of already-")
print("  established torsionverse quantities (Vol(S^2) is generic, CS_(1,2)")
print("  is real but built for Q, a DIFFERENT purpose) that hits 1/pi -- and")
print("  found one, algebraically exactly, for (1,2) specifically. That is")
print("  the DEFINITION of not independently motivated: no physical argument")
print("  says a self-energy correction SHOULD be Vol(S^2)/CS_(p,q) before")
print("  checking that it numerically works. A genuine derivation would need")
print("  its own independent construction (e.g. a real Chern-Simons-style")
print("  integral over the Hopf torus computing self-energy specifically,")
print("  the way gap3_chern_simons.py computes Q specifically) -- NOT a")
print("  reverse-engineered ratio of two numbers that happen to already exist.")
print("  CONCLUSION: the S^3 connection is a real, interesting, UNTESTED lead")
print("  (confirmed via git log -S and notes/ search: never attempted anywhere")
print("  in this repo's history) -- but this script does not close it. A real")
print("  attempt needs an independently-motivated integral, not a fitted ratio.")


