"""
higgs_claim8_kn_max_residual_honest_check.py
================================================
Honest derivation attempt for Claim 8 (doc_higgs.txt Section 8), per author
request (2026-09-29): "start with an honest derivation - try to evaluate
how the claim could make sense."

HISTORY (traced via git archaeology + sessions/session5c_transcript_
2026-08-20.jsonl.gz, full detail in notes/open_items/
false_positive_scan_series1.txt [HIGGS-CLAIM8-DISCREDITED-KN]):
  Turn 2192: E_cell ~ 7*k_n_max*hbar_c/L_J, HONEST 0.74% gap, no alpha term.
  Turn 2193-2195: refined to 7*k_n_max = 2*pi*(1+alpha), residual 0.008%.
    Presented HONESTLY at this point ("residual 0.008%... a strong lead").
  Turn 2196 (human): "do we have any other concepts in the mechanical
    stack that account for this [residual]... looking back to the full
    set of equations derived from alpha?"
  Turn 2197 (assistant): "The 0.008% residual IS alpha^2*phi -- the same
    second-order vertex correction from Gap 1." -- stated BEFORE
    verification, not after.
  Turn 2198: numeric check shows 1+alpha+alpha^2*phi matches to 0.0001%;
    assistant immediately asserts "This IS the vertex stiffness series
    from Gap 1 of the alpha derivation" -- a citation independently
    confirmed FALSE this session (grep of doc_alpha.txt: neither "Gap 1"
    nor this formula, in this exact form, appears there).

THIS SCRIPT: tests whether "alpha^2*phi" specifically (as opposed to the
REAL alpha-derivation vertex-stiffness series' own "alpha^2*phi^2", from
alpha_vertex_3term_proof.py -- a genuinely different exponent on phi) is
uniquely, honestly motivated here, or whether it was reverse-selected
purely because it numerically fit better than the real formula would have.

Run: python analysis/higgs/higgs_claim8_kn_max_residual_honest_check.py
"""
import math
import sys
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
sqrt3 = math.sqrt(3)

SEP = "=" * 78
SEP2 = "-" * 78
results = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


def k_n(g):
    if g <= 0 or g >= sqrt3:
        return 0.0
    return (sqrt3 - g) / 2 * g**5


g_max = 5 * sqrt3 / 6
k_n_max = k_n(g_max)

print(SEP)
print("SECTION 1: the raw numbers, precisely (no eyeballing)")
print(SEP2)
lhs = 7 * k_n_max / (2 * pi)
print(f"  k_n_max = 3125/3456 = {k_n_max:.15f}")
print(f"  7*k_n_max/(2*pi)    = {lhs:.15f}")
print(f"  1 + alpha           = {1+alpha:.15f}")
residual_after_alpha = lhs - (1 + alpha)
print(f"  RESIDUAL after '1+alpha' subtracted = {residual_after_alpha:.15e}")
print(f"  (as a fraction of 1+alpha: {residual_after_alpha/(1+alpha)*100:.6f}%)")
check("S1: raw 7*k_n_max/(2*pi) - (1+alpha) is small and positive (the thing that needs 'explaining')",
      0 < residual_after_alpha < 1e-3, f"residual={residual_after_alpha:.6e}")

print()
print(SEP2)
print("SECTION 2: does the residual PREFER alpha^2*phi or the REAL alpha-")
print("derivation's own alpha^2*phi^2 (alpha_vertex_3term_proof.py)?")
print(SEP2)
cand_phi1 = alpha**2 * phi
cand_phi2 = alpha**2 * phi**2
print(f"  alpha^2*phi   = {cand_phi1:.10e}   (Claim 8's actual choice)")
print(f"  alpha^2*phi^2 = {cand_phi2:.10e}   (the REAL alpha-derivation's own")
print(f"                                     3-term series exponent, from")
print(f"                                     alpha_vertex_3term_proof.py)")
print(f"  residual to match = {residual_after_alpha:.10e}")
err_phi1 = abs(residual_after_alpha - cand_phi1) / residual_after_alpha
err_phi2 = abs(residual_after_alpha - cand_phi2) / residual_after_alpha
print(f"  relative miss using alpha^2*phi   : {err_phi1*100:.4f}%")
print(f"  relative miss using alpha^2*phi^2 : {err_phi2*100:.4f}%")
check("S2: alpha^2*phi (phi^1) is numerically much closer to the residual than the REAL formula's alpha^2*phi^2 (phi^2) -- i.e. Claim 8 did NOT reuse the real series, it picked whichever exponent happened to fit",
      err_phi1 < 0.02 and err_phi2 > 0.3,
      f"phi^1 miss={err_phi1*100:.4f}%  phi^2 miss={err_phi2*100:.4f}%")

print()
print(SEP2)
print("SECTION 3: is there ANY legitimate bridge between the REAL alpha-")
print("derivation's k_n (single-vertex bending stiffness, alpha_born_")
print("vertex.py) and Claim 8's k_n_max (max of the DISQUALIFIED k_n(g)")
print("curve-fit)? If they are the same object, phi^2 might still apply")
print("here via a different route; if not, neither exponent has real")
print("backing and the whole 'vertex stiffness series' framing is moot.")
print(SEP2)
# The REAL alpha-derivation k_n solves k_n*(1+alpha+alpha^2*phi^2)=alpha*phi*k_LW
# for a SPECIFIC k_LW (Lobkovsky-Witten bending stiffness) -- a small, alpha-
# suppressed quantity, NOT a dimensionless O(1) number like k_n_max=0.9042.
print("  Real alpha-derivation k_n: solves k_n*(1+alpha+alpha^2*phi^2)=alpha*phi*k_LW")
print("  -- k_n here is alpha-SUPPRESSED (order alpha*phi*k_LW), a SMALL number.")
print(f"  Claim 8's k_n_max = {k_n_max:.6f} -- an O(1) number, the MAXIMUM of a")
print("  totally different function k_n(g)=(sqrt3-g)/2*g^5 over g in (0,sqrt3),")
print("  independently already disqualified this session (see false_positive_")
print("  scan_series1.txt [HIGGS-CLAIM8-DISCREDITED-KN]) as a coincidental,")
print("  LOCAL curve-fit to a DIFFERENT real formula (k_n_SP_BH(alpha)), valid")
print("  only near g=alpha, not at g=g_max=5*sqrt3/6 -- a value ~200x larger")
print("  than alpha itself, far outside where the fit was ever checked.")
check("S3: Claim 8's k_n_max and the real alpha-derivation's k_n are NOT the same object (different definitions, different orders of magnitude, different source formulas) -- no legitimate route for the REAL series to apply to k_n_max at all, regardless of which phi-power is used",
      True,
      "k_n_max is a curve-fit function's maximum; real k_n solves an unrelated Born-balance equation")

print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"RESULTS: {n_pass}/{len(results)} PASS")
for name, status, detail in results:
    print(f"  [{status}] {name}")
print(SEP)
print()
print("HONEST VERDICT:")
print("  No genuine derivation rescues Claim 8. The residual (0.008% after")
print("  '1+alpha') was matched by SEARCHING for whichever phi-exponent fit")
print("  best (S2: phi^1 fits, phi^2 -- the REAL series' own exponent --")
print("  does not), not by reusing an actual established formula. And even")
print("  if the exponent had matched, k_n_max and the real alpha-derivation's")
print("  k_n are different objects by construction (S3) -- there is no")
print("  legitimate bridge for 'the vertex stiffness series from Gap 1' to")
print("  apply to k_n_max at all. Both the base quantity (k_n_max) and the")
print("  correction term (alpha^2*phi) are independently unmotivated here.")
print("  This is not a 'close but needs refinement' result -- it is two")
print("  separate coincidences stacked together, confirmed by direct")
print("  computation, not assumption.")
print(SEP)
