"""
higgs_ag_channel_born_vertex_check.py
========================================
Multi-model subagent review of item 6 (2026-09-23): Grok 4.7 re-derived the
already-known "2-vertex"/Z_2 algebra (delta_m = 2*alpha*(hbar_c/L_J), pi
cancels via E_cell's own 2*pi) and correctly flagged the SAME still-open
sub-question already sitting in judgment_calls.txt: does each Z_2-identified
phase's contact contribute exactly BARE alpha, or alpha*phi (phi = chi(T_1g,
C5), the T_1g-specific vertex enhancement from alpha_born_vertex.py Part A)?

THE ALREADY-TRIED, ALREADY-FAILED ATTEMPT (higgs_c2_vertex_born_amplitude_
test.py, cited in judgment_calls.txt): reused alpha_born_vertex.py Part A's
T_1g-specific formula at the C2 class, computing chi(T_1g,C5)+|chi(T_1g,C2)|
= phi+1 = phi^2 = 2.618 -- a genuine 31% miss vs the needed value of 2. This
used T_1g's OWN characters throughout, for BOTH classes.

THE GAP THIS SCRIPT CHECKS (NOT the same calculation, a materially different
one): the Higgs/scalar mode is A_g, not T_1g -- and A_g is, by definition,
the TRIVIAL 1D irrep (doc_higgs.txt Section 3a: "the unique totally-
symmetric mode invariant under all I_h operations"), meaning chi(A_g, g) = 1
for EVERY group element g, not just C5. The earlier C2 test never actually
plugged A_g's OWN (necessarily trivial) character into the SAME Born-
projection formula structure -- it stayed within T_1g's own character
throughout. This script does that specific, not-yet-tried substitution and
reports the result honestly, including an explicit epistemic caveat about
whether it is a genuine new physical fact or a near-tautological consequence
of projecting onto a trivial representation.

Run: python analysis/higgs/higgs_ag_channel_born_vertex_check.py
"""

import math

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3

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
print("A_g-CHANNEL BORN VERTEX CHECK (2026-09-23)")
print("Grok 4.7 subagent review: does the SCALAR (A_g) channel's own vertex")
print("coupling carry the T_1g-specific phi factor, or bare alpha?")
print(SEP)

# =============================================================================
print()
print(SEP2)
print("STEP 1: T_1g's OWN character at C5 (reused verbatim from")
print("        alpha_born_vertex.py Part A -- NOT recomputed differently)")
print(SEP2)
chi_T1g_C5 = 1 + 2 * math.cos(2 * pi / 5)
print(f"  chi(T_1g, C5) = 1 + 2*cos(72 deg) = {chi_T1g_C5:.10f}")
check("AG1: chi(T_1g,C5) = phi (reproduces alpha_born_vertex.py Part A exactly)",
      abs(chi_T1g_C5 - phi) < 1e-10, f"chi={chi_T1g_C5:.10f}, phi={phi:.10f}")

# =============================================================================
print()
print(SEP2)
print("STEP 2: A_g's OWN character -- by DEFINITION (trivial 1D irrep),")
print("        not by any matrix-trace computation")
print(SEP2)
print("  A_g is the totally-symmetric, trivial 1D representation of I_h --")
print("  by definition, chi(A_g, g) = 1 for EVERY group element g, including")
print("  BOTH C5 and C2 (this is not a coincidence needing verification for a")
print("  specific class -- it is what 'trivial representation' MEANS).")
chi_Ag_C5 = 1.0
chi_Ag_C2 = 1.0
check("AG2: chi(A_g, C5) = 1 and chi(A_g, C2) = 1 (definitional, trivial rep)",
      chi_Ag_C5 == 1.0 and chi_Ag_C2 == 1.0, "by definition of the trivial irrep")

# =============================================================================
print()
print(SEP2)
print("STEP 3: apply the SAME Born-projection formula structure")
print("        (M(Gamma <- C5 vertex) = chi(Gamma,C5) * alpha * k_LW) to A_g")
print("        instead of T_1g -- the substitution the earlier C2 test never made")
print(SEP2)
M_T1g_C5 = chi_T1g_C5  # * alpha * k_LW, k_LW factored out (ratio quantity)
M_Ag_C5 = chi_Ag_C5    # * alpha * k_LW
print(f"  M(T_1g <- C5 vertex) / (alpha*k_LW) = chi(T_1g,C5) = {M_T1g_C5:.6f}  [the W/Z channel, alpha_born_vertex.py's actual subject]")
print(f"  M(A_g  <- C5 vertex) / (alpha*k_LW) = chi(A_g, C5) = {M_Ag_C5:.6f}  [the Higgs/scalar channel -- NOT what alpha_born_vertex.py actually computed]")
check("AG3: M(A_g<-vertex) carries NO phi enhancement -- exactly bare alpha*k_LW, unlike T_1g's alpha*phi*k_LW",
      abs(M_Ag_C5 - 1.0) < 1e-12, f"M_Ag/M_T1g = {M_Ag_C5/M_T1g_C5:.6f} (1/phi = {1/phi:.6f})")

# =============================================================================
print()
print(SEP2)
print("STEP 4: does this differ from (i.e. is it a DIFFERENT calculation than)")
print("        the already-failed C2 test's chi(T_1g,C5)+|chi(T_1g,C2)| = phi^2?")
print(SEP2)
chi_T1g_C2 = 1 + 2 * math.cos(pi)  # C2 = pi-rotation
prior_c2_test_value = chi_T1g_C5 + abs(chi_T1g_C2)
print(f"  Prior (already-failed) C2 test: chi(T_1g,C5) + |chi(T_1g,C2)| = {prior_c2_test_value:.6f} (= phi^2 = {phi**2:.6f}, needed 2, 31% miss)")
print(f"  This check: chi(A_g,C5) [and chi(A_g,C2), both = 1] -- a DIFFERENT")
print(f"  irrep's own character used throughout, not a 2-class sum within T_1g.")
check("AG4: this check uses a genuinely DIFFERENT quantity (A_g's own character) than the prior C2 test (T_1g's own 2-class sum) -- not a restatement of already-failed work",
      abs(prior_c2_test_value - M_Ag_C5) > 0.5,
      f"prior={prior_c2_test_value:.4f} vs this check={M_Ag_C5:.4f}")

# =============================================================================
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"RESULTS: {n_pass}/{len(results)} PASS")
for name, status, detail in results:
    print(f"  [{status}] {name}")
print(SEP)
print()
print("HONEST VERDICT: using A_g's OWN (necessarily trivial) character, the SAME")
print("Born-projection formula structure that gives T_1g's vertex coupling a phi")
print("enhancement gives the A_g/scalar channel NO such enhancement -- exactly")
print("bare alpha. This DOES directly answer Grok's flagged open question the")
print("way needed for the 'each Z_2 phase contributes bare alpha' picture to work,")
print("and it is NOT a restatement of the already-failed C2 test (AG4) -- that")
print("test stayed within T_1g's own characters throughout; this one uses A_g's.")
print()
print("IMPORTANT EPISTEMIC CAVEAT (do not oversell this): chi(A_g, g)=1 for ALL g")
print("is true by DEFINITION of a trivial representation, for ANY group and ANY")
print("class, regardless of anything specific to icosahedral symmetry or this")
print("framework's vertex-contact mechanism. This means 'projecting any vertex-")
print("localized perturbation onto A_g gives back its own bare strength with no")
print("enhancement' is close to a GENERIC, near-tautological fact about trivial")
print("representations, not a substantive NEW physical derivation specific to")
print("Jobson cells. Its real value here is narrower and more honest: it explains")
print("WHY the earlier C2 test's phi^2 miss was likely the wrong quantity to have")
print("computed for THIS question (it used T_1g characters for a question about")
print("the A_g channel) -- not that a deep new mechanism has been found. The")
print("bare-alpha-per-phase assumption in the Z_2/2-vertex picture is consistent")
print("with, but not independently PROVEN by, this check.")
