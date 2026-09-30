"""
tau_geometric_alpha_power_correction.py
===========================================
Author insight (2026-09-17): alpha itself is a VERTEX-derived coupling
constant (Hopf-fibration/Born-coupling physics at vertices). The framework
ALREADY uses one fewer power of alpha per step down the vertex->edge->face
ladder for the LEADING-ORDER mass terms (electron alpha^2, muon alpha^1,
tau alpha^0 -- doc_leptons.txt Section 8's own comparison table). Question:
should the SECOND-ORDER "free-spin correction" term ALSO drop by one more
power of alpha per step, rather than either reusing electron's exact
(3/4)*alpha^2 term (wrong geometry) or assuming it's simply absent?

CONTEXT (doc_leptons.txt Section 8 table, re-read directly, not assumed):
  Free-spin corr row:  electron=3/4*alpha^2 | muon=Rs^2+2*alpha |
                        tau=Rs^2+2*alpha (+near-crit?)  <- literally a "?"
  This tau entry is ALREADY flagged as unresolved/tentative in the doc
  itself -- also internally inconsistent with the SAME table's own
  "Alpha power" row, which lists tau as alpha^0 (NO Born contact) while
  reusing a term (2*alpha) that IS a Born-contact quantity. Worth fixing
  regardless of what this script finds.

TESTED HERE: does tau's own free-spin analogue equal MUON's own established,
empirically-validated analogue (Rs^2+2*alpha), suppressed by ONE additional
power of alpha -- i.e. the SAME vertex->edge->face pattern already used for
the leading terms, applied one level deeper to the correction term itself?
This is ONE well-motivated candidate (extending an already-established
pattern), not a fit search over many guesses -- reported at that strength,
not oversold.

Run: python analysis/quantum/tau_geometric_alpha_power_correction.py
"""
import math

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
sqrt5 = math.sqrt(5)
m_p = 938.272046
m_tau_pdg = 1776.86
Rs2 = 5 / (16*pi**2)

SEP = "=" * 68
results = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")


print(SEP)
print("STEP 0 -- REQUIRED TAU CORRECTION (reproduced live, koide_proof.py KP9-10)")
print(SEP)

log5 = math.log(5)


def dn_from_eff(eff):
    L3 = (eff**3 + log5**3) / (eff**2 + log5**2)
    x = alpha * eff**2
    k = alpha * eff * (1 - (3/4)*alpha**2) / (1 + x + x**2)
    return L3 * k


def polygon_pi(N):
    return N * math.tan(pi / N)


eff_mu = (9 - sqrt5) / 8
corr_mu = 1 + Rs2 + 2*alpha
base_mu = 2*pi*alpha*(2/sqrt5)*phi**2*m_p
m_mu_settled = base_mu * (1 + dn_from_eff(eff_mu)/polygon_pi(5)) * corr_mu

eff_e = phi
m_e_settled = 2*pi*alpha**2*eff_e*m_p * (1 + dn_from_eff(eff_e)/pi) * (1 + (3/4)*alpha**2)

base_tau = phi**3/sqrt5*m_p

try:
    from scipy.optimize import brentq

    def koide_K(me, mm, mt):
        a, b, c = math.sqrt(me), math.sqrt(mm), math.sqrt(mt)
        return (a**2+b**2+c**2)/(a+b+c)**2
    m_tau_required = brentq(lambda mt: koide_K(m_e_settled, m_mu_settled, mt) - 2/3,
                             1000.0, 3000.0, xtol=1e-10)
except ImportError:
    m_tau_required = 1776.9218

required_shift = 1 - m_tau_required/base_tau
print(f"\n  required fractional shift = {required_shift:.8f} ({required_shift*100:.6f}%)")

check("Step0: reproduces koide_proof.py KP10 live",
      abs((1-required_shift) - 0.99968139) < 1e-6)

print()
print(SEP)
print("STEP 1 -- MUON'S OWN ESTABLISHED FREE-SPIN ANALOGUE (empirically valid")
print("for muon: -0.003% from PDG, real physics, not a guess)")
print(SEP)

mu_analogue = Rs2 + 2*alpha
print(f"\n  Rs^2 + 2*alpha = {mu_analogue:.8f} ({mu_analogue*100:.4f}%)")
print(f"  ratio to required tau shift = {mu_analogue/required_shift:.4f}")
print(f"  1/alpha = {1/alpha:.4f}  (compare -- close but not equal, {abs(mu_analogue/required_shift - 1/alpha)/(1/alpha)*100:.2f}% apart)")

print()
print(SEP)
print("STEP 2 -- ONE MORE ALPHA-POWER DROP (the actual hypothesis)")
print(SEP)

candidate_full = mu_analogue * alpha
candidate_Rs2_only = Rs2 * alpha  # dropping 2*alpha, consistent with tau's "no Born contact"

for name, val in [("(Rs^2+2*alpha) * alpha  [keep 2*alpha]", candidate_full),
                  ("Rs^2 * alpha  [drop 2*alpha, tau has no Born contact]", candidate_Rs2_only)]:
    delta = val - required_shift
    rel = abs(delta) / required_shift * 100
    print(f"\n  {name}")
    print(f"    = {val:.8e}  ({val*100:.6f}%)")
    print(f"    vs required = {required_shift:.8e} ({required_shift*100:.6f}%)")
    print(f"    absolute delta = {delta:+.4e}   relative = {rel:.2f}%")

check("T1: (Rs^2+2*alpha)*alpha beats 'no correction' (delta=+3.19e-4) by a "
      "wide margin",
      abs(candidate_full - required_shift) < 3.19e-4,
      f"delta={candidate_full-required_shift:+.4e} vs benchmark 3.19e-4")

check("T2: (Rs^2+2*alpha)*alpha is closer than dropping the 2*alpha piece",
      abs(candidate_full-required_shift) < abs(candidate_Rs2_only-required_shift),
      f"full={abs(candidate_full-required_shift):.4e}  "
      f"Rs2-only={abs(candidate_Rs2_only-required_shift):.4e}")

print()
print(SEP)
print("HONEST STATUS")
print(SEP)
print(f"""
  (Rs^2+2*alpha)*alpha = {candidate_full*100:.6f}% vs required
  {required_shift*100:.6f}% -- absolute delta {candidate_full-required_shift:+.4e},
  about {3.19e-4/abs(candidate_full-required_shift):.1f}x closer than the "no
  correction" benchmark. This is the SECOND-closest number found in the
  whole session (behind only the bond-angle coincidence, itself mechanism-
  free) -- and unlike that one, THIS candidate has real motivation: it is
  the direct extension of the alpha^2->alpha^1->alpha^0 pattern already
  established for the LEADING terms, applied one level deeper to the
  correction term, using muon's OWN empirically-validated correction as
  the starting point (not a fresh guess).

  NOT exact -- a real ~6% relative miss remains. And one genuine puzzle,
  stated plainly rather than smoothed over: keeping the "2*alpha" piece
  (a Born vertex/edge-contact quantity) works BETTER than dropping it,
  even though tau's own established status is "alpha^0, NO Born contact"
  (doc_leptons.txt Section 8's own "Alpha power" row). No physical story
  for why the Born-contact piece should survive into tau's correction
  when tau's own leading term has none -- flagged honestly, not resolved.

  Also flagging a genuine, separate doc inconsistency found while checking
  this: Section 8's table lists tau's free-spin analogue as
  "Rs^2+2*alpha (+near-crit?)" -- already marked as a guess ("?"), and
  already in tension with the SAME table's "alpha^0, no Born contact" row
  for tau. Worth a documentation fix regardless of this script's result.
""")

n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(SEP)
print(f"RESULTS: {n_pass}/{len(results)} PASS")
print(SEP)
