"""
bound_neutron_exact_distance_solve.py
======================================
Author-directed re-investigation (2026-10-01) of neutron_g_factor.py's bound-
neutron g-factor claim, correcting an earlier same-day framing that wrongly
implied "the neutron sitting at/inside r_p" was an unjustified or implausible
assumption. Author: "first you should research all the code and transcripts
deeply for any suggestion of why the original distance was used. then you
should be reporting the exact distance that replicates the magnetic moment
exactly and THEN you might report the difference in effect at the other
distances." Then, separately: "can you review the falloff code just to
double check that formula is well founded."

WHY THE ORIGINAL "FLAT REUSE AT r_p" HAPPENED (git + transcript archaeology):
  Commit 67e42e1 (2026-08-21) created neutron_g_factor.py. Its own comment
  calls the bound-neutron number a "Rough estimate" and tags it
  "[OPEN: quantitative treatment of proton Zone 3 pressure on bound neutron]"
  -- i.e. the ORIGINAL AUTHOR ALREADY FLAGGED this as not rigorously derived.
  The session10 transcript (2026-08-21T21:50Z, same day) shows the author's
  own working description: "bound neutron g_n differs from free neutron g_n
  by approximately the proton Zone 3 contribution AT THE PROTON-NEUTRON
  SEPARATION DISTANCE" -- i.e. distance-dependence was explicitly anticipated
  from the start. The code simply never implemented that dependence; it used
  a flat full-strength subtraction (mu_Zone3_proton's complete value) as a
  placeholder. No specific distance was ever modeled, chosen, or claimed --
  r_p only appears in the surrounding PROSE as "nuclear separation" context
  (r_grind to r_0), never as an actual input to the formula.

WHICH FALLOFF EXPONENT IS ACTUALLY WELL-FOUNDED (the "double check" requested):
  Two DIFFERENT physical quantities are in play, easy to conflate:
  (a) doc_nucleus.txt Section 1.1: "v(r) = Rs*c*(r_p/r)^2" is the INTERNAL
      co-rotation LINEAR velocity of a Zone 3 cell at radius r (used to
      build mu_Zone3_proton itself, by integrating over lambda_p<r<r_p).
      This IS consistent with real GR Lense-Thirring: the ANGULAR precession
      rate Omega_LT(r) falls as 1/r^3 (standard GR result); a co-rotating
      point's LINEAR/tangential velocity v(r)=Omega_LT(r)*r then falls as
      1/r^2, matching this ansatz exactly. Not a contradiction -- confirmed
      well-founded for what it actually describes (an internal velocity
      profile), once angular vs. linear velocity are not conflated (an
      earlier same-day edit to entanglement_phase_lock.py's comment wrongly
      called this ansatz inconsistent with GR; corrected back).
  (b) analysis/quantum/entanglement_phase_lock.py's actual, tested, USED
      formula for the Zone 3 COUPLING ENERGY between two separate particles
      is E_Z3(r) = E_0*(r_p/r)^3 -- built from a DIFFERENT physical analogy
      ("the Hopf winding IS a topological magnetic dipole -> field ~1/r^3"),
      matching standard E&M dipole-FIELD falloff, not derived by reusing v(r)
      directly (energy != velocity; an extra power of 1/r separates them).

  THE QUESTION THIS SCRIPT ANSWERS IS CASE (b)'s TYPE, NOT CASE (a)'s: "how
  does the proton's OWN Zone-3-generating structure influence a SEPARATE,
  EXTERNAL neutron at distance r?" -- an inter-particle coupling question,
  physically analogous to a dipole FIELD felt at a distance, not the
  internal velocity profile within the proton's own shell. The dipole-field
  exponent (n=3, matching entanglement_phase_lock.py's own tested formula)
  is therefore the better-justified choice for mu_induced(r); n=2 (the
  internal-velocity law, reused naively in an earlier version of this
  script) describes a different quantity and is shown below only for
  comparison/transparency, not as equally well-justified for THIS question.

THIS SCRIPT: builds the r-dependent model using BOTH exponents for
transparency, leading with n=3 (dipole field, better-justified for this
external-coupling question):

  mu_induced(r) = mu_Zone3_proton * (r_p/r)^n       [r = proton-neutron
  g_n(bound; r) = g_n(free) - mu_induced(r)          center separation]

Both reduce EXACTLY to the original flat-reuse number at r = r_p (the
"onset" radius) regardless of n -- so the old calculation is not a separate,
competing assumption, it is simply the r=r_p special case of this family of
laws, for any n.

APPROACH (per author's instruction): solve for the EXACT r that reproduces
the measured g_n = -1.9130 mu_N precisely (primary result, both n=2 and
n=3), THEN separately report g_n(bound; r) at other already-established,
physically-motivated separations for comparison (r_p itself, r_grind, the
pion force range r_0, and the real external Reid (1968) NN potential's
attractive-well minimum).

CAVEAT (author-flagged, not resolved here): real experiments (EMC effect /
"proton spin puzzle", e.g. polarized deep-inelastic scattering measuring the
g1 structure function) DO constrain proton spin structure -- but that is a
QCD quark+gluon+orbital spin decomposition, a different accounting than this
framework's geometric medium co-rotation picture. No established mapping
between the two exists in this repo; this script does NOT attempt one, and
the result below does not rely on or claim support from that measurement.

Run: python sandbox/nuclear/bound_neutron_exact_distance_solve.py
Reference: analysis/nuclear/neutron_g_factor.py (g_n_free, mu_Zone3_proton),
  docs/series1/doc_nucleus.txt Section 1.1 (the (r_p/r)^2 internal-velocity
  law), analysis/quantum/entanglement_phase_lock.py (the (r_p/r)^3
  dipole-field law, actually used/tested there for inter-particle coupling),
  analysis/nuclear/nn_potential_external_comparison.py (Reid 1968 potential)
"""

import sys, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

# ── Established inputs (unchanged from neutron_g_factor.py / nucleus_doc.py) ─
g_n_free        = -1.5666   # mu_N, neutron_g_factor.py SECTION 4 total
mu_Zone3_proton = 0.3602    # mu_N, proton's own Zone 3 self-contribution (Section 4)
g_n_meas        = -1.9130   # mu_N, PDG (the SAME target already used in the
                            # existing calculation -- not re-litigated here)
r_p_fm          = 0.8414    # fm, proton charge radius / Zone 3 outer edge
lambda_p_fm     = 0.2103    # fm, Zone 1/2 boundary
r_grind_fm      = 2 * lambda_p_fm   # 0.4206 fm, Zone 2-boundary contact distance
r_pion_fm       = 1.41445           # fm, pion force range (r_0, doc Section 5.4)

print(SEP)
print("bound_neutron_exact_distance_solve.py")
print("Exact separation reproducing g_n(bound) = PDG precisely, for BOTH")
print("candidate falloff exponents, plus comparison at other distances")
print(SEP)

# ── Section 1: the r-dependent model, two candidate exponents ───────────────
print()
print(SEP2)
print("SECTION 1: two candidate falloff exponents (see docstring for why)")
print(SEP2)

def mu_induced(r_fm, n):
    """Induced external Zone 3 contribution at proton-neutron center separation r."""
    return mu_Zone3_proton * (r_p_fm / r_fm) ** n

def g_n_bound(r_fm, n):
    return g_n_free - mu_induced(r_fm, n)

print("  n=2: doc_nucleus.txt Sec 1.1 internal co-rotation LINEAR velocity law")
print("       (consistent with real GR Lense-Thirring Omega~1/r^3 via v=Omega*r,")
print("       but describes an INTERNAL Zone 3 velocity profile, not an")
print("       external inter-particle coupling -- shown for comparison only)")
print("  n=3: entanglement_phase_lock.py's actual tested inter-particle")
print("       COUPLING ENERGY law (topological-dipole FIELD falloff) --")
print("       better physical match for 'proton's structure influencing a")
print("       separate, external neutron' -- used as PRIMARY below")
print()
for n in (2, 3):
    print(f"  At r = r_p exactly (n={n}): mu_induced = {mu_induced(r_p_fm, n):.4f} mu_N, "
          f"g_n(bound) = {g_n_bound(r_p_fm, n):.4f} mu_N")

check("BD1 Both exponents reproduce the ORIGINAL flat-reuse number exactly at "
      "r=r_p (confirms r_p is the r=r_p special case of either law, not a "
      "separate/competing assumption)",
      all(abs(g_n_bound(r_p_fm, n) - (g_n_free - mu_Zone3_proton)) < 1e-9 for n in (2, 3)),
      f"n=2: {g_n_bound(r_p_fm, 2):.4f}  n=3: {g_n_bound(r_p_fm, 3):.4f}  "
      f"vs flat reuse {g_n_free - mu_Zone3_proton:.4f}")

# ── Section 2: PRIMARY RESULT -- exact distance reproducing PDG precisely ────
print()
print(SEP2)
print("SECTION 2: PRIMARY RESULT -- exact r reproducing g_n = -1.9130 mu_N")
print(SEP2)

needed_induced = g_n_free - g_n_meas
ratio = needed_induced / mu_Zone3_proton

r_exact = {}
for n in (2, 3):
    r_exact[n] = r_p_fm / ratio ** (1.0 / n)

print(f"  Needed induced correction: g_n(free) - PDG = {g_n_free:.4f} - ({g_n_meas:.4f}) "
      f"= {needed_induced:.4f} mu_N")
print(f"  (r_p/r)^n = {needed_induced:.4f} / {mu_Zone3_proton:.4f} = {ratio:.6f}")
print()
for n in (2, 3):
    tag = "PRIMARY (dipole field)" if n == 3 else "secondary (internal velocity)"
    print(f"  n={n} [{tag}]:")
    print(f"    r_exact = r_p / ({ratio:.6f})^(1/{n}) = {r_exact[n]:.4f} fm")
    print(f"    r_exact / r_p = {r_exact[n]/r_p_fm:.4f}  "
          f"({100*(r_exact[n]/r_p_fm - 1):+.2f}% beyond r_p)")
    print(f"    Verification: g_n(bound; r_exact) = {g_n_bound(r_exact[n], n):.6f} mu_N  "
          f"(target {g_n_meas:.4f})")

check("BD2 Both exponents' r_exact solve g_n(bound; r) = PDG to high precision",
      all(abs(g_n_bound(r_exact[n], n) - g_n_meas) < 1e-6 for n in (2, 3)),
      f"n=2: {g_n_bound(r_exact[2], 2):.6f}  n=3: {g_n_bound(r_exact[3], 3):.6f}")
check("BD3 Both exponents' r_exact fall INSIDE their falloff law's own stated "
      "valid domain (r > r_p) -- neither is an extrapolation",
      all(r_exact[n] > r_p_fm for n in (2, 3)),
      f"n=2: r_exact={r_exact[2]:.4f} fm, n=3: r_exact={r_exact[3]:.4f} fm, "
      f"vs r_p={r_p_fm:.4f} fm")

print()
print(f"  HEADLINE (n=3, primary): exact-match separation = {r_exact[3]:.4f} fm, only")
print(f"  {100*(r_exact[3]/r_p_fm - 1):.1f}% beyond the proton's own r_p -- essentially simple")
print(f"  nucleon-nucleon CONTACT at the proton's own Zone 3 edge. CONCLUSION IS")
print(f"  ROBUST to the exponent choice: n=2 gives {r_exact[2]:.4f} fm ({100*(r_exact[2]/r_p_fm-1):.1f}%")
print(f"  beyond r_p) -- both land in the same 'ordinary near-contact' regime,")
print(f"  not a special or implausible assumption, and both are in-domain.")

# ── Section 3: SECONDARY -- effect at other established separations ─────────
print()
print(SEP2)
print("SECTION 3: g_n(bound; r) at other already-established separations (n=3)")
print(SEP2)
print("  (secondary context only, per author's request -- primary result is")
print("   Section 2 above. r < r_p evaluations are OUTSIDE either falloff")
print("   law's own stated domain [r > r_p] and are flagged as extrapolations.)")
print()

def report(label, r_fm, in_domain, n=3):
    g = g_n_bound(r_fm, n)
    err = 100 * (g - g_n_meas) / g_n_meas
    dom = "in-domain (r>r_p)" if in_domain else "EXTRAPOLATION (r<r_p)"
    print(f"  {label:34s} r={r_fm:.4f} fm  g_n(bound)={g:+.4f} mu_N  "
          f"err={err:+.2f}%  [{dom}, n={n}]")

report("r_p (original flat-reuse value)", r_p_fm, True)
report("r_exact (Section 2 primary, n=3)", r_exact[3], True)
report("r_grind = 2*lambda_p (Zone 2 contact)", r_grind_fm, False)
report("r_0 = pion force range", r_pion_fm, True)

# Reid (1968) NN potential attractive-well minimum (reuses the SAME formula
# already committed in nn_potential_external_comparison.py, not a new fit)
def V_reid(r_fm):
    mu = 0.7
    x = mu * r_fm
    return (-10.463 * math.exp(-x) / x
            - 1650.6 * math.exp(-4 * x) / x
            + 6484.2 * math.exp(-7 * x) / x)

# scan finely for the attractive well minimum
best_r, best_V = None, 1e9
r = 0.5
while r <= 2.0:
    V = V_reid(r)
    if V < best_V:
        best_V, best_r = V, r
    r += 0.0005

report(f"Reid(1968) well minimum (r={best_r:.4f} fm)", best_r, best_r > r_p_fm)

print()
print(f"  Reid well minimum location: r = {best_r:.4f} fm, V_Reid = {best_V:.2f} MeV")
print(f"  (real, external, peer-reviewed NN scattering fit -- Reid 1968, already")
print(f"  used in nn_potential_external_comparison.py; not refit or tuned here)")
print(f"  Note: r_exact(n=3) ({r_exact[3]:.4f} fm), Reid well minimum ({best_r:.4f} fm),")
print(f"  and r_p ({r_p_fm:.4f} fm) are all within ~2% of each other -- a real,")
print(f"  reproducible numerical closeness worth noting, not necessarily causal.")

check("BD4 Reid well minimum also falls in-domain (r > r_p), consistent with "
      "r_exact also falling in-domain",
      best_r > r_p_fm,
      f"Reid well minimum r = {best_r:.4f} fm vs r_p = {r_p_fm:.4f} fm")

# ── Section 4: proton-spin-measurement caveat (flagged, not used) ────────────
print()
print(SEP2)
print("SECTION 4: proton spin structure -- caveat, not incorporated")
print(SEP2)
print("  Real experiments (polarized deep-inelastic scattering: EMC, COMPASS,")
print("  HERMES; the 'proton spin puzzle') measure proton spin structure via")
print("  the g1 structure function, finding quark spin contributes only ~30%")
print("  of the proton's total spin (rest: gluon spin + orbital angular")
print("  momentum). This is a REAL, measured constraint on proton spin -- but")
print("  it is a QCD quark+gluon+orbital decomposition, a DIFFERENT accounting")
print("  than this framework's geometric Zone 3 medium co-rotation picture.")
print("  No established mapping between the two exists in this repo. This")
print("  script does NOT attempt one; the r_exact result above does not rely")
print("  on, or claim support from, that measurement. Flagged for future work.")
check("BD5 (informational, always True) proton-spin-measurement caveat stated, "
      "not silently ignored",
      True)

# ── Summary ───────────────────────────────────────────────────────────────
print()
print(SEP)
total = len(results)
passed = sum(1 for _, s, _ in results if s == "PASS")
print(f"  Total: {total}  PASS: {passed}  FAIL: {total-passed}")
print(SEP)
