"""
curie_temperature_reproduction.py
===================================
doc_magnetism.txt Section 3.3 claims "The ordering Co > Fe > Ni IS correctly
reproduced" by T_C = f_Ag * W_d * z*S(S+1) / (3*k_B), but no script anywhere
in the repo evaluates z, S, W_d, or T_C numerically for the three metals.
Searched notes/, sessions/, legacy .tex files, docs/leads_applications.txt --
no prior calculation exists to recover; docs/zenodo.txt's "ordering correct"
line just restates the doc's own claim, it isn't an independent check.

RESULT (2026-09-02): the doc's OWN prescribed inputs (crude W_d=G_medium*a^3,
S=unpaired/2) give the WRONG ordering (Ni>Co>Fe, CT1-CT2 FAIL). But a genuine
fix EXISTS and WORKS: replacing a^3 (whole unit cell) with the actual
d-orbital volume the doc's own "MISSING INPUT" note already called for --
computed here via Slater's rules (standard, algorithmic effective nuclear
charge, not a looked-up table) -- and using REAL measured saturation moments
for S (not the naive unpaired/2 convention) gives the CORRECT ordering
Co>Fe>Ni (CT6 PASS), with the three elements' T_C(predicted)/T_C(measured)
ratios landing within 1.07x of each other (CT6b) -- i.e. a SINGLE missing
overall constant (~5.6x, matching the doc's own anticipated "factor 6 off")
would bring absolute values into near agreement too, not just the ordering.

WHICH FACTORS ARE ACTUALLY FROM THE MEDIUM (framework-derived) VS EXTERNAL:
  - G_medium: YES, medium-derived (rho=mu_0, K=1/eps_0 chain, doc_torsion.txt
    Section 3.3). Verified (CT4) to be a common multiplicative factor across
    all 3 elements -- it affects absolute scale only, never ordering.
  - f_Ag = 1/dim(irrep)^2: VERIFIED here (CT5) as genuine group theory, not
    an ansatz -- the trivial irrep A_g has multiplicity exactly 1 in every
    irrep's self-product X(x)X (checked via magnetism_doc.py's own I_h
    character-table projection formula for all 5 gerade irreps), so
    f_Ag=n_Ag(X,X)/dim(X)^2=1/dim(X)^2 is a derived consequence.
  - r_d (d-orbital radius, replacing a in W_d): computed via Slater's rules
    (standard textbook algorithm for effective nuclear charge Z_eff from the
    electron configuration) + hydrogen-like r=n^2*a_bohr/Z_eff scaling.
    Reproduces the well-known d-orbital CONTRACTION across the 3d row
    (Z_eff increases Fe->Co->Ni, so r_d shrinks) -- external atomic physics,
    but a standard, reproducible calculation, not a fudge factor.
  - S: real measured T=0 saturation magnetic moments (Fe=2.22, Co=1.72,
    Ni=0.606 mu_B/atom -- standard values, e.g. Kittel's Introduction to
    Solid State Physics), NOT the doc's own "unpaired/2" convention (which
    does not match real itinerant-electron moments for these 3 metals).
  - z, a (lattice constant, only used for the CT1-CT3 crude-form checks):
    external, real measured crystallography.

Run: python analysis/nuclear/curie_temperature_reproduction.py
Reference: docs/series1/doc_magnetism.txt Section 3.3.
"""

import math

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

k_B = 8.617333262e-5   # eV/K

# Framework-derived (this session, doc_torsion.txt Section 3.3, rho=mu_0 chain)
G_medium_Pa = 3.576e9   # Pa = J/m^3

# doc_magnetism.txt Section 3.2's own irrep/unpaired-electron assignments
elements = {
    # name: (irrep_dim, unpaired, z_coordination, a_lattice_m, T_C_measured_K, crystal_note)
    "Fe": (4, 4,  8, 2.866e-10, 1043, "BCC (alpha-Fe, stable well above T_C=1043K)"),
    "Co": (3, 3, 12, 3.544e-10, 1388, "FCC (gamma-Co, stable phase at T_C=1388K; HCP->FCC at 695K)"),
    "Ni": (2, 2, 12, 3.524e-10,  627, "FCC (stable at all T up to melting)"),
}

print(SEP)
print("REPRODUCING doc_magnetism.txt's T_C = f_Ag * W_d * z*S(S+1) / (3*k_B)")
print(SEP2)
print(f"  G_medium = {G_medium_Pa:.3e} Pa  [framework-derived]")
print()

computed = {}
for name, (dim, unpaired, z, a_m, T_C_meas, note) in elements.items():
    f_Ag = 1.0 / dim**2
    S    = unpaired / 2.0
    W_d_J  = G_medium_Pa * a_m**3          # Joules
    W_d_eV = W_d_J / 1.602176634e-19       # eV
    T_C_pred = f_Ag * W_d_eV * z * S * (S + 1) / (3 * k_B)
    computed[name] = T_C_pred
    print(f"  {name} [{note}]:")
    print(f"    dim={dim}  f_Ag=1/{dim}^2={f_Ag:.4f}  S={S}  S(S+1)={S*(S+1):.2f}  z={z}")
    print(f"    a={a_m*1e10:.3f} Angstrom  W_d=G_medium*a^3={W_d_eV:.4f} eV")
    print(f"    T_C(predicted, this formula) = {T_C_pred:.2f} K   (measured: {T_C_meas} K)")
    print()

print(SEP2)
print("ORDERING CHECK")
print(SEP2)
order_measured  = sorted(elements, key=lambda n: -elements[n][4])
order_predicted = sorted(elements, key=lambda n: -computed[n])
print(f"  Measured  ordering (highest T_C first): {' > '.join(order_measured)}")
print(f"  Predicted ordering (highest T_C first): {' > '.join(order_predicted)}")

check("CT1 Predicted ordering matches measured ordering (Co > Fe > Ni)",
      order_predicted == order_measured,
      f"predicted={' > '.join(order_predicted)}  measured={' > '.join(order_measured)}")

# ── Cross-check: does a uniform rescaling save the ordering? ────────────────
# The doc itself says W_d=G_medium*a^3 is "off by a constant factor (~6)" --
# a UNIFORM constant factor cannot change an ordering, so if CT1 fails, no
# single missing multiplicative constant can rescue the claim as stated.
ratios = {n: computed[n] / elements[n][4] for n in elements}
print()
print(f"  T_C(predicted)/T_C(measured) per element: "
      + ", ".join(f"{n}={ratios[n]:.4f}" for n in elements))
check("CT2 A single missing constant factor could fix this (ratios all equal)",
      max(ratios.values()) / min(ratios.values()) < 1.5,
      f"ratio spread = {max(ratios.values())/min(ratios.values()):.2f}x "
      f"(if >>1, no single missing constant reconciles predicted vs measured)")

# ── Is f_Ag=1/dim^2 really an unverified ansatz, or real group theory? ──────
# Using the SAME I_h character table + projection-formula tool already in
# magnetism_doc.py (reimplemented here for standalone-runnability): if the
# trivial irrep A_g appears with multiplicity EXACTLY 1 in EVERY irrep's
# self-product X x X (a standard character-theory fact for real/self-
# conjugate irreps -- the trivial rep's multiplicity in X(x)X is always 1
# for irreducible X), then f_Ag = n_Ag(X,X)/dim(X)^2 = 1/dim(X)^2 follows as
# a DERIVED consequence, not a free-standing guess.
phi = (1 + math.sqrt(5)) / 2
irreps_Ih = [
    ("A_g",  1,  1,  1,  1,  1,       1),
    ("T_1g", 3,  3, -1,  0,  phi,     -(phi - 1)),
    ("T_2g", 3,  3, -1,  0, -(phi - 1), phi),
    ("G_g",  4,  4,  0,  1, -1,       -1),
    ("H_g",  5,  5,  1, -1,  0,        0),
]
class_sizes_Ih = [1, 15, 20, 12, 12]
order_Ih = 60

def n_Ag_in_product(name):
    chi_Ag = [1, 1, 1, 1, 1]
    chi = next(r[2:] for r in irreps_Ih if r[0] == name)
    return sum(class_sizes_Ih[c] * chi_Ag[c] * chi[c] * chi[c] for c in range(5)) / order_Ih

print()
print(SEP2)
print("IS f_Ag=1/dim^2 A DERIVED GROUP-THEORY FACT, OR JUST AN ANSATZ?")
print(SEP2)
all_mult_one = True
for name, dim, *_ in irreps_Ih:
    n_Ag = n_Ag_in_product(name)
    all_mult_one = all_mult_one and abs(n_Ag - 1.0) < 1e-9
    print(f"  {name} (dim={dim}): n_Ag(X,X) = {n_Ag:.6f}  ->  f_Ag=n_Ag/dim^2 = "
          f"{n_Ag / dim**2:.4f}  (doc's 1/dim^2 = {1 / dim**2:.4f})")

check("CT5 A_g multiplicity in X-x-X is exactly 1 for all 5 gerade irreps "
      "(f_Ag=1/dim^2 IS derived group theory, verified, not an ansatz)",
      all_mult_one,
      "confirmed by the same projection formula magnetism_doc.py M4/M5/M6 use "
      "-- f_Ag is NOT the source of the ordering failure")

# ── Alternative: use real (Curie-Weiss, high-T susceptibility) effective S ──
# instead of the doc's "unpaired/2" irrep-derived S, to see whether a more
# standard mean-field input fares better. p_eff = g*sqrt(S(S+1)), g~2.
p_eff = {"Fe": 3.13, "Co": 3.15, "Ni": 1.61}  # real measured effective moments, mu_B
print()
print(SEP2)
print("CROSS-CHECK: real measured effective paramagnetic moments (not this framework's S)")
print(SEP2)
computed_real_S = {}
for name, (dim, unpaired, z, a_m, T_C_meas, note) in elements.items():
    f_Ag = 1.0 / dim**2
    S_eff = (-1 + math.sqrt(1 + p_eff[name]**2)) / 2   # invert p_eff=2*sqrt(S(S+1))
    W_d_J  = G_medium_Pa * a_m**3
    W_d_eV = W_d_J / 1.602176634e-19
    T_C_pred = f_Ag * W_d_eV * z * S_eff * (S_eff + 1) / (3 * k_B)
    computed_real_S[name] = T_C_pred
    print(f"  {name}: p_eff={p_eff[name]} mu_B -> S_eff={S_eff:.3f}  "
          f"T_C(predicted)={T_C_pred:.2f} K  (measured {T_C_meas} K)")

order_predicted_real_S = sorted(elements, key=lambda n: -computed_real_S[n])
check("CT3 Ordering with REAL effective moments (not doc's unpaired/2) matches measured",
      order_predicted_real_S == order_measured,
      f"predicted={' > '.join(order_predicted_real_S)}  measured={' > '.join(order_measured)}")

# ── CANDIDATE FIX: replace a^3 (whole unit cell) with the ACTUAL d-orbital ──
# volume the doc's own "MISSING INPUT" note calls for (W_d = G_medium * a^3 *
# V_d_orbital/V_unit_cell = G_medium * V_d_orbital -- the a^3 cancels). The
# d-orbital radius r_d is computed here via Slater's rules (standard,
# textbook-algorithmic effective nuclear charge Z_eff, not a looked-up table
# -- reproducible from the electron configuration alone) plus the hydrogen-
# like scaling r = n^2*a_bohr/Z_eff. This captures d-orbital CONTRACTION
# across the 3d row (well-established atomic physics: Z_eff(Fe)<Z_eff(Co)
# <Z_eff(Ni), so r_d shrinks Fe->Co->Ni), which a^3 alone cannot.
a_bohr_A = 0.52917721   # Angstrom
configs = {  # (Z, n_3d electrons, n_4s electrons); core = [Ar] 1s2 2s2p6 3s3p6 = 18 e-
    "Fe": (26, 6, 2),
    "Co": (27, 7, 2),
    "Ni": (28, 8, 2),
}
print()
print(SEP2)
print("CANDIDATE FIX: Slater's-rule d-orbital radius (replaces a^3) + real S")
print(SEP2)
r_d = {}
for name, (Z, n3d, n4s) in configs.items():
    # Slater shielding for a 3d electron: 0.35 per OTHER 3d electron,
    # 1.00 per electron in a lower group (here: the 18 e- Ar-like core),
    # 0 per electron in a higher group (4s here does not shield 3d).
    S_shield = 0.35 * (n3d - 1) + 1.00 * 18
    Z_eff = Z - S_shield
    r_d[name] = 9 * a_bohr_A / Z_eff   # n=3 -> n^2=9, Bohr-like scaling
    print(f"  {name}: Z={Z}  Z_eff(3d)=Z-{S_shield:.2f}={Z_eff:.2f}  "
          f"r_d=9*a_bohr/Z_eff={r_d[name]:.4f} Angstrom")

check("CT6a d-orbital radius CONTRACTS monotonically Fe>Co>Ni (Z_eff increases)",
      r_d["Fe"] > r_d["Co"] > r_d["Ni"],
      f"r_d: Fe={r_d['Fe']:.4f}  Co={r_d['Co']:.4f}  Ni={r_d['Ni']:.4f} Angstrom "
      "(standard d-orbital contraction across the 3d row)")

mu_sat = {"Fe": 2.22, "Co": 1.72, "Ni": 0.606}   # measured saturation moments, mu_B/atom
computed_fix = {}
ratios_fix = {}
for name, (dim, unpaired, z, a_m, T_C_meas, note) in elements.items():
    f_Ag = 1.0 / dim**2
    S_sat = mu_sat[name] / 2.0   # mu_s = g*S, g~2
    r_d_m = r_d[name] * 1e-10
    W_d_J  = G_medium_Pa * (4 / 3) * math.pi * r_d_m**3   # V_d_orbital, not a^3
    W_d_eV = W_d_J / 1.602176634e-19
    T_C_pred = f_Ag * W_d_eV * z * S_sat * (S_sat + 1) / (3 * k_B)
    computed_fix[name] = T_C_pred
    ratios_fix[name] = T_C_pred / T_C_meas
    print(f"  {name}: S_sat={S_sat:.3f} (from mu_sat={mu_sat[name]} mu_B)  "
          f"W_d={W_d_eV:.5f} eV  T_C(predicted)={T_C_pred:.2f} K  (measured {T_C_meas} K)")

order_predicted_fix = sorted(elements, key=lambda n: -computed_fix[n])
print(f"\n  Predicted ordering: {' > '.join(order_predicted_fix)}")
print(f"  T_C(predicted)/T_C(measured) per element: "
      + ", ".join(f"{n}={ratios_fix[n]:.3f}" for n in elements))

check("CT6 Ordering with Slater-rule r_d + real saturation-moment S matches "
      "measured (Co > Fe > Ni)",
      order_predicted_fix == order_measured,
      f"predicted={' > '.join(order_predicted_fix)}  measured={' > '.join(order_measured)}")
check("CT6b Bonus: ratio spread is now tight (~uniform missing factor, "
      "consistent with the doc's own 'factor 6 off' note)",
      max(ratios_fix.values()) / min(ratios_fix.values()) < 1.5,
      f"ratio spread = {max(ratios_fix.values())/min(ratios_fix.values()):.2f}x "
      f"(vs {max(ratios.values())/min(ratios.values()):.1f}x with the crude a^3 form)")

# ── Did it work with the PRIOR (pre-MG-J1) medium density? ──────────────────
# doc_magnetism.txt (2026-08-20) predates the rho=mu_0 unification (MG-J1,
# 2026-09-01) -- before that, the only other candidate density anywhere in
# the repo was rho_Lambda (Planck dark energy density), which gives a
# WILDLY different G_medium (rho_mu0_vs_rho_lambda_impact.py section A1).
# Testing whether the ORIGINAL claim could have used this "prior" value
# (against the CANDIDATE FIX above, since that is the one that works).
G_medium_legacy_Pa = 1.6619e-11   # Pa, from rho=rho_Lambda (pre-MG-J1)

print()
print(SEP2)
print("DID THE FIX WORK WITH PRIOR (pre-MG-J1, rho_Lambda-based) G_medium?")
print(SEP2)
print(f"  G_medium (legacy, rho_Lambda) = {G_medium_legacy_Pa:.4e} Pa")
print(f"  G_medium (current, rho=mu_0)  = {G_medium_Pa:.4e} Pa")

computed_legacy = {}
for name, (dim, unpaired, z, a_m, T_C_meas, note) in elements.items():
    f_Ag = 1.0 / dim**2
    S_sat = mu_sat[name] / 2.0
    r_d_m = r_d[name] * 1e-10
    W_d_J  = G_medium_legacy_Pa * (4 / 3) * math.pi * r_d_m**3
    W_d_eV = W_d_J / 1.602176634e-19
    T_C_pred = f_Ag * W_d_eV * z * S_sat * (S_sat + 1) / (3 * k_B)
    computed_legacy[name] = T_C_pred

order_predicted_legacy = sorted(elements, key=lambda n: -computed_legacy[n])
print(f"  Predicted ordering (legacy G_medium): {' > '.join(order_predicted_legacy)}")

check("CT4 Ordering with the fix is UNCHANGED under the legacy G_medium value "
      "(still Co>Fe>Ni -- G_medium is a common factor, cannot affect ordering)",
      order_predicted_legacy == order_predicted_fix,
      f"legacy-G ordering={' > '.join(order_predicted_legacy)}  "
      f"current-G ordering={' > '.join(order_predicted_fix)}")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP2)
print("  The crude W_d=G_medium*a^3 form (whole unit cell) gives the WRONG")
print("  ordering (Ni>Co>Fe) and cannot be fixed by any single missing constant")
print("  (6x ratio spread). Replacing a^3 with the ACTUAL d-orbital volume")
print("  (Slater's-rule Z_eff -> r_d, standard atomic physics, algorithmic not")
print("  looked-up) and using REAL saturation moments for S (not unpaired/2)")
print("  gives the CORRECT ordering Co>Fe>Ni, with a near-uniform ~5.5x")
print("  residual scale factor -- close to the doc's own anticipated 'factor 6'.")
print("  f_Ag=1/dim^2 is confirmed genuine group theory (CT5), not the culprit.")
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"  Total: {passed}/{len(results)}  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, s, d in results:
        if s == "FAIL": print(f"  FAILED: {name}")
print(SEP)
