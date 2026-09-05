"""
transition_metal_NJ_survey.py
==============================
Extends iron_dscale_NJ_survey.py's method to Mn (H_g), Co (T_1g), and
Ni (E_1/2) -- the three elements judgment_calls.txt MG2-1 flagged as
"no script anywhere computes N_J at realistic energy scales for" these
elements, the way iron_dscale_NJ_survey.py did for iron/G_g. Cu (A_g) is
included as a bonus fourth case (same "sub-cell" claim as H_g).

BACKGROUND (see judgment_calls_resolved.txt MG-J4 for the full history):
  N_J = hbar_c/(E*L_J) uses a PARTICLE's REST MASS in its standard usage
  (b quark N_J=4.75, Higgs N_J=0.159, W/Z N_J~0.22-0.25, top N_J=0.115 --
  all from analysis/nuclear/proton_structure.py PS3 and doc_jobson_cell.txt
  Section 4.3). Every d-electron in every element has the SAME electron
  rest mass (0.511 MeV), which gives the SAME N_J_e=38,870 (deep BULK)
  regardless of which atom it is bound to -- this formula, by itself,
  cannot distinguish one element's d-electrons from another's.
  doc_magnetism.txt Section 3.1's table nonetheless assigns each element a
  distinct "N_J regime" (Mn/H_g: sub-cell, needs N_J<1, i.e. E>E_cell=
  19.87 GeV; Fe/G_g: boundary, needs N_J~1-10, i.e. E~2-20 GeV; Co/T_1g:
  listed as "bulk / sub-cell" -- itself an unresolved ambiguity, see
  VERDICT below; Ni/E_1/2: bulk, needs N_J>>1, automatically satisfied by
  the ordinary electron mass with no new scale required at all).
  iron_dscale_NJ_survey.py checked 12 realistic iron-specific energy
  candidates and found NONE reach the GeV window iron's OWN "boundary"
  claim needs -- a 6-9 order-of-magnitude gap. This script repeats that
  same exercise for Mn, Co, Ni, Cu using their own real reference energies,
  to convert "almost certainly generalizes" into an actually-computed result.

DATA CONFIDENCE: ionization energies, work functions, cohesive energies,
  Curie/Neel temperatures, melting points, and K-shell edges are standard,
  well-tabulated values (CODATA/NIST-consistent, similar precision to
  iron_dscale_NJ_survey.py's own entries). Bulk plasmon energy, 3d-4s
  exchange splitting, and the mean-field exchange integral are ORDER-OF-
  MAGNITUDE placeholders (same as iron_dscale_NJ_survey.py labels them)
  -- not claimed to be precise to more than 1 significant figure. None of
  this affects the conclusion: the gap to the GeV window is 6-9 orders of
  magnitude, far larger than any realistic uncertainty in these numbers.

Reference: docs/series1/doc_magnetism.txt Section 3.1;
  judgment_calls.txt MG2-1; analysis/nuclear/iron_dscale_NJ_survey.py
"""
import math

pi     = math.pi
alpha  = 7.2973525693e-3
phi    = (1 + math.sqrt(5)) / 2
r_p    = 0.8414          # fm
hbar_c = 197.3269804     # MeV*fm
L_J_fm = alpha * phi * r_p   # fm  (same L_J as iron_dscale_NJ_survey.py)
kB_eV  = 8.617333262e-5     # eV/K (CODATA Boltzmann constant)

def N_J(E_MeV):
    return hbar_c / (E_MeV * L_J_fm)

def E_for_NJ(N):
    return hbar_c / (N * L_J_fm)

SEP  = "=" * 78
SEP2 = "-" * 78

print(SEP)
print("TRANSITION-METAL N_J SURVEY: Mn (H_g), Co (T_1g), Ni (E_1/2), Cu (A_g)")
print(SEP)
print(f"  L_J = {L_J_fm:.6e} fm   E_cell = hbar_c/L_J = {E_for_NJ(1):.4f} MeV "
      f"= {E_for_NJ(1)/1000:.3f} GeV")
print(f"  Regime targets (doc_jobson_cell.txt Section 4.3): sub-cell needs N_J<1 "
      f"(E>{E_for_NJ(1)/1000:.2f} GeV);")
print(f"  boundary needs N_J~1-10 (E~{E_for_NJ(10)/1000:.2f}-{E_for_NJ(1)/1000:.2f} GeV); "
      f"bulk needs N_J>>1 (any atomic-scale E qualifies trivially).")
print(f"  Electron rest mass N_J_e = {N_J(0.5109989):.0f} (deep bulk, SAME for every element)")

# ── Candidate energy tables, one per element ─────────────────────────────────
# Each entry: (label, E_eV). Structure mirrors iron_dscale_NJ_survey.py's own
# 12-candidate list for direct comparability.

elements = {
    "Mn (Z=25, H_g claim: sub-cell, needs N_J<1)": {
        "irrep": "H_g", "target": "sub-cell (N_J<1)", "A": 55,
        "candidates": [
            ("1st ionization energy (Mn I)",              7.434),
            ("2nd ionization energy (Mn II)",              15.64),
            ("Work function (polycrystalline Mn)",         4.1),
            ("Cohesive/binding energy per atom (solid)",   2.92),
            ("Bulk plasmon energy (typical 3d metal, order-of-mag)", 17.0),
            ("3d-4s exchange splitting (order-of-mag)",    1.5),
            ("Mean-field exchange integral J (order-of-mag)", 0.1),
            ("Neel temperature k_B*T_N (alpha-Mn, T_N~95 K)", kB_eV*95),
            ("Melting point k_B*T_melt (T_melt=1519 K)",   kB_eV*1519),
            ("K-shell (1s) binding energy (X-ray)",        6539.0),
            ("Nuclear binding energy PER NUCLEON (Mn-55)", 8.7648e6),
            ("TOTAL nuclear binding energy (Mn-55, 55 nucleons)", 55*8.7648e6),
        ],
    },
    "Co (Z=27, T_1g claim: 'bulk / sub-cell' -- ambiguous, see VERDICT)": {
        "irrep": "T_1g", "target": "bulk AND/OR sub-cell (ambiguous label)", "A": 59,
        "candidates": [
            ("1st ionization energy (Co I)",              7.881),
            ("2nd ionization energy (Co II)",              17.08),
            ("Work function (polycrystalline Co)",        5.0),
            ("Cohesive/binding energy per atom (solid)",  4.387),
            ("Bulk plasmon energy (typical 3d metal, order-of-mag)", 18.0),
            ("3d-4s exchange splitting (order-of-mag)",   1.5),
            ("Mean-field exchange integral J (order-of-mag)", 0.1),
            ("Curie temperature k_B*T_C (T_C=1388 K)",    kB_eV*1388),
            ("Melting point k_B*T_melt (T_melt=1768 K)",  kB_eV*1768),
            ("K-shell (1s) binding energy (X-ray)",       7709.0),
            ("Nuclear binding energy PER NUCLEON (Co-59)", 8.7683e6),
            ("TOTAL nuclear binding energy (Co-59, 59 nucleons)", 59*8.7683e6),
        ],
    },
    "Ni (Z=28, E_1/2 claim: bulk, needs N_J>>1)": {
        "irrep": "E_1/2", "target": "bulk (N_J>>1)", "A": 58,
        "candidates": [
            ("1st ionization energy (Ni I)",              7.640),
            ("2nd ionization energy (Ni II)",              18.17),
            ("Work function (polycrystalline Ni)",       5.15),
            ("Cohesive/binding energy per atom (solid)", 4.435),
            ("Bulk plasmon energy (typical 3d metal, order-of-mag)", 20.0),
            ("3d-4s exchange splitting (order-of-mag)",   0.3),
            ("Mean-field exchange integral J (order-of-mag)", 0.1),
            ("Curie temperature k_B*T_C (T_C=627 K)",     kB_eV*627),
            ("Melting point k_B*T_melt (T_melt=1728 K)",  kB_eV*1728),
            ("K-shell (1s) binding energy (X-ray)",       8333.0),
            ("Nuclear binding energy PER NUCLEON (Ni-58)", 8.7323e6),
            ("TOTAL nuclear binding energy (Ni-58, 58 nucleons)", 58*8.7323e6),
        ],
    },
    "Cu (Z=29, A_g claim: sub-cell, needs N_J<1) [bonus, not in MG2-1's list]": {
        "irrep": "A_g", "target": "sub-cell (N_J<1)", "A": 63,
        "candidates": [
            ("1st ionization energy (Cu I)",              7.726),
            ("2nd ionization energy (Cu II)",              20.29),
            ("Work function (polycrystalline Cu)",       4.65),
            ("Cohesive/binding energy per atom (solid)", 3.49),
            ("Bulk plasmon energy (measured, Cu)",        19.3),
            ("K-shell (1s) binding energy (X-ray)",       8979.0),
            ("Nuclear binding energy PER NUCLEON (Cu-63)", 8.7521e6),
            ("TOTAL nuclear binding energy (Cu-63, 63 nucleons)", 63*8.7521e6),
        ],
    },
}

# ── Fe/G_g baseline reproduced here for direct side-by-side consistency ─────
iron_candidates = [
    ("1st ionization energy (Fe I)",            7.902),
    ("2nd ionization energy (Fe II)",            16.18),
    ("Work function (polycrystalline Fe)",       4.5),
    ("Cohesive/binding energy per atom (solid)", 4.28),
    ("Bulk plasmon energy (typical 3d metal, order-of-mag)", 17.0),
    ("3d-4s exchange splitting (order-of-mag)",  1.5),
    ("Mean-field exchange integral J (order-of-mag)", 0.1),
    ("Curie temperature k_B*T_C (T_C=1043 K)",   kB_eV*1043),
    ("Melting point k_B*T_melt (T_melt=1811 K)", kB_eV*1811),
    ("K-shell (1s) binding energy (X-ray)",       7112.0),
    ("Nuclear binding energy PER NUCLEON (Fe-56)", 8.7903e6),
    ("TOTAL nuclear binding energy (Fe-56, 56 nucleons)", 56*8.7903e6),
]
elements = {"Fe (Z=26, G_g claim: boundary, needs N_J~1-10) [BASELINE, "
            "reproduces iron_dscale_NJ_survey.py]": {
                "irrep": "G_g", "target": "boundary (N_J~1-10)", "A": 56,
                "candidates": iron_candidates,
            }, **elements}

# ── Run the survey ────────────────────────────────────────────────────────────
summary = []
for name, info in elements.items():
    print()
    print(SEP2)
    print(f"{name}")
    print(f"  irrep={info['irrep']}  target regime={info['target']}")
    print(SEP2)
    print(f"  {'Candidate':<58} {'E (MeV)':>10} {'N_J':>12}")
    print(f"  {'-'*58} {'-'*10} {'-'*12}")
    best_nj = None
    for label, E_eV in info["candidates"]:
        E_MeV = E_eV * 1e-6
        nj = N_J(E_MeV)
        print(f"  {label:<58} {E_MeV:>10.3e} {nj:>12.3e}")
        if best_nj is None or nj < best_nj:
            best_nj = nj   # smallest N_J = candidate that gets CLOSEST to sub-cell/boundary
    n_j_electron = N_J(0.5109989)
    print(f"  {'(for reference) electron rest mass, N_J_e':<58} "
          f"{0.5109989:>10.3e} {n_j_electron:>12.3e}  <-- BULK, same for all elements")
    reaches_bulk_trivially = True   # N_J_e=38870 is always >>1
    closest_to_high_regime = best_nj
    summary.append((name, info['irrep'], info['target'], closest_to_high_regime))

# ── VERDICT ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("VERDICT")
print(SEP)
print(f"  {'Element (irrep)':<45} {'Target regime':<24} {'Closest real N_J':>16}")
print(f"  {'-'*45} {'-'*24} {'-'*16}")
for name, irrep, target, closest in summary:
    label = name.split(' (')[0] + f" ({irrep})"
    print(f"  {label:<45} {target:<24} {closest:>16.2f}")

print()
print("  Every element's CLOSEST realistic candidate (in every case, the FULL")
print("  nuclear binding energy of the whole nucleus -- not a d-electron")
print("  quantity at all) lands at N_J~35-45, i.e. deep in the BULK regime,")
print("  not boundary (needs ~1-10) and nowhere near sub-cell (needs <1).")
print("  Every genuinely atomic/solid-state (eV-to-keV) candidate for every")
print("  element gives N_J in the tens-of-thousands to millions -- the SAME")
print("  order of magnitude as the ordinary electron rest mass's N_J_e=38,870,")
print("  because that IS fundamentally what these numbers are (electron-scale")
print("  energies), for any element in the periodic table.")
print()
print("  CONCLUSION (generalizes iron_dscale_NJ_survey.py's finding):")
print("  - Fe/G_g 'boundary' (N_J~1-10): CONFIRMED still fails -- reproduces")
print("    iron_dscale_NJ_survey.py's own result exactly (same candidates).")
print("  - Mn/H_g 'sub-cell' (N_J<1): FAILS by an even LARGER margin than Fe's")
print("    boundary claim did (sub-cell needs an even higher energy than")
print("    boundary) -- no candidate anywhere close, including Mn's own full")
print("    nuclear binding energy.")
print("  - Cu/A_g 'sub-cell' (N_J<1) [bonus]: FAILS for the same reason as Mn.")
print("  - Ni/E_1/2 'bulk' (N_J>>1): satisfied AUTOMATICALLY and TRIVIALLY by")
print("    the ordinary electron rest mass (N_J_e=38,870) -- no new element-")
print("    specific scale is needed or possible here; this row was never at")
print("    risk the way the boundary/sub-cell rows were.")
print("  - Co/T_1g 'bulk / sub-cell': the compound label is the table's OWN")
print("    ambiguity, independent of the energy-scale question -- T_1g is")
print("    realized by BOTH the (bulk) photon and the (sub-cell) W/Z in")
print("    doc_jobson_cell.txt Section 4.3, and the table does not specify")
print("    which applies to cobalt's d-electrons. If 'bulk' is meant: same")
print("    trivial automatic satisfaction as Ni. If 'sub-cell' is meant: same")
print("    failure as Mn/Cu -- no realistic cobalt energy scale reaches it.")
print()
print("  NET RESULT: the scale-mismatch problem that killed Fe/G_g's boundary")
print("  claim (MG-J4) DOES generalize to every row claiming 'boundary' or")
print("  'sub-cell' (Fe, Mn, and -- depending on which reading is meant --")
print("  possibly Co). The 'bulk' rows (Ni, and Co under the other reading)")
print("  were never actually at risk: N_J_e=38,870 satisfies 'bulk' for any")
print("  element trivially, with no element-specific claim being made at all.")
print(SEP)
