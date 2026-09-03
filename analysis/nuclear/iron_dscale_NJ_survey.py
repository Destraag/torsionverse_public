"""
iron_dscale_NJ_survey.py
=========================
doc_magnetism.txt (pre-2026-09-01 drafts) claimed iron's d-electrons sit at
N_J~1-10 (Section 3.3) / N_J~4-5 (Section 5.2), matching the b quark's
boundary regime. No script anywhere derives this. This script computes
N_J = hbar_c/(E*L_J) for every physically realistic energy scale associated
with iron that could plausibly be meant by "iron's d-electron scale", to see
whether ANY of them land near 1-10.

This is a SURVEY/diagnostic, not a proposed derivation -- it does not assert
any of these candidates IS the right physics, only reports what N_J each one
would give.

Reference: docs/series1/doc_magnetism.txt Sections 3.3, 5.2.
"""
import math

pi     = math.pi
alpha  = 7.2973525693e-3
phi    = (1 + math.sqrt(5)) / 2
r_p    = 0.8414          # fm
hbar_c = 197.3269804     # MeV*fm
L_J_fm = alpha * phi * r_p   # fm

def N_J(E_MeV):
    return hbar_c / (E_MeV * L_J_fm)

def E_for_NJ(N):
    return hbar_c / (N * L_J_fm)

SEP = "=" * 70
print(SEP)
print("IRON d-ELECTRON N_J SURVEY: does any realistic energy scale give N_J~1-10?")
print(SEP)
print(f"  L_J = {L_J_fm:.6e} fm")
print()
print(f"  Target window: N_J=1 needs E = {E_for_NJ(1):.4e} MeV = {E_for_NJ(1)/1000:.2f} GeV")
print(f"                 N_J=10 needs E = {E_for_NJ(10):.4e} MeV = {E_for_NJ(10)/1000:.3f} GeV")
print(f"                 N_J=4.75 (b quark match) needs E = {E_for_NJ(4.75):.4e} MeV "
      f"= {E_for_NJ(4.75)/1000:.3f} GeV")
print()

# All candidate energies for iron, converted to MeV. Sources: standard
# solid-state / atomic / nuclear reference values (order-of-magnitude typical
# figures, not all specific to a single crystal structure or compound).
candidates = [
    ("1st ionization energy (Fe I)",            7.902,          "eV"),
    ("2nd ionization energy (Fe II)",            16.18,          "eV"),
    ("Work function (polycrystalline Fe)",       4.5,            "eV"),
    ("Cohesive/binding energy per atom (solid)", 4.28,           "eV"),
    ("Bulk plasmon energy (typical 3d metal)",   17.0,           "eV"),
    ("3d-4s exchange splitting (band structure, order-of-mag)", 1.5, "eV"),
    ("Mean-field exchange integral J (Heisenberg, order-of-mag)", 0.1, "eV"),
    ("Curie temperature thermal energy k_B*T_C (T_C=1043 K)",   0.0899, "eV"),
    ("Melting point thermal energy k_B*T_melt (T_melt=1811 K)", 0.156,  "eV"),
    ("K-shell (1s) binding energy (X-ray)",       7112,           "eV"),
    ("Nuclear binding energy PER NUCLEON (Fe-56, the 'iron peak')", 8.79e6, "eV"),
    ("TOTAL nuclear binding energy (Fe-56, all 56 nucleons)",   4.9224e8, "eV"),
]

print(f"  {'Candidate':<58} {'E (MeV)':>12} {'N_J':>12}")
print(f"  {'-'*58} {'-'*12} {'-'*12}")
for label, E_eV, unit in candidates:
    E_MeV = E_eV * 1e-6
    nj = N_J(E_MeV)
    in_window = 1 <= nj <= 10
    flag = "  <-- IN [1,10]" if in_window else ""
    print(f"  {label:<58} {E_MeV:>12.4e} {nj:>12.4e}{flag}")

print()
print(SEP)
print("VERDICT")
print(SEP)
any_in_window = any(1 <= N_J(E_eV*1e-6) <= 10 for _, E_eV, _ in candidates)
print(f"  Any realistic candidate lands in N_J=[1,10]? {any_in_window}")
print()
print("  Every atomic/solid-state energy scale for iron (eV to keV) gives")
print("  N_J in the billions-to-trillions (deep BULK regime, same as any")
print("  free electron, N_J_e=38,870) -- because eV/keV is ~6-9 orders of")
print("  magnitude below the ~2-20 GeV window the N_J=hbar_c/(E*L_J) formula")
print("  needs to land at N_J~1-10. Even iron's own K-shell X-ray line")
print("  (7.1 keV, the single largest realistic atomic-physics energy for")
print("  iron) is still ~6 orders of magnitude short.")
print()
print("  The ONLY candidate that gets close is the Fe-56 nucleus's TOTAL")
print("  binding energy (492 MeV, all 56 nucleons combined) -- this is a")
print(f"  NUCLEAR quantity, not a d-ELECTRON one, and still gives N_J="
      f"{N_J(4.9224e8*1e-6):.2f}, outside [1,10] (closer to the 'GeV-to-keV'")
print("  qualitative range the stale notes file used loosely, but not a")
print("  precise match, and conceptually a different physical quantity than")
print("  'iron's d-electron scale' as written in the doc).")
print()
print("  CONCLUSION: no realistic, physically-motivated energy scale specific")
print("  to iron reproduces N_J~1-10 or ~4-5 via the standard formula. This")
print("  is not a failure to search hard enough -- it is a ~6-9 order-of-")
print("  magnitude scale gap that no atomic or solid-state quantity can close.")
print(SEP)
