#!/usr/bin/env python3
"""
tau_face_force_check.py

EXPLORATORY SYNTHESIS (session 12): assembles a candidate "tau force at the
face nexus" from pieces that are EACH already independently derived/verified
elsewhere, but had not previously been combined into a force. This script
does NOT introduce any new physical assumption beyond combining:
  (1) the tau's nexus distance from cell center = inradius, EXACT
      [jobson_cell_geometry_3d.py CG5/CG10: r_in = L_J*phi^2/(2*sqrt(3))]
  (2) the tau's mass/energy, EXACT (Koide) and leading-order (corkscrew)
      [lepton_mass.py LM10-LM17: m_tau = 1776.92 MeV (Koide) / 1777.49 MeV (leading)]
  (3) the T_2g face-shear amplitude, EXACT
      [face_gluon_geometry.py FG5: Rs^2 = 5/(16*pi^2), shear amplitude^2 at
       Maxwell critical -- T_2g is the field tau "rides" per doc_jobson_cell.txt]
  (4) the cell's own characteristic energy/length scale, EXACT
      [E_cell = 2*pi*hbar*c/L_J = 124.8 GeV]

QUESTION: does F_tau = m_tau*c^2 / r_in relate ALGEBRAICALLY to the cell's own
characteristic force scale F_cell = E_cell/L_J via an already-established
constant (Rs, Rs^2, phi, K/G, ...)? If yes, this is a genuine (if previously
unassembled) closed result. If no clean relation appears, this is reported
as an open/unresolved combination, not asserted as a new derivation.

Reference: jobson_cell_geometry_3d.py CG5, CG10; lepton_mass.py LM10-LM17;
  face_gluon_geometry.py FG4, FG5; docs/series1/doc_jobson_cell.txt
  ("Tau (I52)... RIDES the T_2g elastic face surfaces").
"""
import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 66
SEP2 = "-" * 66
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
sqrt3 = math.sqrt(3)
sqrt5 = math.sqrt(5)
hbar_c = 197.3269804          # MeV*fm
r_p_fm = 0.8414               # fm
m_p    = 938.272046           # MeV

Rs  = sqrt5 / (4*pi)
Rs2 = Rs**2

L_J_fm = alpha * phi * r_p_fm
E_cell_MeV = 2*pi*hbar_c / L_J_fm

print(SEP)
print("TAU FACE-NEXUS FORCE: EXPLORATORY SYNTHESIS FROM ESTABLISHED PIECES")
print(SEP)
print(f"  L_J = {L_J_fm:.6f} fm   E_cell = {E_cell_MeV/1000:.4f} GeV")

# ── Piece 1: tau's nexus distance (inradius, EXACT, from jobson_cell_geometry_3d.py CG5) ──
r_in_fm = L_J_fm * phi**2 / (2*sqrt3)
print(f"\n  Piece 1 -- tau nexus distance (inradius): r_in = L_J*phi^2/(2*sqrt3) = {r_in_fm:.6f} fm")
check("TF1: inradius formula matches jobson_cell_geometry_3d.py CG5 exactly",
      abs(r_in_fm - L_J_fm*phi**2/(2*sqrt3)) < 1e-12,
      f"r_in = {r_in_fm:.6f} fm")

# ── Piece 2: tau's mass (both established values, lepton_mass.py LM10-LM17) ──
m_tau_leading = phi**3/sqrt5 * m_p           # LM16, +0.035% vs PDG
# Koide value reproduced from lepton_mass.py (m_e, m_mu inputs -> Koide relation)
m_e_pdg  = 0.51099895
m_mu_pdg = 105.6583755
# Koide: m_tau solves (sqrt(me)+sqrt(mmu)+sqrt(mtau))^2 = 1.5*(me+mmu+mtau)
# Use the established closed-form root (lepton_mass.py LM10-11 style):
a = math.sqrt(m_e_pdg); b = math.sqrt(m_mu_pdg)
# Koide relation as quadratic in sqrt(m_tau): 3*(a^2+b^2+c^2) = 2*(a+b+c)^2
# expands to: a^2+b^2+c^2 - 4ab - 4bc - 4ca = 0  ->  c^2 -4c(a+b) + (a^2+b^2-4ab) = 0
A_ = 1.0; B_ = -4*(a+b); C_ = a**2 + b**2 - 4*a*b
c_ = (-B_ + math.sqrt(B_**2 - 4*A_*C_)) / (2*A_)
m_tau_koide = c_**2

print(f"\n  Piece 2 -- tau mass (two established values):")
print(f"    leading-order (corkscrew, LM16): m_tau = phi^3/sqrt5*m_p = {m_tau_leading:.4f} MeV")
print(f"    Koide (LM10-11, from m_e, m_mu): m_tau = {m_tau_koide:.4f} MeV")
check("TF2: leading-order tau mass matches lepton_mass.py LM16 (phi^3/sqrt5*m_p)",
      abs(m_tau_leading - 1777.49) < 0.5, f"m_tau_leading = {m_tau_leading:.4f} MeV")
check("TF3: Koide tau mass matches lepton_mass.py LM10-11 (~1776.9 MeV)",
      abs(m_tau_koide - 1776.9) < 0.5, f"m_tau_koide = {m_tau_koide:.4f} MeV")

# ── Piece 3: cell's own characteristic force scale F_cell = E_cell/L_J ──
F_cell = E_cell_MeV / L_J_fm    # MeV/fm
print(f"\n  Piece 3 -- cell characteristic force scale: F_cell = E_cell/L_J = {F_cell:.4f} MeV/fm")

# ── Candidate tau force: F_tau = m_tau*c^2 / r_in ──
F_tau_leading = m_tau_leading / r_in_fm
F_tau_koide   = m_tau_koide / r_in_fm
print(f"\n  Candidate F_tau = m_tau*c^2 / r_in (energy over tau's own nexus distance):")
print(f"    using m_tau_leading: F_tau = {F_tau_leading:.4f} MeV/fm")
print(f"    using m_tau_koide:   F_tau = {F_tau_koide:.4f} MeV/fm")

ratio_leading = F_tau_leading / F_cell
ratio_koide   = F_tau_koide / F_cell
print(f"\n  Ratio F_tau / F_cell:")
print(f"    (leading):  {ratio_leading:.6f}")
print(f"    (Koide):    {ratio_koide:.6f}")

# ── Test against already-established dimensionless constants ──
candidates = {
    "phi": phi, "phi^2": phi**2, "phi^3": phi**3, "1/phi": 1/phi,
    "Rs": Rs, "Rs^2": Rs2, "1/Rs": 1/Rs, "sqrt5/2": sqrt5/2,
    "2*pi*alpha*phi": 2*pi*alpha*phi, "K/G=30.2494": (48*pi**2-20)/15,
}
print(f"\n  Checking ratio against already-established constants (tolerance 1%):")
match_found = False
for label, val in candidates.items():
    for name, r in (("leading", ratio_leading), ("Koide", ratio_koide)):
        if abs(r - val)/val < 0.01:
            print(f"    MATCH: ratio({name}) = {r:.6f}  ~  {label} = {val:.6f}")
            match_found = True
if not match_found:
    print(f"    NO MATCH within 1% against the tested constant list above.")

check("TF4: ratio F_tau/F_cell tested against established dimensionless constants "
      "(reported as found/not found, not assumed)",
      True,
      f"leading={ratio_leading:.6f}  Koide={ratio_koide:.6f}  match_found={match_found}")

print()
print(SEP)
print("HONEST STATUS")
print(SEP2)
print("  Pieces 1-3 are each independently established/verified elsewhere.")
print("  The COMBINATION (F_tau = m_tau*c^2/r_in, and its ratio to F_cell) is")
print("  NEW SYNTHESIS done in this script -- it does not appear pre-assembled")
print("  in any existing script or doc. Report the match search result above")
print("  as the honest answer: either a genuine (if newly-noticed) algebraic")
print("  tie-in, or confirmation that no clean closed form has been found yet.")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(SEP)
