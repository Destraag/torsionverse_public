"""
higgs_self_energy_isotropic_derivation_attempt.py
==================================================
GENUINE DERIVATION ATTEMPT (not an audit): can the Higgs mass's (1+alpha/pi)
correction be built from torsion-medium-native geometry, rather than imported
from standard QED? This is a follow-on to higgs_alpha_pi_import_audit.py,
which confirmed the import and ruled out 4 candidate substitutions (delta_n
family + bare alpha). This script tries NEW candidates motivated by the
mode's own ISOTROPY, rather than reusing the T_1g "contact correction" rule
(which is for a DIRECTIONAL mode, and already failed as Candidate 4).

HONEST FRAMING (stated up front, not after the fact):
  Two already-tried candidates PASSED numerically but were REJECTED as
  independently-derived: Vol(S^2)/CS_(1,2)=1/pi exactly (higgs_alpha_pi_
  import_audit.py Part D) is algebraically exact but reverse-engineered --
  no physical argument said self-energy = Vol(S^2)/CS(p,q) BEFORE checking
  it numerically worked. This script tries to do better: state WHY each
  candidate should be the right formula BEFORE computing it, using
  quantities this framework ALREADY uses for a DIFFERENT, independent
  purpose (not invented to hit 1/pi).

CANDIDATE SOURCE 1 -- Icosahedral discrete geometry (vertices=12, faces=20,
  edges=30, solid angle per face=pi/5 [gap1_icosahedral_hopf.py, established]):
  tested RIGOROUSLY here for the first time (previously only reasoned about
  qualitatively, not computed) whether any discrete sum/ratio of these
  counts reproduces 1/pi.

CANDIDATE SOURCE 2 -- Isotropic Green's-function normalization: this
  framework's OWN Claim 7 (Coulomb's law V=-alpha*hbar*c/r) is ALREADY
  derived from the 3D Poisson Green's function, which carries a physical
  1/(4*pi) normalization from the solid angle of a point source's field
  (Gauss's law) -- an isotropic (A_g-type) source, unlike the T_1g vertex's
  directional contact correction. Testing whether the SAME normalization,
  reused (not invented), explains the Higgs mass correction.

Run: python analysis/higgs/higgs_self_energy_isotropic_derivation_attempt.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

pi = math.pi
Rs = math.sqrt(5) / (4 * pi)

SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "FAIL"
    results.append((name, s, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


m_H_pdg22 = 125.25
m_H_unc = 0.17
m_H_target_ratio = 1 + alpha / pi  # what we're trying to reproduce, 1.002322...

print(SEP)
print("SOURCE 1: RIGOROUS TEST OF DISCRETE ICOSAHEDRAL GEOMETRY (not just reasoning)")
print(SEP2)
print()

V, E, F = 12, 30, 20
solid_angle_per_face = 4 * pi / F  # = pi/5, established (gap1_icosahedral_hopf.py)
solid_angle_per_vertex = 4 * pi / V  # = pi/3
print(f"  V={V}, E={E}, F={F}  (icosahedron, established geometry)")
print(f"  Solid angle per face = 4*pi/{F} = {solid_angle_per_face:.6f} = pi/5")
print(f"  Solid angle per vertex = 4*pi/{V} = {solid_angle_per_vertex:.6f} = pi/3")
print()

# A uniform sum over ALL faces/vertices is trivially 1 (complete tiling) --
# confirmed here explicitly rather than just asserted:
uniform_face_sum = F * (solid_angle_per_face / (4 * pi))
uniform_vertex_sum = V * (solid_angle_per_vertex / (4 * pi))
print(f"  CONFIRMED (not just reasoned): uniform sum over all faces = "
      f"{uniform_face_sum:.10f} (trivial, complete tiling)")
print(f"  CONFIRMED: uniform sum over all vertices = {uniform_vertex_sum:.10f} (trivial)")
print()

# Only genuinely non-trivial discrete ratios: single-unit quantities compared
# to pi itself (not to 4*pi, which is guaranteed trivial by completeness).
candidates_discrete = [
    ("1/F = 1/20",                    1 / F),
    ("1/V = 1/12",                    1 / V),
    ("1/E = 1/30",                    1 / E),
    ("solid_angle_per_face / 4 = (pi/5)/4",   solid_angle_per_face / 4),
    ("solid_angle_per_vertex / 3 = (pi/3)/3", solid_angle_per_vertex / 3),
    ("F/(V*E) = 20/360",              F / (V * E)),
    ("V/(F*E)",                       V / (F * E)),
]
print("  Testing discrete ratios against the target 1/pi = "
      f"{1/pi:.8f}:")
for name, val in candidates_discrete:
    print(f"    {name:32s} = {val:.8f}   (ratio to 1/pi: {val*pi:.4f})")

check("S1 no bare discrete icosahedral ratio (V,E,F alone) equals 1/pi",
      not any(abs(val - 1 / pi) < 1e-6 for _, val in candidates_discrete),
      "confirmed: none of the tested bare discrete ratios hit 1/pi")

print()
print("  RESULT: pi does not appear from bare vertex/edge/face COUNTS alone --")
print("  expected, since finite combinatorial ratios of integers are rational,")
print("  and pi is transcendental. Two of the ABOVE rows (pi/5)/4 and (pi/3)/3")
print("  DO contain a pi (since solid-angle-per-face/vertex already has one),")
print("  but dividing by an arbitrary integer (4, 3) to try to cancel it to")
print("  1/pi requires pi/5/4 = 1/pi <=> pi^2=20 (false) -- confirmed FAIL")
print("  numerically above, not just algebraically asserted.")

# =============================================================================
print()
print(SEP2)
print("SOURCE 2: ISOTROPIC GREEN'S-FUNCTION NORMALIZATION (from this")
print("framework's OWN Claim 7, reused not invented)")
print(SEP2)
print()
print("  doc_higgs.txt Claim 7: Coulomb's law V=-alpha*hbar*c/r is derived from")
print("  the torsion medium's 3D Poisson Green's function -- an ISOTROPIC point")
print("  source, carrying the standard 1/(4*pi) solid-angle normalization from")
print("  Gauss's law (G(r) = 1/(4*pi*r) for the 3D Poisson equation). This is")
print("  the ONE already-established isotropic (not directional-contact)")
print("  self-coupling normalization in this framework -- reused here, not")
print("  invented for this purpose.")
print()

candidate_isotropic = alpha / (4 * pi)
m_H_isotropic = E_cell_GeV * (1 + candidate_isotropic)
sigma_isotropic = abs(m_H_isotropic - m_H_pdg22) / m_H_unc
print(f"  CANDIDATE: alpha/(4*pi) [pure Green's-function isotropic normalization]")
print(f"    = {candidate_isotropic:.8f}   vs needed alpha/pi = {alpha/pi:.8f}")
print(f"    E_cell*(1+alpha/(4*pi)) = {m_H_isotropic:.4f} GeV  ({sigma_isotropic:.2f} sigma)")
print()

check("S2 bare alpha/(4*pi) does NOT match m_H as well as alpha/pi",
      sigma_isotropic > abs(E_cell_GeV * (1 + alpha / pi) - m_H_pdg22) / m_H_unc,
      f"{sigma_isotropic:.2f} sigma vs 0.95 sigma")

print()
print("  FAILS ALONE by a factor of exactly 4 (alpha/(4*pi) vs alpha/pi) --")
print("  needs a multiplicity of 4 to close. doc_higgs.txt Section 4 ALREADY")
print("  uses '4' for an independent reason: 'the number of real scalar")
print("  components in the SU(2) Higgs doublet' (lambda=(1-nu)/4). Testing")
print("  whether reusing THIS SAME 4 (not inventing a new one) bridges the gap:")
print()

candidate_4x_isotropic = 4 * alpha / (4 * pi)
m_H_4x = E_cell_GeV * (1 + candidate_4x_isotropic)
sigma_4x = abs(m_H_4x - m_H_pdg22) / m_H_unc
print(f"  CANDIDATE: 4 * alpha/(4*pi) = alpha/pi exactly (algebraic identity)")
print(f"    E_cell*(1+4*alpha/(4*pi)) = {m_H_4x:.4f} GeV  ({sigma_4x:.2f} sigma)")

check("S3 4*alpha/(4*pi) reduces ALGEBRAICALLY to alpha/pi (not a new number)",
      abs(candidate_4x_isotropic - alpha / pi) < 1e-15,
      f"{candidate_4x_isotropic:.15f} vs {alpha/pi:.15f}")

print()
print("  HONEST VERDICT ON S2/S3: this reproduces the target EXACTLY, but I")
print("  must flag the SAME weakness as the Vol(S^2)/CS_(1,2) case: '4 real")
print("  scalar components' (Section 4's OWN reason for its own '4', in a")
print("  QUARTIC self-coupling |H|^4 context) has NO independent physical link")
print("  to 'how many isotropic contact corrections sum for a MASS term' --")
print("  the match (4/(4*pi)=1/pi) is ALGEBRAICALLY guaranteed for ANY '4' I")
print("  might reach for, not evidence the mechanism is right. This is a")
print("  second instance of the reverse-engineering pattern already flagged,")
print("  not an independent confirmation of it -- finding a candidate whose")
print("  arithmetic works is not the same as finding one whose PHYSICS forces")
print("  that arithmetic. I am NOT claiming S2/S3 as a recovery.")

print()
print(SEP)
print("FINAL SUMMARY")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass}/{len(results)} checks PASS")
print()
print("  Discrete icosahedral vertex/edge/face counts (rigorously tested here,")
print("  not just reasoned about): confirmed to NOT produce 1/pi -- expected,")
print("  since finite integer ratios are rational and pi is transcendental.")
print("  The isotropic Green's-function normalization (1/(4*pi), reused from")
print("  this framework's own Claim 7) gets within a clean factor of 4 of the")
print("  target, and that same '4' already exists in this doc for an UNRELATED")
print("  reason (SU(2) doublet component count) -- but no independent physical")
print("  argument connects the two uses. STATUS: still no honestly-recovered")
print("  derivation. The genuinely untried avenue remains a real, from-scratch")
print("  Chern-Simons-style integral (analogous to gap3_chern_simons.py) with")
print("  its own physical motivation stated before computation -- not")
print("  attempted in this script, which only tested combinations of ALREADY-")
print("  established quantities.")
