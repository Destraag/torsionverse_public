"""
kim_eraser_creation_nexus_test.py
==================================
Tests H6, the "creation-nexus rule" (notes/open_items/
kim_eraser_explanation_proposal.txt, STATUS UPDATE 2026-10-04): does
reusing double_slit_field_interference.py's ALREADY-VERIFIED corpuscle-
definite/field-spans-both-slits mechanism a SECOND time -- at the SPDC
crystal and the idler's own beamsplitters -- reproduce Kim et al. 2000's
own reported R01/R02/R03/R04 structure, with no G32 chain, no winding,
no new ontology?

THE ONE ASSUMPTION BEING TESTED (flagged honestly, NOT claimed derived
from Jobson-cell medium mechanics -- this is the same "imported from
standard SPDC momentum conservation" gap doc_qm.txt Section 5.2/6
already flags): spontaneous down-conversion couples to the FIELD at the
crystal, not to the corpuscle's already-committed slit. Per Section
4.1's own mechanism (same one double_slit_field_interference.py checks),
that field already coherently spans BOTH slits at the crystal, regardless
of which slit the real corpuscle used. If down-conversion inherits this
two-path coherence -- phase SUM locked to the pump by energy/momentum
conservation, phase DIFFERENCE left free to carry the slit-A-vs-slit-B
information -- then the signal's conditional field amplitude at D0,
GIVEN which SLIT the idler is correlated with, carries the SAME
x-dependent relative phase an ordinary double slit would produce. This
one assumption is the entire content of H6; everything downstream of it
(beamsplitter unitarity, Huygens recombination) is machinery already
established and verified elsewhere in this framework.

SIGNAL/IDLER MODEL (matches the paper's own apparatus, Section 2 of
kim_eraser_explanation_proposal.txt):
  PSI_A(x), PSI_B(x) = SIGNAL's conditional field amplitude at D0 position
  x, given the idler is slit-A-correlated / slit-B-correlated. Same
  envelope for both (already-verified normal-incidence result,
  double_slit_field_interference.py DS4a: each slit's own far-field
  envelope is centered at y=0 regardless of lateral slit position), equal
  and opposite relative phase beta(x) = pi*d*x/(lambda*f) (identical
  parameter, identical formula, to that script's "theta_term_d").
    PSI_A(x) = sqrt(env(x)) * exp(+i*beta(x))
    PSI_B(x) = sqrt(env(x)) * exp(-i*beta(x))
  IDLER ROUTING (apparatus, not assumed -- this is the real Kim et al.
  beam path): each of PSI_A, PSI_B first meets its own 50/50 beamsplitter
  (BSa, BSb). Transmitted portion -> directly to D3 (from A) / D4 (from
  B): which-path information PRESERVED. Reflected portion -> mirror (Ma/
  Mb) -> shared final beamsplitter BSc, which recombines the two
  reflected portions: which-path information ERASED. BSc's own
  transmission/reflection unitary transformation is what the paper
  itself cites for the +/- sign between its D1 and D2 formulas -- this
  script does not add anything new there, it is standard 50/50
  beamsplitter unitarity: out_D1=(a+b)/sqrt(2), out_D2=(a-b)/sqrt(2).

Checks:
  CN1  R01(x) = |out_D1|^2 shows high-visibility cos^2-type fringes
  CN2  R02(x) = |out_D2|^2 shows high-visibility sin^2-type fringes,
       exactly in antiphase with R01 (peaks of one are troughs of other)
  CN3  R03(x) = R04(x) exactly (same envelope, no fringes -- which-path
       preserved on both direct ports, by the same symmetry DS4a already
       established)
  CN4  Unitarity/probability conservation: R01+R02+R03+R04 = total
       incoming probability |PSI_A|^2+|PSI_B|^2 at every x (nothing lost
       or gained by routing through beamsplitters)
  CN5  R01+R02 alone (erasure port, NOT sorted into D1 vs D2 individually)
       washes out to a smooth envelope -- fringes require sorting all the
       way down to a SPECIFIC erasure-port detector, not just "erasure
       happened"
  CN6  CROSS-VALIDATION against double_slit_field_interference.py: R01(x)
       here is IDENTICAL in form to that script's own case-A (coherent)
       result -- same mechanism, run a second time, not new machinery
  CN7  R01(x), R02(x) match the paper's own cited functional form,
       cos^2(x*pi*d/(lambda*f)) and sin^2(x*pi*d/(lambda*f))

WHAT THIS DOES NOT CLOSE (flagged explicitly, not overclaimed): the
phase-sum-lock rule itself is still imported from standard SPDC
momentum/energy conservation, not derived from Jobson-cell mechanics --
this narrows, but does not fully close, doc_qm.txt Section 6's existing
"NOT YET DERIVED" item on the down-conversion process itself.

Run: python analysis/quantum/kim_eraser_creation_nexus_test.py
Reference: notes/open_items/kim_eraser_explanation_proposal.txt H6;
           analysis/quantum/double_slit_field_interference.py (reused
           mechanism, same parameter conventions)
"""
import math
import numpy as np

SEP  = "=" * 72
SEP2 = "-" * 72
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail:
        print(f"         {detail}")


pi = math.pi

# ── Same electron-beam / slit-geometry parameters as double_slit_field_
# interference.py, reused verbatim for direct comparability (CN6) ───────────
hbar = 1.054571817e-34
m_e  = 9.1093837015e-31
v    = 1.0e6
lambda_dB = 2 * pi * hbar / (m_e * v)
d = 1.0e-6
a = 2.0e-8
L = 1.0

print(SEP)
print("kim_eraser_creation_nexus_test.py -- H6 creation-nexus rule check")
print(SEP)
print(f"  Reusing double_slit_field_interference.py's own parameters:")
print(f"  lambda_dB = {lambda_dB*1e12:.4f} pm, d = {d*1e6:.2f} um, a = {a*1e9:.0f} nm, L = {L} m")
print()


def sinc2(x):
    return np.where(np.abs(x) < 1e-12, 1.0, (np.sin(x) / np.where(x == 0, 1, x)) ** 2)


y = np.linspace(-4e-2, 4e-2, 800001)
theta_term_a = pi * a * y / (lambda_dB * L)
theta_term_d = pi * d * y / (lambda_dB * L)          # this IS beta(x) from the docstring
env = sinc2(theta_term_a)                             # single-slit envelope, shared by both paths

# ── PSI_A, PSI_B: signal's conditional field amplitude at D0, given the
# idler is slit-A-correlated / slit-B-correlated (the H6 assumption) ───────
PSI_A = np.sqrt(env) * np.exp(1j * theta_term_d)
PSI_B = np.sqrt(env) * np.exp(-1j * theta_term_d)

# ── BSa/BSb: 50/50 split, transmitted -> D3/D4 (direct), reflected -> BSc ──
to_D3_direct = PSI_A / math.sqrt(2)
to_D4_direct = PSI_B / math.sqrt(2)
to_BSc_from_A = PSI_A / math.sqrt(2)
to_BSc_from_B = PSI_B / math.sqrt(2)

# ── BSc: standard 50/50 unitary recombination (paper's own cited +/- sign) ─
out_D1 = (to_BSc_from_A + to_BSc_from_B) / math.sqrt(2)
out_D2 = (to_BSc_from_A - to_BSc_from_B) / math.sqrt(2)

R01 = np.abs(out_D1) ** 2
R02 = np.abs(out_D2) ** 2
R03 = np.abs(to_D3_direct) ** 2
R04 = np.abs(to_D4_direct) ** 2

fringe_spacing_nominal = lambda_dB * L / d
narrow = np.abs(y) < 3 * fringe_spacing_nominal


def visibility(I):
    Imax, Imin = I.max(), I.min()
    return (Imax - Imin) / (Imax + Imin) if (Imax + Imin) > 0 else 0.0


V_01 = visibility(R01[narrow])
V_02 = visibility(R02[narrow])

print(SEP2)
print("CN1/CN2 -- erasure ports (D1, D2): fringe visibility and antiphase")
print(SEP2)
print(f"  V(R01) = {V_01:.4f}   V(R02) = {V_02:.4f}")
check("CN1: R01(x) shows high-visibility fringes (V > 0.95)", V_01 > 0.95, f"V_01={V_01:.4f}")
check("CN2: R02(x) shows high-visibility fringes (V > 0.95)", V_02 > 0.95, f"V_02={V_02:.4f}")

antiphase_err = np.max(np.abs((R01[narrow] + R02[narrow]) - env[narrow]))
rel_antiphase_err = antiphase_err / env[narrow].max()
check("CN2b: R01 and R02 are in exact antiphase (R01+R02 == env(x) itself, "
      "cos^2+sin^2=1 identically -- no residual fringe-frequency oscillation "
      "beyond the slowly-varying envelope, which is NOT perfectly flat on "
      "its own and should not be compared to one)",
      rel_antiphase_err < 1e-9, f"max|R01+R02-env| / env_max = {rel_antiphase_err:.3e}")

print()
print(SEP2)
print("CN3 -- direct ports (D3, D4): identical envelopes, no fringes")
print(SEP2)
max_diff_R3R4 = np.max(np.abs(R03 - R04))
check("CN3: R03(x) == R04(x) exactly (same envelope, which-path preserved)",
      max_diff_R3R4 < 1e-12, f"max|R03-R04| = {max_diff_R3R4:.3e}")
V_03 = visibility(R03[narrow])
check("CN3b: R03(x) shows NO fringes (V < 0.01)", V_03 < 0.01, f"V_03={V_03:.6f}")

print()
print(SEP2)
print("CN4 -- unitarity: nothing lost/gained routing through beamsplitters")
print(SEP2)
total_out = R01 + R02 + R03 + R04
total_in  = np.abs(PSI_A) ** 2 + np.abs(PSI_B) ** 2
max_rel_err = np.max(np.abs(total_out - total_in) / np.where(total_in > 1e-15, total_in, 1))
check("CN4: R01+R02+R03+R04 == |PSI_A|^2+|PSI_B|^2 at every x (<1e-9 rel err)",
      max_rel_err < 1e-9, f"max relative error = {max_rel_err:.3e}")

print()
print(SEP2)
print("CN5 -- erasure port WITHOUT sorting into D1 vs D2: fringes vanish again")
print(SEP2)
V_erasure_unsorted = visibility((R01 + R02)[narrow])
check("CN5: R01+R02 (erasure happened, but D1-vs-D2 not distinguished) shows "
      "NO fringes (V < 0.01) -- sorting must go all the way to a specific "
      "erasure-port detector, matching the paper's own coincidence requirement",
      V_erasure_unsorted < 0.01, f"V(R01+R02) = {V_erasure_unsorted:.6f}")

print()
print(SEP2)
print("CN6 -- cross-validation vs double_slit_field_interference.py's case A")
print(SEP2)
I_A_reference = env * np.cos(theta_term_d) ** 2     # that script's own case-A formula, reproduced here
shape_diff = np.max(np.abs(R01[narrow] / R01[narrow].max() - I_A_reference[narrow] / I_A_reference[narrow].max()))
check("CN6: R01(x) has the SAME normalized shape as double_slit_field_"
      "interference.py's own coherent case A (shape_diff < 1e-9) -- same "
      "mechanism, run a second time, no new machinery",
      shape_diff < 1e-9, f"shape_diff = {shape_diff:.3e}")

print()
print(SEP2)
print("CN7 -- match paper's own cited functional form")
print(SEP2)
cos2_form = env * np.cos(theta_term_d) ** 2
sin2_form = env * np.sin(theta_term_d) ** 2
err_cos2 = np.max(np.abs(R01 - cos2_form))
err_sin2 = np.max(np.abs(R02 - sin2_form))
check("CN7a: R01(x) == env(x)*cos^2(pi*d*x/(lambda*f)) exactly",
      err_cos2 < 1e-12, f"max abs error = {err_cos2:.3e}")
check("CN7b: R02(x) == env(x)*sin^2(pi*d*x/(lambda*f)) exactly",
      err_sin2 < 1e-12, f"max abs error = {err_sin2:.3e}")

# ── Summary ──────────────────────────────────────────────────────────────
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
n_fail = len(results) - n_pass
print(f"SUMMARY: {len(results)} checks, {n_pass} PASS, {n_fail} FAIL")
print()
if n_fail == 0:
    print("  ALL CHECKS PASSED. The creation-nexus/phase-lock rule (H6), GIVEN")
    print("  as an assumption, reproduces Kim et al.'s own R01/R02/R03/R04")
    print("  structure (cos^2 / sin^2 / matched single-peak envelopes / exact")
    print("  unitarity) using ONLY: (a) the already-verified double-slit field")
    print("  mechanism (double_slit_field_interference.py), applied a second")
    print("  time at the idler's own beamsplitters, and (b) standard 50/50")
    print("  beamsplitter unitarity, which the paper itself cites for the same")
    print("  +/- sign. NO G32 chain, NO winding, NO new ontology was used.")
print()
print("  WHAT THIS DOES NOT CLOSE: the phase-sum-lock rule itself (PSI_A,")
print("  PSI_B inheriting the double slit's x-dependent relative phase)")
print("  remains IMPORTED from standard SPDC energy/momentum conservation,")
print("  not derived from Jobson-cell medium mechanics. This narrows, but")
print("  does not fully close, doc_qm.txt Section 6's existing NOT-YET-")
print("  DERIVED item on the down-conversion process itself.")
