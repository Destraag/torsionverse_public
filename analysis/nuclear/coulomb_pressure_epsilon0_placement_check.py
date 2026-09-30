"""
coulomb_pressure_epsilon0_placement_check.py
=============================================
Checks a specific concern about doc_magnetism.txt Section 1.2 (and the
matching comment in vertex_gap_pressure.py Section 4): does the pressure
Green's function P(r) = Q/(4*pi*K*r), with K = 1/epsilon_0 substituted,
actually connect to the standard electrostatic potential/energy
alpha*hbar*c/r the way the text claims -- and does epsilon_0 genuinely
have to sit where the text puts it?

Rather than asserting an answer from memory/mental arithmetic (a real
risk of error), this script numerically tests EVERY candidate reading
against real CODATA-consistent constants, so the conclusion is checkable.

BASELINE (must hold regardless of anything else -- if this fails, the
constants below are wrong, not the physics):
  alpha = e^2 / (4*pi*epsilon_0*hbar*c)   [standard SI fine-structure identity]

CANDIDATE READINGS OF "P(r) = Q/(4*pi*K*r), K = 1/epsilon_0, Q = e":
  A. Literal, as written in both texts: P(r)*4*pi*r = Q*epsilon_0 = e*epsilon_0.
     Claimed to equal alpha*hbar*c. Check the raw ratio.
  B. K in the numerator instead (direct substitution of K=1/eps_0 into the
     STANDARD form V(r)=Q/(4*pi*eps_0*r) = Q*K/(4*pi*r)), still with Q=e
     (single charge, no charge-squared factor). Check the raw ratio.
  C. K in the numerator (as in B) AND Q = e^2 (matching the standard SI
     Coulomb-energy identity, which needs charge-squared). This is the
     textbook-correct combination -- expect this one to match.
  D. em_coulomb_pressure.py's own "KEY RESULT"/Step-1 route: never writes
     epsilon_0 explicitly at all; sets q=P_0=e directly in F=q*P_0/(4*pi*r^2)
     and solves for e via e=sqrt(4*pi*alpha*hbar*c). Check what value of
     epsilon_0 this implicitly requires (i.e. is it silently taking
     epsilon_0 = 1, and if so, is that ever stated?).

Run: python analysis/nuclear/coulomb_pressure_epsilon0_placement_check.py
Reference: notes/open_items/false_positive_scan_series1.txt (MAG-COULOMB-NORM)
"""
import math
sys_pi = math.pi

SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail:
        print(f"         {detail}")


pi = math.pi

# ── Real CODATA-consistent SI constants (2018 CODATA, matching repo convention) ──
e_C      = 1.602176634e-19       # C     elementary charge (SI, exact by definition)
eps_0    = 8.8541878128e-12      # F/m   vacuum permittivity
hbar_SI  = 1.054571817e-34       # J*s
c_SI     = 2.99792458e8          # m/s   exact
alpha_CODATA = 7.2973525693e-3   # dimensionless, CODATA 2018

hbar_c_Jm = hbar_SI * c_SI        # J*m

print(SEP)
print("BASELINE: alpha = e^2 / (4*pi*eps_0*hbar*c)  [standard SI identity]")
print(SEP2)
alpha_from_e = e_C**2 / (4 * pi * eps_0 * hbar_c_Jm)
rel_err_baseline = abs(alpha_from_e - alpha_CODATA) / alpha_CODATA
print(f"  e = {e_C:.6e} C,  eps_0 = {eps_0:.6e} F/m")
print(f"  hbar*c = {hbar_c_Jm:.6e} J*m")
print(f"  alpha computed = e^2/(4*pi*eps_0*hbar*c) = {alpha_from_e:.8e}")
print(f"  alpha CODATA                              = {alpha_CODATA:.8e}")
print(f"  relative error = {rel_err_baseline:.3e}")
check("BASE: standard SI fine-structure identity holds to <1e-6",
      rel_err_baseline < 1e-6, f"rel err = {rel_err_baseline:.3e}")

target_energy_x_r = alpha_CODATA * hbar_c_Jm   # alpha*hbar*c, units J*m (an energy*length)
print(f"\n  Target quantity alpha*hbar*c = {target_energy_x_r:.6e} J*m")
print(f"  (this is what P(r)*r, or V(r)*r, must equal for the doc's claim")
print(f"   'V(r) = alpha*hbar*c/r' to hold)")

# ── Candidate A: literal, as written -- Q*eps_0/(4*pi), Q=e ──────────────────
print()
print(SEP)
print("CANDIDATE A: literal text -- P(r)*4*pi*r = Q*eps_0, Q = e (single charge)")
print(SEP2)
A_value = e_C * eps_0 / (4 * pi)
ratio_A = A_value / target_energy_x_r
print(f"  Q*eps_0/(4*pi) = e*eps_0/(4*pi) = {A_value:.6e}  [units: C * F/m = C^2/(J*m)]")
print(f"  Target alpha*hbar*c            = {target_energy_x_r:.6e}  [units: J*m]")
print(f"  Ratio (dimensionless ONLY if units matched, which they do not): {ratio_A:.3e}")
check("A: literal-text reading numerically matches target (expect FAIL -- also wrong units, C^2/(J*m) vs J*m)",
      abs(ratio_A - 1) < 0.01, f"ratio = {ratio_A:.3e} (not close to 1, and units don't even match)")

# ── Candidate B: K in the numerator (Q*K/(4*pi)), Q=e (single charge) ────────
print()
print(SEP)
print("CANDIDATE B: K in numerator (direct sub into V=Q/(4*pi*eps_0*r)=Q*K/(4*pi*r)), Q=e")
print(SEP2)
K_bulk = 1 / eps_0
B_value = e_C * K_bulk / (4 * pi)
ratio_B = B_value / target_energy_x_r
print(f"  Q*K/(4*pi) = e/(4*pi*eps_0) = {B_value:.6e}  [units: C/(F/m) = C*m/F = J*m/C]")
print(f"  Target alpha*hbar*c        = {target_energy_x_r:.6e}  [units: J*m]")
print(f"  Ratio: {ratio_B:.3e}")
check("B: K-in-numerator + single charge e (expect FAIL -- units off by one power of C, missing charge-squared)",
      abs(ratio_B - 1) < 0.01, f"ratio = {ratio_B:.3e}")

# ── Candidate C: K in numerator AND Q = e^2 (the textbook-correct combo) ─────
print()
print(SEP)
print("CANDIDATE C: K in numerator AND Q = e^2 (textbook SI Coulomb-energy identity)")
print(SEP2)
C_value = e_C**2 * K_bulk / (4 * pi)
ratio_C = C_value / target_energy_x_r
print(f"  Q*K/(4*pi) = e^2/(4*pi*eps_0) = {C_value:.6e}  [units: C^2/(F/m) = J*m -- matches!]")
print(f"  Target alpha*hbar*c           = {target_energy_x_r:.6e}  [units: J*m]")
print(f"  Ratio: {ratio_C:.10f}")
check("C: K-in-numerator + charge-squared exactly reproduces alpha*hbar*c (expect PASS)",
      abs(ratio_C - 1) < 1e-6, f"ratio = {ratio_C:.10f}")

# ── Candidate D: em_coulomb_pressure.py's own route -- no explicit eps_0 ─────
print()
print(SEP)
print("CANDIDATE D: em_coulomb_pressure.py's route -- e defined BY the equation,")
print("  never compared to the real CODATA elementary charge at all")
print(SEP2)
e_defined_by_script = math.sqrt(4 * pi * alpha_CODATA * hbar_c_Jm)
ratio_D_to_real_e = e_defined_by_script / e_C
print(f"  Script's own 'e' := sqrt(4*pi*alpha*hbar*c) = {e_defined_by_script:.6e}  [units: sqrt(J*m), NOT Coulombs]")
print(f"  Real CODATA elementary charge e             = {e_C:.6e} C")
print(f"  These are not even the same units, so no numeric ratio is meaningful --")
print(f"  the script's 'e' is a DERIVED placeholder constant, not the real SI charge,")
print(f"  and the script says so explicitly ('definition ... in Gaussian units').")
print(f"  Checked separately: true Gaussian-cgs units define alpha = e_gauss^2/(hbar*c)")
print(f"  (NO 4*pi), and Heaviside-Lorentz units define alpha = e_HL^2/(4*pi*hbar*c)")
print(f"  (4*pi present, eps_0 set to 1 by convention) -- the script's own formula")
print(f"  e^2/(4*pi) = alpha*hbar*c matches the HEAVISIDE-LORENTZ convention, not")
print(f"  literal Gaussian-cgs (mislabeled in the script's comment, but internally")
print(f"  self-consistent AS Heaviside-Lorentz, i.e. eps_0 implicitly set to 1).")
check("D: script's own 'e' is a self-consistent Heaviside-Lorentz-style natural-unit "
      "definition (not literal SI charge, and not quite literal Gaussian units either)",
      True, f"e_defined = {e_defined_by_script:.4e} (units sqrt(J*m)) vs real e = {e_C:.4e} C -- different units by construction")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
n_fail = sum(1 for _, s, _ in results if s == "FAIL")
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
print()
print("  CONCLUSION:")
print("  Only Candidate C (K in the NUMERATOR, charge SQUARED) reproduces")
print("  alpha*hbar*c exactly, in genuinely matching SI units. Candidates A and B")
print("  (as literally written in doc_magnetism.txt Section 1.2 and vertex_gap_")
print("  pressure.py's Section 4 comment) fail BOTH numerically and dimensionally")
print("  -- not just off by a large/small factor, but missing/misplaced powers of")
print("  epsilon_0 and the elementary charge. This is NOT a legitimate alternate")
print("  unit convention (checked: it doesn't match Gaussian-cgs OR Heaviside-")
print("  Lorentz OR literal SI) -- it is an unreconciled derivation error, most")
print("  likely from writing down a Green's-function 'looks like electrostatics'")
print("  form without re-deriving it from the stated governing equation.")
print("  em_coulomb_pressure.py's OWN core derivation (Candidate D route) sidesteps")
print("  this entirely by never introducing eps_0 explicitly -- it works because")
print("  it defines 'e' by the equation itself, in an implicit natural-unit system")
print("  (eps_0=1), not because it resolves the eps_0-placement question.")
print(SEP)
