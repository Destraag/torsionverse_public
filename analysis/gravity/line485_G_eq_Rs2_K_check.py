"""
line485_G_eq_Rs2_K_check.py
=============================
doc_orbit_pressure.txt line 485 (Lense-Thirring/Einstein-Cartan note)
originally stated "Einstein-Cartan (torsion gravity) = our framework with
shear wave (G = Rs^2*K)". A 2026-09-05 session pass wrongly "corrected" this
to "G = Rs^2*K is only a leading-order approximation of the precise
K/G=30.25 relation (doc_torsion T3.2)", treating the two formulas as
zeroth-order vs. precise versions of the SAME relation. Re-examination
(prompted by the question "doc_orbit_pressure came after doc_torsion -- is
torsion's formula really more precise?") shows that framing was wrong.

doc_torsion.txt Section 3.3 ("NOTE ON K -- TWO REGIMES, NOT ONE VALUE",
2026-09-04) and analysis/medium/fluid_vs_jammed_K_resolution.py (FJ1-FJ6,
already committed) establish that this medium has TWO distinct, both-valid
bulk moduli, not one true value and one approximation:
  K_fluid  = 1/eps_0                        (static/Coulomb, unjammed state)
  K_jammed = rho*(v_p^2 - 4/3*v_s^2)        (propagating-wave, jammed state)
  G        = rho*v_s^2                       (shear modulus -- only exists
                                               jammed; G=0 in the fluid state)
These checks establish which K pairs correctly with G in the Einstein-Cartan/
shear-wave note, and whether either relation supplies post-Newtonian content.

Reference: doc_orbit_pressure.txt Lense-Thirring/Einstein-Cartan note;
  doc_torsion.txt Section 3.1/3.3 (T3.1-T3.2); fluid_vs_jammed_K_resolution.py;
  judgment_calls_resolved.txt EC-1
"""
import math

pi = math.pi
eps_0 = 8.8541878128e-12          # F/m, CODATA
mu_0 = 1.25663706212e-6           # torsion-medium density identity
c = 2.99792458e8                  # m/s, exact by SI definition
Rs = math.sqrt(5) / (4 * pi)
Rs2 = Rs ** 2
rho = mu_0

SEP = "=" * 72
results = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status))
    print(f"  [{'PASS' if cond else '*** FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


K_fluid = 1.0 / eps_0
K_jammed = rho * (c ** 2 - (4.0 / 3.0) * Rs2 * c ** 2)
G = rho * Rs2 * c ** 2

print(SEP)
print("Check 1: is 'G = Rs^2*K_fluid' an approximation, or exact?")
print(SEP)
print(f"  K_fluid  = 1/eps_0                    = {K_fluid:.6e} Pa")
print(f"  K_jammed = rho*(c^2 - 4/3*Rs^2*c^2)   = {K_jammed:.6e} Pa")
print(f"  G        = rho*Rs^2*c^2               = {G:.6e} Pa")
print()

# LC1: G = Rs^2*K_fluid is EXACT, not approximate -- because K_fluid=1/eps_0
# is defined via c=sqrt(K_fluid/rho), i.e. K_fluid=rho*c^2 identically, and
# G=rho*Rs^2*c^2=Rs^2*(rho*c^2)=Rs^2*K_fluid follows by trivial substitution.
check("LC1: G = Rs^2*K_fluid holds to float precision (EXACT identity, "
      "not an approximation)",
      abs(Rs2 * K_fluid - G) / G < 1e-9,
      f"Rs^2*K_fluid = {Rs2*K_fluid:.6e} Pa vs G = {G:.6e} Pa, "
      f"ratio = {(Rs2*K_fluid)/G:.12f}")

# LC2: K_jammed/G = 30.25 (doc_torsion T3.2), ALSO exact -- for the OTHER K.
KG_jammed = K_jammed / G
check("LC2: K_jammed/G = 30.25 also holds exactly (doc_torsion T3.2) -- for "
      "the JAMMED K, not the fluid K",
      abs(KG_jammed - 30.2494) < 1e-3,
      f"K_jammed/G = {KG_jammed:.4f}")

# LC3: K_fluid and K_jammed are genuinely different numbers (the already-
# established 4.22% gap from fluid_vs_jammed_K_resolution.py FJ2), so LC1
# and LC2 are two DIFFERENT relations, not the same relation at two
# precisions.
gap = 1 - K_jammed / K_fluid
check("LC3: K_fluid != K_jammed (differ by the established ~4.22% gap) -- "
      "confirms LC1 and LC2 are two distinct exact relations, not one "
      "relation stated imprecisely",
      abs(gap - 4.0/3.0*Rs2) < 1e-9,
      f"1 - K_jammed/K_fluid = {gap*100:.3f}%, 4/3*Rs^2 = {4/3*Rs2*100:.3f}%")

print()
print("  CORRECTED CONCLUSION: 'G = Rs^2*K' is NOT a zeroth-order approximation")
print("  of 'K/G=30.25' -- both are exact, but for two DIFFERENT, deliberately")
print("  distinct bulk moduli (doc_torsion's own fluid/jammed resolution, not")
print("  a precision issue). The 2026-09-05 same-session edit that called")
print("  G=Rs^2*K 'stale' and 'corrected' it to the 30.25 relation was WRONG:")
print("  it treated two non-competing quantities as one quantity measured at")
print("  two precisions. Reverted.")

print()
print(SEP)
print("Check 2: WHICH pairing is the physically apt one for a SHEAR-WAVE /")
print("Einstein-Cartan note specifically (not a precision question -- a")
print("regime-matching question)?")
print(SEP)
print("""  G (shear modulus) exists ONLY in the jammed state -- G=0 in the fluid
  ground state (doc_torsion Section 3.1: the medium "flows freely, G -> 0"
  below the shear threshold). K_fluid is specifically the UNJAMMED, static/
  Coulomb-regime modulus (a non-propagating point charge never shears the
  medium). So "G = Rs^2*K_fluid" is an exact identity, but it PAIRS a
  jammed-only quantity (G) with a fluid-only quantity (K_fluid) that never
  co-occur in the same physical state of the medium -- it is true, but not
  a same-regime physical relation.

  K_jammed and G, by contrast, are BOTH jammed/wave-state quantities --
  doc_torsion.txt Section 3.3 lists them together in the same "jammed/wave
  state" block, both solved from the same v_p=c, v_s=Rs*c inputs. K_jammed/G
  =30.25 is the same-regime relation: it describes one coherent mechanical
  state (the medium as actively sheared by a passing wave), matching the
  note's own subject (Einstein-Cartan TORSION, a dynamical/wave phenomenon,
  not a static point charge).

  CONCLUSION: for a shear-wave/Einstein-Cartan context specifically, the
  K_jammed-based K/G=30.25 relation is the more physically apt one to cite
  -- not because G=Rs^2*K_fluid is imprecise (it is exact) but because it
  cross-pairs quantities from two different, non-co-occurring regimes.""")

print()
print(SEP)
print("Check 3: does EITHER relation supply POST-NEWTONIAN/dynamical content?")
print(SEP)
print("""  Both K/G relations (however correctly paired) are STATIC ratios between
  two LINEAR elastic constants -- exactly the same kind of relation as the
  bulk/shear modulus ratio of an ordinary Hookean solid (steel, rock). They
  describe the medium's response to INFINITESIMAL, LINEAR strain only.

  A post-Newtonian derivation (Mercury precession, Tier 4) requires
  EXPANDING the medium's equations of motion to higher (nonlinear) order
  in the field strength/velocity -- e.g. third-order elastic constants, or
  a relativistic stress-energy tensor and field equations for the medium,
  analogous to the Parameterized Post-Newtonian (PPN) formalism's higher-
  order metric terms in GR. Nothing of that kind is present in either K/G
  ratio, or anywhere else already committed for this medium: no nonlinear
  stress-strain terms, no medium Lagrangian, no field equations beyond the
  LINEAR pressure_isobar_orbit.py EOM (F=ma, Newtonian order only, PO1-PO9).

  CONCLUSION: K/G (regardless of which K, however precisely stated or
  correctly regime-paired) is necessary but nowhere near SUFFICIENT to
  derive post-Newtonian corrections. It is a single static number, not a
  dynamical field theory. This part of the original finding is UNCHANGED
  by the LC1-LC3 correction above.""")

print()
print(SEP)
print("SUMMARY")
print(SEP)
n_pass = sum(1 for _, s in results if s == "PASS")
for name, status in results:
    print(f"  [{status}] {name}")
print(f"\n  Total: {len(results)}  PASS: {n_pass}  FAIL: {len(results) - n_pass}")
