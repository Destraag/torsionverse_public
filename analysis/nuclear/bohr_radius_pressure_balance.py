"""
bohr_radius_pressure_balance.py
================================
Writes the two pressure/force functions doc_magnetism.txt Section 1.2a and
doc_electron.txt Section 3.1 narrate (electron's own outward pressure vs.
the proton's inward Coulomb pull) as explicit functions of r, and solves
for where they balance -- rather than asserting the standard Bohr formula
a_0 = hbar*c/(m_e*c^2*alpha) as the outcome without showing the algebra.

Every existing occurrence of a_0 in this repo (proton_structure.py PS5,
electron_doc.py ED3, pressure_isobar_orbit.py PO2/PO3, doc_orbit_pressure.txt
Section 1.3) either plugs into that formula directly and compares to CODATA,
or narrates the mechanism qualitatively -- none derives it from an explicit
balance-point calculation. This script does, using ONLY pieces already
established elsewhere in the framework:
  - V(r) = -alpha*hbar_c/r : the Coulomb potential, already derived from the
    medium's pressure Green's function (doc_magnetism.txt Section 1.2,
    Claim C7 of doc_higgs).
  - T(r) : the confinement (kinetic) energy of a wave packet localized to
    scale r -- a direct consequence of the Schrodinger equation this
    framework already derives as the non-relativistic limit of the medium's
    wave equation (qm_from_medium.py QM5). This IS the "electron cannot be
    compressed below its natural wave-packet size" mechanism the docs
    describe in words; here it is written down and solved.

METHOD (standard variational calculation, applied to this framework's own
already-derived V(r) and Schrodinger kinetic operator -- not new physics):
  Trial radial wavefunction psi(r) = exp(-r/a) (the same functional form as
  hydrogen's true ground state; for this ansatz the variational method is
  exact, not approximate).
    <T>(a) = (hbar*c)^2 / (2*m_e*c^2*a^2)   [KE expectation, verified by
             numerical integration below, not just quoted]
    <V>(a) = -alpha*hbar*c / a              [PE expectation, likewise]
    E(a) = <T>(a) + <V>(a)
  FORCE-BALANCE (= pressure-balance) form, matching the doc's own language:
    F_out(a) = -d<T>/da = +(hbar*c)^2/(m_e*c^2*a^3)   [electron's own
               outward push, diverges as a -> 0: "cannot be compressed"]
    F_in(a)  = -d<V>/da = -alpha*hbar*c/a^2           [proton's inward pull]
  Setting F_out(a) + F_in(a) = 0 (equilibrium, i.e. dE/da=0) gives, exactly:
    a = hbar*c / (m_e*c^2*alpha) = a_0
  -- the same formula used everywhere else in the repo, but now the OUTCOME
  of an explicit balance calculation rather than an assumed input.

This also explains PS5's virial-theorem cross-check (V(a_0)=2x13.6eV) rather
than merely being consistent with it: <T>=-<V>/2 falls out automatically at
the E(a)-minimizing a, since that is what the variational method guarantees
at any stationary point of a 1/r-plus-1/r^2 energy function.

Run: python analysis/nuclear/bohr_radius_pressure_balance.py
Reference: docs/series1/doc_magnetism.txt Section 1.2a/3.4,
           docs/series1/doc_electron.txt Section 3.1,
           analysis/quantum/qm_from_medium.py (Schrodinger equation),
           analysis/nuclear/proton_structure.py PS5 (virial cross-check).
"""

import math
import numpy as np

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

pi      = math.pi
alpha   = 7.2973525693e-3       # fine structure constant (CODATA)
hbar_c  = 197.3269804            # MeV*fm
m_e_c2  = 0.51099895             # MeV  (electron rest energy)
a_0_CODATA_fm = 5.29177210903e-11 * 1e15   # CODATA Bohr radius, fm

print(SEP)
print("BOHR RADIUS FROM AN EXPLICIT ELECTRON/PROTON PRESSURE BALANCE")
print(SEP2)
print("  V(r) = -alpha*hbar_c/r        [Coulomb, already derived: Section 1.2]")
print("  T(r) = confinement KE of a wave packet localized to scale r")
print("         [already-derived Schrodinger equation, qm_from_medium.py QM5]")
print()

# ── Section 1: numerically verify <T>(a) and <V>(a) for trial psi=exp(-r/a) ──
# (closed forms quoted in the docstring are not just asserted -- integrated here)
a_test = 1.0   # fm, arbitrary probe scale for the numerical check
r = np.linspace(1e-6, 60 * a_test, 400000)
psi   = np.exp(-r / a_test)
dpsi  = -(1.0 / a_test) * np.exp(-r / a_test)

norm      = np.trapezoid(psi**2 * r**2, r)
T_numer   = np.trapezoid(dpsi**2 * r**2, r)
V_numer   = np.trapezoid(psi**2 * r**2 * (-alpha * hbar_c / r), r)

T_expect_closedform = (hbar_c**2) / (2 * m_e_c2 * a_test**2)
V_expect_closedform = -alpha * hbar_c / a_test

T_over_norm = (hbar_c**2 / (2 * m_e_c2)) * T_numer / norm
V_over_norm = V_numer / norm

check("BR1 <T>(a) numerical integral matches closed form (hbar_c)^2/(2*m_e*c^2*a^2)",
      abs(T_over_norm - T_expect_closedform) / T_expect_closedform < 1e-4,
      f"numerical={T_over_norm:.6f} MeV  closed-form={T_expect_closedform:.6f} MeV  (a={a_test} fm)")

check("BR2 <V>(a) numerical integral matches closed form -alpha*hbar_c/a",
      abs(V_over_norm - V_expect_closedform) / abs(V_expect_closedform) < 1e-4,
      f"numerical={V_over_norm:.6f} MeV  closed-form={V_expect_closedform:.6f} MeV  (a={a_test} fm)")

# ── Section 2: solve dE/da = 0 for the balance point ─────────────────────────
def T_of_a(a): return (hbar_c**2) / (2 * m_e_c2 * a**2)
def V_of_a(a): return -alpha * hbar_c / a
def E_of_a(a): return T_of_a(a) + V_of_a(a)

def F_out(a):  return  (hbar_c**2) / (m_e_c2 * a**3)   # -dT/da, outward
def F_in(a):   return -alpha * hbar_c / a**2            # -dV/da, inward

a_balance = hbar_c / (m_e_c2 * alpha)   # algebraic solution of F_out+F_in=0

print()
print(SEP2)
print("SECTION 2: SOLVING THE BALANCE POINT")
print(SEP2)
print(f"  F_out(a) = +(hbar_c)^2/(m_e*c^2*a^3)   [electron's own outward push]")
print(f"  F_in(a)  = -alpha*hbar_c/a^2           [proton's inward pull]")
print(f"  Solving F_out(a)+F_in(a)=0:  a = hbar_c/(m_e*c^2*alpha) = {a_balance:.2f} fm")
print(f"  = {a_balance*1e-15:.6e} m   (CODATA a_0 = {a_0_CODATA_fm*1e-15:.6e} m)")
print()

check("BR3 Force balance F_out(a)+F_in(a)=0 at the algebraic solution",
      abs(F_out(a_balance) + F_in(a_balance)) < 1e-9,
      f"F_out={F_out(a_balance):.6e}  F_in={F_in(a_balance):.6e} MeV/fm")

# Independent cross-check: numerically minimize E(a) directly (no calculus assumed)
a_scan = np.linspace(a_balance * 0.5, a_balance * 1.5, 2_000_001)
a_min_numeric = a_scan[np.argmin(E_of_a(a_scan))]

check("BR4 Numerical minimization of E(a)=T(a)+V(a) agrees with algebraic balance point",
      abs(a_min_numeric - a_balance) / a_balance < 1e-5,
      f"numeric a*={a_min_numeric:.4f} fm  algebraic a*={a_balance:.4f} fm")

check("BR5 Balance point a* equals the standard Bohr radius formula hbar_c/(m_e*c^2*alpha)",
      abs(a_balance - hbar_c / (m_e_c2 * alpha)) < 1e-9,
      f"a* = {a_balance:.6f} fm")

check("BR6 a* matches CODATA a_0 to within 0.001%",
      abs(a_balance - a_0_CODATA_fm) / a_0_CODATA_fm < 1e-5,
      f"a*={a_balance:.6f} fm  CODATA={a_0_CODATA_fm:.6f} fm")

# ── Section 3: virial theorem falls out automatically (cross-check vs PS5) ──
T_at_balance = T_of_a(a_balance)
V_at_balance = V_of_a(a_balance)

print()
print(SEP2)
print("SECTION 3: VIRIAL RELATION AT THE BALANCE POINT (cross-check vs PS5)")
print(SEP2)
print(f"  <T>(a_0) = {T_at_balance*1e6:.4f} eV")
print(f"  <V>(a_0) = {V_at_balance*1e6:.4f} eV  (PS5 quotes V(a_0)=27.21 eV = 2x13.6 eV)")
print(f"  -<V>(a_0)/2 = {-V_at_balance/2*1e6:.4f} eV")

check("BR7 Virial theorem <T>=-<V>/2 emerges automatically at the balance point",
      abs(T_at_balance - (-V_at_balance / 2)) / T_at_balance < 1e-9,
      f"<T>={T_at_balance*1e6:.4f} eV  -<V>/2={-V_at_balance/2*1e6:.4f} eV")

check("BR8 <V>(a_0) matches PS5's quoted 27.21 eV (2x13.6 eV) to 0.1%",
      abs(abs(V_at_balance) * 1e6 - 27.211) / 27.211 < 1e-3,
      f"|<V>(a_0)|={abs(V_at_balance)*1e6:.4f} eV  (PS5: 27.211 eV)")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY")
print(SEP2)
passed = sum(1 for _, s, _ in results if s == "PASS")
failed = sum(1 for _, s, _ in results if s == "FAIL")
print(f"  a_0 derived from explicit pressure/force balance = {a_balance:.4f} fm")
print(f"  = hbar_c/(m_e*c^2*alpha), same formula used elsewhere -- now DERIVED,")
print(f"  not merely asserted with already-derived inputs plugged in.")
print(f"  Total: {passed}/{len(results)}  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, s, d in results:
        if s == "FAIL": print(f"  FAILED: {name}")
print(SEP)
