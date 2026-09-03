"""
jobson_depletion_rho_mu0_impact.py
===================================
Quantifies the impact of switching rho_medium from rho_Lambda to mu_0 on
analysis/molecular/jobson_depletion.py
(JD1-JD12), the ONLY other place found in the repo where an ABSOLUTE (non-
cancelling) value of G_shear = rho*(Rs*c)^2 feeds a published check, rather
than a ratio of two rho-proportional quantities (which cancels rho and is
therefore unaffected -- confirmed separately for every other rho_Lambda
consumer in analysis/alpha, analysis/gravity, analysis/medium).

jobson_depletion.py is the companion script for docs/series2/doc_chemistry.txt
Mechanism 3 (T_2g toroidal vortex medium depletion) -- a mechanism the doc's
OWN text already calls "OPEN" (depletion depth "not yet derived"). This script
does not decide anything; it only reports what changes.

Reference: analysis/molecular/jobson_depletion.py JD1-JD6,
           docs/series2/doc_chemistry.txt Section 3b Mechanism 3.
"""
import math

pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
Rs    = math.sqrt(5) / (4 * pi)
c_ms  = 2.99792458e8
tau_relax = 1.0e8            # s (lower bound, wave_dispersion.py)
L_grain   = 9.9347e-18       # m
lam_m     = 2.70e-6          # m (O-H resonance)
I_center  = 1.31e10          # W/cm^2 (BD12 working point)
I_IRMPD   = 1.47e9           # W/cm^2

SEP = "=" * 65
results = []
def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status))
    print(f"  [{'PASS' if cond else '*** FAIL'}] {name}")
    if detail: print(f"         {detail}")

def jd_block(rho, label):
    G_shear = rho * (Rs * c_ms)**2
    I_sat_Wcm2 = G_shear * c_ms * 1e-4
    P_rad = (I_center * 1e4) / c_ms
    v_creep = P_rad * L_grain / (G_shear * tau_relax)
    t_deplete_phi2 = (lam_m / v_creep) / phi**2
    F_crit_cm2 = (lam_m * G_shear * tau_relax * c_ms / L_grain) * 1e-4
    print(f"  [{label}] rho={rho:.4e} kg/m^3  G_shear={G_shear:.4e} Pa")
    print(f"    I_sat      = {I_sat_Wcm2:.4e} W/cm^2   (JD1 threshold: < 1e-5)")
    print(f"    v_creep    = {v_creep:.4e} m/s")
    print(f"    t_deplete  = {t_deplete_phi2:.4e} s  ({t_deplete_phi2/60:.2f} min; JD4 threshold: < 3600 s)")
    print(f"    F_crit     = {F_crit_cm2:.4e} J/cm^2")
    print(f"    Ordering I_sat << I_IRMPD << I_center? "
          f"{I_sat_Wcm2 < I_IRMPD < I_center}  "
          f"(I_sat={I_sat_Wcm2:.2e}, I_IRMPD={I_IRMPD:.2e}, I_center={I_center:.2e})")
    return I_sat_Wcm2, v_creep, t_deplete_phi2, F_crit_cm2

print(SEP)
print("IMPACT OF rho=mu_0 vs rho=rho_Lambda ON jobson_depletion.py (JD1-JD12)")
print(SEP)
print()
rho_Lambda = 5.8424e-27
rho_mu0    = 4 * pi * 1e-7
I_sat_L, v_creep_L, t_dep_L, F_crit_L = jd_block(rho_Lambda, "as-published, rho_Lambda")
print()
I_sat_M, v_creep_M, t_dep_M, F_crit_M = jd_block(rho_mu0, "rho=mu_0")
print()

check("JD1 (as published, rho_Lambda): I_sat < 1e-5 W/cm^2",
      I_sat_L < 1e-5, f"I_sat={I_sat_L:.3e} W/cm^2")
check("JD1 (rho=mu_0): I_sat < 1e-5 W/cm^2 -- does the ORIGINAL claim survive?",
      I_sat_M < 1e-5, f"I_sat={I_sat_M:.3e} W/cm^2")
check("JD4 (as published, rho_Lambda): t_deplete < 3600 s",
      t_dep_L < 3600, f"t_deplete={t_dep_L:.3e} s")
check("JD4 (rho=mu_0): t_deplete < 3600 s -- does the ORIGINAL claim survive?",
      t_dep_M < 3600, f"t_deplete={t_dep_M:.3e} s")
check("Operating window I_sat << I_IRMPD << I_center holds (rho_Lambda, as published)",
      I_sat_L < I_IRMPD < I_center, "")
check("Operating window I_sat << I_IRMPD << I_center holds (rho=mu_0)",
      I_sat_M < I_IRMPD < I_center, "")

print()
print(SEP)
n_pass = sum(1 for _, s in results if s == "PASS")
print(f"SUMMARY: {n_pass}/{len(results)} PASS, {len(results)-n_pass} FAIL")
print("(FAILs here are EXPECTED/DIAGNOSTIC, not bugs -- this script exists to")
print(" show exactly which claims break under rho=mu_0, not to enforce PASS.)")
print(SEP)
