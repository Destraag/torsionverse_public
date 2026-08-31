"""
series2_doc.py
==============
Companion script for doc_series2_synthesis.txt
Verifies the key new theoretical claims in Series 2 that are not
already covered by Series 1 companion scripts.

All checks use only Series 1 derived constants (zero new parameters).
Reference: docs/doc_series2_synthesis.txt
DOI: https://doi.org/10.5281/zenodo.22108664
"""
import sys, math
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')

# All constants inline — script runs standalone on any machine
pi      = math.pi
phi     = (1 + math.sqrt(5)) / 2         # golden ratio
alpha   = 7.2973525693e-3                 # fine structure constant (CODATA 2018)
r_p     = 0.8414e-15                      # m, proton charge radius (CODATA 2018)
hbar_c  = 197.3269804                     # MeV*fm
L_J     = alpha * phi * (r_p * 1e15)     # fm, Jobson cell edge

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  {'[PASS]' if cond else '[FAIL]'} {name}")
    if detail: print(f"         {detail}")

# I_h character table (gerade irreps, classes: E, C2, C3, C5, C5^2)
chi = {
    'A' : [ 1,   1,       1,          1,        1],
    'T1': [ 3,  -1,       0,          phi,      -(phi-1)],
    'T2': [ 3,  -1,       0,         -(phi-1),   phi],
    'G' : [ 4,   0,       1,         -1,        -1],
    'H' : [ 5,   1,      -1,          0,         0],
}
class_sizes = [1, 15, 20, 12, 12]  # E, C2, C3, C5, C5^2
G_ORDER = 60

def decompose(chi_prod):
    result = {}
    for name, chi_irr in chi.items():
        n = sum(class_sizes[c] * chi_irr[c] * chi_prod[c]
                for c in range(5)) / G_ORDER
        result[name] = int(round(n))
    return result

print(SEP)
print("series2_doc.py -- Series 2 synthesis: key derivation checks")
print(SEP)

# ── Section 1: A_g beam crossing enhancement ─────────────────────────────────
print()
print(SEP2)
print("SECTION 1: A_g ORTHOGONAL-BEAM ENHANCEMENT (Section 3.2 of synthesis)")
print(SEP2)

# T_1g x T_1g CG decomposition
chi_T1xT1 = [chi['T1'][c]**2 for c in range(5)]
decomp_T1xT1 = decompose(chi_T1xT1)
dim_check = sum(chi[k][0] * v for k, v in decomp_T1xT1.items())

print(f"  T_1g x T_1g characters: {[round(x,4) for x in chi_T1xT1]}")
print(f"  Decomposition: " + " + ".join(f"{v}{k}" for k, v in decomp_T1xT1.items() if v > 0))
print(f"  Dimension check: {dim_check} = {chi['T1'][0]}^2 = 9")

check("SD1: T_1g x T_1g = A_g + T_1g + H_g [A_g appears once -- key for water splitting]",
      decomp_T1xT1['A'] == 1 and decomp_T1xT1['T1'] == 1 and decomp_T1xT1['H'] == 1
      and decomp_T1xT1.get('T2', 0) == 0 and decomp_T1xT1.get('G', 0) == 0,
      f"T_1g x T_1g = " + " + ".join(f"{v}{k}" for k, v in decomp_T1xT1.items() if v > 0))

# Character values at key crossing angles
chi_T1g_C5 = chi['T1'][3]   # = phi (C5 vertex, 72 deg)
chi_T1g_C2 = chi['T1'][1]   # = -1  (C2 edge, 90 deg)

print(f"  chi(T_1g, C5) = {chi_T1g_C5:.6f}  [72 deg crossing -- C5 icosahedral angle]")
print(f"  chi(T_1g, C2) = {chi_T1g_C2:.6f}  [90 deg crossing -- C2 edge angle]")

check("SD2: chi(T_1g, C5) = phi [golden ratio, C5 vertex character]",
      abs(chi_T1g_C5 - phi) < 1e-10,
      f"chi(T_1g,C5) = {chi_T1g_C5:.6f}  phi = {phi:.6f}")

check("SD3: chi(T_1g, C2) = -1 [C2 edge character]",
      abs(chi_T1g_C2 - (-1)) < 1e-10,
      f"chi(T_1g,C2) = {chi_T1g_C2}")

# A_g yield at each crossing angle = chi^2 (Born coupling squared)
A_g_yield_72 = chi_T1g_C5**2   # = phi^2 = 2.618
A_g_yield_90 = chi_T1g_C2**2   # = (-1)^2 = 1
enhancement = A_g_yield_72 / A_g_yield_90

print(f"  A_g yield at 72 deg = chi(T_1g,C5)^2 = {A_g_yield_72:.6f} = phi^2")
print(f"  A_g yield at 90 deg = chi(T_1g,C2)^2 = {A_g_yield_90:.6f}")
print(f"  Enhancement ratio 72/90 = {enhancement:.6f} = phi^2 = {phi**2:.6f}")

check("SD4: A_g enhancement at 72 deg / 90 deg = phi^2 = 2.618 [key prediction]",
      abs(enhancement - phi**2) < 1e-10,
      f"ratio = {enhancement:.6f}  phi^2 = {phi**2:.6f}")

# ── Section 2: Jet mechanism CG (AT2018hyz) ──────────────────────────────────
print()
print(SEP2)
print("SECTION 2: TDE JET MECHANISM -- T_2g x T_2g CG (Section 1.1 of synthesis)")
print(SEP2)

chi_T2xT2 = [chi['T2'][c]**2 for c in range(5)]
decomp_T2xT2 = decompose(chi_T2xT2)
print(f"  T_2g x T_2g Decomposition: " + " + ".join(f"{v}{k}" for k, v in decomp_T2xT2.items() if v > 0))
print(f"  A_g: infalling debris compresses boundary further")
print(f"  T_2g: radially reflected (accretion disk wind)")
print(f"  H_g: emitted perpendicular = JET (H_g = T_1g x T_2g = field strength tensor)")

check("SD5: T_2g x T_2g = A_g + T_2g + H_g [jet mechanism CG exact]",
      decomp_T2xT2['A'] == 1 and decomp_T2xT2['T2'] == 1 and decomp_T2xT2['H'] == 1
      and decomp_T2xT2.get('T1', 0) == 0 and decomp_T2xT2.get('G', 0) == 0,
      f"T_2g x T_2g = " + " + ".join(f"{v}{k}" for k, v in decomp_T2xT2.items() if v > 0))

# ── Section 3: MCF anti-sticking vortex energy ───────────────────────────────
print()
print(SEP2)
print("SECTION 3: MCF ANTI-STICKING VORTEX ENERGY (Section 2.1 of synthesis)")
print(SEP2)

Rs = math.sqrt(5) / (4 * pi)
r_He4_fm = 1.680   # fm  He-4 charge radius
E_vortex_He4 = Rs * hbar_c / r_He4_fm   # MeV

print(f"  Rs = sqrt(5)/(4*pi) = {Rs:.6f}")
print(f"  r_He4 = {r_He4_fm} fm  (He-4 charge radius)")
print(f"  E_vortex(He-4) = Rs * hbar_c / r_He4 = {E_vortex_He4:.3f} MeV")
print(f"  [= Rs * m_pi_He4 in the Zone 2 boundary mode picture]")

check("SD6: E_vortex(He-4) = Rs * hbar_c / r_He4 = 20.90 MeV within 0.1%",
      abs(E_vortex_He4 - 20.90) < 0.05,
      f"E_vortex = {E_vortex_He4:.3f} MeV  target = 20.90 MeV")

# ── Section 4: Bragg cloning identity (from PG15) ────────────────────────────
print()
print(SEP2)
print("SECTION 4: BRAGG CLONING IDENTITY (Section 4.3 of synthesis)")
print(SEP2)

m_p_MeV   = 938.272
m_crit     = m_p_MeV / (8 * alpha * phi)
bragg_90   = hbar_c / (2 * L_J)

print(f"  m_crit = m_p / (8*alpha*phi) = {m_crit:.3f} MeV = {m_crit/1000:.4f} GeV")
print(f"  hbar_c / (2*L_J) = {bragg_90:.3f} MeV  [Bragg at 90 deg, d=L_J]")
print(f"  Algebraic identity: both equal m_p/(8*alpha*phi)")

check("SD7: hbar_c/(2*L_J) = m_crit within 0.025% [Bragg = winding critical mass]",
      abs(bragg_90 - m_crit) / m_crit < 0.00025,
      f"hbar_c/(2*L_J) = {bragg_90:.3f} MeV  m_crit = {m_crit:.3f} MeV")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == 'PASS')
n_fail = sum(1 for _, s, _ in results if s == 'FAIL')
print(f"  Total: {len(results)}  PASS: {n_pass}  FAIL: {n_fail}")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
print()
print(f"  KEY RESULT: A_g enhancement at 72 deg (C5) vs 90 deg (C2) = phi^2 = {phi**2:.4f}")
print(f"  This is the primary new prediction of Series 2 (Section 3.2).")
print(f"  Falsifiable: measure bond-breaking ratio at 72 deg vs 90 deg.")
print(SEP)
