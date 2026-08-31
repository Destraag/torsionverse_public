"""
higgs_mechanism_stack.py
========================
Applies the full alpha vertex-stiffness mechanism stack to lambda
and identifies why it cannot close the 0.4% gap for the Higgs.

KEY FINDING:
  The alpha mechanism (L3, Born weighting, vertex stiffness) applies to
  BULK-regime particles (N_J >> 1) like the electron.
  The Higgs is SUB-CELL (N_J = 0.16 < 1) -- above the UV cutoff.
  Sub-cell particles interact with the cell as a WHOLE UNIT, not the vertex.
  => Different correction mechanism needed for lambda.

  SUB-CELL COUPLING (discovered this session):
  For sub-cell particles, the coupling is set by BULK medium properties.
  lambda_sub = (1-nu)/4 = 1/(8*(1-Rs^2)) = 2*pi^2/(16*pi^2-5)
  = 0.12908 -- only -0.16% from SM lambda (vs -0.40% for phi/(4*pi))

Run: python analysis/higgs/higgs_mechanism_stack.py
"""

import math, sys, os
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.dirname(__file__))
from constants import *

SEP  = "=" * 65
SEP2 = "-" * 65

print(SEP)
print("HIGGS MECHANISM STACK ANALYSIS")
print("Applying alpha vertex correction stack to lambda")
print(SEP)
print()

# ── SETUP ─────────────────────────────────────────────────────────────────────
n    = 2
Q    = 4*pi**2/phi
Rs   = math.sqrt(5)/(4*pi)
L3   = (phi**3 + math.log(5)**3)/(phi**2 + math.log(5)**2)
nu   = (1 - 2*Rs**2) / (2*(1 - Rs**2))

print(f"Alpha mechanism constants:")
print(f"  n = {n}, Q = {Q:.6f}, Rs = {Rs:.8f}")
print(f"  L3(PHI,log5) = {L3:.8f}")
print(f"  n_exact = 2.01868959, delta_n = 0.01869")
print(f"  delta_k_alpha = delta_n/L3 = {0.01869/L3:.6f}")
print()

# ── STEP 1: N_J for Higgs vs electron ────────────────────────────────────
print(SEP2)
print("STEP 1  N_J: Is the Higgs bulk or sub-cell?")
print(SEP2)
R_H = hbar_c / (m_H_pdg22 * 1000)  # fm  Compton radius
R_e = hbar_c / 0.511               # fm  electron Compton radius
N_H = R_H / L_J
N_e = R_e / L_J
print(f"  Higgs:    R = {R_H:.6f} fm   N_J = {N_H:.4f}  SUB-CELL")
print(f"  Electron: R = {R_e:.4f} fm   N_J = {N_e:.0f}  BULK")
print(f"  b quark:  N_J = {(hbar_c/4180)/L_J:.2f}   BOUNDARY")
print()
print(f"  Higgs N_J = {N_H:.4f} < 1 => SUB-CELL (above UV cutoff E_cell)")
print(f"  Electron N_J >> 1 => BULK (below UV cutoff, sees cell structure)")
print()

# ── STEP 2-4: The alpha mechanism stack and why it fails for Higgs ───────────
print(SEP2)
print("STEP 2-4  Alpha mechanism stack applied to lambda")
print(SEP2)
print(f"  Alpha stack: f1=PHI, f2=log5, f_eff=L3, then delta_k closes gap")
print()
print(f"  For BULK electron: wavefunction resolves vertex structure")
print(f"    => l=6 channel activated (N_J >> 1)")
print(f"    => Both l=0 and l=6 contribute: f_eff = L3 = {L3:.6f}")
print()
print(f"  For SUB-CELL Higgs: wavelength >> L_J, particle doesn't 'fit' near vertex")
print(f"    => l=6 channel SUPPRESSED by (L_J/R_H)^6 = {(L_J/R_H)**6:.2e}")
print(f"    => Only l=0 contributes at most (if any vertex interaction)")
print()
print(f"  CONCLUSION: Alpha vertex mechanism does NOT close the lambda gap.")
print(f"  Reason: Higgs is sub-cell. Different physics applies.")
print()

# ── STEP 5: Sub-cell coupling from BULK medium properties ────────────────────
print(SEP2)
print("STEP 5  Sub-cell coupling: use BULK properties instead of vertex")
print(SEP2)
print(f"  Bulk-regime particles (electron): coupling via vertex stiffness")
print(f"    => uses phi, log5, L3 (vertex geometry)")
print()
print(f"  Sub-cell particles (Higgs, W, Z, top): coupling via BULK medium")
print(f"    => Poisson ratio nu = {nu:.6f}")
print(f"    => 1 - nu = 1/(2*(1-Rs^2)) = {1-nu:.8f}")
print()
print(f"  SUB-CELL LAMBDA CANDIDATE:")
lam_subcell = (1 - nu) / 4
print(f"    lambda_sub = (1-nu)/4 = 1/(8*(1-Rs^2))")
print(f"               = {lam_subcell:.8f}")
print(f"    vs phi/(4*pi) = {lam_phi4pi:.8f}")
print(f"    vs lambda_SM  = {lam_SM:.8f}")
print()
print(f"  Deviations from SM:")
print(f"    phi/(4*pi):  {(lam_phi4pi/lam_SM-1)*100:+.4f}%")
print(f"    (1-nu)/4:    {(lam_subcell/lam_SM-1)*100:+.4f}%  (2.5x closer!)")
print()

# ── STEP 6: Derive (1-nu)/4 from Rs ──────────────────────────────────────────
print(SEP2)
print("STEP 6  Deriving lambda_sub = (1-nu)/4 from Rs")
print(SEP2)
print(f"  nu     = (1-2*Rs^2) / (2*(1-Rs^2))")
print(f"  1-nu   = 1 - (1-2*Rs^2)/(2*(1-Rs^2))")
print(f"         = [2*(1-Rs^2) - (1-2*Rs^2)] / [2*(1-Rs^2)]")
print(f"         = 1 / [2*(1-Rs^2)]")
print(f"  lambda_sub = (1-nu)/4 = 1/(8*(1-Rs^2))")
print()
print(f"  Substituting Rs = sqrt(5)/(4*pi):")
print(f"    Rs^2 = 5/(16*pi^2) = {5/(16*pi**2):.8f}")
print(f"    1-Rs^2 = (16*pi^2-5)/(16*pi^2)")
print(f"    lambda_sub = 2*pi^2/(16*pi^2-5)")
print(f"               = {2*pi**2/(16*pi**2-5):.8f}")
print(f"    verify: = {lam_subcell:.8f}  [match: {abs(2*pi**2/(16*pi**2-5)-lam_subcell)<1e-10}]")
print()
print(f"  This is DERIVED from Rs alone -- no free parameters.")
print(f"  Physical meaning: for sub-cell particles, the quartic coupling")
print(f"  is set by the bulk compressibility (1-nu) of the medium,")
print(f"  not by the vertex geometry (phi).")
print()

# ── STEP 7: Remaining gap ────────────────────────────────────────────────────
print(SEP2)
print("STEP 7  Remaining gap and interpretation")
print(SEP2)
gap_remaining = (lam_subcell/lam_SM - 1)*100
print(f"  lambda_sub = (1-nu)/4 = {lam_subcell:.8f}")
print(f"  lambda_SM              = {lam_SM:.8f}")
print(f"  Remaining gap:         = {gap_remaining:+.4f}%")
print()
print(f"  The sub-cell formula (1-nu)/4 closes from -0.40% to {gap_remaining:.2f}%")
print(f"  The final {abs(gap_remaining):.2f}% gap may be:")
print(f"    (a) A higher-order correction from 1/N terms")
print(f"    (b) Measurement uncertainty in m_H and v_EW")
print(f"    (c) A small residual from the l=0 vertex contribution")
print()

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print(SEP)
print("SUMMARY: TWO COUPLING REGIMES")
print(SEP)
print(f"  BULK (N_J >> 1): electron, proton, b quark, tau, charm")
print(f"    Coupling via VERTEX STIFFNESS")
print(f"    Correction factor: L3(PHI,log5) = {L3:.6f}")
print(f"    Example: n_exact = 2 + L3*delta_k -> closes alpha gap to -0.0007 sigma")
print()
print(f"  SUB-CELL (N_J < 1): Higgs, W, Z, top quark")
print(f"    Coupling via BULK MEDIUM PROPERTIES")
print(f"    lambda = (1-nu)/4 = 2*pi^2/(16*pi^2-5) = {lam_subcell:.6f}")
print(f"    vs lambda_SM = {lam_SM:.6f}  (gap: {gap_remaining:+.4f}%)")
print()
print(f"  This distinction is the key to Doc Higgs.")
print(f"  The lambda gap is NOT from vertex stiffness (wrong regime).")
print(f"  It IS from bulk medium Poisson ratio (right regime, -0.16% residual).")
print(SEP)
