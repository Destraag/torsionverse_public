"""
higgs_g_su2_from_alpha_weinberg.py

FEASIBILITY CHECK: can the SU(2) gauge coupling g be derived from cell
geometry, closing the gap flagged in doc_jobson_cell.txt ("g is not yet
derivable from the cell geometry", higgs_w_vertex.py) and doc_higgs?

IDEA (not yet tried elsewhere in the repo -- higgs_w_vertex.py searched for
NEW I_h mixing-angle mechanisms for theta_W itself and for direct T_1g
vertex corrections to m_W; it did not try combining the two quantities
that are ALREADY independently closed elsewhere):

  Standard EW relation (external, standard-model, not I_h-specific):
    e = g*sin(theta_W)          =>   g = e / sin(theta_W)
    alpha_em = e^2/(4*pi)       =>   e = sqrt(4*pi*alpha_em)
  Combining:  g = sqrt(4*pi*alpha_em) / sin(theta_W) = sqrt(4*pi*alpha_em/sin^2(theta_W))

  Torsionverse ALREADY independently derives both inputs:
    alpha         -- doc_alpha.txt, Born balance, closed to 0.00000022%
    sin^2(theta_W)* -- jobson_cell_doc.py J19 (GAP-C + 2*alpha^2*phi^2),
                     closed to 4.6e-6 [formula copied verbatim below]

  If g follows just from combining these two ALREADY-closed numbers via the
  standard (external) EW relation, g would ALSO become a zero-free-parameter
  geometric prediction. This script checks whether it actually does, and
  flags the one place a SEPARATE external input (RG running of alpha with
  energy scale) may be unavoidably needed.

Run: python analysis/higgs/higgs_g_su2_from_alpha_weinberg.py
"""
import math

pi = math.pi
phi = (1 + math.sqrt(5)) / 2
alpha_CODATA = 7.2973525693e-3     # low-energy (Thomson-limit) alpha, torsionverse's usual input
v_EW = 246.22                       # GeV, PDG (used throughout analysis/higgs/constants.py)
m_W_pdg = 80.377                    # GeV
m_Z_pdg = 91.1880                   # GeV

SEP = "=" * 70
print(SEP)
print("STEP 1: sin^2(theta_W)* -- torsionverse's own closed formula (J19)")
print(SEP)
# Copied verbatim from analysis/demos/jobson_cell_doc.py (uses alpha_CODATA)
sin2_tw = (1 - (math.sqrt(phi/math.sqrt(5))*(1+5*alpha_CODATA))**2) + 2*alpha_CODATA**2*phi**2
sin2_tw_pdg = 0.22290
print(f"  sin^2(theta_W)* = {sin2_tw:.8f}  (PDG {sin2_tw_pdg}, gap {sin2_tw-sin2_tw_pdg:.2e})")
print()

print(SEP)
print("STEP 2: g = sqrt(4*pi*alpha/sin^2(theta_W)) using LOW-ENERGY alpha")
print(SEP)
g_lowE = math.sqrt(4*pi*alpha_CODATA/sin2_tw)
g_measured = 2*m_W_pdg/v_EW
print(f"  alpha (CODATA, low-energy)   = {alpha_CODATA:.8e}")
print(f"  g_predicted (low-E alpha)    = {g_lowE:.6f}")
print(f"  g_measured (= 2*m_W_pdg/v_EW)= {g_measured:.6f}")
print(f"  error = {(g_lowE/g_measured - 1)*100:+.4f}%")
m_W_from_g_lowE = g_lowE * v_EW / 2
print(f"  => m_W predicted = g_lowE*v_EW/2 = {m_W_from_g_lowE:.4f} GeV"
      f"  (PDG {m_W_pdg} GeV, {(m_W_from_g_lowE/m_W_pdg-1)*100:+.4f}%)")
print()

print(SEP)
print("STEP 3: same, using RG-RUN alpha_em AT THE Z SCALE (external QED input)")
print(SEP)
# Standard literature value for the running EM coupling evaluated at m_Z
# (vacuum polarization from all charged fermions between Q~0 and Q~m_Z).
# This is a well-known external QED fact, NOT derived anywhere in torsionverse.
alpha_em_mZ_external = 1/127.955
g_running = math.sqrt(4*pi*alpha_em_mZ_external/sin2_tw)
print(f"  alpha_em(m_Z) [external, PDG RG-running value] = {alpha_em_mZ_external:.8e}")
print(f"  g_predicted (running alpha)  = {g_running:.6f}")
print(f"  g_measured                   = {g_measured:.6f}")
print(f"  error = {(g_running/g_measured - 1)*100:+.4f}%")
m_W_from_g_running = g_running * v_EW / 2
print(f"  => m_W predicted = {m_W_from_g_running:.4f} GeV"
      f"  (PDG {m_W_pdg} GeV, {(m_W_from_g_running/m_W_pdg-1)*100:+.4f}%)")
print()

print(SEP)
print("STEP 4: cross-check against m_Z = m_W/cos(theta_W)")
print(SEP)
cos_tw = math.sqrt(1 - sin2_tw)
m_Z_from_running = m_W_from_g_running / cos_tw
m_Z_from_lowE = m_W_from_g_lowE / cos_tw
print(f"  cos(theta_W)* = {cos_tw:.8f}")
print(f"  m_Z from low-E g  = {m_Z_from_lowE:.4f} GeV  (PDG {m_Z_pdg}, {(m_Z_from_lowE/m_Z_pdg-1)*100:+.4f}%)")
print(f"  m_Z from running g= {m_Z_from_running:.4f} GeV  (PDG {m_Z_pdg}, {(m_Z_from_running/m_Z_pdg-1)*100:+.4f}%)")
print()

print(SEP)
print("VERDICT")
print(SEP)
print("  If STEP 2 (low-energy alpha) is close: g is essentially already closed")
print("  by combining two ALREADY-derived torsionverse numbers (alpha, sin^2thetaW)")
print("  via the standard (external, not I_h-specific) e=g*sin(theta_W) relation.")
print("  If STEP 3 (running alpha) is needed instead: closing g fully requires")
print("  torsionverse to ALSO derive (or import) the RG running of alpha with")
print("  energy scale, which is not attempted anywhere in this repo currently --")
print("  a genuine, separate, well-defined open item, not a vague gap.")
