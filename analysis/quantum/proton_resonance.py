import math
phi=(1+math.sqrt(5))/2; alpha=7.2973525693e-3; hbar_c=197.3269804
pi=math.pi

# Chi table including H_g with chi(C2)=+1
chi = {
    'A_g':  ( 1.0,  1.0,  1,  1, 1),
    'T_1g': ( phi, -1/phi, 0, -1, 3),
    'T_2g': (-1/phi, phi,  0, -1, 3),
    'G_g':  (-1.0, -1.0,  1,  0, 4),
    'H_g':  ( 0.0,  0.0, -1,  1, 5),
    'E+':   ( phi, -1/phi, 1,  0, 2),
    'E-':   (-1/phi, phi,  1,  0, 2),
    'G32':  ( 1.0,  1.0,  1,  0, 4),
    'I52':  (-1.0, -1.0,  0,  0, 6),
}

def pc(r1,r2): return tuple(chi[r1][i]*chi[r2][i] for i in range(4))+(chi[r1][4]*chi[r2][4],)
def matches(p,t): return all(abs(p[i]-chi[t][i])<1e-8 for i in range(5))

def sig_scale(m_MeV, chi_sq):
    """Peak cross-section: alpha * chi_sq * (hbar_c/m)^2 in cm^2."""
    return alpha * chi_sq * (hbar_c/m_MeV)**2 * 1e-26

# Verify T_1g x T_1g = A_g + T_1g + H_g  (neutron Zone 2 contains singlet A_g)
p_T1g_T1g = pc('T_1g','T_1g')
p_sum = tuple(chi['A_g'][i]+chi['T_1g'][i]+chi['H_g'][i] for i in range(4))+(9,)
assert all(abs(p_T1g_T1g[i]-p_sum[i])<1e-8 for i in range(5)), "T_1g x T_1g != A_g+T_1g+H_g"

E_cell = 124799.  # MeV
m_p=938.272; m_e=0.511; m_mu=105.66; m_tau=1776.86

# chi_sq = |chi(T_1g x target, C5)|^2 = (phi * chi(target,C5))^2 consistently
#   Electron   (E+,  chi_target=+phi):   product chi = phi*phi = phi^2     -> chi_sq = phi^4
#   Muon nu    (G32, chi_target=+1):     product chi = phi*1   = phi       -> chi_sq = phi^2

targets = [
    # (name,              irrep,  m_peak, cg_result,          chi_sq,   note)
    ('Free cell (A_g)',   'A_g',  0,      'T_1g (unchanged)',  0,        'passes through'),
    ('Proton Z2 (T_2g)',  'T_2g', m_p,    'G_g + H_g',        1.0,      'no A_g singlet'),
    ('Neutron Z2 (T_1g)', 'T_1g', m_p,    'A_g+T_1g+H_g',    phi**4,   'A_g: singlet allowed!'),
    ('Freed nu_e (E-)',   'E-',   m_tau,  'I52 (tau)',         1.0,      'tau resonance (WI9)'),
    ('Electron (E+)',     'E+',   m_e,    'E++G32 mix',       phi**4,   'Thomson; chi(C5)=phi^2'),
    ('Muon nu (G32)',     'G32',  m_mu,   'dim=12, chi=phi',  phi**2,   'chi(C5)=phi; sub-tau'),
]

print(f'PHOTON / FREED LEPTON RESONANCE RANGES vs PROTON ZONE 2 SHELL')
print(f'  Mode                   Resonance    Regime at E_nu=5 MeV   Regime at E_nu=1 GeV')
print(f'  {"-"*85}')
modes_res = [
    ('Free cell (A_g)',     E_cell,  'far sub-resonant',  'sub-resonant'),
    ('Proton Z2 (T_2g)',   m_p,     'sub-resonant',      'OVER-resonant (falling)'),
    ('Neutron Z2 (T_1g)',  m_p,     'sub-resonant',      'OVER-resonant (falling)'),
    ('Freed nu_e (E-)',    m_tau,   'sub-resonant',      'sub-resonant'),
    ('Electron (E+)',      m_e,     'OVER-resonant',     'OVER-resonant'),
]
for name, E_res, at5, at1000 in modes_res:
    print(f'  {name:<22}  {E_res:<12.1f}  {at5:<22}  {at1000}')
print()
print(f'  KEY: Proton Zone 2 resonance = m_p = 938 MeV  <<  Free cell = E_cell = 124,799 MeV')
print(f'  The proton shell resonates LOWER because N_J=21 locked cells act as heavier mass.')
print(f'  A neutrino at E_nu=5 MeV is sub-resonant to Zone 2 but MUCH less sub-resonant than')
print(f'  to the free cell: sigma ratio (proton/free_cell) = (E_cell/m_p)^4 = {(E_cell/m_p)**4:.2e}')
print()
print(f'  IN THE RANGE m_p < E_nu < E_cell  (0.94 GeV to 124.8 GeV):')
print(f'    Proton Zone 2: OVER-RESONANT (sigma falling as (m_p/E)^2)')
print(f'    Free cell:     still sub-resonant (sigma rising as (E/E_cell)^2)')
print(f'    Crossover at E_nu = E_cell = {E_cell:.0f} MeV -- they become equal there.')
print()

# Sigma vs energy table
print(f'  COUPLING sigma (cm^2) AT DIFFERENT NEUTRINO ENERGIES:')
print(f'  {"E_nu":<12} {"Proton Z2":<20} {"Free cell":<20} {"Ratio P/C":<15} Proton regime')
print(f'  {"-"*80}')
for E_nu, label in [(1,  '1 MeV'), (5, '5 MeV'), (30, '30 MeV (SN)'),
                    (938.3,'938 MeV (peak)'), (5000,'5 GeV'), (E_cell,'E_cell')]:
    # Proton Zone 2 (chi_sq=1, peak at m_p=938)
    if E_nu <= m_p:
        s_p = sig_scale(m_p, 1.0) * (E_nu/m_p)**2  # sub-resonant rising
        regime = 'sub-resonant'
    else:
        s_p = sig_scale(m_p, 1.0) * (m_p/E_nu)**2  # above-resonant falling
        regime = 'OVER-resonant'
    # Free cell (chi_sq=1, peak at E_cell)
    s_c = sig_scale(E_cell, 1.0) * (E_nu/E_cell)**2  # always sub-resonant
    ratio = s_p / s_c if s_c > 0 else float('inf')
    print(f'  {label:<12} {s_p:<20.2e} {s_c:<20.2e} {ratio:<15.2e} {regime}')

print()
print(f'  E_cell = {E_cell:.0f} MeV = absolute coupling ceiling')
print(f'  Sigma at E above peak falls as (E_peak/E)^2')
print(f'  Sigma at E below peak rises as (E/E_peak)^2')
print()
# Ratio: proton Zone 2 coupling vs free cell at E=1 GeV
E_GeV = 1000.  # MeV = 1 GeV
r_proton = (m_p/E_GeV)**2 * (hbar_c/m_p)**2  # fm^2, above proton resonance
r_cell   = (E_GeV/E_cell)**2 * (hbar_c/E_cell)**2  # fm^2, below cell resonance
print(f'  At E=1 GeV, above proton resonance and below cell resonance:')
print(f'    sigma(proton Z2) ~ {r_proton*1e-26:.2e} cm^2  (above-resonant, falling)')
print(f'    sigma(free cell) ~ {r_cell*1e-26:.2e} cm^2  (sub-resonant, rising)')
print(f'    Ratio:  sigma(proton)/sigma(free cell) = {r_proton/r_cell:.2e}')
print(f'    => Proton Zone 2 dominates by factor {r_proton/r_cell:.2e} at 1 GeV')