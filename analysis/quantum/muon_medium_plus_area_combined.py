"""
muon_medium_plus_area_combined.py

Precise combined test, using the ACTUAL functional form from open_items.txt's
G-ALTERATION-MEDIUM mechanism (additive: alpha_eff ~ alpha_topological +
vertex_correction*(rho/rho_Earth)) rather than folding density into the
N inside a 1/(1+N*alpha) denominator (which is a different, not-yet-
validated functional form of my own invention).

STRUCTURE:
  eff_mu = eff_topological + delta_vertex * density_ratio

  eff_topological: a clean, group-theory candidate (phi or 1/phi -- both
    already independently motivated from Section 4.1's own G32|D5 roots).
  delta_vertex: the SAME kind of small "vertex correction" alpha produces
    for itself (order alpha, analogous to alpha's own +0.000560% shift).
  density_ratio: the muon's local medium compression relative to the
    electron's, from its larger displaced volume -- tested at both the
    linear mass ratio (207) and the area-scaling power (207^(2/3)=35),
    since "did you combine both" specifically asks about area scaling.

Run: python analysis/quantum/muon_medium_plus_area_combined.py
"""
import math

phi = (1 + math.sqrt(5)) / 2
alpha = 7.2973525693e-3
log5 = math.log(5)
Rs2 = (math.sqrt(5) / (4 * math.pi)) ** 2
m_p = 938.272046
m_e_pdg = 0.51099895
m_mu_pdg = 105.6583755
poly = 5 * math.tan(math.pi / 5)
mass_ratio = m_mu_pdg / m_e_pdg

def mass_formula(eff_mu):
    L3 = (eff_mu**3 + log5**3) / (eff_mu**2 + log5**2)
    x = alpha * eff_mu**2
    k = alpha * eff_mu * (1 - 0.75 * alpha**2) / (1 + x + x**2)
    dn = L3 * k
    base = 2 * math.pi * alpha * (2 / math.sqrt(5)) * phi**2 * m_p
    corr = 1 + Rs2 + 2 * alpha
    return base * (1 + dn / poly) * corr

print(f"mass_ratio = {mass_ratio:.4f}, ^(2/3) = {mass_ratio**(2/3):.4f}")
print()

eff_topologicals = {"phi": phi, "1/phi": 1/phi}
density_ratios = {
    "linear mass ratio (207)": mass_ratio,
    "area scaling (207^2/3=35)": mass_ratio**(2/3),
}
# delta_vertex tested across a small range since its exact size isn't derived
delta_scales = [alpha, alpha*phi, alpha**2, alpha/phi]

print(f"  {'eff_topo':10s} {'density basis':28s} {'delta_vertex':14s} {'eff_mu':>10s} {'m_mu (MeV)':>12s} {'error %':>10s}")
print(f"  {'-'*10} {'-'*28} {'-'*14} {'-'*10} {'-'*12} {'-'*10}")
best = None
for topo_name, eff_topo in eff_topologicals.items():
    for dens_name, dens in density_ratios.items():
        for delta in delta_scales:
            eff = eff_topo + delta * dens
            if eff <= 0:
                continue
            m_pred = mass_formula(eff)
            err = (m_pred - m_mu_pdg) / m_mu_pdg * 100
            print(f"  {topo_name:10s} {dens_name:28s} {delta:14.6e} {eff:10.6f} {m_pred:12.4f} {err:+10.4f}")
            if best is None or abs(err) < abs(best[4]):
                best = (topo_name, dens_name, delta, eff, err)
    print()

print(f"Best: eff_topo={best[0]}, density={best[1]}, delta={best[2]:.4e}, eff={best[3]:.6f}, error={best[4]:+.4f}%")
print(f"Reference bipyramid: error={(mass_formula((9-math.sqrt(5))/8)-m_mu_pdg)/m_mu_pdg*100:+.4f}%")
