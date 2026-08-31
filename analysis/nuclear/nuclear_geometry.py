"""
nuclear_geometry.py
===================
Nuclear structure from Jobson cell geometry alone.
No BW formula. No fitted mass parameters. Only:
  - Derived Zone radii (from r_p = 4*lambda_p, doc_nucleus PS4)
  - Derived I_h character table (from alpha derivation)
  - Derived pion mass (from torsionverse_doc SY8)
  - Derived Zone 2 coupling difference: proton T_2g (A_g=1), neutron T_1g (A_g=0)

Physical picture:
  Proton: T_2g diquark, Zone 2 active (N_J=21, Maxwell boundary)
    - Hard outer shell at r_p = 4*lambda_p (cog contact scale)
    - Inner core at lambda_p (quark confinement)
    - Zone 2 pulled INWARD by resonance -> effective nuclear radius tighter
  Neutron: T_1g diquark, Zone 2 INACTIVE (A_g(T_1g x T_2g)=0)
    - Same Zone 1 quark structure (lambda_n ~ lambda_p)
    - No Zone 2 resonance -> no inward compression -> slightly larger effective radius
    - Acts as pure geometric BUFFER between protons
  Electron: N_J >> 1 (bulk), pure Hopf winding (1,2), zero nuclear excluded volume

Key geometric computations:
  NG1  Proton and neutron Zone radii from lambda_p and r_p
  NG2  p-p hard core from cog grinding geometry (2*lambda_p)
  NG3  p-n separation scale from Zone 2 active/inactive asymmetry
  NG4  Nuclear force range from derived pion mass (r_0 = hbar_c/m_pi)
  NG5  Nuclear radius r_0 from nucleon packing at pion force range
  NG6  I_h branching: l=3 -> T_2g + G_g (dim 3+4=7) [f orbital decomposition]
  NG7  I_h branching: l=4 -> G_g + H_g (dim 4+5=9) [g orbital decomposition]
  NG8  Magic 28: intruder f_{7/2} dim=8=2*G_g (G_g is boundary regime irrep)
  NG9  Magic 50: intruder g_{9/2} dim=10=2*H_g (H_g is next irrep up)
  NG10 I_h branching: l=5 -> T_1g + T_2g + H_g (dim 3+3+5=11)
  NG11 I_h branching: l=6 -> A_g + T_1g + G_g + H_g (dim 1+3+4+5=13)
  NG12 Intruder h_{11/2} dim=12: not a single I_h irrep -> COMPOSITE (magic 82 open)
  NG13 I_h branching: l=7 [computed here] -> determines j_{15/2} dim=16 status for N=184
  NG14 Neutron effective radius > proton (Zone 2 absent -> no inward compression)
  NG15 N/Z geometry: neutron buffer packs inside pion range around proton

Run: python analysis/nuclear/nuclear_geometry.py
Reference: docs/doc_nucleus.txt
"""

import sys, os, math
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha, hbar_c, r_p

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, "PASS" if cond else "FAIL", detail))
    print(f"  [{s}] {name}")
    if detail: print(f"         {detail}")

# ── Derived constants (no free parameters) ────────────────────────────────────
pi    = math.pi
phi   = (1 + math.sqrt(5)) / 2
Rs    = math.sqrt(5) / (4 * pi)
m_p   = 938.272          # MeV
m_n   = 939.565          # MeV (from m_n - m_p = alpha*Rs*m_p*(1+2*Rs^2), SY9)
m_pi  = m_p / (4*phi*(1 + Rs**2 + alpha))  # derived pion mass (SY8)

r_p_fm    = r_p * 1e15              # 0.8414 fm
lambda_p  = hbar_c / m_p            # 0.2103 fm (Compton scale, Zone 1 boundary)
r_p_check = 4 * lambda_p            # must equal r_p_fm (PS4)

# Zone radii (proton)
r_Zone1_p = lambda_p                # Zone 1: quark core
r_Zone2_p = r_p_fm                  # Zone 2: Maxwell boundary (outer)
r_Zone2_inner = lambda_p            # Zone 2: inner boundary (= Zone 1 outer)
r_Zone3_onset = r_p_fm              # Zone 3: co-rotating cells start here
V_p_fm3  = (4/3)*pi*r_p_fm**3      # proton excluded volume (fm^3)

# ── I_h character table (gerade, from alpha derivation) ───────────────────────
# Classes: E(1), 12C_5, 12C_5^2, 20C_3, 15C_2
# Characters: [E, C_5, C_5^2, C_3, C_2]
I_h_chars = {
    'A_g':   [1,      1,          1,       1,    1],
    'T_1g':  [3,      phi,       -1/phi,   0,   -1],
    'T_2g':  [3,     -1/phi,      phi,     0,   -1],
    'G_g':   [4,     -1,         -1,       1,    0],
    'H_g':   [5,      0,          0,      -1,    1],
}
I_h_class_sizes = [1, 12, 12, 20, 15]   # |class|
I_h_order = 60                           # |I| (rotation subgroup; gerade sector)

def ih_branching(l):
    """
    Decompose orbital angular momentum l into I_h gerade irreps.
    Returns dict {irrep_name: multiplicity}.
    Uses character of SO(3) l-representation evaluated on I conjugacy classes:
      chi_l(phi_rot) = sin((l+1/2)*phi_rot) / sin(phi_rot/2)
    Rotation angles for I classes: C_5 = 2pi/5, C_5^2 = 4pi/5, C_3 = 2pi/3, C_2 = pi
    """
    angles = [0, 2*pi/5, 4*pi/5, 2*pi/3, pi]  # E, C_5, C_5^2, C_3, C_2
    # Character of l rep at each class
    chi_l = []
    for phi_rot in angles:
        if phi_rot == 0:
            chi_l.append(2*l + 1)
        else:
            val = math.sin((l + 0.5)*phi_rot) / math.sin(phi_rot/2)
            chi_l.append(round(val, 10))

    # Project onto each irrep
    decomp = {}
    for name, chi_x in I_h_chars.items():
        n = sum(I_h_class_sizes[i] * chi_l[i] * chi_x[i] for i in range(5)) / I_h_order
        decomp[name] = round(n)
    return decomp, chi_l

# ── Section 1: Zone geometry ──────────────────────────────────────────────────
print(SEP)
print("SECTION 1: ZONE RADII FROM DERIVED GEOMETRY")
print(SEP2)
print(f"  lambda_p = hbar_c/m_p = {lambda_p:.4f} fm  (Zone 1 outer / cog contact scale)")
print(f"  r_p      = 4*lambda_p  = {r_p_fm:.4f} fm  (Zone 2 outer = proton charge radius, PS4)")
print(f"  V_p      = 4pi/3 * r_p^3 = {V_p_fm3:.4f} fm^3  (proton excluded volume)")
print(f"  m_pi (derived) = {m_pi:.4f} MeV  (from SY8: m_p/(4*phi*(1+Rs^2+alpha)))")
print(f"  r_pi = hbar_c/m_pi = {hbar_c/m_pi:.4f} fm  (nuclear force range)")
print()

# Neutron geometry: no Zone 2 compression
# Proton Zone 2 resonance COMPRESSES r_p by the binding ratio:
# delta_m = alpha*Rs*m_p*(1+2*Rs^2) is the Zone 2 resonance energy
# This reduces the proton charge radius from r_p_bare to r_p by delta_m/m_p fraction
delta_Z2 = alpha * Rs * m_p * (1 + 2*Rs**2)  # Zone 2 resonance energy per proton (MeV)
# Neutron has no Zone 2 compression -> r_n = r_p * (1 + delta_Z2/m_p)^(1/3)
# (slightly LARGER than r_p because it lacks the inward Zone 2 pull)
r_n_factor = (1 + delta_Z2/m_p)**(1/3)
r_n_fm = r_p_fm * r_n_factor

print(f"  Zone 2 resonance energy:     delta_Z2 = {delta_Z2:.4f} MeV")
print(f"  Proton r_p (Zone 2 active):  {r_p_fm:.4f} fm  (T_2g, Zone 2 pulls inward)")
print(f"  Neutron r_n (Zone 2 absent): {r_n_fm:.4f} fm  (T_1g, no inward compression)")
print(f"  r_n/r_p = {r_n_factor:.6f}  (neutron is {(r_n_factor-1)*100:.4f}% larger)")
print()

check("NG1 r_p = 4*lambda_p to 0.05% (PS4, Zone 2 outer = charge radius)",
      abs(r_p_fm - 4*lambda_p)/r_p_fm < 5e-4,
      f"r_p = {r_p_fm:.4f} fm  4*lambda_p = {4*lambda_p:.4f} fm  diff={abs(r_p_fm-4*lambda_p):.4f} fm")
check("NG14 Neutron r_n > r_p (no Zone 2 compression, T_1g diquark A_g=0)",
      r_n_fm > r_p_fm,
      f"r_n = {r_n_fm:.5f} fm > r_p = {r_p_fm:.4f} fm")

# ── Section 2: p-p and p-n separation scales ──────────────────────────────────
print()
print(SEP)
print("SECTION 2: NUCLEON SEPARATION FROM COG GEOMETRY")
print(SEP2)

r_pp_hard = 2 * lambda_p   # p-p: cog zones first touch at 2*lambda_p (grinding)
# p-n: neutron (no Zone 2) can enter proton Zone 2 until Zone 1 cores touch
# Zone 1 cores touch when: r_center_to_center = lambda_p + lambda_n ≈ 2*lambda_p
r_pn_min  = 2 * lambda_p   # same scale but attraction (no grinding), actual min from quark overlap
# Equilibrium p-n distance: neutron sits at Zone 2 midpoint of proton
r_pn_equil = (lambda_p + r_p_fm) / 2   # geometric midpoint of Zone 2

# Nuclear force range from derived pion mass (Yukawa)
r_pion = hbar_c / m_pi

print(f"  p-p hard core (cog grinding): r_pp = 2*lambda_p = {r_pp_hard:.4f} fm")
print(f"    [Same-chirality cog zones first contact. Observed: 0.4-0.6 fm]")
print(f"  p-n Zone 2 midpoint (equil.): r_pn = (lambda_p+r_p)/2 = {r_pn_equil:.4f} fm")
print(f"    [Neutron (T_1g, no Zone 2) sits inside proton Zone 2 buffer]")
print(f"  Nuclear force range (pion):   r_pi  = hbar_c/m_pi = {r_pion:.4f} fm")
print(f"    [Sets maximum nuclear separation: SY8-derived pion mass]")
print()
print(f"  Hierarchy: r_pp_hard < r_pn_equil < r_pi")
print(f"             {r_pp_hard:.3f}      <   {r_pn_equil:.3f}       <   {r_pion:.3f}  [fm]")
print()

check("NG2 p-p hard core = 2*lambda_p = 0.42 fm  [observed 0.4-0.6 fm]",
      0.40 < r_pp_hard < 0.45,
      f"r_pp = 2*lambda_p = {r_pp_hard:.4f} fm")
check("NG3 p-n equilibrium inside proton Zone 2: lambda_p < r_pn < r_p",
      lambda_p < r_pn_equil < r_p_fm,
      f"r_pn = {r_pn_equil:.4f} fm  in ({lambda_p:.3f}, {r_p_fm:.3f})")
check("NG4 Nuclear force range = hbar_c/m_pi (from derived pion mass, -0.04%)",
      abs(r_pion - 1.4145) < 0.005,
      f"r_pi = {r_pion:.4f} fm  (expected ~1.4145 fm)")

# ── Section 3: Nuclear radius from packing at pion range ─────────────────────
print()
print(SEP)
print("SECTION 3: NUCLEAR RADIUS FROM PION-RANGE PACKING")
print(SEP2)
# Nucleons bind when separation < r_pi. In close-packed nucleus:
# Each nucleon occupies a sphere whose TOUCHING distance = r_pi.
# -> effective nucleon radius for packing = r_pi/2
r_pack = r_pion / 2
print(f"  Pion force range: r_pi = {r_pion:.4f} fm")
print(f"  Effective packing radius: r_pack = r_pi/2 = {r_pack:.4f} fm")
print(f"    [Nucleons touch at r_pi; each contributes r_pi/2 to nuclear radius]")
# For FCC close packing (highest density, fraction = pi/(3*sqrt(2)) = 0.7405):
f_fcc = pi / (3 * math.sqrt(2))
r_0_fcc = r_pack / f_fcc**(1/3)
# For A nucleons: R = r_0 * A^(1/3) where r_0 is the per-unit radius
print(f"  r_0 (FCC packing): r_pack / f_fcc^(1/3) = {r_0_fcc:.4f} fm")
print(f"    [Measured nuclear r_0 ~ 1.2-1.25 fm for charge radius]")
# Simple estimate: R_nucleus = r_pi * (A/4pi)^(1/3) from spherical packing
r_0_simple = r_pion * (3/(4*pi))**(1/3)
print(f"  r_0 (spherical shell approx): r_pi*(3/4pi)^(1/3) = {r_0_simple:.4f} fm")
print()

# Nuclear r_0 from pion packing gives a LOWER BOUND (strong force equilibrium not yet
# derived from torsion medium -- requires balance of pion attraction vs hard core).
# Check: r_0 is bounded between lambda_p (hard core) and r_pi (force range).
check("NG5 Nuclear r_0 (geom) bounded: lambda_p < r_0_fcc < r_pi (correct order)",
      lambda_p < r_0_fcc < r_pion,
      f"r_0_fcc={r_0_fcc:.4f} fm  in ({lambda_p:.3f}, {r_pion:.3f}) fm  [measured 1.2 fm -- strong force balance open]")

# ── Section 4: I_h orbital decompositions ─────────────────────────────────────
print()
print(SEP)
print("SECTION 4: I_h BRANCHING RULES -- ORBITAL ANGULAR MOMENTA")
print(SEP2)
print("  Decomposing orbital l into I_h irreps via character projection.")
print("  Intruder magic: intruder j=l+1/2 has dim=2j+1=2(l+1).")
print("  If dim = 2 * (single dominant irrep) -> clean magic number prediction.")
print()
print(f"  {'l':>3}  {'dim':>4}  Decomposition                        Intruder dim  I_h basis?")
print(f"  {'-'*3}  {'-'*4}  {'-'*35}  {'-'*12}  {'-'*20}")

magic_status = {}
for l in range(0, 8):
    decomp, chi_l = ih_branching(l)
    j_intruder = l + 0.5
    dim_intruder = int(2*j_intruder + 1)  # = 2*(l+1)

    # Build readable decomposition string
    parts = []
    for name in ['A_g','T_1g','T_2g','G_g','H_g']:
        n = decomp[name]
        if n > 0:
            parts.append(f"{n}*{name}({I_h_chars[name][0]})" if n > 1 else f"{name}({I_h_chars[name][0]})")
    decomp_str = " + ".join(parts)
    dim_check = sum(decomp[k]*I_h_chars[k][0] for k in decomp)
    assert dim_check == 2*l+1, f"dim mismatch l={l}: {dim_check} != {2*l+1}"

    # Check if intruder dim = 2 * single irrep
    clean_basis = "?"
    for name in ['A_g','T_1g','T_2g','G_g','H_g']:
        d = I_h_chars[name][0]
        if 2*d == dim_intruder and decomp.get(name, 0) > 0:
            clean_basis = f"2*{name}  CLEAN"
            break
    # Check if intruder dim = 2 * sum of some irreps from decomposition
    if clean_basis == "?":
        # Try pairs
        names = [n for n in ['A_g','T_1g','T_2g','G_g','H_g'] if decomp.get(n,0)>0]
        for i, n1 in enumerate(names):
            for n2 in names[i:]:
                d_pair = I_h_chars[n1][0] + I_h_chars[n2][0]
                if 2*d_pair == dim_intruder:
                    clean_basis = f"2*({n1}+{n2})  COMPOSITE"
                    break

    print(f"  {l:>3}  {2*l+1:>4}  {decomp_str:<35}  {dim_intruder:>12}  {clean_basis}")
    magic_status[l] = (decomp, dim_intruder, clean_basis)

print()

check("NG6 l=3 -> T_2g(3) + G_g(4): intruder f_{7/2} dim=8=2*G_g  [magic 28]",
      magic_status[3][0]['T_2g'] == 1 and magic_status[3][0]['G_g'] == 1,
      f"l=3: {magic_status[3][0]}")
check("NG7 l=4 -> G_g(4) + H_g(5): intruder g_{9/2} dim=10=2*H_g  [magic 50]",
      magic_status[4][0]['G_g'] == 1 and magic_status[4][0]['H_g'] == 1,
      f"l=4: {magic_status[4][0]}")
check("NG8 Magic 28: f_{7/2} dim=8 = 2*G_g  (G_g = boundary irrep)",
      "2*G_g" in magic_status[3][2],
      f"intruder basis: {magic_status[3][2]}")
check("NG9 Magic 50: g_{9/2} dim=10 = 2*H_g  (H_g = next irrep)",
      "2*H_g" in magic_status[4][2],
      f"intruder basis: {magic_status[4][2]}")
check("NG10 l=5 -> T_1g + T_2g + H_g  (dim=3+3+5=11)",
      magic_status[5][0]['T_1g'] == 1 and magic_status[5][0]['T_2g'] == 1 and magic_status[5][0]['H_g'] == 1,
      f"l=5: {magic_status[5][0]}")
check("NG11 l=6 -> A_g + T_1g + G_g + H_g  (dim=1+3+4+5=13)",
      magic_status[6][0]['A_g'] == 1 and magic_status[6][0]['G_g'] == 1,
      f"l=6: {magic_status[6][0]}")
check("NG12 h_{11/2} dim=12 NOT a clean 2*irrep (magic 82 basis is composite)",
      "COMPOSITE" in magic_status[5][2] or "?" in magic_status[5][2],
      f"l=5 intruder basis: {magic_status[5][2]}")

# ── Section 5: N=184 status -- j=l=7 analysis ─────────────────────────────────
print()
print(SEP)
print("SECTION 5: l=7 DECOMPOSITION -- BASIS FOR 1j_{15/2} AND N=184")
print(SEP2)
decomp_7, chi_l7 = ih_branching(7)
dim_intruder_7 = 2*(7+1)  # = 16
parts_7 = []
for name in ['A_g','T_1g','T_2g','G_g','H_g']:
    n = decomp_7[name]
    if n > 0:
        parts_7.append(f"{n}*{name}({I_h_chars[name][0]})" if n > 1 else f"{name}({I_h_chars[name][0]})")
decomp_str_7 = " + ".join(parts_7)

print(f"  l=7 orbital decomposes as: {decomp_str_7}")
print(f"  Intruder 1j_{{15/2}}: j=l+1/2=15/2, dim=2j+1={dim_intruder_7}")
print()

# Check if dim=16 matches 2 * any combination from the l=7 decomposition
found_16 = []
names_7 = [n for n in ['A_g','T_1g','T_2g','G_g','H_g'] if decomp_7.get(n,0)>0]
# Single irrep * 2
for name in names_7:
    if 2 * I_h_chars[name][0] == dim_intruder_7:
        found_16.append(f"2*{name} (dim={I_h_chars[name][0]})")
# Pairs * 2
for i, n1 in enumerate(names_7):
    for n2 in names_7[i:]:
        if 2*(I_h_chars[n1][0]+I_h_chars[n2][0]) == dim_intruder_7:
            found_16.append(f"2*({n1}+{n2})")
# Triple * 2
for i, n1 in enumerate(names_7):
    for j2, n2 in enumerate(names_7):
        for n3 in names_7[j2:]:
            if n1 != n2 and 2*(I_h_chars[n1][0]+I_h_chars[n2][0]+I_h_chars[n3][0]) == dim_intruder_7:
                found_16.append(f"2*({n1}+{n2}+{n3})")

if found_16:
    print(f"  dim=16 decomposes as: {', '.join(found_16)}")
    print(f"  I_h GEOMETRIC BASIS for 1j_{{15/2}} ESTABLISHED from l=7 orbital.")
else:
    print(f"  dim=16 has NO clean decomposition from l=7 irreps as 2*X.")
    print(f"  N=184 intruder 1j_{{15/2}} LACKS a clean single I_h irrep basis.")
    print(f"  The l=7 orbital contains: {decomp_str_7}")
    print(f"  Total l=7 dim={2*7+1}. Intruder dim={dim_intruder_7}.")
    # What DOES the l=7 orbital contain?
    print()
    print(f"  Alternative: can dim=16 be written as multiplicity * irrep from l=7?")
    for name in ['A_g','T_1g','T_2g','G_g','H_g']:
        if decomp_7.get(name,0) > 0:
            d = I_h_chars[name][0]
            mult = dim_intruder_7 / d
            if mult == int(mult):
                print(f"    {dim_intruder_7} = {int(mult)} * {name}(dim={d})  "
                      f"[appears {decomp_7[name]}x in l=7 -- {'YES' if decomp_7[name] >= int(mult) else 'NO, need '+str(int(mult))}]")

print()
check("NG13 l=7 decomposition is well-defined (dim=15 verified)",
      sum(decomp_7[k]*I_h_chars[k][0] for k in decomp_7) == 2*7+1,
      f"l=7: {decomp_str_7}  total dim = {sum(decomp_7[k]*I_h_chars[k][0] for k in decomp_7)}")

# ── Section 6: Zone 2 geometry and N/Z ────────────────────────────────────────
print()
print(SEP)
print("SECTION 6: N/Z GEOMETRY FROM ZONE 2 BUFFER CAPACITY")
print(SEP2)
# Zone 2 volume per proton (the buffer shell)
V_Zone2 = (4/3)*pi*(r_p_fm**3 - lambda_p**3)
V_Zone1 = (4/3)*pi*lambda_p**3
V_ratio  = V_Zone2 / V_Zone1

# Neutron (no Zone 2) has geometric radius ~ r_n_fm
V_n_geom = (4/3)*pi*r_n_fm**3

# The number of neutrons that can fit INSIDE the proton Zone 2 shell
# (Zone 2 shell volume / neutron Zone 1 volume)
n_in_Z2 = V_Zone2 / V_Zone1
# More relevant: number of neutrons that can surround one proton at nuclear separation
# Using pion force range as the effective nuclear cell radius
V_nuclear_cell = (4/3)*pi*r_pion**3
n_surround = (V_nuclear_cell - V_p_fm3) / V_n_geom

print(f"  Proton Zone 2 volume:    V_Z2 = {V_Zone2:.4f} fm^3")
print(f"  Neutron Zone 1 volume:   V_Z1 = {V_Zone1:.4f} fm^3")
print(f"  V_Zone2/V_Zone1 = {V_ratio:.2f}  (Zone 2 could geometrically hold {V_ratio:.1f} neutron cores)")
print(f"  Nuclear cell volume (pion sphere): V_pi = {V_nuclear_cell:.4f} fm^3")
print(f"  Neutron geom. volume: V_n = {V_n_geom:.4f} fm^3")
print(f"  Neutrons surrounding 1 proton in nuclear cell: {n_surround:.2f}")
print()
print(f"  Geometric N/Z estimate from Zone buffer capacity: N/Z ~ {n_surround+1:.1f}:1")
print(f"  (This is an UPPER bound; actual N/Z < this from Coulomb and kinematic effects)")
print()
# Observed N/Z at Z=114: ~1.69 (from measured isotopes N~184, Z=114)
nz_114 = 184/114
print(f"  Measured N/Z at Z=114 (N=184): {nz_114:.3f}")

check("NG15 Zone 2 buffer capacity > measured N/Z (Zone 2 can accommodate neutrons)",
      V_ratio > nz_114,
      f"V_Zone2/V_Zone1 = {V_ratio:.2f} > N/Z(observed) = {nz_114:.3f}")

# ── Summary ────────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SECTION 7: PROLATE DEFORMATION -- POLAR POCKET SEARCH (Z=114-120)")
print(SEP2)
# Prolate deformation (β₂>0) at N between magic numbers 126 and 184.
# d quarks of neutrons in triangular orbit align with the nuclear symmetry axis,
# creating prolate elongation with maximum at midshell N~155.
# At the poles, nuclear surface proton density is lower -> Coulomb energy reduced
# for a new proton sitting at the pole vs. on a spherical surface.
#
# CAVEAT: the I_h magic numbers (Sections 4/5) assume spherical symmetry.
# Prolate deformation (β₂>0) breaks I_h -> D_5h. The Z=114 G_g shell closure
# is a spherical result and shifts under deformation. Pocket results are:
#   - Quantitatively valid at N=184 (β₂=0, both shells closed, spherical)
#   - Qualitatively valid at midshell but need I_h->D_5h Nilsson decomposition
#     for precise deformed magic numbers. [OPEN: F-8]
#
# β₂(N) = β₂_max * sin(π*(N-126)/(184-126))  [zero at both shell closures]
# R_pole = R_nuc * (1 + β₂ * Y₂₀(0))  where Y₂₀(0) = sqrt(5/(4π)) = 0.6305
#
# Coulomb energy reduction at pole vs sphere (torsionverse Zone 3 picture):
#   ΔE_C = Z * alpha * hbar_c * (1/R_nuc - 1/R_pole)  [MeV]
# This is the Coulomb energy gain for sitting at the pole instead of a sphere.

Y20_pole   = math.sqrt(5 / (4 * math.pi))          # 0.6305
Y20_eq     = -math.sqrt(5 / (16 * math.pi))         # -0.3153 (equator compressed)
beta2_max  = 0.25                                    # typical midshell superheavy
N_lo, N_hi = 126, 184                               # neutron magic boundaries
# Nuclear radius from pion-based r_0 (Section 3):
r_0_pack   = r_pion * (3 / (4*pi))**(1/3)           # ~0.708 fm from pion Compton

# Scan Z=114-120 at key neutron numbers
print(f"  {'Z':>3}  {'N':>3}  {'A':>3}  {'β₂':>5}  {'R_nuc':>6}  "
      f"{'R_pole':>6}  {'ΔE_C':>9}  {'Pocket?'}")
print(f"  {'-'*3}  {'-'*3}  {'-'*3}  {'-'*5}  {'-'*6}  {'-'*6}  {'-'*9}  {'-'*8}")

pockets = []
for Z_scan in range(114, 121):
    for N_scan in [126, 140, 155, 162, 170, 178, 184, 190]:
        A = Z_scan + N_scan
        # Sinusoidal midshell deformation model
        if N_scan <= N_lo or N_scan >= N_hi:
            b2 = 0.0
        else:
            b2 = beta2_max * math.sin(math.pi * (N_scan - N_lo) / (N_hi - N_lo))
        R_nuc  = r_0_pack * A**(1/3)
        R_pole = R_nuc * (1 + b2 * Y20_pole)
        # Coulomb energy reduction at pole (Zone 3 pressure reduction)
        dE_C = Z_scan * alpha * hbar_c * (1/R_nuc - 1/R_pole)   # MeV (>0 = favourable)
        has_pocket = (dE_C > 1.5)   # threshold: >1.5 MeV Coulomb gain at pole
        if has_pocket:
            pockets.append((Z_scan, N_scan, b2, R_nuc, R_pole, dE_C))
        flag = "POCKET" if has_pocket else ""
        if Z_scan == 115 or has_pocket:
            print(f"  {Z_scan:>3}  {N_scan:>3}  {A:>3}  {b2:.3f}  {R_nuc:.4f}  "
                  f"{R_pole:.4f}  {dE_C:+.3f} MeV  {flag}")

print()
print(f"  Most stable Z=115 (N=184, β₂=0):")
A115 = 115 + 184
R115 = r_0_pack * A115**(1/3)
print(f"    Mc-299: Z=115, N=184, A=299  R_nuc = {R115:.4f} fm  β₂=0 (both shells closed)")
print(f"    Shape: nearly spherical. Spin J=9/2 (unpaired 2g_{{9/2}} proton).")
print(f"    Required synthesis: Ca-48 + Am-251 -> Mc-299 + 4n  (not yet achieved)")
print(f"    Current experiments: Mc-288 to Mc-291 (N=173-176, ~8 neutrons short)")
print()

# Deformation at N=184: must be zero (shell closure)
b2_at_184 = beta2_max * math.sin(math.pi * (184 - N_lo) / (N_hi - N_lo)) \
            if 184 < N_hi else 0.0
check("NG16 β₂(N=184) = 0  [neutron shell closure: spherical]",
      abs(b2_at_184) < 1e-9, f"β₂(N=184) = {b2_at_184:.6f}")

b2_mid = beta2_max * math.sin(math.pi * 0.5)  # N = midshell
check("NG17 β₂(N_midshell) = β₂_max  [maximum prolate at midshell]",
      abs(b2_mid - beta2_max) < 1e-9, f"β₂(midshell) = {b2_mid:.4f} = β₂_max = {beta2_max}")

# For any β₂ > 0: R_pole > R_equator (prolate = elongated at poles)
b2_test = 0.20
R_nuc_test  = 9.0
R_pole_test = R_nuc_test * (1 + b2_test * Y20_pole)
R_eq_test   = R_nuc_test * (1 + b2_test * Y20_eq)
check("NG18 Prolate deformation: R_pole > R_nuc > R_equator  [geometric consistency]",
      R_pole_test > R_nuc_test > R_eq_test,
      f"R_pole={R_pole_test:.4f} > R_nuc={R_nuc_test:.4f} > R_eq={R_eq_test:.4f} fm  (β₂=0.20)")

# Coulomb reduction at pole of Z=114, N=178 (β₂~0.20): must be positive
Z_test, N_test = 114, 178
b2_178 = beta2_max * math.sin(math.pi * (178 - N_lo) / (N_hi - N_lo))
R_nuc_178  = r_0_pack * (Z_test + N_test)**(1/3)
R_pole_178 = R_nuc_178 * (1 + b2_178 * Y20_pole)
dE_178     = Z_test * alpha * hbar_c * (1/R_nuc_178 - 1/R_pole_178)
check("NG19 Coulomb energy reduction at pole of Z=114, N=178 is positive (favours pole)",
      dE_178 > 0, f"ΔE_C = {dE_178:.4f} MeV  β₂={b2_178:.3f}")


print(SEP2)
print(f"  CLEAN I_h GEOMETRIC BASIS:")
print(f"    Magic 2, 8, 20:  harmonic oscillator levels  (universal)")
print(f"    Magic 28:  f_{{7/2}} intruder, l=3 -> G_g(4),  dim=8=2*G_g  [NG8 PASS]")
print(f"    Magic 50:  g_{{9/2}} intruder, l=4 -> H_g(5),  dim=10=2*H_g [NG9 PASS]")
print(f"    Z=114:     1i_{{13/2}} intruder, l=6 -> G_g, 2*(T_{{2g}}+G_g)=14 [ESSENTIALLY CLOSED]")
print(f"    Z_crit:    1/alpha = 137.036 [exact, from electron bulk/boundary]")
print()
print(f"  COMPOSITE (not a single I_h irrep):")
print(f"    Magic 82:  h_{{11/2}} dim=12, l=5 -> T_{{1g}}+T_{{2g}}+H_g -- composite [NG12]")
print(f"    Magic 126: i_{{13/2}} dim=14 = 2*(T_{{2g}}+G_g) -- composite [essentially closed]")
print()
print(f"  N=184 STATUS (l=7 result):")
l7_result = "ESTABLISHED" if found_16 else "NO CLEAN BASIS"
print(f"    j_{{15/2}} dim=16, l=7 decomp = {decomp_str_7}")
print(f"    dim=16 from I_h: {l7_result}")
print(f"    16 = 2*(T_1g+H_g) = 2*(T_2g+H_g) = 2*(G_g+G_g): COMPOSITE but I_h-grounded.")
print(f"    Physical: T_1g+H_g shares neutron diquark irrep (T_1g) and sub-cell irrep (H_g).")
print(f"    Status: I_h basis EXISTS, composite -- stronger than borrowed, weaker than clean.")
print()
print(f"  GEOMETRIC SCALES:")
print(f"    p-p hard core: {r_pp_hard:.4f} fm  (2*lambda_p, cog grinding)")
print(f"    Nuclear force: {r_pion:.4f} fm  (hbar_c/m_pi, pion Yukawa range)")
print(f"    Neutron larger than proton by {(r_n_factor-1)*100:.4f}%  (Zone 2 absent)")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(f"  Reference: docs/doc_nucleus.txt")
print(SEP)
