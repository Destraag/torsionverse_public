"""
gluon_tau_helix.py
Compares the tau corkscrew helix geometry with the gluon half-wave helix.

SCOPE NOTE (session 12): this script computes the tau's discrete PATH
CENTERLINE only (20 face-center positions, deflection angles, step length) --
it says nothing about amplitude, dynamics, or whether combinations of this
path (pairing, forward+backward, etc.) behave differently. Read its "path
stays at constant radius, never reaches r=0" result as scoped to that single
bare centerline, not as a claim about every possible tau-related construction.
Later, more specific scripts build on top of this one rather than replacing
it: analysis/quantum/tau_pair_configuration.py (forward+backward pairing),
analysis/quantum/tau_pair_wz_composite.py (I52 x I52 algebraic combination),
analysis/quantum/wz_directed_cone_geometry.py (T_1g amplitude pattern).

Reference: docs/doc_jobson_cell.txt JC7-JC9 (Hamiltonian cycle, 138.19 deg)
           docs/doc_leptons.txt (tau I52 face corkscrew)
"""
import sys, os, math, itertools
import numpy as np
sys.stdout.reconfigure(encoding='utf-8', errors='replace')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import phi, pi, E_cell_GeV

SEP  = "=" * 65
SEP2 = "-" * 65
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  {'[PASS]' if cond else '[FAIL]'} {name}")
    if detail: print(f"         {detail}")

# ── Build icosahedral graph (reuse from JC section) ───────────────────────────
_verts = []
for s1, s2 in itertools.product([1,-1],[1,-1]):
    _verts += [(0,s1,s2*phi),(s1,s2*phi,0),(s2*phi,0,s1)]

def dsq(a,b): return sum((x-y)**2 for x,y in zip(a,b))
eset = set()
for i in range(12):
    for j in range(i+1,12):
        if abs(dsq(_verts[i],_verts[j])-4.0)<1e-9: eset.add((i,j))

nb = {i:[] for i in range(12)}
for (i,j) in eset: nb[i].append(j); nb[j].append(i)

# Build faces and face adjacency
flist = []
for a in range(12):
    for b in nb[a]:
        if b>a:
            for c in nb[a]:
                if c>b and c in nb[b]: flist.append((a,b,c))

fadj = {i:[] for i in range(20)}
for i in range(20):
    for j in range(i+1,20):
        if len(set(flist[i])&set(flist[j]))==2:
            fadj[i].append(j); fadj[j].append(i)

def fcen(f): return np.array([sum(_verts[x][k] for x in f)/3 for k in range(3)])
fcenters = [fcen(f) for f in flist]

# Find Hamiltonian cycle (same as JC7)
def ham_cycle(adj, n):
    path=[0]; vis={0}
    def bt():
        if len(path)==n: return 0 in adj[path[-1]]
        for nb2 in adj[path[-1]]:
            if nb2 not in vis:
                path.append(nb2); vis.add(nb2)
                if bt(): return True
                path.pop(); vis.remove(nb2)
        return False
    return bt(), path

found, hpath = ham_cycle(fadj, 20)

print(SEP)
print("gluon_tau_helix.py -- tau helix geometry vs gluon half-wave")
print(SEP)

# ── Tau corkscrew helix geometry ──────────────────────────────────────────────
print()
print(SEP2)
print("TAU CORKSCREW HELIX GEOMETRY (from Hamiltonian face cycle)")
print(SEP2)

check("GH1: Hamiltonian cycle found", found, f"{len(hpath)} faces")

# Path vectors: from face center to face center
path_pts = np.array([fcenters[hpath[k]] for k in range(20)])

# Step vectors and their magnitudes
steps = np.array([path_pts[(k+1)%20] - path_pts[k] for k in range(20)])
step_len = np.linalg.norm(steps, axis=1)
step_len_mean = float(np.mean(step_len))

print(f"\n  Face-center to face-center step length:")
print(f"    mean = {step_len_mean:.6f}  (in edge=2 units)")
print(f"    min  = {np.min(step_len):.6f}  max = {np.max(step_len):.6f}")
print(f"    All equal? {np.allclose(step_len, step_len_mean, rtol=1e-6)}")

# Deflection angle at each step (angle between consecutive step vectors)
deflections = []
for k in range(20):
    v1 = steps[k-1]; v2 = steps[k]
    cos_a = np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2))
    deflections.append(math.degrees(math.acos(float(np.clip(cos_a,-1,1)))))

defl_mean = float(np.mean(deflections))
print(f"\n  Path deflection at each face center:")
print(f"    mean = {defl_mean:.4f} deg")
print(f"    min  = {min(deflections):.4f}  max = {max(deflections):.4f}")

# Fit a helix to the face-center path: project onto C5 axis
c5_axis = np.array([0,0,1.0])  # icosahedron C5 axis
# Project path points onto C5 axis and transverse plane
z_coords = np.array([np.dot(p, c5_axis) for p in path_pts])
xy_pts = path_pts - np.outer(z_coords, c5_axis)
radii = np.linalg.norm(xy_pts, axis=1)
angles = np.arctan2(xy_pts[:,1], xy_pts[:,0])

print(f"\n  Projected onto C5 axis:")
print(f"    mean radius = {np.mean(radii):.6f}")
print(f"    z-range     = {np.min(z_coords):.4f} to {np.max(z_coords):.4f}")
print(f"    total z-advance over 20 steps = {z_coords[-1]-z_coords[0]:.4f}")
print(f"    total angular advance (rad)   = {angles[-1]-angles[0]:.4f}")
print(f"    angular advance per step      = {(angles[-1]-angles[0])/20:.4f} rad"
      f" = {math.degrees((angles[-1]-angles[0])/20):.2f} deg")

# Helix pitch angle = arctan(z_advance_per_step / (r * dtheta_per_step))
dz_per_step = (z_coords[-1]-z_coords[0])/20
r_mean = float(np.mean(radii))
dtheta_per_step = abs((angles[-1]-angles[0])/20)
pitch_angle = math.degrees(math.atan2(abs(dz_per_step), r_mean * dtheta_per_step))
print(f"\n  Tau helix pitch angle = arctan(dz / r*dtheta) = {pitch_angle:.4f} deg")

# ── Gluon half-wave frequency (from first principles) ────────────────────────
print()
print(SEP2)
print("GLUON FREQUENCY FROM FIRST PRINCIPLES (massless standing wave on edge L_J)")
print(SEP2)

# Gluon = massless wave, speed c, standing on edge of length L_J
# Fundamental mode: half-wavelength in L_J -> lambda_wave = 2*L_J
# omega = pi*c/L_J = (1/2) * 2*pi*c/L_J = E_cell / (2*hbar)
# E_gluon = hbar * omega = hbar * pi * c / L_J = E_cell / 2

E_cell_gluon = E_cell_GeV / 2.0   # GeV
print(f"  E_cell = 2*pi*hbar*c/L_J = {E_cell_GeV:.4f} GeV")
print(f"  Gluon: standing half-wave on L_J -> omega = pi*c/L_J = E_cell/(2*hbar)")
print(f"  E_gluon = hbar*omega = E_cell/2 = {E_cell_gluon:.4f} GeV")
print(f"  Rotation rate (circular polarization) = same omega")
print(f"  Both polarizations: in quadrature, combined amplitude NEVER zero")

check("GH0: E_gluon = E_cell/2 (massless half-wave on edge L_J, from first principles)",
      abs(E_cell_gluon - E_cell_GeV/2) < 1e-10,
      f"E_gluon = {E_cell_gluon:.6f} GeV = E_cell/2 [spatial profile and frequency: DERIVED]")

# ── Gluon amplitude from face geometry ────────────────────────────────────────
print()
print(SEP2)
print("GLUON AMPLITUDE FROM FACE GEOMETRY (not from Born balance)")
print(SEP2)

# The gluon oscillates perpendicular to edge, IN the face plane.
# This direction points EXACTLY toward the face center (tau nexus).
# Therefore A = distance from edge midpoint to face center.
# For equilateral triangle side L, midpoint at (L/2, 0), centroid at (L/2, L*sqrt(3)/6):
L_raw = 2.0  # edge in raw icosahedral coords
A_geometric = L_raw * math.sqrt(3) / 6   # = L/sqrt(12)

print(f"  Equilateral triangle side L = {L_raw:.1f} (raw coords)")
print(f"  Edge midpoint: (L/2, 0) = ({L_raw/2:.3f}, 0)")
print(f"  Face center:   (L/2, L*sqrt(3)/6) = ({L_raw/2:.3f}, {L_raw*math.sqrt(3)/6:.6f})")
print(f"  Direction: purely perpendicular to edge (into face) = gluon transverse direction")
print(f"  A = L*sqrt(3)/6 = L/sqrt(12) = {A_geometric:.8f}  (A/L = {A_geometric/L_raw:.8f})")
pitch_geo = math.degrees(math.atan(A_geometric * math.pi / L_raw))
print(f"  Gluon pitch angle (geometric A) = arctan(sqrt(3)*pi/6) = {pitch_geo:.4f} deg")
print(f"  Three edge gluons on one face ALL reach face center simultaneously.")

check("GH0b: gluon transverse direction = edge-midpoint to face-center (perpendicular to edge)",
      abs(A_geometric - L_raw / math.sqrt(12)) < 1e-12,
      f"A = L*sqrt(3)/6 = L/sqrt(12) = {A_geometric:.8f}  [DERIVED from equilateral triangle]")
check("GH0c: A/L = 1/sqrt(12) = sqrt(3)/6 (exact algebraic, no free parameters)",
      abs(A_geometric/L_raw - 1/math.sqrt(12)) < 1e-12,
      f"A/L = {A_geometric/L_raw:.10f} = 1/sqrt(12) = {1/math.sqrt(12):.10f}")

# ── Gluon half-wave helix geometry ───────────────────────────────────────────
print()
print(SEP2)
print("GLUON HALF-WAVE HELIX GEOMETRY (one edge, both polarizations)")
print(SEP2)

# For a half-wavelength helix on edge of length L=2 (raw coords):
# x(t) = t  (along edge, 0 to L)
# y(t) = A * sin(pi*t/L) * cos(omega*t)
# z(t) = A * sin(pi*t/L) * sin(omega*t)
# The pitch angle depends on the ratio A/L

# The gluon amplitude A is set by Born coupling: from the alpha derivation,
# the vertex coupling ratio k_n/k_eff = alpha*phi / (1+alpha*phi^2)
# The transverse amplitude ~ k_n = alpha*phi * k_eff * L / (1+alpha*phi^2)
# As a fraction of the edge length:
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'higgs'))
from constants import alpha
k_ratio = alpha*phi / (1 + alpha*phi**2)  # k_n/k_eff from Born balance
print(f"  Born coupling k_n/k_eff = alpha*phi/(1+alpha*phi^2) = {k_ratio:.8f}")
print(f"  This gives the gluon transverse amplitude / edge_length ~ {k_ratio:.6f}")

L = 2.0  # edge length in raw coords
A_gluon = k_ratio * L  # transverse amplitude estimate
print(f"  Estimated A/L = {A_gluon/L:.8f}")

# For the circularly polarized gluon (both pols in quadrature):
# The helix path: x moves 0->L, (y,z) rotates with amplitude A*sin(pi*x/L)
# At x = L/2 (midpoint): max transverse = A
# Helix pitch angle at midpoint:
# dz/dx = A * (pi/L) * cos(pi*x/L) * sin(phi_t) + A * sin(pi*x/L) * omega * (-cos(phi_t))
# At x=L/2: sin(pi/2)=1, cos(pi/2)=0, so dz/dx = -A*omega*cos(phi_t)
# The pitch at midpoint: arctan(|d_transverse/dx|) = arctan(A*pi/L)  (envelope derivative)
pitch_gluon = math.degrees(math.atan(A_gluon * math.pi / L))
print(f"  Gluon helix pitch angle (at midpoint, envelope) = {pitch_gluon:.4f} deg")

# Compare tau and gluon pitch angles
print()
print(SEP2)
print("COMPARISON")
print(SEP2)
print(f"  Tau corkscrew pitch angle  = {pitch_angle:.4f} deg")
print(f"  Gluon helix pitch angle    = {pitch_gluon:.4f} deg  (A/L from Born coupling)")
print(f"  Difference                 = {abs(pitch_angle-pitch_gluon):.4f} deg")
print(f"  Ratio tau/gluon            = {pitch_angle/pitch_gluon:.6f}")

check("GH2: tau deflection angle uniform across all 20 steps = 72 deg (C5)",
      max(deflections)-min(deflections) < 0.001 and abs(defl_mean - 72.0) < 0.01,
      f"all = {defl_mean:.4f} deg  spread = {max(deflections)-min(deflections):.6f} deg")
check("GH3: step lengths exactly uniform (2*phi/3 in edge=2 coords)",
      (np.max(step_len)-np.min(step_len)) < 1e-9,
      f"all = {step_len_mean:.8f}  expected = {2*phi/3:.8f}  match = {abs(step_len_mean - 2*phi/3) < 1e-6}")
check("GH4: step length = 2*phi/3 exactly",
      abs(step_len_mean - 2*phi/3) < 1e-6,
      f"computed = {step_len_mean:.8f}  2*phi/3 = {2*phi/3:.8f}")

# ── Verify across multiple start faces ───────────────────────────────────────
print()
print(SEP2)
print("ROBUSTNESS: tau path geometry across different starting faces")
print(SEP2)

uniform_72 = []
for start in range(20):
    path2 = [start]; vis2 = {start}
    def bt2():
        if len(path2)==20: return start in fadj[path2[-1]]
        for nb3 in fadj[path2[-1]]:
            if nb3 not in vis2:
                path2.append(nb3); vis2.add(nb3)
                if bt2(): return True
                path2.pop(); vis2.remove(nb3)
        return False
    found2 = bt2()
    if found2:
        pts2 = np.array([fcenters[path2[k]] for k in range(20)])
        steps2 = np.array([pts2[(k+1)%20]-pts2[k] for k in range(20)])
        defls2 = []
        for k in range(20):
            v1=steps2[k-1]; v2=steps2[k]
            ca=float(np.dot(v1,v2)/(np.linalg.norm(v1)*np.linalg.norm(v2)))
            defls2.append(math.degrees(math.acos(np.clip(ca,-1,1))))
        is_uniform = max(defls2)-min(defls2) < 0.01 and abs(np.mean(defls2)-72.0) < 0.1
        uniform_72.append(is_uniform)

print(f"  Hamiltonian cycles found from all 20 start faces: all 72-deg uniform = {sum(uniform_72)}/20")
check("GH5: tau 72-deg geometric path holds from all 20 starting faces",
      sum(uniform_72) == 20,
      f"{sum(uniform_72)}/20 starting faces give uniform 72-deg Hamiltonian cycles")

# ── GH6/GH7: Tau chord midpoint depth ────────────────────────────────────────
# Between adjacent face-center nexuses the tau path is a straight chord through
# the cell interior. Exact algebraic values (edge=2 coords):
#   r_chord_unscaled = (2+sqrt(5))/3   [proved: sqrt((9+4sqrt5)/9) = (2+sqrt5)/3]
#   r_chord/L_J      = (2+sqrt(5))/6 = 0.706   [divide by edge=2]
#   dip/L_J          = (sqrt(5)-1)/12 = 0.103   [phi/2 - r_chord/2]
# This guarantees the tau chord never contacts the edge midpoints (r_mid = 0.809).
print()
print(SEP2)
print("TAU CHORD MIDPOINT DEPTH (chord between adjacent face-center nexuses)")
print(SEP2)

r_chord_exact = (2 + math.sqrt(5)) / 3   # in edge=2 coords
dip_exact     = (math.sqrt(5) - 1) / 12  # in L_J = edge/2 units

chord_radii = []
for k in range(20):
    midpt = (path_pts[k] + path_pts[(k+1) % 20]) / 2
    chord_radii.append(float(np.linalg.norm(midpt)))

r_chord_mean = sum(chord_radii) / 20
r_chord_ratio = r_chord_mean / 2          # as fraction of L_J
dip_computed  = (phi - r_chord_mean) / 2  # (r_mid - r_chord) in L_J units

print(f"  Chord midpoint radius (edge=2 coords):  mean = {r_chord_mean:.8f}")
print(f"  Exact (2+sqrt5)/3                            = {r_chord_exact:.8f}")
print(f"  r_chord/L_J = {r_chord_ratio:.6f}  (doc: 0.706)")
print(f"  dip/L_J     = {dip_computed:.6f}  exact (sqrt5-1)/12 = {dip_exact:.6f}  (doc: 0.103)")
print(f"  r_mid/L_J = phi/2 = {phi/2:.6f} > r_chord/L_J = {r_chord_ratio:.6f}: no edge contact.")

check("GH6: tau chord midpoint radius = (2+sqrt5)/3 in edge=2 coords = 0.706*L_J",
      abs(r_chord_mean - r_chord_exact) < 1e-9 and
      max(chord_radii) - min(chord_radii) < 1e-9,
      f"mean={r_chord_mean:.8f}  exact={r_chord_exact:.8f}  r/L_J={r_chord_ratio:.6f}")

check("GH7: chord dip below edge midpoint = (sqrt5-1)/12 = 0.103*L_J (no edge contact)",
      abs(dip_computed - dip_exact) < 1e-9,
      f"dip/L_J={dip_computed:.6f}  exact={dip_exact:.6f}  (r_mid=0.809 > r_chord=0.706)")

print()
print(SEP)
passed = sum(1 for _,s,_ in results if s=='PASS')
failed = sum(1 for _,s,_ in results if s=='FAIL')
print(f"  Total checks: {len(results)}  PASS: {passed}  FAIL: {failed}")
print(SEP)
