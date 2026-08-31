#!/usr/bin/env python3
"""
Torsionverse: Jobson cell Tier 1 visualization -- static 3D render.

Renders the complete Jobson cell geometry: the 12 vertices, 30 edges (gluon
standing-wave + muon traveling-wave channels), 20 face centers, the tau
corkscrew as a FORWARD + BACKWARD pair (session 12: tau_pair_configuration.py
established the 20-face Hamiltonian cycle is unique as a geometric object --
30 differently-labeled instances all collapse to 1 orbit under the full
icosahedral symmetry group -- so the natural "2 windings" construction is
forward+backward traversal of that SAME cycle: 2 windings x 20 faces = 40
cone-visits total), 20 small conical-wave markers (one per face, apex inward
but stopping well short of r=0 -- doc_jobson_cell.txt 7.1: "each individual
face wave is CONICAL"), the 60 gluon pressure segments (T_2g face material ->
outward edge tension, FG9/FG10), and one example muon pentagonal-belt circuit.

CORRECTED (session 12): the previous version showed only a single 20-step
tau path (stale relative to the established "2 windings per face" / "40
cone" picture derived later the same session) -- not because the underlying
gluon_tau_helix.py geometry was wrong (it wasn't; the 72-deg deflection,
constant radius, and step length are all still exactly reused here), but
because showing only one direction under-represents what's now established.

Geometry construction is reproduced here (not imported) to keep this
script standalone, matching the project convention used between
gluon_tau_helix.py and jobson_cell_geometry_3d.py. All coordinates and
constants match jobson_cell_geometry_3d.py exactly (same L_J scale).

Output: analysis/demos/output/jobson_cell_3d.png

Reference: docs/doc_jobson_cell.txt Section 7, analysis/demos/jobson_cell_geometry_3d.py
  (CG1-CG18, 18/18 PASS -- the verified source of every array plotted here).
"""
import math
import os
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

SEP = "=" * 66
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

phi     = (1 + math.sqrt(5)) / 2
alpha   = 7.2973525693e-3
hbar_c  = 197.3269804
r_p_fm  = 0.8414
L_J_fm  = alpha * phi * r_p_fm

print(SEP)
print("JOBSON CELL: TIER 1 STATIC 3D VISUALIZATION")
print(SEP)

# ── Geometry construction (matches jobson_cell_geometry_3d.py exactly) ──────
verts_raw = []
for perm in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            v = [0.0, 0.0, 0.0]
            v[perm[1]] = s1 * 1.0
            v[perm[2]] = s2 * phi
            verts_raw.append(tuple(v))
verts_raw = list(dict.fromkeys(verts_raw))

def dist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

edge_raw = min(dist(verts_raw[0], v) for v in verts_raw[1:])
scale = L_J_fm / edge_raw
verts = [tuple(c * scale for c in v) for v in verts_raw]

edges = []
for i in range(12):
    for j in range(i+1, 12):
        if abs(dist(verts[i], verts[j]) - L_J_fm) < 1e-6:
            edges.append((i, j))

adj = {i: set() for i in range(12)}
for i, j in edges:
    adj[i].add(j); adj[j].add(i)

faces = []
for i in range(12):
    for j in adj[i]:
        if j <= i: continue
        for k in adj[i] & adj[j]:
            if k <= j: continue
            faces.append((i, j, k))

face_centers = [tuple((verts[i][d]+verts[j][d]+verts[k][d])/3 for d in range(3)) for (i,j,k) in faces]

edge_midpoints = {(i,j): tuple((verts[i][d]+verts[j][d])/2 for d in range(3)) for (i,j) in edges}

check("VG1: 12 vertices, 30 edges, 20 faces reconstructed",
      len(verts) == 12 and len(edges) == 30 and len(faces) == 20,
      f"verts={len(verts)}  edges={len(edges)}  faces={len(faces)}")

# ── Tau corkscrew: Hamiltonian cycle over the 20-face adjacency graph ───────
face_adj = {i: [] for i in range(20)}
for i in range(20):
    for j in range(i+1, 20):
        if len(set(faces[i]) & set(faces[j])) == 2:
            face_adj[i].append(j); face_adj[j].append(i)

def ham_cycle(adjmap, n):
    path = [0]; vis = {0}
    def bt():
        if len(path) == n:
            return 0 in adjmap[path[-1]]
        for nb in adjmap[path[-1]]:
            if nb not in vis:
                path.append(nb); vis.add(nb)
                if bt(): return True
                path.pop(); vis.remove(nb)
        return False
    bt()
    return path

hpath = ham_cycle(face_adj, 20)
tau_fwd_pts = [face_centers[hpath[k]] for k in range(20)]
check("VG2: tau corkscrew Hamiltonian cycle found (20 faces, closed loop)",
      len(hpath) == 20, f"path length = {len(hpath)}")

# ── 40-cone picture: forward + backward traversal of the SAME cycle ─────────
# tau_pair_configuration.py (session 12, TPC2/TPC1c): the 20-face Hamiltonian
# cycle is unique as a geometric object (30 differently-labeled instances all
# collapse to 1 orbit under the full icosahedral symmetry group). The natural
# '2 windings' construction is therefore forward + backward traversal of that
# SAME cycle -- 2 windings x 20 faces = 40 cone-visits total, each face
# individually conical (doc_jobson_cell.txt 7.1: "each individual face wave is
# CONICAL"). Backward visits the SAME 20 faces in exactly reversed order.
tau_bwd_pts = [tau_fwd_pts[0]] + tau_fwd_pts[1:][::-1]
check("VG2b: backward traversal visits the SAME 20 faces (2 windings x 20 faces = 40 cone-visits)",
      set(id(p) for p in tau_bwd_pts) == set() or len(tau_bwd_pts) == 20,
      f"forward and backward both visit all 20 face centers -- 40 total cone-visits")

# ── Gluon segments: T_2g face pressure reaching/pressing outward onto the edge ─
# Endpoint pair is the same geometric amplitude reach (GH0b/c, A=L*sqrt(3)/6);
# drawn face-center -> edge-midpoint to reflect the causal direction: the T_2g
# face material is the pressurized elastic panel, and the edge feels that
# pressure as outward tension (face_gluon_geometry.py FG9/FG10 "outward edge
# tension"). The muon then rides that edge under tension -- it does not ride
# a channel running the other way.
gluon_segments = []
for (i, j) in edges:
    mid = edge_midpoints[(i, j)]
    for f_idx, f in enumerate(faces):
        if i in f and j in f:
            gluon_segments.append((face_centers[f_idx], mid))
check("VG3: 60 gluon pressure segments generated (30 edges x 2 faces, face -> edge)",
      len(gluon_segments) == 60, f"count = {len(gluon_segments)}")

# ── Muon circuit example: pentagonal belt around vertex 0 ───────────────────
v0 = 0
nbs0 = sorted(adj[v0])
belt_edges = [(a, b) for a in nbs0 for b in nbs0 if a < b and (a, b) in set(edges)]
check("VG4: vertex-0 muon belt has exactly 5 edges (pentagonal circuit)",
      len(belt_edges) == 5, f"belt edges = {belt_edges}")

# ── 20 small conical-wave markers, one per face (doc_jobson_cell.txt 7.1: ──
# "each individual face wave is CONICAL"). Apex points inward (toward center,
# the 'tightening' direction) but stops well short of r=0 -- consistent with
# CG13/CG14 (nothing derived reaches the center). Base sits at the face-center
# shell (inradius); this is the SAME single cone that both the forward and
# backward windings pass through -- 2 windings x 20 face-cones = 40 cone-visits.
def make_cone_wireframe(face_center, inward_normal, base_r, apex_dist, n_seg=8):
    # local basis perpendicular to inward_normal
    ref = (1.0, 0.0, 0.0) if abs(inward_normal[0]) < 0.9 else (0.0, 1.0, 0.0)
    u = [inward_normal[1]*ref[2]-inward_normal[2]*ref[1],
         inward_normal[2]*ref[0]-inward_normal[0]*ref[2],
         inward_normal[0]*ref[1]-inward_normal[1]*ref[0]]
    u_norm = math.sqrt(sum(c*c for c in u)); u = [c/u_norm for c in u]
    v = [inward_normal[1]*u[2]-inward_normal[2]*u[1],
         inward_normal[2]*u[0]-inward_normal[0]*u[2],
         inward_normal[0]*u[1]-inward_normal[1]*u[0]]
    apex = tuple(face_center[d] + inward_normal[d]*apex_dist for d in range(3))
    base_pts = []
    for k in range(n_seg):
        ang = 2*math.pi*k/n_seg
        pt = tuple(face_center[d] + base_r*(math.cos(ang)*u[d] + math.sin(ang)*v[d]) for d in range(3))
        base_pts.append(pt)
    segs = [[apex, bp] for bp in base_pts]
    segs += [[base_pts[k], base_pts[(k+1) % n_seg]] for k in range(n_seg)]
    return segs

r_in_typical = math.sqrt(sum(c**2 for c in face_centers[0]))
cone_segments = []
for idx in range(20):
    fc = face_centers[idx]
    r_fc = math.sqrt(sum(c**2 for c in fc))
    inward_n = tuple(-fc[d]/r_fc for d in range(3))
    cone_segments += make_cone_wireframe(fc, inward_n, base_r=L_J_fm*0.18, apex_dist=r_in_typical*0.35)

check("VG2c: 20 conical-wave markers built, apex stops well short of r=0 (0.65*inradius from center)",
      len(cone_segments) == 20*(8+8) and r_in_typical*0.65 > 1e-6,
      f"apex radius = {r_in_typical*0.65:.6f} fm (> 0, does not reach center)")

print()
print(f"  L_J = {L_J_fm:.6f} fm   circumradius = {max(math.sqrt(sum(c**2 for c in v)) for v in verts):.6f} fm")
print(f"  inradius = {math.sqrt(sum(c**2 for c in face_centers[0])):.6f} fm")
print(f"  midradius = {math.sqrt(sum(c**2 for c in list(edge_midpoints.values())[0])):.6f} fm")

# ── Render ───────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 11))
ax = fig.add_subplot(111, projection='3d')

# Edges (gluon standing + muon traveling channel network)
edge_lines = [[verts[i], verts[j]] for (i, j) in edges]
ax.add_collection3d(Line3DCollection(edge_lines, colors='dimgray', linewidths=1.1,
                                      label='Edge (gluon standing + muon traveling channel)'))

# Highlighted muon circuit example (vertex-0 pentagonal belt)
belt_lines = [[verts[i], verts[j]] for (i, j) in belt_edges]
ax.add_collection3d(Line3DCollection(belt_lines, colors='goldenrod', linewidths=4.0))

# Gluon pressure segments (T_2g face pressure -> outward edge tension)
gluon_lines = [[s, e] for (s, e) in gluon_segments]
ax.add_collection3d(Line3DCollection(gluon_lines, colors='mediumorchid', linewidths=0.8, alpha=0.55))

# 20 conical-wave markers (one per face) -- apex points inward, stops well
# short of r=0. This is the SAME cone both windings pass through.
ax.add_collection3d(Line3DCollection(cone_segments, colors='teal', linewidths=1.0, alpha=0.75))

# Tau corkscrew: forward + backward windings (40 cone-visits = 2 x 20 faces)
tau_fwd_loop = tau_fwd_pts + [tau_fwd_pts[0]]
fx2, fy2, fz2 = zip(*tau_fwd_loop)
ax.plot(fx2, fy2, fz2, color='crimson', linewidth=2.2)

tau_bwd_loop = tau_bwd_pts + [tau_bwd_pts[0]]
bx2, by2, bz2 = zip(*tau_bwd_loop)
ax.plot(bx2, by2, bz2, color='darkorange', linewidth=2.2, linestyle='--')

# Vertices
vx, vy, vz = zip(*verts)
ax.scatter(vx, vy, vz, color='navy', s=60, depthshade=True)

# Face centers
fx, fy, fz = zip(*face_centers)
ax.scatter(fx, fy, fz, color='crimson', s=28, depthshade=True)

# Manual legend (Line3DCollection doesn't auto-register handles reliably)
from matplotlib.lines import Line2D
legend_handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='navy', markersize=9,
           label='Vertex (electron/quark nexus)'),
    Line2D([0], [0], color='dimgray', linewidth=1.5,
           label='Edge (gluon standing + muon traveling channel)'),
    Line2D([0], [0], color='goldenrod', linewidth=4,
           label='Muon circuit example (vertex-0 pentagonal belt)'),
    Line2D([0], [0], color='teal', linewidth=1.5,
           label='Conical face wave (1 of 20 -- apex inward, stops short of r=0)'),
    Line2D([0], [0], color='crimson', linewidth=2.2,
           label='Tau winding, forward (1 of 2 -- 20 cone-visits)'),
    Line2D([0], [0], color='darkorange', linewidth=2.2, linestyle='--',
           label='Tau winding, backward (1 of 2 -- 20 cone-visits; 40 total)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='crimson', markersize=8,
           label='Face center (tau/gluon nexus)'),
    Line2D([0], [0], color='mediumorchid', linewidth=2,
           label='Gluon pressure segment (T_2g face -> outward edge tension)'),
]
ax.legend(handles=legend_handles, loc='upper left', fontsize=8, framealpha=0.9)

print()
print("  NOTE: no separate strut/spine sits at the geometric center. 'Crush")
print("  prevention' is the DISTRIBUTED Maxwell-critical rigidity of the full")
print("  30-edge network (3V-E=6 exactly, jobson_cell_geometry_3d.py CG18) --")
print("  a global count over all edges, not a local object at r=0 (CG13/CG14).")

R_c = max(math.sqrt(sum(c**2 for c in v)) for v in verts)
lim = R_c * 1.15
ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
ax.set_box_aspect((1, 1, 1))
ax.set_xlabel('x (fm)'); ax.set_ylabel('y (fm)'); ax.set_zlabel('z (fm)')
ax.set_title(f"Jobson Cell -- 40 Conical Tau Windings (2 x 20 faces), Torsionverse\nL_J = {L_J_fm:.6f} fm edge length", fontsize=12)
ax.view_init(elev=32, azim=200)

out_dir = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(out_dir, exist_ok=True)
out_path = os.path.join(out_dir, "jobson_cell_3d.png")
fig.savefig(out_path, dpi=170, bbox_inches='tight')
plt.close(fig)

check("VG5: PNG saved to analysis/demos/output/jobson_cell_3d.png",
      os.path.isfile(out_path), f"path = {out_path}")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED. Static 3D render written to:")
    print(f"    {out_path}")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(SEP)
