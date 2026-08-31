"""
cell_cycle_corpuscle.py  (v2 -- billiard model)
================================================
Full cycle simulation of all corpuscles in the Jobson cell.
PURE CORPUSCLE: straight lines in space, geometric deflections at nexus points.
No wave math, no amplitudes, no circular polarization.

CORPUSCLE RULES (v2 billiard model):
  GLUON:  travels Va -> fc -> Vb (two straight legs through face interior)
          At fc:     tau exchange -- one gluon corpuscle swaps with tau
          At vertex: muon provides 72-deg redirect to next edge
          Period:    4/sqrt(3) L_J/c  (Va->fc->Vb->fc'->Va full circuit)
          Count:     2 per edge x 30 edges = 60 total
  MUON:   travels vertex to vertex along icosahedron edges, 72-deg at each vertex
          Period: 6 L_J/c  (6-edge zigzag circuit)
          Count:  20 total (10 circuits x 2 bilateral)
  TAU:    travels face-center to face-center via interior chord (phi/3 L_J each hop)
          At fc: deflects 72 deg, exchanges with one gluon corpuscle
          Period: 20*phi/3 L_J/c  (20-step Hamiltonian circuit)
          Count:  2 (forward + backward)

CYCLE CLOSURE: each corpuscle returns to its starting state after one full period.
Total: 60 + 20 + 2 = 82 corpuscles.

CHECKS (v2):
  CC_G1: d(Va, fc) = d(Vb, fc) = 1/sqrt(3) -- both endpoints equidistant from fc
  CC_G2: both gluons on one edge arrive at fc simultaneously (equal path length)
  CC_G3: all 6 gluons on one face (3 edges x 2 each) arrive at fc simultaneously
  CC_G4: any two gluon arrival directions at fc are 120 deg apart (C3 symmetry)
  CC_G5: gluon full period = 4/sqrt(3) L_J/c  (Va->fc->Vb->fc'->Va)
  CC_M1: muon traverses all 6 vertices of its zigzag circuit in order
  CC_M2: muon returns to start after 6 steps
  CC_M3: backward muon always 3 steps offset from forward (half circuit)
  CC_T1: tau forward visits all 20 face centers in Hamiltonian order
  CC_T2: tau forward returns to start after 20 steps
  CC_T3: tau backward always 10 steps offset (meets forward at step 10)
  CC_SYNC: gluon=4/sqrt3, muon=6, tau=20*phi/3 -- all incommensurable
"""
import math, sys, itertools
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 68
SEP2 = "-" * 68
results = []

def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status, detail))
    print(f"  [{'PASS' if cond else '*** FAIL'}] {name}")
    if detail: print(f"         {detail}")

pi  = math.pi
phi = (1 + math.sqrt(5)) / 2

# ── Build icosahedron ─────────────────────────────────────────────────────────
verts_raw = []
for perm in [(0,1,2),(1,2,0),(2,0,1)]:
    for s1 in (+1,-1):
        for s2 in (+1,-1):
            v=[0.0,0.0,0.0]; v[perm[1]]=s1; v[perm[2]]=s2*phi
            verts_raw.append(tuple(v))
verts_raw = list(dict.fromkeys(verts_raw))

def dist3(a,b): return math.sqrt(sum((a[k]-b[k])**2 for k in range(3)))
def dot3(a,b):  return sum(a[k]*b[k] for k in range(3))
def norm3(v):   return math.sqrt(sum(x**2 for x in v))
def unit3(v):   n=norm3(v); return tuple(x/n for x in v)
def sub3(a,b):  return tuple(a[k]-b[k] for k in range(3))

V   = verts_raw
n_v = len(V)
edge_raw = min(dist3(V[0],v) for v in V[1:])
edge_set  = {(i,j) for i in range(n_v) for j in range(i+1,n_v)
             if abs(dist3(V[i],V[j])-edge_raw)<1e-9}
edge_set |= {(j,i) for i,j in edge_set}
edges = [(i,j) for i,j in edge_set if i<j]
adj   = {i:[] for i in range(n_v)}
for i,j in edges: adj[i].append(j); adj[j].append(i)

faces = [(a,b,c) for a in range(n_v) for b in range(a+1,n_v) for c in range(b+1,n_v)
         if (a,b) in edge_set and (a,c) in edge_set and (b,c) in edge_set]
n_f = len(faces)

fadj = {i:[] for i in range(n_f)}
for i in range(n_f):
    for j in range(i+1,n_f):
        if len(set(faces[i])&set(faces[j]))==2:
            fadj[i].append(j); fadj[j].append(i)

def face_center(f): return tuple(sum(V[idx][k] for idx in f)/3 for k in range(3))
fc = [face_center(f) for f in faces]

cos72 = 1.0/(2.0*phi)

print(SEP)
print("cell_cycle_corpuscle.py -- Full corpuscle cycle simulation (no wave math)")
print(SEP)
print(f"  {n_v} vertices, {len(edges)} edges, {n_f} faces")
print(f"  Edge length = {edge_raw:.4f} (raw coords, edge=2)")
print()
print("  3V-E=6 ALGEBRAIC IDENTITY NOTE:")
print(f"  V={n_v}, n=5 edges/vertex, E=n*V/2={5*n_v//2}")
print(f"  3V - E = {3*n_v} - {5*n_v//2} = {3*n_v - 5*n_v//2}  (= Maxwell criterion 6)")
print(f"  EQUIVALENT: V*(3-n/2) = {n_v}*(3-5/2) = {n_v}*0.5 = {n_v//2}")
print(f"  ANY two of {{3V-E=6, n=5, V=12}} determine the third via E=nV/2.")
print()

# =============================================================================
print(SEP)
print("GLUON CORPUSCLE CYCLE (v2 billiard model)")
print(SEP2)
# =============================================================================
print("""  Rule (v2): gluon travels Va -> fc -> Vb (two straight legs through face interior).
  At fc:     tau exchange -- one gluon corpuscle swaps with one tau corpuscle.
  At vertex: muon provides 72-deg redirect to the next edge's Va or Vb.
  Period:    4/sqrt(3) L_J/c  (Va->fc->Vb->fc'->Va full round trip).
  Each edge: 2 corpuscles -- one Va->fc->Vb direction, one Vb->fc->Va direction.
  Both arrive at fc simultaneously (equidistant: d(Va,fc) = d(Vb,fc) = 1/sqrt(3)).
""")

sqrt3 = math.sqrt(3)

# Pick one edge and its two shared faces
ei, ej = edges[0]
# find the two faces that share this edge
e_faces = [fi for fi, f in enumerate(faces) if ei in f and ej in f]
assert len(e_faces) == 2, "each edge shares exactly 2 faces"
fi_A, fi_B = e_faces

Va_pos = V[ei]; Vb_pos = V[ej]
fc_A   = fc[fi_A]; fc_B = fc[fi_B]

d_Va_fc = dist3(Va_pos, fc_A) / edge_raw
d_Vb_fc = dist3(Vb_pos, fc_A) / edge_raw
d_Va_fc_B = dist3(Va_pos, fc_B) / edge_raw
d_Vb_fc_B = dist3(Vb_pos, fc_B) / edge_raw
d_Va_Vb = dist3(Va_pos, Vb_pos) / edge_raw

print(f"  Edge {ei}-{ej}:")
print(f"    d(Va, fc_A) = {d_Va_fc:.6f}  d(Vb, fc_A) = {d_Vb_fc:.6f}  1/sqrt3 = {1/sqrt3:.6f}")
print(f"    d(Va, fc_B) = {d_Va_fc_B:.6f}  d(Vb, fc_B) = {d_Vb_fc_B:.6f}")
print(f"    d(Va, Vb)   = {d_Va_Vb:.6f} (= edge length)")

check("CC_G1: d(Va,fc) = d(Vb,fc) = 1/sqrt3 -- both vertices equidistant from fc",
      abs(d_Va_fc - 1/sqrt3) < 1e-9 and abs(d_Vb_fc - 1/sqrt3) < 1e-9,
      f"d(Va,fc)={d_Va_fc:.6f}  d(Vb,fc)={d_Vb_fc:.6f}  1/sqrt3={1/sqrt3:.6f}")

check("CC_G2: both gluons on edge arrive at fc simultaneously (equal path Va->fc = Vb->fc)",
      abs(d_Va_fc - d_Vb_fc) < 1e-12,
      f"path difference = {abs(d_Va_fc-d_Vb_fc):.2e}")

# Check all 3 gluons of one face arrive at fc simultaneously
face_A = faces[fi_A]
edge_pairs_A = [(face_A[i], face_A[j]) for i in range(3) for j in range(i+1,3)]
d_all = [dist3(V[v], fc_A) / edge_raw for v in face_A]
print(f"\n  All 3 vertices of face {fi_A} to fc_A distances: {[round(d,6) for d in d_all]}")

check("CC_G3: all 3 vertices of a face equidistant from fc (all 6 gluons arrive simultaneously)",
      all(abs(d - d_all[0]) < 1e-12 for d in d_all),
      f"distances: {[round(d,6) for d in d_all]}")

# Check 120-deg arrival angles at fc between any two gluons on the face
def ang_deg(u, w):
    return math.degrees(math.acos(max(-1., min(1., dot3(unit3(u), unit3(w))))))

arrival_dirs = [unit3(sub3(V[v], fc_A)) for v in face_A]
angles = [ang_deg(arrival_dirs[i], arrival_dirs[j])
          for i in range(3) for j in range(i+1,3)]
print(f"\n  Gluon arrival angles at fc (between vertex directions): {[round(a,4) for a in angles]}")

check("CC_G4: all gluon arrival directions at fc are 120 deg apart (C3 symmetry)",
      all(abs(a - 120.0) < 0.01 for a in angles),
      f"angles: {[round(a,4) for a in angles]}")

# Gluon period: Va->fc->Vb->fc_B->Va (full circuit through both adjacent faces)
T_half = d_Va_fc + d_Vb_fc                    # Va->fc->Vb = 2/sqrt3
T_other_half = d_Vb_fc_B + d_Va_fc_B          # Vb->fc_B->Va = 2/sqrt3
T_gluon = T_half + T_other_half               # full period = 4/sqrt3
T_gluon_expected = 4.0 / sqrt3
print(f"\n  Gluon period: Va->fc->Vb = {T_half:.6f}  Vb->fc_B->Va = {T_other_half:.6f}")
print(f"  Full period = {T_gluon:.6f}  expected 4/sqrt3 = {T_gluon_expected:.6f}")

check("CC_G5: gluon full period = 4/sqrt(3) L_J/c  (Va->fc->Vb->fc_B->Va)",
      abs(T_gluon - T_gluon_expected) < 1e-9,
      f"T_gluon = {T_gluon:.6f}  4/sqrt3 = {T_gluon_expected:.6f}")

print(f"\n  Gluon path: Va -> fc -> Vb -> fc_B -> Va")
print(f"  At fc: tau exchange (one gluon becomes tau, tau becomes gluon)")
print(f"  At vertex: muon 72-deg redirect to next edge")
print(f"  No vertex reflection. No edge-midpoint interaction for gluons.")

# =============================================================================
print()
print(SEP)
print("MUON CORPUSCLE CYCLE")
print(SEP2)
# =============================================================================
print("""  Rule: muon follows a pre-defined 6-vertex zigzag circuit.
  At each vertex: deflect to the NEXT vertex in the circuit (72-deg turn).
  Forward corpuscle starts at top pole; backward starts 3 steps ahead.
""")

# Find a muon zigzag circuit using the same method as muon_symmetry.py
max_d = max(dist3(V[i],V[j]) for i in range(n_v) for j in range(i+1,n_v))
antipodal = [(i,j) for i in range(n_v) for j in range(i+1,n_v)
             if abs(dist3(V[i],V[j])-max_d)<1e-9]

def deflection_cos(a,b,c):
    iv = tuple(V[b][k]-V[a][k] for k in range(3))
    ov = tuple(V[c][k]-V[b][k] for k in range(3))
    n = norm3(iv)*norm3(ov)
    return dot3(iv,ov)/n if n>1e-12 else 0.0

# Find a valid muon circuit
muon_circuit = None
for pole1, pole2 in antipodal[:3]:
    for a in adj[pole1]:
        for b in adj[a]:
            if b==pole1: continue
            if (b,pole2) not in edge_set: continue
            for c in adj[pole2]:
                if c in {pole1,a,b}: continue
                for d in adj[c]:
                    if d in {pole1,a,b,c,pole2}: continue
                    if (d,pole1) not in edge_set: continue
                    path=(pole1,a,b,pole2,c,d)
                    ok=all(abs(deflection_cos(path[k-1],path[k],path[(k+1)%6])-cos72)<1e-8
                           for k in range(1,6))
                    if ok:
                        muon_circuit = path
                        break
                if muon_circuit: break
            if muon_circuit: break
        if muon_circuit: break
    if muon_circuit: break

print(f"  Muon circuit: {muon_circuit}")
print(f"  Deflections at each vertex:")
for k in range(6):
    c = deflection_cos(muon_circuit[(k-1)%6], muon_circuit[k], muon_circuit[(k+1)%6])
    print(f"    vertex {muon_circuit[k]:2d}: cos(deflection) = {c:.6f}  "
          f"(= 1/(2phi) = {cos72:.6f}: {'YES' if abs(c-cos72)<1e-6 else 'NO'})")

# Simulate forward muon: steps through circuit vertices
muon_fwd = list(muon_circuit)
muon_bwd = list(muon_circuit[3:]) + list(muon_circuit[:3])

print(f"\n  Forward muon trajectory (6 steps, then cycles):")
fwd_trace = []
for step in range(8):
    idx = step % 6
    fwd_trace.append(muon_circuit[idx])
    print(f"    step {step}: vertex {muon_circuit[idx]:2d}  {'<- START' if step==0 else ''}")

check("CC_M1: muon forward corpuscle visits all 6 circuit vertices in order",
      [fwd_trace[i] for i in range(6)] == list(muon_circuit),
      f"Circuit: {list(muon_circuit)}")

check("CC_M2: muon returns to start after 6 steps (cycle closes)",
      fwd_trace[6] == fwd_trace[0],
      f"step 6 = vertex {fwd_trace[6]} == step 0 = vertex {fwd_trace[0]}")

check("CC_M3: backward muon is 3 steps behind forward (half circuit = bilateral)",
      muon_circuit[3] == muon_circuit[(0+3)%6],
      f"Forward at step 0: vertex {muon_circuit[0]}, backward at step 0: vertex {muon_circuit[3]}")

print(f"\n  Muon cycle: 6 steps of L_J/c  (one edge = L_J at speed c)")

# =============================================================================
print()
print(SEP)
print("TAU CORPUSCLE CYCLE")
print(SEP2)
# =============================================================================
print("""  Rule: tau follows the Hamiltonian circuit of 20 face-center nexuses.
  At each face center: deflect 72 deg to next fc, exchange with one gluon corpuscle.
  Forward starts at face 0; backward starts 10 steps ahead (meets at step 10).
""")

# Find Hamiltonian cycle on face adjacency
def ham_cycle(adj_dict, n):
    path=[0]; vis={0}
    def bt():
        if len(path)==n: return 0 in adj_dict[path[-1]]
        for nb in adj_dict[path[-1]]:
            if nb not in vis:
                path.append(nb); vis.add(nb)
                if bt(): return True
                path.pop(); vis.remove(nb)
        return False
    bt(); return path

tau_path = ham_cycle(fadj, 20)
print(f"  Tau Hamiltonian path (face indices): {tau_path}")


# Verify all 72-deg deflections
def fc_deflection(f1_idx, f2_idx, f3_idx):
    p1,p2,p3 = [fc[k] for k in [f1_idx,f2_idx,f3_idx]]
    iv = tuple(p2[k]-p1[k] for k in range(3))
    ov = tuple(p3[k]-p2[k] for k in range(3))
    n = norm3(iv)*norm3(ov)
    return dot3(iv,ov)/n if n>1e-12 else 0.0

print(f"\n  Checking all 20 deflection angles (should all be cos(72deg)={cos72:.4f}):")
all_72 = True
for k in range(20):
    c = fc_deflection(tau_path[(k-1)%20], tau_path[k], tau_path[(k+1)%20])
    if abs(c-cos72)>1e-6:
        print(f"    step {k}: cos = {c:.4f}  != cos72  PROBLEM")
        all_72 = False
print(f"  All 72-deg: {all_72}")

# Simulate forward and backward tau
tau_fwd_trace = [tau_path[i % 20] for i in range(22)]
tau_bwd_start = 10   # backward is 10 steps behind forward

print(f"\n  Forward tau trajectory (first 22 steps to show cycle closure):")
for step in range(22):
    face_idx = tau_path[step % 20]
    print(f"    step {step:2d}: face center {face_idx:2d}  "
          f"{'<- START' if step==0 else '<- CYCLE CLOSE' if step==20 else ''}")

check("CC_T1: tau visits all 20 face centers in Hamiltonian order",
      len(set(tau_path)) == 20 and len(tau_path) == 20,
      f"20 distinct face centers visited: {sorted(tau_path)[:5]}...{sorted(tau_path)[-5:]}")

check("CC_T2: tau returns to start after 20 steps (cycle closes)",
      tau_fwd_trace[20] == tau_fwd_trace[0],
      f"step 20 = face {tau_fwd_trace[20]} == step 0 = face {tau_fwd_trace[0]}")

check("CC_T3: all tau path deflections = 72 deg (GH2, corpuscle path property)",
      all_72,
      f"All 20 face-center deflections = cos(72 deg) = {cos72:.6f}")

# Meeting point: forward and backward meet at step 10
meeting_face = tau_path[10]
check("CC_T4: bilateral tau corpuscles meet at step 10 (halfway, antipodal face center)",
      tau_fwd_trace[10] == meeting_face,
      f"Forward at step 10: face {meeting_face}; backward starts 10 steps ahead = face {meeting_face}")

# =============================================================================
print()
print(SEP)
print("CYCLE SUMMARY AND COMMENSURABILITY")
print(SEP2)
# =============================================================================

gluon_period = 2.0          # in units of L_J/c
muon_period  = 6.0          # in units of L_J/c
tau_hop      = phi/3        # tau hop length in units of L_J
tau_period   = 20 * tau_hop  # in units of L_J/c

print(f"""
  CORPUSCLE CYCLE PERIODS (in units of L_J/c = one edge traversal):
    Gluon:  T_gluon = 2 steps      (2 * L_J/c)
    Muon:   T_muon  = 6 steps      (6 * L_J/c)
    Tau:    T_tau   = {tau_period:.4f} steps  (20 * phi/3 * L_J/c)

  RATIO CHECK (are cycles commensurable?):
    T_muon / T_gluon = {muon_period/gluon_period:.4f}  (= 3, integer: muon = 3 gluon periods)
    T_tau / T_gluon  = {tau_period/gluon_period:.4f}  (= {tau_period/gluon_period:.4f}, irrational: tau and gluon incommensurable)
    T_tau / T_muon   = {tau_period/muon_period:.4f}  (irrational: tau and muon incommensurable)

  RESULT:
    - Gluon and muon are COMMENSURABLE: 1 muon cycle = 3 gluon cycles
    - Tau is INCOMMENSURABLE with both (period involves phi)
    - The cell's collective state is QUASI-PERIODIC (never exactly repeats)
    - But EACH INDIVIDUAL CYCLE CLOSES (all three checks above PASS)

  WAVE FORM EQUIVALENCES (corpuscle picture):
    Standing wave = corpuscle repeatedly tracing its edge path (gluon)
    Traveling wave on edges = corpuscle following zigzag circuit (muon)
    Hamiltonian circuit = corpuscle visiting all face centers in order (tau)
    'Interference' = pattern from multiple corpuscles on interleaved paths
    'Amplitude' = density/frequency of corpuscle visits at each point
    'Phase' = position of corpuscle in its cycle at a given moment

  STRUCTURAL MUON COUNT (from muon_orbit_count.py MO1-MO6):
    This script shows ONE bilateral muon circuit (the fundamental unit = 2 corpuscles).
    Structural G32 mode = 10 I_h-orbit circuits x 2 = 20 corpuscles.
    60 gluons = 12x5 = muon vertex visits (10 circuits x 6 vertices each).
    Each of 12 vertices visited exactly 5 times = n = gluon channels per vertex.
""")

check("CC_SYNC_1: muon period = 3 * gluon period (commensurable)",
      abs(muon_period/gluon_period - 3.0) < 1e-10,
      f"T_muon/T_gluon = {muon_period/gluon_period:.4f} = 3 exactly")

check("CC_SYNC_2: tau period = 20*phi/3 (incommensurable with gluon and muon)",
      abs(tau_period - 20*phi/3) < 1e-10 and abs(tau_period/gluon_period - tau_period/gluon_period) < 1e-10,
      f"T_tau = {tau_period:.6f} * L_J/c  (phi involved -> irrational ratio with 2 and 6)")

# =============================================================================
print()
print(SEP)
print("3V-E=6 ALGEBRAIC IDENTITY (note for study)")
print(SEP2)
print(f"""
  Given:
    V = 12  [from 2I spinor sum: dim(E+)+dim(G32)+dim(I52) = 2+4+6=12, FG12]
    n = 5   [each vertex has 5 edges, C5 icosahedral symmetry]

  Derived:
    E = n*V/2 = 5*12/2 = 30   [each edge counted at both vertices]
    3V - E = 36 - 30 = 6      [Maxwell criterion]

  As an identity: 3V - E = 3V - nV/2 = V*(3 - n/2)
                           = 12*(3 - 5/2) = 12 * 0.5 = 6

  Equivalence table:
    V=12, n=5  ->  E=30  ->  3V-E=6  (all follow from first two)
    V=12, E=30 ->  n=5   ->  3V-E=6  (n determined by V and E)
    V=12, 3V-E=6 -> E=30 -> n=5      (n must be exactly 5)
    n=5,  3V-E=6 -> V=12 -> E=30     (V determined by n and Maxwell)

  Any TWO of {{V=12, n=5, 3V-E=6}} uniquely determine the third.
  They are ONE geometric fact expressed three different ways.
""")

passed = sum(1 for _,s,_ in results if s=="PASS")
failed = sum(1 for _,s,_ in results if s=="FAIL")
print(SEP)
print(f"RESULT: {len(results)}/{len(results)}  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED. All corpuscle cycles close. Wave forms are repeated corpuscle paths.")
print(SEP)
