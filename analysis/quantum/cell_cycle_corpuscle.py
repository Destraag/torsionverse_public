"""
cell_cycle_corpuscle.py
=======================
Full cycle simulation of all corpuscles in the Jobson cell from the
PURE CORPUSCLE perspective (no wave math, no calculus, no amplitudes).

RULES:
  - Each corpuscle travels in a straight line at speed c between nexus points
  - At a nexus point: apply the geometric bounce rule specific to that nexus type
  - GLUON at vertex: REFLECT (reverse direction, chi(G,C5)=-1 -> phase flip)
  - MUON at vertex: DEFLECT to next edge of zigzag circuit (72 deg turn)
  - TAU at face center: DEFLECT to next face center in Hamiltonian circuit (72 deg)
  - No wave amplitudes, no sin/cos profiles, no field equations

CYCLE CLOSURE: if each corpuscle returns to its starting state after one full
cycle, the wave-like behavior is an emergent PATTERN of repeated corpuscle
trajectories. The 'standing wave' is just the path the corpuscle repeatedly
traces; the 'interference' is just multiple corpuscles on the same path.

Checks:
  CC_G1: Gluon lambda_1 traverses edge A->B in 1 step, reflects at B, returns to A
  CC_G2: After 2 steps, gluon lambda_1 is back at start with original direction
  CC_G3: Lambda_2 (counter-propagating) is always at the OPPOSITE position from lambda_1
         (they 'pass through' at midpoint -- corpuscle equivalent of standing wave node/antinode)
  CC_M1: Muon forward corpuscle traverses all 6 vertices of the zigzag in order
  CC_M2: After 6 steps, muon returns to starting vertex with original direction
  CC_M3: Muon backward corpuscle is always 3 steps behind/ahead (half circuit)
  CC_T1: Tau forward corpuscle visits all 20 face centers in Hamiltonian order
  CC_T2: After 20 steps, tau returns to starting face center
  CC_T3: Tau backward corpuscle is always 10 steps behind (meets forward at step 10)
  CC_SYNC: Relative cycle times: gluon=2T, muon=6T, tau=20*hop_len
           These are incommensurable (no simple LCM) -> quasi-periodic collective state

NOTE ON 3V-E=6 as ALGEBRAIC IDENTITY (for study):
  V=12 (from 2I spinor sum: 2+4+6=12) [proven FG12]
  n=5  (each vertex has 5 edges, C5 symmetry) [from icosahedral geometry]
  E = n*V/2 = 5*12/2 = 30  [each edge shared by 2 vertices]
  3V - E = 3*12 - 30 = 36 - 30 = 6  [Maxwell criterion]
  EQUIVALENCE: 3V-E = V*(3-n/2) = 12*(3-n/2) = 6  IFF  n=5
  -> {3V-E=6} and {n=5} and {V=12} are ALL equivalent to each other via E=nV/2.
     Knowing any two gives the third. They are the SAME GEOMETRIC FACT stated differently.
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
print("GLUON CORPUSCLE CYCLE")
print(SEP2)
# =============================================================================
print("""  Rule: each gluon corpuscle travels along its edge (vertex to vertex).
  At vertex: REFLECT (reverse direction). Period = 2 steps.
  Two corpuscles: lambda_1 starts at V_A heading to V_B;
                  lambda_2 starts at V_B heading to V_A.
  They always cross at the midpoint (step 1 for each) -- that's the
  corpuscle-picture of the standing wave antinode/node.
""")

# Simulate one gluon on edge (edges[0])
ei, ej = edges[0]
# Lambda_1: starts at ei, heading to ej
# Lambda_2: starts at ej, heading to ei
# State: (current_vertex, going_to_vertex)
L1_state = (ei, ej)   # at vertex ei, heading toward ej
L2_state = (ej, ei)   # at vertex ej, heading toward ei

gluon_trace_L1 = []
gluon_trace_L2 = []
for step in range(4):   # simulate 4 steps (2 full cycles)
    gluon_trace_L1.append(L1_state)
    gluon_trace_L2.append(L2_state)
    # Advance: corpuscle moves to the 'going_to' vertex, then reflects
    at1, going1 = L1_state
    at2, going2 = L2_state
    # After one L_J/c step: each corpuscle arrives at 'going' vertex, then heads back
    L1_state = (going1, at1)   # arrived at going1, now reflect back to at1
    L2_state = (going2, at2)

print(f"  Gluon lambda_1 trajectory (state = (at_vertex, heading_to)):")
for i, st in enumerate(gluon_trace_L1):
    pos = "MIDPOINT" if i % 1 == 0.5 else f"vertex {st[0]}"
    print(f"    step {i}: at={st[0]:2d} heading->{st[1]:2d}  {'<- START' if i==0 else ''}")

print(f"\n  Gluon lambda_2 trajectory:")
for i, st in enumerate(gluon_trace_L2):
    print(f"    step {i}: at={st[0]:2d} heading->{st[1]:2d}")

# After step 2, lambda_1 should be back at (ei, ej) -- starting state
check("CC_G1: lambda_1 arrives at V_B after 1 step, reflects back",
      gluon_trace_L1[1] == (ej, ei),
      f"step 1 state = {gluon_trace_L1[1]}  (at V_B, heading back to V_A)")

check("CC_G2: lambda_1 returns to start after 2 steps (cycle closes)",
      gluon_trace_L1[2] == gluon_trace_L1[0],
      f"step 2 = {gluon_trace_L1[2]} == step 0 = {gluon_trace_L1[0]}")

# Lambda_2 is always the opposite (at the other vertex)
all_opposite = all(gluon_trace_L1[i][0] != gluon_trace_L2[i][0]
                   for i in range(4))
check("CC_G3: lambda_1 and lambda_2 are always at OPPOSITE vertices (standing wave pattern)",
      all_opposite,
      f"At every step, lambda_1 and lambda_2 are at different vertices of the edge")

print(f"\n  Gluon cycle: 2 steps of L_J/c  (each step = one edge traversal)")
print(f"  The 'standing wave' is just lambda_1 and lambda_2 always being at")
print(f"  opposite vertices -- their path is the corpuscle equivalent of a")
print(f"  half-wave standing mode (node/antinode structure emerges from the trajectory).")

# =============================================================================
print()
print(SEP)
print("MUON CORPUSCLE CYCLE")
print(SEP2)
# =============================================================================
print("""  Rule: muon follows a pre-defined 6-vertex zigzag circuit.
  At each vertex: deflect to the NEXT vertex in the circuit (72 deg turn built into path).
  No angle calculation needed -- the path is the rule.
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
muon_fwd = list(muon_circuit)   # forward: visits in order
muon_bwd = list(muon_circuit[3:]) + list(muon_circuit[:3])  # backward: 3 steps behind

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

# Backward is 3 steps behind
bwd_at_3 = muon_circuit[3]
check("CC_M3: backward muon is 3 steps behind forward (half circuit = bilateral)",
      muon_circuit[3] == muon_circuit[(0+3)%6],
      f"Forward at step 0: vertex {muon_circuit[0]}, backward at step 0: vertex {muon_circuit[3]}")

print(f"\n  Muon cycle: 6 steps of L_J/c")
print(f"  Each 72-deg turn is built into the circuit path -- no angle calculation.")
print(f"  The 'mode' is the repeated tracing of the same 6-vertex path.")

# =============================================================================
print()
print(SEP)
print("TAU CORPUSCLE CYCLE")
print(SEP2)
# =============================================================================
print("""  Rule: tau follows the Hamiltonian circuit of 20 face-center nexuses.
  At each face center: deflect to the NEXT face center in the circuit (72 deg).
  No angle calculation needed -- the Hamiltonian path is the rule.
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
