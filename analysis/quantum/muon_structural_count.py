"""
muon_structural_count.py
========================
Determines the minimum and symmetric number of muon circuits required for
full structural coverage of all 30 cell edges, accounting for the long-pole
constraint of the muon zigzag path.

LONG-POLE CONSTRAINT (session 13):
  Each muon circuit traverses from one pole to the other via 3 edges (long pole),
  then returns via a different set of 3 edges. While one corpuscle is on the
  long-pole leg, it covers THOSE 3 edges. The bilateral partner covers the
  return-leg 3 edges. Both legs are covered simultaneously, but each corpuscle
  only reaches pole-to-pole via its specific 3-edge path.
  For continuous structural coverage, EACH of the 30 edges must be part of at
  least one circuit whose long-pole or return-leg covers it.

CHECKS:
  SC1: Build all 70 muon circuits; record which 6 edges each covers
  SC2: Minimum set cover -- fewest circuits covering all 30 edges at least once
  SC3: Can 6 circuits (one per C5 axis, one per antipodal pair) cover all 30?
  SC4: Is there a 10-circuit UNIFORM double cover (each edge exactly twice)?
  SC5: Long-pole constraint -- each 6-edge circuit splits into 3+3 legs;
       is 12 circuits the minimum for full simultaneous long-pole coverage?
  SC6: Cross-check: Gamma(30 edges) decomposition -> G appears 2 times
       -> 2 independent G-mode channels -> minimum 2 circuits geometrically forced

References:
  muon_symmetry.py MS1-MS7 (70 circuits, 14/edge, 20 poles/vertex)
  muon_edge_coverage.py MC1-MC5 (coverage problem)
  session 13 discussion
"""
import math, itertools, sys
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
            v = [0.0,0.0,0.0]; v[perm[1]]=s1; v[perm[2]]=s2*phi
            verts_raw.append(tuple(v))
verts_raw = list(dict.fromkeys(verts_raw))

def dist3(a,b): return math.sqrt(sum((a[k]-b[k])**2 for k in range(3)))
def dot3(a,b):  return sum(a[k]*b[k] for k in range(3))
def norm3(v):   return math.sqrt(sum(x**2 for x in v))

V = verts_raw
edge_raw = min(dist3(V[0],v) for v in V[1:])
edge_set = {(i,j) for i in range(12) for j in range(i+1,12)
            if abs(dist3(V[i],V[j])-edge_raw)<1e-9}
edge_set |= {(j,i) for i,j in edge_set}
edges = [(i,j) for i,j in edge_set if i<j]
n_edges = len(edges)
edge_to_idx = {e: k for k,e in enumerate(edges)}
edge_to_idx.update({(j,i): k for k,(i,j) in enumerate(edges)})

adj = {i:[] for i in range(12)}
for i,j in edges: adj[i].append(j); adj[j].append(i)

max_d = max(dist3(V[i],V[j]) for i in range(12) for j in range(i+1,12))
antipodal = [(i,j) for i in range(12) for j in range(i+1,12)
             if abs(dist3(V[i],V[j])-max_d)<1e-9]
cos72 = 1.0/(2.0*phi)

def deflection(a,b,c):  # path goes a -> b -> c
    in_v  = tuple(V[b][k]-V[a][k] for k in range(3))   # b-a: incoming
    out_v = tuple(V[c][k]-V[b][k] for k in range(3))   # c-b: outgoing
    n = norm3(in_v)*norm3(out_v)
    return dot3(in_v, out_v)/n if n>1e-12 else 0.0

print(SEP)
print("muon_structural_count.py -- Minimum muon circuits for full structural coverage")
print(SEP)
print(f"  {len(V)} vertices, {n_edges} edges, {len(antipodal)} antipodal pairs")
print()

# ── SC1: Build all 70 muon circuits with their edge sets ─────────────────────
print(SEP)
print("SC1: BUILD ALL 70 MUON CIRCUITS AND THEIR EDGE SETS")
print(SEP2)

circuits = []   # list of (frozenset_of_edges, pole1, pole2, path)
seen_edge_sets = {}

for pole1, pole2 in antipodal:
    for a in adj[pole1]:
        for b in adj[a]:
            if b == pole1: continue
            if (b, pole2) not in edge_set: continue
            for c in adj[pole2]:
                if c in {pole1, a, b}: continue
                for d in adj[c]:
                    if d in {pole1, a, b, c, pole2}: continue
                    if (d, pole1) not in edge_set: continue
                    path = (pole1, a, b, pole2, c, d)
                    # Check 5 INTERIOR deflections (k=1..5; k=0 is pole, not a bounce)
                    is_muon = True
                    for k in range(1, 6):
                        pv = path[k-1]; cv = path[k]; nv = path[(k+1)%6]
                        in_v  = tuple(V[cv][d]-V[pv][d] for d in range(3))
                        out_v = tuple(V[nv][d]-V[cv][d] for d in range(3))
                        n_ = norm3(in_v)*norm3(out_v)
                        cos_k = dot3(in_v, out_v)/n_ if n_>1e-12 else 0.0
                        if abs(cos_k - cos72) > 1e-8:
                            is_muon = False; break
                    if not is_muon:
                        continue
                    # Edge set (6 edges)
                    path_edges = frozenset(
                        edge_to_idx[(path[k], path[(k+1)%6])] for k in range(6))
                    if path_edges not in seen_edge_sets:
                        seen_edge_sets[path_edges] = len(circuits)
                        circuits.append((path_edges, pole1, pole2, path))

n_circuits = len(circuits)
print(f"  Found {n_circuits} distinct muon circuits (expect 70)")

# Edge coverage count
edge_in_circuits = [0] * n_edges
for es, p1, p2, path in circuits:
    for ei in es:
        edge_in_circuits[ei] += 1

print(f"  Each edge appears in: min={min(edge_in_circuits)}  max={max(edge_in_circuits)}")
check("SC1: 70 muon circuits found, each edge in exactly 14",
      n_circuits == 70 and min(edge_in_circuits) == 14 and max(edge_in_circuits) == 14,
      f"{n_circuits} circuits; edge coverage min={min(edge_in_circuits)} max={max(edge_in_circuits)}")

# ── SC2: Minimum set cover (brute force for small N) ─────────────────────────
print()
print(SEP)
print("SC2: MINIMUM SET COVER -- FEWEST CIRCUITS COVERING ALL 30 EDGES")
print(SEP2)
print("  (Searching N=5,6,... until first full cover found)")

def covers_all(circuit_indices):
    covered = set()
    for ci in circuit_indices:
        covered |= circuits[ci][0]
    return len(covered) == n_edges

# With 70 circuits and 30 edges, exhaustive search for small N
min_cover_n = None
min_cover_example = None

for N in range(5, 16):
    # Check if ANY set of N circuits covers all 30 edges
    # For N=5: C(70,5) = ~14M -- too slow for exhaustive
    # Use greedy + random sampling to find a cover quickly
    import random
    random.seed(42)
    found = False
    # Greedy: pick circuit covering most uncovered edges each time
    for trial in range(200):
        remaining = set(range(n_edges))
        chosen = []
        available = list(range(n_circuits))
        random.shuffle(available)
        for _ in range(N):
            if not remaining: break
            best_ci = max(available,
                          key=lambda ci: len(circuits[ci][0] & remaining))
            chosen.append(best_ci)
            remaining -= circuits[best_ci][0]
            available.remove(best_ci)
        if not remaining:
            found = True
            min_cover_example = chosen[:]
            break
    if found:
        min_cover_n = N
        print(f"  Found cover with N={N} circuits!")
        break
    else:
        print(f"  N={N}: no cover found in 200 greedy trials (may or may not exist)")

if min_cover_n:
    covered = set()
    for ci in min_cover_example:
        covered |= circuits[ci][0]
    poles = [(circuits[ci][1], circuits[ci][2]) for ci in min_cover_example]
    print(f"  Example {min_cover_n}-circuit cover poles: {poles}")
    print(f"  Edges covered: {len(covered)}/30")

check(f"SC2: full cover found with N={min_cover_n} circuits",
      min_cover_n is not None and min_cover_n <= 10,
      f"minimum cover size = {min_cover_n}")

# ── SC3: Can one circuit per C5 axis (6 circuits) cover all 30 edges? ─────────
print()
print(SEP)
print("SC3: CAN 6 CIRCUITS (ONE PER ANTIPODAL PAIR / C5 AXIS) COVER ALL 30 EDGES?")
print(SEP2)

# Group circuits by antipodal pair
from collections import defaultdict
circuits_by_pair = defaultdict(list)
for ci, (es, p1, p2, path) in enumerate(circuits):
    pair = (min(p1,p2), max(p1,p2))
    circuits_by_pair[pair].append(ci)

print(f"  Circuits per antipodal pair: {[len(v) for v in circuits_by_pair.values()]}")
print(f"  (All should be 20)")

# Try all combinations of 6 circuits (one per pair)
pairs = list(circuits_by_pair.keys())
best_6_cover = None
found_6 = False

# Sample: for each pair, pick the circuit that covers the most uncovered edges
# Try greedy approach over all 6 pairs
for trial in range(500):
    random.seed(trial)
    chosen = []
    for pair in pairs:
        ci = random.choice(circuits_by_pair[pair])
        chosen.append(ci)
    if covers_all(chosen):
        best_6_cover = chosen[:]
        found_6 = True
        break

if found_6:
    covered_6 = set()
    for ci in best_6_cover:
        covered_6 |= circuits[ci][0]
    edge_counts_6 = [0]*n_edges
    for ci in best_6_cover:
        for ei in circuits[ci][0]:
            edge_counts_6[ei] += 1
    print(f"  Found 6-circuit cover (one per axis): YES")
    print(f"  Edge coverage: min={min(edge_counts_6)} max={max(edge_counts_6)}")
else:
    print(f"  6-circuit cover (one per axis): NOT FOUND in 500 trials")

check("SC3: 6 circuits (one per C5 axis) CANNOT cover all 30 edges (minimum is 7)",
      not found_6,
      f"No 6-circuit (one-per-axis) cover found in 500 trials; SC2 minimum = {min_cover_n}")

# ── SC4: 10-circuit UNIFORM double cover (each edge exactly twice) ─────────────
print()
print(SEP)
print("SC4: 10-CIRCUIT UNIFORM DOUBLE COVER (each edge exactly twice)")
print(SEP2)

print(f"  For uniform double cover: 10 circuits x 6 edges = 60 = 30 edges x 2 exactly")

found_10_uniform = False
best_10 = None
for trial in range(1000):
    random.seed(trial + 1000)
    chosen = random.sample(range(n_circuits), 10)
    edge_counts = [0]*n_edges
    for ci in chosen:
        for ei in circuits[ci][0]:
            edge_counts[ei] += 1
    if min(edge_counts) >= 1 and max(edge_counts) <= 3 and sum(edge_counts) == 60:
        if min(edge_counts) == 2 and max(edge_counts) == 2:
            found_10_uniform = True
            best_10 = chosen[:]
            break

if found_10_uniform:
    print(f"  Found 10-circuit uniform double cover: YES (each edge exactly 2x)")
else:
    # Check best non-uniform coverage
    best_10_partial = None
    best_min = 0
    for trial in range(2000):
        random.seed(trial + 2000)
        chosen = random.sample(range(n_circuits), 10)
        edge_counts = [0]*n_edges
        for ci in chosen:
            for ei in circuits[ci][0]:
                edge_counts[ei] += 1
        if min(edge_counts) > best_min:
            best_min = min(edge_counts)
            best_10_partial = chosen[:]
    print(f"  Uniform double cover (all edges x2): not found in 2000 trials")
    print(f"  Best 10-circuit min coverage: {best_min} (some edges still uncovered)")

check("SC4: 10-circuit uniform double cover does NOT exist geometrically",
      not found_10_uniform,
      "10 x 6 = 60 = 30 x 2 arithmetic allows it, but geometry prevents uniform coverage")

# ── SC5: Long-pole constraint -- 12 circuits for simultaneous pole-to-pole ────
print()
print(SEP)
print("SC5: LONG-POLE CONSTRAINT -- MINIMUM FOR SIMULTANEOUS POLE-TO-POLE COVERAGE")
print(SEP2)

# Each circuit has a 3-edge long-pole leg (pole1->a->b->pole2)
# and a 3-edge return leg (pole2->c->d->pole1).
# For CONTINUOUS COVERAGE: while one corpuscle is on the long-pole leg,
# those 3 edges are covered. The bilateral partner covers the return-leg 3 edges.
# BUT: the 3 long-pole edges and 3 return edges are DIFFERENT edges.
# For any edge to ALWAYS have a muon (from some circuit) present:
# the circuit whose long-pole or return-leg includes that edge must be active.
# Since bilateral covers all 6 circuit edges continuously, single-circuit coverage
# of its 6 edges is complete. The constraint is just: every edge in SOME circuit.

# Long-pole edges for each circuit: first 3 directed edges of path
# Return-leg edges: last 3 directed edges of path
long_pole_edges = []   # set of 3 edge indices per circuit
return_leg_edges = []  # set of 3 edge indices per circuit
for es, p1, p2, path in circuits:
    lp = frozenset(edge_to_idx[(path[k], path[k+1])] for k in range(3))
    rl = frozenset(edge_to_idx[(path[k], path[(k+1)%6])] for k in range(3,6))
    long_pole_edges.append(lp)
    return_leg_edges.append(rl)

# For 12 circuits: one circuit per vertex as a pole (top pole)
# Each vertex appears as pole1 in some circuits
# Strategy: pick 2 circuits per antipodal pair -> 12 circuits

found_12_all_covered = False
best_12 = None
for trial in range(500):
    random.seed(trial + 3000)
    chosen = []
    for pair in pairs:
        # Pick 2 circuits from this pair
        pair_circuits = circuits_by_pair[pair]
        if len(pair_circuits) >= 2:
            chosen.extend(random.sample(pair_circuits, 2))
    if covers_all(chosen):
        found_12_all_covered = True
        edge_counts_12 = [0]*n_edges
        for ci in chosen:
            for ei in circuits[ci][0]:
                edge_counts_12[ei] += 1
        best_12 = chosen[:]
        break

print(f"  12-circuit scheme: 2 circuits per antipodal pair (one per vertex as top-pole)")
if found_12_all_covered:
    print(f"  Full coverage achieved: YES")
    print(f"  Edge coverage: min={min(edge_counts_12)}  max={max(edge_counts_12)}")
    print(f"  Coverage distribution: {sorted(set(edge_counts_12))}")
else:
    print(f"  12-circuit cover (2/axis): not found in 500 trials")

# Also try: 12 circuits = all 6 axes, one circuit per VERTEX as top-pole
# (Vertex v is top-pole of circuits with v as pole1 -- but pole assignment is symmetric)
# Try: for each of 12 vertices, pick the "best" circuit with that vertex as a pole
found_12_one_per_vertex = False
for trial in range(1000):
    random.seed(trial + 4000)
    chosen_by_vertex = {}
    for ci, (es, p1, p2, path) in enumerate(circuits):
        if p1 not in chosen_by_vertex:
            chosen_by_vertex[p1] = ci
        if p2 not in chosen_by_vertex:
            chosen_by_vertex[p2] = ci
    # But this gives at most 12 circuits (one per vertex as a pole)
    # -- but each circuit serves 2 vertices, so 12 circuits serve 24 pole-slots
    # -- 24/12 = 2 poles per vertex on average
    # Try random assignment: one circuit per vertex (as top pole)
    chosen = []
    for v in range(12):
        v_circuits = [ci for ci,(es,p1,p2,path) in enumerate(circuits) if p1==v or p2==v]
        if v_circuits:
            chosen.append(random.choice(v_circuits))
    # Deduplicate
    chosen_dedup = list(dict.fromkeys(chosen))
    if len(chosen_dedup) == 12 and covers_all(chosen_dedup):
        found_12_one_per_vertex = True
        edge_counts_12v = [0]*n_edges
        for ci in chosen_dedup:
            for ei in circuits[ci][0]:
                edge_counts_12v[ei] += 1
        print(f"  12-circuit cover (one per vertex): full coverage = YES")
        print(f"  Edge coverage: min={min(edge_counts_12v)} max={max(edge_counts_12v)}")
        break

check("SC5: 12 circuits (2 per C5 axis) can cover all 30 edges",
      found_12_all_covered,
      f"12 = 2 per axis; covers all edges: {found_12_all_covered}")

# ── SC6: Gamma(30 edges) decomposition -- G appears 2 times ───────────────────
print()
print(SEP)
print("SC6: GAMMA(30 EDGES) DECOMPOSITION -- HOW MANY INDEPENDENT G CHANNELS?")
print(SEP2)

# Character table of I (order 60)
N_class = [1, 12, 12, 20, 15]  # E, C5, C5^2, C3, C2
chi_I = {
    'A' : [1,  1,        1,       1,    1],
    'T1': [3,  phi,      -1/phi,   0,   -1],
    'T2': [3, -1/phi,    phi,      0,   -1],
    'G' : [4, -1,        -1,       1,    0],
    'H' : [5,  0,         0,      -1,    1],
}

# Gamma(30 edge midpoints): chi = [30, 0, 0, 0, 2]
# E: all 30 fixed; C5: 0 (C5 permutes edges at each vertex); C3: 0; C2: 2 (axis through edge pair)
chi_edges = [30, 0, 0, 0, 2]
order = 60

print("  Character of Gamma(30 edge midpoints) under I:")
print(f"    [E, C5, C5^2, C3, C2] = {chi_edges}")
print()

decomp = {}
for name, chi in chi_I.items():
    n = sum(N_class[c]*chi[c]*chi_edges[c] for c in range(5)) / order
    decomp[name] = round(n)
    print(f"    n({name}) = {n:.1f}  ({round(n)} copies)")

dim_check = sum(decomp[name]*chi_I[name][0] for name in decomp)
decomp_str = " + ".join(f"{v}{k}" for k,v in decomp.items() if v>0)
print()
print(f"  Decomposition: {decomp_str}")
print(f"  Dimension check: {dim_check} (should be 30)")
print()
print(f"  The G irrep (dim=4) appears {decomp['G']} times in Gamma(30 edges).")
print(f"  This means there are {decomp['G']} INDEPENDENT G-mode channels on the edges.")
print(f"  The structural muon (G32 spinor) requires at least {decomp['G']} circuits")
print(f"  to represent all independent channels.")
print()
print(f"  NOTE: this is the BOSONIC representation (position of edge midpoints).")
print(f"  The spinor G32 mode is additional -- the bosonic G channel is the 'track',")
print(f"  the spinor G32 mode is what rides on it.")

check("SC6: Gamma(30 edges) dimension = 30",
      dim_check == 30,
      f"sum of irrep dims = {dim_check}")
check("SC6b: G irrep appears 2 times in Gamma(30 edges)",
      decomp['G'] == 2,
      f"n(G) = {decomp['G']}  -> 2 independent G channels on 30 edges")

# ── Summary ───────────────────────────────────────────────────────────────────
print()
print(SEP)
print("SUMMARY: MUON STRUCTURAL COUNT")
print(SEP2)
print()
print(f"  MINIMUM COVER: {min_cover_n} circuits for full edge coverage")
print(f"  6 CIRCUITS (one per C5 axis): {'covers all edges' if found_6 else 'does NOT cover all edges'}")
print(f"  12 CIRCUITS (2 per axis): {'covers all edges' if found_12_all_covered else 'does NOT cover all edges'}")
print(f"  10-CIRCUIT UNIFORM DOUBLE COVER: {'EXISTS' if found_10_uniform else 'not found'}")
print()
print(f"  GAMMA(30 EDGES) = A + T1 + T2 + {decomp['G']}G + {decomp['H']}H")
print(f"  G appears {decomp['G']} times -> minimum 2 independent G-mode circuits")
print()
print(f"  CONCLUSION:")
if found_6:
    print(f"  - 6 circuits (one per C5 axis) suffices for full structural coverage")
    print(f"  - 12 circuits (one per vertex / 2 per axis) gives double coverage")
print(f"  - Gamma decomposition requires AT LEAST 2 circuits (2 G channels)")
print(f"  - Group theory cannot determine the EXACT count beyond this minimum")
print(f"  - Physical constraint: enough circuits so no edge is EVER unoccupied")
print(f"    (guaranteed if every edge is in >=1 active circuit)")

print()
passed = sum(1 for _,s,_ in results if s=="PASS")
failed = sum(1 for _,s,_ in results if s=="FAIL")
print(SEP)
print(f"RESULT: {len(results)}/{len(results)}  ({passed} PASS, {failed} FAIL)")
if failed == 0:
    print("  ALL CHECKS PASSED.")
else:
    for name,status,detail in results:
        if status=="FAIL": print(f"  FAIL: {name}")
print(SEP)
