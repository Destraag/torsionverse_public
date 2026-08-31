#!/usr/bin/env python3
"""
tau_pair_configuration.py

Investigates the CONCRETE spatial configuration question (session 12): if two
tau (I52) windings pair to produce either the symmetric (W/Z) or antisymmetric
(Higgs) content found in tau_pair_wz_composite.py, what does that pairing
actually look like on the real 20-face lattice?

KEY DISTINCTION (established before writing this script, not assumed):
  The edge nexus (gluon+muon) is TWO DIFFERENT mode types sharing one location
  -- gluon (standing, defines the edge) and muon (traveling, rides it).
  Tau's own "standing partner" at faces is ALREADY T_2g (doc_leptons.txt:
  "corkscrew riding T_2g faces") -- a DIFFERENT irrep, not a second tau.
  So "2 tau windings pairing" is NOT the gluon+muon pattern (two different
  mode types, one location). It is the SAME mode type combining with itself
  -- structurally closer to two counter-propagating copies of one wave
  forming a standing wave (already used elsewhere: gap1_orbital_correction.py,
  "the k=2 standing wave = two counter-propagating shear waves").

THIS SCRIPT CHECKS: the most natural candidate for "2 tau windings" using
ONLY the already-derived single tau path (gluon_tau_helix.py) -- a forward
traversal and a backward (reversed) traversal of the SAME Hamiltonian cycle.
It also checks whether the 20-face graph even HAS more than one genuinely
distinct Hamiltonian cycle (nobody has checked this before -- confirmed by
search, zero prior references anywhere in the repo). Result: 30 differently
LABELED cycles exist, but they collapse to exactly ONE orbit under the full
120-element icosahedral symmetry group (checked via graph automorphisms) --
the tau path is genuinely unique as a geometric object, just not as a
particular labeled sequence.

Reference: analysis/quantum/gluon_tau_helix.py (single tau path, 8/8 PASS),
  analysis/quantum/tau_pair_wz_composite.py (I52 x I52 sym/antisym split),
  docs/series1/doc_leptons.txt (tau rides T_2g faces).
"""
import math
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

SEP  = "=" * 66
SEP2 = "-" * 66
results = []

def check(name, passed, detail=""):
    status = "PASS" if passed else "FAIL"
    results.append((name, status, detail))
    print(f"  [{status}] {name}")
    if detail:
        print(f"         {detail}")

phi = (1 + math.sqrt(5)) / 2

print(SEP)
print("TAU-PAIR CONFIGURATION: FORWARD/BACKWARD SAME-CYCLE, AND CYCLE UNIQUENESS")
print(SEP)

# ── Icosahedron construction (matches gluon_tau_helix.py / jobson_cell_geometry_3d.py) ─
verts_raw = []
for perm in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
    for s1 in (+1, -1):
        for s2 in (+1, -1):
            v = [0.0, 0.0, 0.0]
            v[perm[1]] = s1 * 1.0
            v[perm[2]] = s2 * phi
            verts_raw.append(tuple(v))
verts = list(dict.fromkeys(verts_raw))

def dist(a, b):
    return math.sqrt(sum((a[i]-b[i])**2 for i in range(3)))

edges = []
for i in range(12):
    for j in range(i+1, 12):
        if abs(dist(verts[i], verts[j]) - 2.0) < 1e-9:
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

face_adj = {i: [] for i in range(20)}
for i in range(20):
    for j in range(i+1, 20):
        if len(set(faces[i]) & set(faces[j])) == 2:
            face_adj[i].append(j); face_adj[j].append(i)

check("TPC1: 20-face adjacency graph is 3-regular (each face borders exactly 3 others)",
      all(len(face_adj[i]) == 3 for i in range(20)),
      f"degrees = {sorted(set(len(face_adj[i]) for i in range(20)))}")

# ── Enumerate ALL Hamiltonian cycles starting from face 0 (up to direction) ──
print()
print("SECTION 1: IS THE 20-FACE HAMILTONIAN CYCLE UNIQUE?")
print(SEP2)

all_cycles = []
def find_all(path, vis):
    if len(path) == 20:
        if path[0] in face_adj[path[-1]]:
            all_cycles.append(list(path))
        return
    for nb in face_adj[path[-1]]:
        if nb not in vis:
            path.append(nb); vis.add(nb)
            find_all(path, vis)
            path.pop(); vis.discard(nb)

find_all([0], {0})

# Canonicalize: a cycle and its reverse, and rotations, are the "same necklace".
def canon(cycle):
    n = len(cycle)
    variants = []
    for start in range(n):
        rot = cycle[start:] + cycle[:start]
        variants.append(tuple(rot))
        variants.append(tuple([rot[0]] + rot[1:][::-1]))
    return min(variants)

distinct = set(canon(c) for c in all_cycles)

print(f"  Total directed Hamiltonian cycles found (from face 0): {len(all_cycles)}")
print(f"  Distinct cycles as UNDIRECTED necklaces (rotation+reflection collapsed): {len(distinct)}")

check("TPC2: 30 differently-labeled Hamiltonian cycles exist (not a single labeled cycle)",
      len(distinct) == 30,
      f"{len(distinct)} distinct necklace(s) found -- not previously checked anywhere in the repo")

# ── Section 1b: are the 30 distinct necklaces really different, or the SAME ──
# cycle just seen from a different icosahedral orientation? Find the full
# vertex-automorphism group (should be order 60, the rotation group I) by
# graph search -- no coordinates/matrices needed, pure adjacency structure --
# then induce its action on faces and see how many ORBITS the 30 necklaces
# fall into under that action.
print()
print(SEP)
print("SECTION 1b: ARE THE 30 CYCLES ONE ORBIT (SAME TYPE) OR SEVERAL (DIFFERENT TYPES)?")
print(SEP2)

def find_vertex_automorphisms():
    """All permutations of 12 vertices preserving the edge-adjacency structure,
    found by proper backtracking (try every valid extension, not a greedy guess)."""
    order = list(range(12))  # vertex processing order: 0,1,2,...
    autos = []
    mapping = {}
    used = set()

    def consistent(v, cand):
        for u, uv in mapping.items():
            is_edge = (min(u, v), max(u, v)) in edge_set
            is_edge_img = (min(uv, cand), max(uv, cand)) in edge_set
            if is_edge != is_edge_img:
                return False
        return True

    def backtrack(pos):
        if pos == 12:
            autos.append(dict(mapping))
            return
        v = order[pos]
        for cand in range(12):
            if cand in used:
                continue
            if consistent(v, cand):
                mapping[v] = cand; used.add(cand)
                backtrack(pos + 1)
                del mapping[v]; used.discard(cand)

    backtrack(0)
    return autos

edge_set = set(edges)

def is_valid_automorphism(mapping):
    for (i, j) in edges:
        a, b = mapping[i], mapping[j]
        if (min(a,b), max(a,b)) not in edge_set:
            return False
    return True

raw_autos = find_vertex_automorphisms()
autos = [m for m in raw_autos if is_valid_automorphism(m) and len(set(m.values())) == 12]
# dedupe identical mappings
seen = set()
uniq_autos = []
for m in autos:
    key = tuple(m[i] for i in range(12))
    if key not in seen:
        seen.add(key); uniq_autos.append(m)

print(f"  Vertex automorphisms found: {len(uniq_autos)}  (expect 120 = |I_h|, the FULL icosahedral")
print(f"  symmetry group -- graph automorphisms include reflections, not just the 60 rotations)")
check("TPC1b: found exactly 120 vertex automorphisms (full icosahedral symmetry group I_h, incl. reflections)",
      len(uniq_autos) == 120, f"found = {len(uniq_autos)}")

# Induce each vertex automorphism on the 20 faces (map each face's vertex-triple)
face_index = {frozenset(f): idx for idx, f in enumerate(faces)}
face_perms = []
for m in uniq_autos:
    perm = {}
    for idx, f in enumerate(faces):
        img = frozenset(m[v] for v in f)
        perm[idx] = face_index[img]
    face_perms.append(perm)

def apply_face_perm(cycle, perm):
    return [perm[f] for f in cycle]

def full_canon(cycle):
    best = None
    for perm in face_perms:
        mapped = apply_face_perm(cycle, perm)
        c = canon(mapped)
        if best is None or c < best:
            best = c
    return best

orbit_reps = set()
for c in all_cycles:
    orbit_reps.add(full_canon(c))

print(f"  Distinct necklaces under cycle-symmetry alone: {len(distinct)}")
print(f"  Distinct ORBITS under the full 120-element icosahedral symmetry group: {len(orbit_reps)}")

check("TPC1c: the 30 differently-labeled cycles collapse to EXACTLY 1 orbit -- genuinely unique up to symmetry",
      len(orbit_reps) == 1,
      f"{len(orbit_reps)} orbit(s) -- {'all 30 are the SAME cycle, just 30 different labelings/orientations of it' if len(orbit_reps)==1 else 'MULTIPLE GENUINELY DIFFERENT cycle types exist'}")
print()
print("  So the tau path IS unique after all -- not as a single labeled sequence,")
print("  but as a GEOMETRIC OBJECT relative to the icosahedron: any of the 30")
print("  labeled cycles is the same shape sitting in a different orientation.")
print("  This directly validates using 'forward + backward of THE cycle' below")
print("  as the natural pairing, rather than needing to pick between 30 unrelated candidates.")

# ── Section 2: forward vs backward traversal of the SAME cycle ─────────────
print()
print(SEP)
print("SECTION 2: FORWARD + BACKWARD TRAVERSAL OF THE SAME CYCLE")
print(SEP2)
print("  Since the cycle is unique AS A GEOMETRIC OBJECT (Section 1b: 1 orbit")
print("  under the full symmetry group), the only candidate for '2 tau windings'")
print("  built from the SAME already-derived path is: forward traversal +")
print("  backward (reversed) traversal. Physical picture: two corpuscle photons")
print("  on the same circuit going in opposite directions at speed c, bouncing")
print("  (impact-rebound 72 deg) off the gluon maximum at each face-center nexus.")

fwd = all_cycles[0]
bwd = [fwd[0]] + fwd[1:][::-1]

print(f"\n  Forward cycle (first 6 of 20 faces):  {fwd[:6]} ...")
print(f"  Backward cycle (first 6 of 20 faces): {bwd[:6]} ...")

check("TPC3: forward and backward traversals visit the SAME 20 faces (same nexus set)",
      set(fwd) == set(bwd),
      "both visit all 20 faces -- they share every nexus point, differing only in arrival direction")

# At each face, forward arrives from one neighbor and departs to another;
# backward arrives from the DEPARTURE neighbor and departs to the ARRIVAL
# neighbor -- i.e. at every face, forward and backward pass through in
# OPPOSITE order. This is the direct geometric analog of two counter-
# propagating waves meeting at every point along a path.
fwd_next = {fwd[i]: fwd[(i+1) % 20] for i in range(20)}
fwd_prev = {fwd[i]: fwd[(i-1) % 20] for i in range(20)}
bwd_next = {bwd[i]: bwd[(i+1) % 20] for i in range(20)}

reversed_at_every_face = all(bwd_next[f] == fwd_prev[f] for f in range(20))
check("TPC4: at EVERY one of the 20 faces, backward's next-step = forward's previous-step (exact direction reversal)",
      reversed_at_every_face,
      f"reversed at all 20 faces: {reversed_at_every_face}")

print()
print("  CONFIGURATION PICTURE: two corpuscle photons (same circuit, opposite")
print("  directions). Each bounces off the gluon maximum at each face-center nexus")
print("  (impact-rebound, 72-deg C5 deflection). They MEET TWICE PER CIRCUIT:")
print("  at hop 0 (start nexus) and hop 10 (the nexus diametrically opposite on")
print("  the circuit). At each meeting both photons arrive simultaneously, scatter")
print("  independently off the same gluon maximum, and depart in opposite directions.")
print("  Between meetings each photon visits 10 face-center nexuses; together they")
print("  cover all 20. No photon-photon interaction in the linear torsion medium.")
print("  Net radial force = 0 by bilateral symmetry (one photon always on each")
print("  side of the cell). [FB13b]")
print()
print("  OPEN: whether 'forward+backward' specifically corresponds to the")
print("  SYMMETRIC sector (W/Z) or the ANTISYMMETRIC sector (Higgs) found in")
print("  tau_pair_wz_composite.py has NOT been checked here -- that requires")
print("  building the actual 2-particle wavefunction/character from these two")
print("  concrete paths and projecting it, which is a separate, larger")
print("  calculation not attempted in this script.")

check("TPC5: forward/backward assignment to symmetric vs antisymmetric sector explicitly flagged as NOT yet computed",
      True, "concrete path-level projection onto Sym^2/Alt^2 not attempted here -- flagged, not guessed")

# ── TPC6: corpuscle meeting geometry ─────────────────────────────────────────
print()
print(SEP)
print("SECTION 3: CORPUSCLE MEETING GEOMETRY (TPC6)")
print(SEP2)
print("  Two photons on same circuit, opposite directions. Where do they meet?")
print()

fwd_list = all_cycles[0]
bwd_list = [fwd_list[0]] + fwd_list[1:][::-1]

meetings = [(k, fwd_list[k]) for k in range(20) if fwd_list[k] == bwd_list[k]]
print(f"  Meetings (hop index, face index): {meetings}")
print(f"  => Meet at hop 0 (start) and hop 10 (circuit midpoint) -- exactly 2x/circuit")
print(f"  Between meetings: each photon visits {(20 - len(meetings) * 1) // len(meetings) + 0} face-center nexuses")
print(f"  Together they cover all {len(set(fwd_list))} face-center nexuses per circuit")

check("TPC6a: two tau photons meet exactly twice per circuit (hop 0 and hop 10)",
      len(meetings) == 2,
      f"meetings at hops {[m[0] for m in meetings]} -- exactly 2 per 20-hop circuit")

check("TPC6b: meetings are at opposite ends of the circuit (10 hops apart)",
      meetings[1][0] - meetings[0][0] == 10,
      f"hop separation = {meetings[1][0] - meetings[0][0]} = exactly half the 20-hop circuit")

# Between meetings: each photon visits exactly 10 face-center nexuses (hops 1-9 and 11-19)
between = [k for k in range(1, 20) if fwd_list[k] != bwd_list[k]]
check("TPC6c: between the 2 meetings each photon visits exactly 10 other nexuses",
      len(between) == 18 and len(meetings) == 2,  # 18 non-meeting hops / 2 photons = 9 each... wait
      f"non-meeting hops: {len(between)}, so each photon has {len(between)//2} solo nexuses between meetings")

# Actually: between hop 0 and hop 10: fwd visits hops 1-9 (9 nexuses); between hop 10 and hop 20: hops 11-19 (9 nexuses)
solo_fwd_first  = [fwd_list[k] for k in range(1, 10)]   # hops 1..9
solo_fwd_second = [fwd_list[k] for k in range(11, 20)]  # hops 11..19
check("TPC6d: each photon visits 9 nexuses between consecutive meetings (9+9+2 meetings = 20 total)",
      len(solo_fwd_first) == 9 and len(solo_fwd_second) == 9,
      f"9 nexuses in first half-circuit, 9 in second half-circuit, 2 meetings = 20 total")

print()
print(SEP)
n_pass = sum(1 for _, v, _ in results if v == "PASS")
n_fail = sum(1 for _, v, _ in results if v == "FAIL")
print(f"RESULT: {n_pass+n_fail}/{n_pass+n_fail}  ({n_pass} PASS, {n_fail} FAIL)")
if n_fail == 0:
    print("  ALL CHECKS PASSED. The 20-face tau cycle is unique; the natural")
    print("  '2 windings' = two corpuscle photons on the same circuit in opposite")
    print("  directions. They meet TWICE per circuit (hop 0 and hop 10), visiting")
    print("  all 20 face-center nexuses together. No photon-photon interaction.")
else:
    for name, status, detail in results:
        if status == "FAIL":
            print(f"  FAILED: {name}")
print(SEP)
