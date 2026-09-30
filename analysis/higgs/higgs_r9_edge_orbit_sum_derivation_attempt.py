"""
higgs_r9_edge_orbit_sum_derivation_attempt.py
================================================
Second, distinct derivation attempt for the R9 coupling coefficient (see
notes/open_items/false_positive_scan_series1.txt [HIGGS-R9-CHANNEL-WEIGHT]).

The first attempt (higgs_r9_coupling_weight_derivation_attempt.py) proved
that NO purely representation-theoretic construction (characters,
multiplicities, projectors) can ever fix a coupling MAGNITUDE. This
script tries the one genuinely different, well-defined, non-representation-
theoretic idea surfaced during that investigation: higgs_lagrangian_h2.py's
own "the Higgs feels all 30 edges" language, taken LITERALLY and actually
computed (rather than asserted) for the first time.

SETUP: treat the T_1g field as living on the 12 vertices (one 3-vector per
vertex, W(v)), transforming as (g.W)(v) = R(g) W(g^-1 v) -- the natural
"field decorating an orbit" representation. Construct the manifestly
A_g-invariant quadratic form Q(W) = sum over the 30 edges (i,j) of
W(v_i).W(v_j), and see what it becomes when restricted to the copy (or
copies) of the abstract T_1g irrep living inside this 36-dim space --
this is the closest honest, computable meaning of the doc's own claim.

This does NOT assume the answer is phi, phi^2, or anything else -- it
computes whatever the edge-sum operator's eigenvalue actually is on the
physical T_1g subspace, and reports it, whether that turns out to match
anything already in the framework or not.

Run: python analysis/higgs/higgs_r9_edge_orbit_sum_derivation_attempt.py
"""
import numpy as np
import math

phi = (1 + 5 ** 0.5) / 2
SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "*** FAIL"
    results.append((name, cond, detail))
    print(f"  [{s}] {name}")
    if detail:
        print(f"         {detail}")


print(SEP)
print("R9 EDGE-ORBIT-SUM DERIVATION ATTEMPT -- 'the Higgs feels all 30")
print("edges', computed literally instead of asserted")
print(SEP)

# ── SECTION 1: geometry + group (validated construction, reused) ────────────
verts_raw = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        verts_raw.append([0, s1, s2 * phi])
        verts_raw.append([s1, s2 * phi, 0])
        verts_raw.append([s2 * phi, 0, s1])
verts = np.array(verts_raw, dtype=float)
N = len(verts)
all_dists = sorted(np.linalg.norm(verts[i] - verts[j])
                    for i in range(N) for j in range(i + 1, N))
r_nn = all_dists[0]
tol = r_nn * 0.05
edges = [(i, j) for i in range(N) for j in range(i + 1, N)
         if abs(np.linalg.norm(verts[i] - verts[j]) - r_nn) < tol]
V, E = N, len(edges)


def rot_matrix(axis, angle):
    axis = axis / np.linalg.norm(axis)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


gen1 = rot_matrix(verts[0], 2 * math.pi / 5)
i0, j0 = edges[0]
gen2 = rot_matrix(verts[i0] + verts[j0], math.pi)


def mat_key(M):
    return tuple(np.round(M, 6).flatten())


group = {mat_key(np.eye(3)): np.eye(3), mat_key(gen1): gen1, mat_key(gen2): gen2}
changed = True
while changed:
    changed = False
    for a in list(group.values()):
        for b in list(group.values()):
            c = a @ b
            k = mat_key(c)
            if k not in group:
                group[k] = c
                changed = True
rotations = list(group.values())
check("V=12, E=30, group order 60", V == 12 and E == 30 and len(rotations) == 60,
      f"V={V}, E={E}, |G|={len(rotations)}")


def vertex_perm(R):
    perm = []
    for i in range(N):
        rv = R @ verts[i]
        dists = [np.linalg.norm(rv - verts[k]) for k in range(N)]
        perm.append(int(np.argmin(dists)))
    return perm


# ── SECTION 2: build the 36x36 "vertex-field" representation ────────────────
print(f"\nSECTION 2: build the 36-dim vertex-field representation D(g)")
print(SEP2)
D = []
for R in rotations:
    perm = vertex_perm(R)
    Dg = np.zeros((36, 36))
    for v in range(N):
        # (g.W)(v) = R @ W(g^-1 v); g^-1 v = the vertex that maps TO v under perm
        src = perm.index(v)
        Dg[3*v:3*v+3, 3*src:3*src+3] = R
    D.append(Dg)

# verify it's a genuine representation: D(g1)D(g2) = D(g1 g2)
g1, g2 = rotations[1], rotations[2]
Dg1 = D[1]
Dg2 = D[2]
prod_mat = g1 @ g2
rot_index = {mat_key(r): idx for idx, r in enumerate(rotations)}
prod_idx = rot_index[mat_key(prod_mat)]
check("D(g1)D(g2) = D(g1*g2) (genuine representation, not just permutation)",
      np.allclose(Dg1 @ Dg2, D[prod_idx], atol=1e-6),
      f"max diff = {np.max(np.abs(Dg1@Dg2 - D[prod_idx])):.2e}")

# ── SECTION 3: decompose the 36-dim representation into irreps ──────────────
print(f"\nSECTION 3: decompose the 36-dim vertex-field rep into I_h irreps")
print(SEP2)


def classify(R):
    tr = np.trace(R)
    theta = math.degrees(math.acos(np.clip((tr - 1) / 2, -1, 1)))
    if abs(theta) < 1: return 'E'
    if abs(theta - 72) < 2: return 'C5'
    if abs(theta - 144) < 2: return 'C5^2'
    if abs(theta - 120) < 2: return 'C3'
    if abs(theta - 180) < 2: return 'C2'
    return f'?{theta:.1f}'


class_sizes = {'E': 1, 'C2': 15, 'C3': 20, 'C5': 12, 'C5^2': 12}
order_I = 60
chi_36 = {}
for R, Dg in zip(rotations, D):
    c = classify(R)
    chi_36.setdefault(c, []).append(np.trace(Dg))
chi_36 = {c: np.mean(v) for c, v in chi_36.items()}
print(f"  character of 36-dim rep by class: "
      + ", ".join(f"{c}={v:.4f}" for c, v in chi_36.items()))

char = {
    'A_g':  {'E': 1, 'C2': 1,  'C3': 1,  'C5': 1,      'C5^2': 1},
    'T_1g': {'E': 3, 'C2': -1, 'C3': 0,  'C5': phi,    'C5^2': -1/phi},
    'T_2g': {'E': 3, 'C2': -1, 'C3': 0,  'C5': -1/phi, 'C5^2': phi},
    'G_g':  {'E': 4, 'C2': 0,  'C3': 1,  'C5': -1,     'C5^2': -1},
    'H_g':  {'E': 5, 'C2': 1,  'C3': -1, 'C5': 0,      'C5^2': 0},
}
classes = ['E', 'C2', 'C3', 'C5', 'C5^2']
mult = {}
for name, ch in char.items():
    n = sum(class_sizes[c] * chi_36[c] * ch[c] for c in classes) / order_I
    mult[name] = round(n)
    print(f"  multiplicity of {name}: {n:.6f} -> {round(n)}")
dim_check = sum(mult[k] * char[k]['E'] for k in char)
check("36-dim decomposition dimension check", dim_check == 36, f"sum={dim_check}")
print(f"  36-dim vertex field = "
      + " + ".join(f"{v}*{k}" if v > 1 else k for k, v in mult.items() if v > 0))

# ── SECTION 4: the edge-sum quadratic form Q(W) = sum_edges W_i . W_j ───────
print(f"\nSECTION 4: build the edge-sum operator M (36x36), Q(W)=W^T M W")
print(SEP2)
M = np.zeros((36, 36))
for (i, j) in edges:
    M[3*i:3*i+3, 3*j:3*j+3] += np.eye(3)
    M[3*j:3*j+3, 3*i:3*i+3] += np.eye(3)

# verify M commutes with every D(g) (i.e. Q is really A_g-invariant)
max_comm = max(np.max(np.abs(Dg @ M @ Dg.T - M)) for Dg in D[:10])
check("Edge-sum operator M commutes with the group action (genuinely invariant)",
      max_comm < 1e-6, f"max commutator (10-element sample) = {max_comm:.2e}")

# ── SECTION 5: restrict M to the T_1g-isotypic subspace ─────────────────────
print(f"\nSECTION 5: restrict the edge-sum operator to the T_1g-isotypic subspace")
print(SEP2)
n_T1g = mult['T_1g']
print(f"  T_1g appears {n_T1g} time(s) in the 36-dim vertex-field space.")

P_T1g = np.zeros((36, 36))
for R, Dg in zip(rotations, D):
    c = classify(R)
    P_T1g += char['T_1g'][c] * Dg
P_T1g *= char['T_1g']['E'] / order_I  # standard projector normalization

rank_T1g = np.linalg.matrix_rank(P_T1g, tol=1e-6)
check(f"rank(P_T1g projector) = {3*n_T1g} (dim T_1g x multiplicity)",
      rank_T1g == 3 * n_T1g, f"rank={rank_T1g}, expected={3*n_T1g}")

# Get an orthonormal basis for the T_1g isotypic subspace
eigvals, eigvecs = np.linalg.eigh(P_T1g)
basis = eigvecs[:, np.abs(eigvals - 1.0) < 1e-6]
print(f"  T_1g isotypic subspace dimension found: {basis.shape[1]}")

# Project M onto this subspace
M_restricted = basis.T @ M @ basis
print(f"  M restricted to T_1g-isotypic subspace ({M_restricted.shape[0]}x"
      f"{M_restricted.shape[1]}):")
print(np.round(M_restricted, 6))

if n_T1g == 1:
    # Schur's lemma: M restricted to a SINGLE irrep copy must be proportional
    # to identity on that copy (since M commutes with the group action).
    eigval_M = np.linalg.eigvalsh(M_restricted)
    check("M|_T1g is proportional to identity (Schur's lemma, single copy)",
          np.allclose(eigval_M, eigval_M[0], atol=1e-6),
          f"eigenvalues: {np.round(eigval_M, 6)}")
    coupling_value = eigval_M[0]
    print(f"\n  *** GENUINE EDGE-SUM COUPLING VALUE: {coupling_value:.8f} ***")
    for name, val in [("phi", phi), ("phi^2", phi**2), ("1", 1.0),
                       ("2*phi", 2*phi), ("5 (num edges/6)", 5.0),
                       ("30 (num edges)", 30.0), ("-1 (chi(T1g,C2))", -1.0),
                       ("6 (edges per vertex... no, 5)", 6.0)]:
        print(f"    vs {name:<20} = {val:.8f}  ratio = {coupling_value/val:.6f}")
else:
    print(f"  T_1g appears {n_T1g} times -- by Schur's lemma, M restricted must")
    print(f"  be of the form A_mult (x) I_3 for some {n_T1g}x{n_T1g} matrix A_mult:")
    print(f"  each eigenvalue of A_mult appears with multiplicity 3 (one per")
    print(f"  T_1g vector component). Diagonalizing confirms/extracts this:")
    eigvals_full = np.linalg.eigvalsh(M_restricted)
    print(f"  Full 9 eigenvalues: {np.round(eigvals_full, 6)}")
    # group into clusters of 3 (degenerate sets)
    distinct = []
    for ev in eigvals_full:
        if not any(abs(ev - d) < 1e-4 for d in distinct):
            distinct.append(ev)
    check("Eigenvalues cluster into 3 distinct values, each degeneracy 3 "
          "(confirms Schur block structure A_mult (x) I_3)",
          len(distinct) == 3 and all(
              sum(1 for ev in eigvals_full if abs(ev - d) < 1e-4) == 3
              for d in distinct),
          f"distinct values found: {np.round(sorted(distinct), 6)}")
    print(f"\n  *** THE {len(distinct)} GEOMETRIC EDGE-SUM COUPLING VALUES ***")
    print(f"  (one of these -- which one requires further physical input to")
    print(f"  select -- would be 'the' W-boson field's own coupling strength)")
    for ev in sorted(distinct):
        print(f"\n  Candidate value: {ev:.8f}")
        for name, val in [("phi", phi), ("phi^2", phi**2), ("1", 1.0),
                           ("2*phi", 2*phi), ("-1 (chi(T1g,C2))", -1.0),
                           ("5", 5.0), ("-1/phi", -1/phi), ("3", 3.0)]:
            print(f"    vs {name:<20} = {val:>10.6f}  ratio = {ev/val if val != 0 else float('nan'):.6f}")

# ── SUMMARY ───────────────────────────────────────────────────────────────────
print(f"\n{SEP}")
print("SUMMARY / VERDICT")
print(SEP)
n_pass = sum(1 for _, c, _ in results if c)
print(f"  {n_pass}/{len(results)} checks passed")
print()
print("  This is a genuinely different, actually-computed (not asserted)")
print("  attempt at 'the Higgs feels all 30 edges' -- whatever numeric value")
print("  came out above is reported honestly, whether or not it matches any")
print("  existing framework constant.")
