"""Closed-book check: can a harmonic/jamming lattice route give F = E_cell*L_J/hbar_c = 2*pi?

Proposed by the "Opus 5.5" arm of the 5-model closed-book review dispatched to
test whether Claim 8's underlying question ("can E_cell emerge from an EW/QCD-
scale jamming argument, no density needed") has ANY genuine answer, now that
Claim 8's own specific formula is independently confirmed as a reverse-fit.
That model could not run this itself (read-only tools that session) -- this
run is the independent verification of its proposed no-go arguments.
"""
import math, itertools
import numpy as np

pi, phi = math.pi, (1 + math.sqrt(5)) / 2
Rs = math.sqrt(5) / (4 * pi)
KG = 1 / Rs**2 - 4 / 3
res = []
def check(n, c, d=""):
    res.append(c); print(f"  [{'PASS' if c else 'FAIL'}] {n}  {d}")

# N1: A_g breathing eigenvalue of unit-spring/unit-mass icosahedron
V = np.array([v for s1 in (1, -1) for s2 in (1, -1)
              for v in ([0, s1, s2*phi], [s1, s2*phi, 0], [s2*phi, 0, s1])], float)
d = np.linalg.norm(V[:, None] - V[None], axis=-1)
a = d[d > 1e-9].min()
edges = [(i, j) for i in range(12) for j in range(i+1, 12) if abs(d[i, j] - a) < 1e-9]
D = np.zeros((36, 36))
for i, j in edges:
    P = np.outer(V[j] - V[i], V[j] - V[i]) / a**2
    for p, q, s in ((i, i, 1), (j, j, 1), (i, j, -1), (j, i, -1)):
        D[3*p:3*p+3, 3*q:3*q+3] += s * P
b = (V / np.linalg.norm(V, axis=1)[:, None]).ravel(); b /= np.linalg.norm(b)
lam = b @ D @ b
check("N1a: 30 edges", len(edges) == 30)
check("N1b: breathing is exact eigenvector", np.linalg.norm(D @ b - lam*b) < 1e-12)
check("N1c: lambda_Ag = 5 - sqrt5", abs(lam - (5 - math.sqrt(5))) < 1e-12, f"{lam:.12f}")
print("  full spectrum:", np.round(np.unique(np.round(np.linalg.eigvalsh(D), 10)), 8))

# N2: monatomic Bravais NN-spring bound F_max = omega_max*a_nn/v_max < 2*pi
def F_max(prim, nn):
    prim, nn = np.array(prim, float), np.array(nn, float)
    a_nn = np.linalg.norm(nn[0]); Rh = nn / np.linalg.norm(nn, axis=1)[:, None]
    PP = np.einsum('ni,nj->nij', Rh, Rh)
    Dq = lambda q: np.einsum('n,nij->ij', 1 - np.cos(nn @ q), PP)
    Dc = lambda q: np.einsum('n,nij->ij', 0.5 * (nn @ q)**2, PP)
    rng = np.random.default_rng(0)
    dirs = rng.normal(size=(4000, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
    v_max = math.sqrt(max(np.linalg.eigvalsh(Dc(n))[-1] for n in dirs))
    B = 2 * pi * np.linalg.inv(prim).T
    g = np.linspace(0, 1, 25, endpoint=False)
    w2 = max(np.linalg.eigvalsh(Dq(np.array(x) @ B))[-1] for x in itertools.product(g, g, g))
    gap = min(np.linalg.eigvalsh(Dc(q) - Dq(q))[0] for q in rng.normal(size=(500, 3)) * 5)
    return math.sqrt(w2) * a_nn / v_max, gap
sc  = F_max(np.eye(3), [s*e for e in np.eye(3) for s in (1, -1)])
bcc = F_max([[-.5, .5, .5], [.5, -.5, .5], [.5, .5, -.5]],
            [[x, y, z] for x in (.5, -.5) for y in (.5, -.5) for z in (.5, -.5)])
fcc = F_max([[0, .5, .5], [.5, 0, .5], [.5, .5, 0]],
            [p for x in (.5, -.5) for y in (.5, -.5)
             for p in ([x, y, 0], [x, 0, y], [0, x, y])])
for n, (F, gap) in (("SC", sc), ("BCC", bcc), ("FCC", fcc)):
    check(f"N2 {n}: D(q) <= D_cont(q) (PSD)", gap > -1e-9, f"min eig gap={gap:.2e}")
    check(f"N2 {n}: F_max < 2*pi", F < 2*pi, f"F_max={F:.6f} vs 2pi={2*pi:.6f}")

# N3: the only pi sources in the inputs are Rs/K/G -- pure rewriting, not mechanism
check("N3a: 2*pi == sqrt5/(2*Rs) (definitional)", abs(math.sqrt(5)/(2*Rs) - 2*pi) < 1e-12)
check("N3b: 2*pi == sqrt(5*(K/G+4/3))/2 (definitional)", abs(math.sqrt(5*(KG + 4/3))/2 - 2*pi) < 1e-12)
print(f"RESULTS: {sum(res)}/{len(res)} PASS")
