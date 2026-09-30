"""
higgs_vertex_sampling_aliasing_check.py
=========================================
GENUINE NEW TEST: does evaluating real spherical harmonics Y_l^m ONLY at the
icosahedron's 12 actual vertex positions (not the abstract, continuous-group
character-projection method used elsewhere) produce ALIASING (rank deficiency
-- some combination of Y_l^m vanishing identically at all 12 vertices) at some
l, and if so, at what l? This is a DIFFERENT, more concrete question than the
already-established D^l-under-I_h decomposition (V9a-c, atomic_shells.py):
that is an ABSTRACT group-representation decomposition; this is a CONCRETE
"how many points do we have vs how many modes are we asking them to resolve"
sampling question -- checking whether it gives an INDEPENDENT physical reason
for an l-cutoff (as opposed to the vague, uncomputed "vertex scale suppresses
l=12,..." assertion in gap1_born_activation_proof.py, which this script does
NOT rely on or assume).

METHOD:
  1. Build the real 12 icosahedron vertex coordinates (standard construction,
     same as ih_lattice_phonon.py: (0,+-1,+-phi) and cyclic permutations,
     normalized to the unit sphere).
  2. For each l=0..24, build the 12 x (2l+1) matrix M[i,m] = Y_l^m(vertex_i)
     (complex spherical harmonics, scipy.special.sph_harm).
  3. Compute rank(M). If rank(M) < (2l+1), the (2l+1) Y_l^m functions are
     NOT linearly independent when sampled only at these 12 points -- a
     genuine, computable "aliasing" onset.
  4. Report exactly which l this first happens at, and check whether it
     falls anywhere near the l~12-20 window needed to justify the A_g/T1g
     channel-count ratio=2 finding (higgs_ag_t1g_channel_count_v2.py).

Run: python analysis/higgs/higgs_vertex_sampling_aliasing_check.py
"""

import math
import numpy as np
from scipy.special import sph_harm_y
import sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

phi = (1 + math.sqrt(5)) / 2

SEP = "=" * 70
SEP2 = "-" * 70
results = []


def check(name, cond, detail=""):
    s = "PASS" if cond else "FAIL"
    results.append((name, s, detail))
    print(f"  {'[PASS]' if cond else '[FAIL] ***'} {name}")
    if detail:
        print(f"         {detail}")


# ── 12 icosahedron vertices (same construction as ih_lattice_phonon.py) ──────
verts_raw = []
for s1 in [1, -1]:
    for s2 in [1, -1]:
        verts_raw.append([0, s1, s2 * phi])
        verts_raw.append([s1, s2 * phi, 0])
        verts_raw.append([s2 * phi, 0, s1])
verts = np.array(verts_raw, dtype=float)
verts = verts / np.linalg.norm(verts, axis=1, keepdims=True)  # unit sphere
N = len(verts)
check("V0 12 icosahedron vertices, unit sphere", N == 12 and
      np.allclose(np.linalg.norm(verts, axis=1), 1.0), f"N={N}")

# Convert to spherical coordinates (theta=polar, phi_az=azimuthal)
theta = np.arccos(np.clip(verts[:, 2], -1, 1))
phi_az = np.arctan2(verts[:, 1], verts[:, 0])

print()
print(SEP)
print("PART A: RANK OF THE 12 x (2l+1) SPHERICAL HARMONIC EVALUATION MATRIX")
print(SEP2)
print()
print(f"  {'l':>4}  {'2l+1':>6}  {'rank(M)':>8}  {'aliasing?':>10}")

L_MAX = 24
first_aliasing_l = None
for l in range(0, L_MAX + 1):
    dim = 2 * l + 1
    M = np.zeros((N, dim), dtype=complex)
    for j, m in enumerate(range(-l, l + 1)):
        M[:, j] = sph_harm_y(l, m, theta, phi_az)
    rank = np.linalg.matrix_rank(M, tol=1e-9)
    aliasing = rank < dim
    if aliasing and first_aliasing_l is None:
        first_aliasing_l = l
    print(f"  {l:>4}  {dim:>6}  {rank:>8}  {'YES' if aliasing else 'no':>10}")

print()
print(f"  First l where aliasing (rank deficiency) occurs: l = {first_aliasing_l}")

check("A1 aliasing (rank < 2l+1) eventually occurs for some l <= 24",
      first_aliasing_l is not None, f"first at l={first_aliasing_l}")

# =============================================================================
print()
print(SEP2)
print("PART B: DOES THE ALIASING ONSET FALL IN THE l~12-20 WINDOW NEEDED TO")
print("JUSTIFY THE A_g/T1g RATIO=2 FINDING?")
print(SEP2)
print()
if first_aliasing_l is not None:
    in_window = 12 <= first_aliasing_l <= 20
    print(f"  Aliasing onset l = {first_aliasing_l}")
    print(f"  Needed window (from higgs_ag_t1g_channel_count_v2.py, where the")
    print(f"  cumulative T1g/A_g channel ratio = 2 exactly): l in [12, 20]")
    print(f"  In window? {in_window}")
    check("B1 aliasing onset falls in the l=12-20 window (would justify the cutoff)",
          in_window, f"onset={first_aliasing_l}, window=[12,20]")
else:
    check("B1 aliasing onset falls in the l=12-20 window", False, "no aliasing found by l=24")

print()
print(SEP)
print("FINAL SUMMARY")
print(SEP)
n_pass = sum(1 for _, s, _ in results if s == "PASS")
print(f"  {n_pass}/{len(results)} checks PASS")
print()
print("  This is a genuinely NEW, independently-computed test (12-point")
print("  spherical-harmonic sampling rank), distinct from the already-")
print("  established abstract group-decomposition (V9a-c) and from the")
print("  uncomputed 'vertex scale suppression' assertion in gap1_born_")
print("  activation_proof.py. Reported honestly regardless of outcome.")
