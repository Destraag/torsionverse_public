# TORSIONVERSE — AI FRAMEWORK REFERENCE (v1, 2026-09-23)

Condensed orientation for any AI session (any model) starting fresh on this
repo. Read this FIRST, before reading the full docs — it should give enough
context to understand a request and know where to look for detail, without
loading the entire `docs/` + `analysis/` corpus.

**This file is an INDEX, not a source of truth.** Every number and claim
below is a pointer to where it's actually derived/computed. If this file
and the cited doc/script ever disagree, THE SCRIPT WINS — see EPISTEMIC
POLICY below.

---

## 0. EPISTEMIC POLICY (read this first)

1. **Math and computed results outrank prose.** This repo's own working
   history (see `docs/series1/judgment_calls.txt`) has repeatedly found
   cases where a doc's narrative sentence was imprecise, stale, or
   overclaiming while the underlying script's actual numbers were fine
   (and occasionally the reverse). When a doc's prose and a script's
   computed output disagree, or when a doc cites a PASS count / value that
   doesn't match what the cited script currently prints: **trust the
   script's live output**, not the doc's sentence describing it. Re-run
   the script if in doubt — every calculation in this repo lives in a
   committed, runnable script (the repo's own ABSOLUTE RULE), never in an
   inline one-off calculation, specifically so this is always possible.
2. **A "PASS" script result is stronger evidence than a confident-sounding
   paragraph.** Group-theory/topological facts (characters, dimensions,
   winding numbers) and closed-form algebra are the most trustworthy tier;
   narrative analogies ("this is like a screw tightening") are the least
   trustworthy tier and should be treated as illustrative, not derived,
   unless a script backs the specific numeric claim.
3. **Keep this file up to date.** If you find this file conflicts with the
   actual state of the docs or scripts (a constant changed, a status
   flipped from open to resolved, a file got renamed/moved), UPDATE this
   file in the same session you find the conflict — don't just work
   around the discrepancy silently. A stale reference file is worse than
   no reference file, because it's trusted at face value more readily
   than a doc buried three directories deep.

---

## 1. THE CORE IDEA, IN ONE PARAGRAPH

Torsionverse proposes that the vacuum is a real, discrete, granular
(contact-force, **not** Hookean/elastic-continuum) medium made of tiny
icosahedral cells ("Jobson cells"), packed together and normally free-
flowing (like a fluid, shear modulus G→0) except where an active
disturbance pushes local cells past a rigidity/jamming threshold
(Maxwell criticality), at which point they lock elastic. Particles are
topological windings (Hopf-fibration-type (p,q) torus knots) threading
through this medium; particle properties (mass, charge, spin) emerge from
the winding's own topology and its interaction with the medium's
icosahedral (I_h point group) symmetry — not from separately-posited
quantum numbers. The framework claims to derive real, measured physical
constants (alpha, particle masses, G, etc.) from this geometry with few or
no free parameters, cross-checked against real experimental data
throughout.

---

## 2. CANONICAL CONSTANTS (single source — currently duplicated ad hoc in
## dozens of scripts; if you add a new script, prefer importing from
## `analysis/higgs/constants.py` or reproducing these exact values)

| Symbol | Value | Definition / origin |
|---|---|---|
| `alpha` | 7.2973525693e-3 | Fine-structure constant. CODATA-matching; ALSO independently derived from (1,2) Hopf-winding topology elsewhere (doc_alpha.txt) — the CODATA value and the derived value agree to ~1e-10 relative, so scripts freely use the CODATA constant without it being a fitted input. |
| `phi` | (1+√5)/2 = 1.6180339887 | Golden ratio. Appears via icosahedral geometry — e.g. `chi(T_1g, C5) = phi` exactly (C5-axis character of the T_1g irrep), `2*r_mid = phi*L_J` (edge-midpoint contact geometry). |
| `hbar_c` | 197.3269804 MeV·fm | Standard CODATA. |
| `m_p` | 938.272 MeV | Proton mass, standard. |
| `r_p` | 0.8414e-15 m (0.8414 fm) | Proton charge radius, CODATA. |
| `lambda_p` | hbar_c / m_p | Proton Compton-type length scale. |
| `L_J` | alpha·phi·r_p ≈ 9.927e-3 fm | **The Jobson cell's own edge length.** Fixed — no work anywhere in this repo supports an individual cell changing geometric size. Do not confuse with the LOCAL CELL *SPACING*, which does vary with local packing density (see doc_redshift.txt's void/filament mechanism, `L_J_local ~ n_cells^(-1/3)`). |
| `R_c` | L_J·√(1+phi²)/2 | Circumradius (center-to-vertex) of one Jobson cell. |
| `r_in` | L_J·phi²/(2√3) | Inradius (center-to-face) of one Jobson cell. |
| `r_mid` | L_J·phi/2 | Center-to-edge-midpoint radius. `2*r_mid = phi*L_J` exactly — this is the real, algebraically-exact **inter-cell center-to-center spacing** when two cells touch at edge midpoints (cell_rotation_propagation.py RP1a; doc_higgs.txt Section 7.1). **This is phi*L_J, NOT L_J** — a distinction that has caused real confusion more than once (see judgment_calls history), because several scripts use a simplified/schematic lattice constant of L_J directly instead. |
| `E_cell` | 2π·hbar_c/L_J ≈ 124.799 GeV | The Jobson cell's characteristic energy scale ("zone-boundary A_g phonon" energy, per doc_higgs.txt — treat this label as narratively imprecise per standard lattice-dispersion convention, but the NUMBER itself is empirically solid: every "standard convention" alternative tried misses the measured Higgs mass by 400+ sigma). |
| `Rs` | √5/(4π) ≈ 0.35682 | Independently-derived geometric ratio, cross-validated at 4 unrelated physical scales (nuclear, hadronic, galactic, flyby) with zero adjusted free parameters — one of the framework's strongest, most load-bearing constants. Sets `v_s = Rs*c` (medium shear/transverse wave speed) and appears throughout (MOND a0, m_n-m_p splitting, etc.). |
| `v_p` | c (exactly, by construction/GW170817 bound) | Medium longitudinal wave speed. |
| `v_s` | Rs·c | Medium transverse/shear wave speed. |
| `K` (EM sector) | 1/eps_0 | Medium bulk modulus, EM sector. **Not the same K as gravity's own derivation** — gravity's own G comes from a purely topological exponent (below), not from K at all; do not casually reuse this K for gravity-sector calculations (a real, previously-committed mistake — see `lj_local_K_modulus_consistency_check.py`). |
| `rho` | mu_0 | Medium mass-density analog, EM sector. |
| `alpha_grav` | (m_p/E_cell)^18 | Gravity's coupling, purely topological/dimensional-counting derivation (NOT via K or V_p·K_eff). `n=18 = dim(T_1g)·[dim(T_1g)+dim(T_2g)]` — the same exponent as the Maxwell-rigidity/isostatic jamming transition of one Jobson cell (independently validated, ~0.27% match to measured G). |
| `N_lock` | 2π/(alpha·phi) ≈ 532.1 | Tube-closure number: how many (1,2)-Hopf-winding arc segments close one toroidal loop. |

---

## 3. CORE GEOMETRIC / GROUP-THEORY PICTURE

- **Jobson cell** = one icosahedron (V=12, E=30, F=20), point group **I_h**
  (order 120) or its rotation-only subgroup **I** (order 60).
- **I_h irreps** (rotation-only labels used most often): `A_g` (trivial,
  1D, χ=1 for every class — the scalar/Higgs mode), `T_1g`/`T_2g` (3D
  vector-type, W/Z and nucleon-diquark-adjacent modes respectively),
  `G_g` (4D), `H_g` (5D, spin-2-adjacent). The `_u` (odd-parity) siblings
  exist too (`A_u`, `T_1u`, `T_2u`, `G_u`, `H_u`) and matter for quark/
  fermionic assignments.
- **Maxwell criticality / isostatic rigidity**: a single Jobson cell's own
  spring network (12 vertices, 30 edges) has exactly 6 zero
  ("floppy"/soft) modes at q=0 (3 translations `T_1u` + 3 rotations
  `T_1g`) — the `3V-E=6` Maxwell count. In 3D bulk, all 6 must engage
  simultaneously → `n=3×6=18`, the exponent reused for both the gravity
  coupling AND (candidate, unconfirmed) Higgs-related channel-counting.
- **Zone structure of a nucleon** (proton/neutron), inside out:
  - **Zone 1** (r < lambda_p): sub-cell, no Jobson cells present, quarks
    behave as free particles (asymptotic freedom).
  - **Zone 2** (r ~ lambda_p): the Maxwell-jamming BOUNDARY — a thin
    transition, cells locked (3V-E=6 exhausted).
  - **Zone 3** (lambda_p < r < r_p): co-rotating cells, spin here, outer
    edge = r_p = 4·lambda_p. Coulomb field = this zone's rolling imprint
    propagating outward.
  - **Zone 4** (r > r_p): ordinary free (unjammed, mobile) Jobson-cell
    medium — same L_J-scale grain as everywhere else in "empty" space,
    just not locally jammed by this particular nucleon.
- **Chirality-as-charge (same-spin repel, opposite-spin attract)**:
  established at the level of a real PyBullet rigid-body simulation
  (same-chirality Zone-2 driver pairs never achieve real contact under
  confinement; opposite-chirality pairs can, at short range) — see
  `analysis/nuclear/mesh_grind/pybullet/`. The finer causal picture ("cells
  between them are dragged in opposing directions... grind" vs "mesh") is
  still an unverified kinematic narrative layered on top of that
  aggregate result — tracked as N2-1 in `docs/series1/judgment_calls.txt`.

---

## 4. SOLID VS. GENUINELY OPEN — QUICK INDEX

**Treat as SOLID (don't re-derive without a specific reason):**
- alpha (topological derivation + CODATA match), Rs (4-scale cross-
  validated), r_p (Zone 2 jamming boundary), G/alpha_grav (topological
  exponent n=18), m_n−m_p splitting, Koide's 2/3 as a pure dimension
  ratio (mass-side closure is NOT solid — see below), the chirality-
  as-charge AGGREGATE contact-gating result (PyBullet sim).

**Genuinely OPEN (don't assume closed, don't re-litigate from scratch —
these have each already had multiple honest, documented attempts fail):**
- Higgs mass alpha/pi coefficient's native origin: the FORM (1+alpha/pi,
  standard QED scalar self-energy) and WHICH mode it applies to (scalar
  vs vector, from topology) are both solid; WHY the magnitude is alpha/pi
  specifically (not alpha/(2pi) or any other value), derived from this
  medium's own granular/contact mechanics rather than imported from QED,
  is not — 14+ independent, well-motivated attempts have all failed.
- Satellite/gravitational time dilation's independent medium derivation:
  GM/(r·c²) is confirmed to be GR's own weak-field formula imported
  wholesale (algebraically exact match, traced to primary session
  sources), not derived from this medium's own pressure/packing
  mechanics — multiple independent derivation attempts have all missed
  by many orders of magnitude, proven independent of which modulus or
  orbital scale is used.
- Koide mass-side algebraic closure: the 2/3 ratio's geometric origin
  (a pure dimension-count) is proven; whether the ACTUAL Born-balance-
  derived lepton masses must satisfy Koide's specific quadratic relation
  is not — no algebraic proof exists after exhaustive search, and the
  needed experimental precision (tau mass ~1000x better) doesn't exist.
- N2-1 (nucleon mesh/grind fine causal mechanism, see Section 3 above):
  the aggregate same-repels/opposite-can-touch result is simulation-
  confirmed; the specific claimed cell-level kinematic picture ("cells
  dragged in opposing directions... grind" vs "mesh") is not.
- Many smaller open derivation gaps exist across nuclear, lepton/quark
  mass, and medium-physics topics — this file covers only the handful
  most likely to come up repeatedly, not an exhaustive list.

**Where the day-to-day working history lives** (why something is hedged,
what's been tried and failed, judgment calls needing author sign-off):
tracked internally in this repo (not part of the public sync) — ask if
you need pointers to the internal tracking files for a specific topic.

---

## 5. REPO MAP (where things live)

- `docs/series1/` — the 11 published (Zenodo) core papers (`doc_alpha.txt`,
  `doc_higgs.txt`, `doc_nucleus.txt`, `doc_magnetism.txt`, etc.) + their
  `revisions/` subfolder (dated changelog per doc).
- `docs/series2/`, `docs/series3/` — applications/implications papers and
  newer, less-settled series respectively.
- `analysis/` — every calculation, organized by physics domain (alpha,
  higgs, nuclear, gravity, quantum, molecular, cosmo, medium, corpuscle,
  demos, etc.). Each `doc_X.txt`'s own companion script(s) are cited in
  that doc's text. `archive/` subfolders hold superseded scripts (kept for
  history, not part of the current best-known approach). Scripts are
  meant to be self-contained (standalone-runnable, comments don't rely on
  external notes/tracking files existing) — a real convention, follow it
  in anything new you write here.
- This repo also maintains internal-only working notes and issue-tracking
  files (investigation histories, judgment-call logs, an open-items
  index) that are NOT part of the public sync — don't cite specific paths
  to them in anything meant to be self-contained (scripts, this file), since
  they won't travel with the code. Summarize the relevant fact inline
  instead.
- `/memories/repo/*.md` — this AI's own persistent, repo-scoped working
  notes (conventions, past mistakes, standing policies) — a DIFFERENT,
  complementary thing from this file (private to one AI tool's memory
  system, not git-tracked, not visible to other models/sessions).

---

## 6. STANDING CONVENTIONS (apply repo-wide)

- **ABSOLUTE RULE**: every calculation lives in a committed, runnable
  script — never an inline one-off calculation.
- Docs distinguish **CLOSED / DERIVED / ESSENTIALLY CLOSED** (settled)
  from **OPEN** (genuine gap) — check a claim's own bracket tag before
  assuming it's settled.
- A hedge in a doc (e.g. "plausible, not yet independently verified") is
  usually there because a real investigation already happened and found a
  genuine gap — read the cited judgment-call/notes entry before assuming
  the hedge is just caution that can be casually strengthened.
- Long-running scripts (simulations, sweeps) are launched in the
  background (async) and checked on later, not blocked on synchronously.

---
*This file is synced to the public repo (torsionverse_public) on
`scripts/sync_public.py` runs — keep anything written here appropriate for
a public audience (no internal-only speculation, no unpublished-draft
specifics beyond what's already public in docs/series1-2, and no session/
personal notes).*
