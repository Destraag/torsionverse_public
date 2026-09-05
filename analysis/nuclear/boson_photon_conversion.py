"""
boson_photon_conversion.py
=============================
Verifies (does not just restate) notes/nuclear_pressure.txt LEAD N.2/N.4:
"Under a strong magnetic field B, the torsion medium is partially rigidified
(shear modulus G increases). This changes the T_1g dispersion relation,
potentially MIXING photon with W/Z modes... Conversion probability ~
(B/B_crit)^2 where B_crit ~ m_W^2/e... Planetary B ~ 1T << B_crit... near
neutron stars (B ~ 1e9 T), conversion could be detectable."

That lead was flagged as UNSCRIPTED (no supporting .py ever existed) and was
an order-of-magnitude, back-of-envelope estimate only. This script computes
the actual numbers properly:
  1. Derive B_crit via the standard Schwinger-critical-field CONSTRUCTION
     (the same one QED uses for e+e- pair production), substituting m_W for
     m_e -- sanity-checked first against the well-known QED electron value
     (~4.41e9 T) before trusting the W-boson version.
  2. Use REALISTIC magnetar surface field strengths (1e10-1e11 T, standard
     astrophysical range) rather than the notes file's own stated "~1e9 T",
     which undershoots the standard literature range by 1-2 orders of
     magnitude.
  3. Compute (B/B_crit)^2 at both the notes file's and the realistic range,
     and assess -- on pure order-of-magnitude grounds ONLY, not a full
     emission/luminosity model, which is out of scope here -- whether this
     could plausibly be large enough to matter for anything, including
     MG-J7 (K vs rho universality across the fluid/jammed transition).

SCOPE NOTE (important): this lead, even fully verified, is about G's
response to field strength (shear rigidification) -- it says nothing about
K. A large result here would not by itself resolve whether K is universal
across the fluid/jammed transition; a small result (the expected outcome)
closes out this specific lead without bearing on MG-J7 either way.

Reference: notes/nuclear_pressure.txt LEAD N.2, N.4; judgment_calls.txt MG-J7
"""
import math

SEP = "=" * 72
results = []


def check(name, cond, detail=""):
    status = "PASS" if cond else "FAIL"
    results.append((name, status))
    print(f"  [{'PASS' if cond else '*** FAIL'}] {name}")
    if detail:
        print(f"         {detail}")


# ---- Physical constants (CODATA / PDG SI) ----
c = 2.99792458e8              # m/s, exact
e_charge = 1.602176634e-19    # C, exact (SI 2019)
hbar = 1.054571817e-34        # J*s
eV_to_J = 1.602176634e-19

m_e_MeV = 0.51099895000       # PDG, electron mass
m_W_GeV = 80.3692              # PDG 2022, W boson mass
m_e_kg = m_e_MeV * 1e6 * eV_to_J / c ** 2
m_W_kg = m_W_GeV * 1e9 * eV_to_J / c ** 2

print(SEP)
print("boson_photon_conversion.py -- verifying LEAD N.2/N.4 (magnetized medium, T_1g/W-Z mixing)")
print(SEP)

# ---- Step 1: Schwinger-style critical field, sanity-checked on the electron ----
# Standard QED Schwinger critical field: B_c = m_e^2 * c^2 / (e * hbar)
B_crit_electron = m_e_kg ** 2 * c ** 2 / (e_charge * hbar)
KNOWN_B_SCHWINGER = 4.41e9  # T, well-established QED literature value

print(f"  m_e = {m_e_MeV} MeV = {m_e_kg:.6e} kg")
print(f"  m_W = {m_W_GeV} GeV = {m_W_kg:.6e} kg")
print(f"  B_crit(electron) = m_e^2*c^2/(e*hbar) = {B_crit_electron:.4e} T")
print(f"  Known QED Schwinger field (literature) = {KNOWN_B_SCHWINGER:.2e} T")

check("BC1: Schwinger-field construction reproduces the known QED electron "
      "critical field to <1% (validates the formula/units before reusing "
      "for m_W)",
      abs(B_crit_electron - KNOWN_B_SCHWINGER) / KNOWN_B_SCHWINGER < 0.01,
      f"computed = {B_crit_electron:.4e} T, known = {KNOWN_B_SCHWINGER:.2e} T, "
      f"ratio = {B_crit_electron/KNOWN_B_SCHWINGER:.4f}")

# ---- Step 2: same construction with m_W (the lead's own proposed formula) ----
B_crit_W = m_W_kg ** 2 * c ** 2 / (e_charge * hbar)
NOTES_FILE_B_CRIT = 5e15   # T, the notes file's own stated estimate

print()
print(f"  B_crit(W) = m_W^2*c^2/(e*hbar) = {B_crit_W:.4e} T")
print(f"  notes/nuclear_pressure.txt LEAD N.2's own stated estimate = {NOTES_FILE_B_CRIT:.1e} T")
print(f"  ratio (computed / notes-file estimate) = {B_crit_W/NOTES_FILE_B_CRIT:.4e}")

check("BC2: notes/nuclear_pressure.txt's own '~5e15 T' B_crit estimate is "
      "OFF by many orders of magnitude from the properly computed value "
      "(the lead's own arithmetic was never verified before this script)",
      abs(math.log10(B_crit_W / NOTES_FILE_B_CRIT)) > 2,
      f"log10(ratio) = {math.log10(B_crit_W/NOTES_FILE_B_CRIT):+.2f} "
      f"({B_crit_W/NOTES_FILE_B_CRIT:.2e}x)")

# ---- Step 3: realistic magnetar field strengths ----
# Standard astrophysical range for magnetar surface fields: 1e14-1e15 Gauss
# = 1e10-1e11 Tesla (1 T = 1e4 G). The notes file's "~1e9 T" undershoots
# this standard range by 1-2 orders of magnitude.
B_magnetar_notes_file = 1e9    # T, notes file's own (low) estimate
B_magnetar_low = 1e10           # T, standard literature low end
B_magnetar_high = 1e11          # T, standard literature high end (most extreme known)

print()
print(SEP)
print("STEP 3: conversion probability ~ (B/B_crit)^2 at various field strengths")
print(SEP)

for label, B in [
    ("notes file's own stated B~1e9 T", B_magnetar_notes_file),
    ("standard magnetar range, low end (1e10 T)", B_magnetar_low),
    ("standard magnetar range, high end (1e11 T)", B_magnetar_high),
]:
    prob = (B / B_crit_W) ** 2
    print(f"  {label}:")
    print(f"    B/B_crit(W)        = {B/B_crit_W:.4e}")
    print(f"    conversion prob ~  = {prob:.4e}")

prob_notes = (B_magnetar_notes_file / B_crit_W) ** 2
prob_high = (B_magnetar_high / B_crit_W) ** 2

check("BC3: even at the highest known magnetar field (1e11 T), conversion "
      "probability (B/B_crit)^2 using the PROPERLY computed B_crit is "
      "utterly negligible (< 1e-15)",
      prob_high < 1e-15,
      f"prob (1e11 T) = {prob_high:.4e}")

check("BC4: using the correctly-computed B_crit (not the notes file's own "
      "erroneous ~5e15 T estimate) makes the effect MUCH smaller than the "
      "notes file itself believed, not larger",
      prob_high < (B_magnetar_high / NOTES_FILE_B_CRIT) ** 2,
      f"prob with correct B_crit = {prob_high:.4e} vs notes-file's own "
      f"(already-tiny) estimate using their B_crit = "
      f"{(B_magnetar_high/NOTES_FILE_B_CRIT)**2:.4e}")

print()
print(SEP)
print("SUMMARY")
print(SEP)
n_pass = sum(1 for _, s in results if s == "PASS")
for name, status in results:
    print(f"  [{status}] {name}")
print(f"\n  Total: {len(results)}  PASS: {n_pass}  FAIL: {len(results) - n_pass}")

print()
print("  CONCLUSION: the Schwinger-field-style B_crit ~ m_W^2*c^2/(e*hbar)")
print("  construction, properly computed and unit-checked (BC1), gives")
print(f"  B_crit(W) = {B_crit_W:.3e} T -- roughly 5 orders of magnitude LARGER")
print("  than notes/nuclear_pressure.txt's own stated ~5e15 T (BC2), meaning")
print("  that file's own order-of-magnitude estimate contained an arithmetic")
print("  error never caught before this script existed. Using the CORRECT")
print("  B_crit, even the most extreme known magnetar field (1e11 T) gives a")
print("  photon<->W/Z conversion probability of order 1e-19 or smaller (BC3) --")
print("  many orders of magnitude below any plausible detection threshold,")
print("  and smaller than the notes file's own (already tiny) estimate, not")
print("  larger (BC4). This lead is CLOSED as a candidate for resolving")
print("  MG-J7 or predicting any observable effect: the mechanism may or may")
print("  not be real, but even if it is, it is far too small to matter here,")
print("  and it only ever bore on G's field-dependence, not on K, so it")
print("  could not have distinguished MG-J7's K-vs-rho fork regardless.")
print(SEP)
