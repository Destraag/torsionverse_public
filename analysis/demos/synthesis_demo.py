"""
synthesis_demo.py
=================
Master verification for docs/doc_torsionverse.txt — all results in one run.
Runs every component script and reports a combined pass/fail summary.
Excludes: analysis/gravity/ih_lattice_phonon.py (separate phonon model).

Total expected: 172 checks across 10 scripts.
  47  alpha_doc        [fine structure constant]
  15  torsion_doc      [medium properties, MOND, E=mc^2]
  12  higgs_doc        [Higgs sector, Weinberg angle]
  28  jobson_cell_doc  [Jobson cell geometry]
   8  magnetism_doc    [EM field physics]
  26  nucleus_doc      [nuclear structure, g_p, g_n, magic numbers]
  16  orbit_doc        [G derivation, orbital mechanics, MOND]
   9  torsionverse_doc [unified coupling, GPS, heat, pion, neutron mass]
   4  molien_n18       [I_h invariant ring, n=18 framework]
   7  neutron_g_factor [neutron g_n free/bound, in-medium form factor]

Usage: python analysis/demos/synthesis_demo.py

Reference: docs/doc_torsionverse.txt
"""

import sys, os, subprocess, re, time

SEP  = "=" * 70
SEP2 = "-" * 70

SCRIPTS = [
    # (display_name, path, expected_pass)
    ("alpha_doc",          "analysis/demos/alpha_doc.py",          47),
    ("torsion_doc",        "analysis/demos/torsion_doc.py",         15),
    ("higgs_doc",          "analysis/demos/higgs_doc.py",           12),
    ("jobson_cell_doc",    "analysis/demos/jobson_cell_doc.py",      28),
    ("magnetism_doc",      "analysis/demos/magnetism_doc.py",        8),
    ("nucleus_doc",        "analysis/demos/nucleus_doc.py",         26),
    ("orbit_doc",          "analysis/demos/orbit_doc.py",           16),
    ("torsionverse_doc",   "analysis/demos/torsionverse_doc.py",     9),
    ("molien_n18",         "analysis/demos/molien_n18.py",           4),
    ("neutron_g_factor",   "analysis/nuclear/neutron_g_factor.py",   7),
]

EXPECTED_TOTAL = sum(e for _, _, e in SCRIPTS)

print(SEP)
print("synthesis_demo.py — Master verification for doc_torsionverse.txt")
print(f"Running {len(SCRIPTS)} scripts, {EXPECTED_TOTAL} checks expected")
print(SEP)

results = []
total_pass = 0
total_fail = 0
t_start = time.time()

for name, path, expected in SCRIPTS:
    t0 = time.time()
    proc = subprocess.run(
        [sys.executable, path],
        capture_output=True, text=True,
        cwd=os.path.join(os.path.dirname(__file__), '..', '..')
    )
    elapsed = time.time() - t0
    out = proc.stdout

    # Parse pass/fail from output
    m_pass = re.search(r'PASS:\s*(\d+)', out)
    m_fail = re.search(r'FAIL:\s*(\d+)', out)
    m_total = re.search(r'Total[: ]+(\d+)', out)
    m_nn    = re.search(r'(\d+)/\d+ PASS', out)

    if m_pass and m_fail:
        passed = int(m_pass.group(1))
        failed = int(m_fail.group(1))
    elif re.search(r'ALL CHECKS PASSED', out):
        passed = expected
        failed = 0
    else:
        passed = 0
        failed = expected

    ok = (failed == 0 and passed >= expected)
    results.append((name, expected, passed, failed, ok, elapsed))
    total_pass += passed
    total_fail += failed

    marker = "[PASS]" if ok else "[FAIL] ***"
    print(f"  {marker} {name:<20}  {passed:>3}/{expected:<3}  ({elapsed:.1f}s)")

t_total = time.time() - t_start
print()
print(SEP2)
print("RESULTS BY PAPER")
print(SEP2)
section_map = {
    "alpha_doc":        "(1) Fine structure constant",
    "torsion_doc":      "(3) Medium properties",
    "higgs_doc":        "(2) Higgs boson",
    "jobson_cell_doc":  "(4) Jobson cell",
    "magnetism_doc":    "(5) Electromagnetism",
    "nucleus_doc":      "(6) Nuclear structure",
    "orbit_doc":        "(7) Gravity and orbits",
    "torsionverse_doc": "(8-13) New results [synthesis]",
    "molien_n18":       "[n=18 algebraic framework]",
    "neutron_g_factor": "(13) Neutron magnetic moment",
}
for name, exp, p, f, ok, _ in results:
    section = section_map.get(name, name)
    status  = "PASS" if ok else "FAIL ***"
    print(f"  {status:<8} {section:<40} {p}/{exp}")

print()
print(SEP)
passed_scripts = sum(1 for *_, ok, __ in results if ok)
failed_scripts = len(results) - passed_scripts
print(f"  Scripts: {passed_scripts}/{len(results)} pass")
print(f"  Checks:  {total_pass}/{total_pass + total_fail} pass")
print(f"  Time:    {t_total:.1f}s")
print()
if total_fail == 0 and total_pass == EXPECTED_TOTAL:
    print(f"  ALL {total_pass} CHECKS PASSED.")
else:
    print(f"  *** {total_fail} CHECKS FAILED or {EXPECTED_TOTAL - total_pass} missing ***")
    for name, exp, p, f, ok, _ in results:
        if not ok:
            print(f"    FAILED: {name}  (pass={p}, fail={f}, expected={exp})")
print(SEP)
print("  Reference: docs/doc_torsionverse.txt")
print("  Phonon model: python analysis/gravity/ih_lattice_phonon.py  [12/12]")
print(SEP)
