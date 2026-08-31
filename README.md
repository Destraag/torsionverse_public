# Torsionverse

A geometric derivation of particle physics and field phenomena from a mechanical
medium with icosahedral (I_h) unit cell geometry (the Jobson cell). All results
are derived from two constants — the fine structure constant α and the proton
charge radius r_p — with zero additional free parameters.

**Author:** Robert Jobson  
**Status:** Working theoretical framework. All numerical predictions verified by
companion scripts. Peer review pending.

---

## Quick Start

```bash
# Any companion script runs standalone — no dependencies beyond Python stdlib
python analysis/demos/electron_doc.py

# Run everything at once (synthesis across all published papers)
python analysis/demos/synthesis_demo.py
```

---

## Published Papers

### Series 1 — Foundation

| # | Paper | Script | DOI |
|---|-------|--------|-----|
| R1 | A Geometric Derivation of the Fine Structure Constant | [alpha_doc.py](analysis/demos/alpha_doc.py) | [10.5281/zenodo.22013651](https://doi.org/10.5281/zenodo.22013651) |
| R2 | Torsion Medium Properties and the MOND Scale | [torsion_doc.py](analysis/demos/torsion_doc.py) | [10.5281/zenodo.22016573](https://doi.org/10.5281/zenodo.22016573) |
| R3 | The Jobson-Higgs Connection | [higgs_doc.py](analysis/demos/higgs_doc.py) | [10.5281/zenodo.22032555](https://doi.org/10.5281/zenodo.22032555) |
| R4 | The Jobson Cell: Icosahedral Unit Cell | [jobson_cell_doc.py](analysis/demos/jobson_cell_doc.py) | [10.5281/zenodo.22032906](https://doi.org/10.5281/zenodo.22032906) |
| R5 | Magnetism, Electricity and Torsion | [magnetism_doc.py](analysis/demos/magnetism_doc.py) | [10.5281/zenodo.22036406](https://doi.org/10.5281/zenodo.22036406) |
| R6 | The Proton as Spinning Icosahedral Engine | [nucleus_doc.py](analysis/demos/nucleus_doc.py) | [10.5281/zenodo.22042337](https://doi.org/10.5281/zenodo.22042337) |
| R7 | Orbit as Torsion Medium Pressure | [orbit_doc.py](analysis/demos/orbit_doc.py) | [10.5281/zenodo.22044378](https://doi.org/10.5281/zenodo.22044378) |
| R8 | Quantum Entanglement as Shared A_g Mode | [entanglement_doc.py](analysis/demos/entanglement_doc.py) | [10.5281/zenodo.22052350](https://doi.org/10.5281/zenodo.22052350) |
| R9 | Quantum Mechanics from the Jobson Cell Medium | [qm_doc.py](analysis/demos/qm_doc.py) | [10.5281/zenodo.22052395](https://doi.org/10.5281/zenodo.22052395) |
| R10 | A Mechanical Universe Derived from Light and Geometry | [torsionverse_doc.py](analysis/demos/torsionverse_doc.py) | [10.5281/zenodo.22052870](https://doi.org/10.5281/zenodo.22052870) |
| R11 | Lepton Mass Spectrum from Icosahedral Group Theory | [leptons_doc.py](analysis/demos/leptons_doc.py) | [10.5281/zenodo.22057421](https://doi.org/10.5281/zenodo.22057421) |
| R12 | Particle Generation Thresholds and Winding Formation | [particle_generation_doc.py](analysis/demos/particle_generation_doc.py) | [10.5281/zenodo.22068557](https://doi.org/10.5281/zenodo.22068557) |
| R13 | The Orbital Electron in the Torsion Medium | [electron_doc.py](analysis/demos/electron_doc.py) | [10.5281/zenodo.22105622](https://doi.org/10.5281/zenodo.22105622) |

### Series 2 — Applications and Predictions

| # | Paper | Script | DOI |
|---|-------|--------|-----|
| R14 | Series 2 Synthesis | [series2_doc.py](analysis/demos/series2_doc.py) | [10.5281/zenodo.22108664](https://doi.org/10.5281/zenodo.22108664) |
| R15 | Acoustic Implications of the Torsion Medium | [acoustics_doc.py](analysis/demos/acoustics_doc.py) | [10.5281/zenodo.22139731](https://doi.org/10.5281/zenodo.22139731) |
| R16 | Chemical Bonds and Molecular Manipulation | [chemistry_doc.py](analysis/demos/chemistry_doc.py) | [10.5281/zenodo.22139803](https://doi.org/10.5281/zenodo.22139803) |
| R17 | Gyroscope Precession and Fluid Vortex Sensing | [gyroscope_doc.py](analysis/demos/gyroscope_doc.py) | [10.5281/zenodo.22139868](https://doi.org/10.5281/zenodo.22139868) |
| R18 | Directed Hadron Manipulation | [hadron_manipulation.py](analysis/demos/hadron_manipulation.py) | [10.5281/zenodo.22140068](https://doi.org/10.5281/zenodo.22140068) |
| R19 | Cosmological Redshift — Candidate Mechanisms | [redshift_doc.py](analysis/demos/redshift_doc.py) | [10.5281/zenodo.22140127](https://doi.org/10.5281/zenodo.22140127) |
| R20 | Counter-Rotating Fluid Bells | [fluid_bells_doc.py](analysis/demos/fluid_bells_doc.py) | [10.5281/zenodo.22140336](https://doi.org/10.5281/zenodo.22140336) |
| R21 | Muon Lubrication in Muon-Catalysed Fusion | [muon_lubrication_doc.py](analysis/demos/muon_lubrication_doc.py) | [10.5281/zenodo.22140481](https://doi.org/10.5281/zenodo.22140481) |
| R22 | Piezoelectricity and the Electron Escalator | [acoustics_doc.py](analysis/demos/acoustics_doc.py) | [10.5281/zenodo.22140517](https://doi.org/10.5281/zenodo.22140517) |
| R23 | Spinning Pyramid: T_1g Chi Coupling | [acoustics_doc.py](analysis/demos/acoustics_doc.py) | [10.5281/zenodo.22140877](https://doi.org/10.5281/zenodo.22140877) |

### Series 3 — Medium Structure

| Paper | Script | DOI |
|-------|--------|-----|
| Jobson Cell Lattice Coherence and Rigidity | [cell_coherence_doc.py](analysis/demos/cell_coherence_doc.py) | [10.5281/zenodo.22151160](https://doi.org/10.5281/zenodo.22151160) |

---

## Repository Structure

```
docs/series1/          published papers — Series 1 foundation (.txt, .pdf)
docs/series2/          published papers — Series 2 applications (.txt, .pdf)
docs/series3/          Series 3: doc_cell_coherence (published); others in progress
analysis/demos/        companion scripts, one per paper (all standalone)
analysis/quantum/      proof scripts cited in papers (group theory, CG tables)
analysis/nuclear/      proof scripts cited in papers (nuclear structure)
analysis/gravity/      proof scripts cited in papers (orbital mechanics, MOND)
analysis/higgs/        proof scripts cited in papers (Higgs/EW sector)
analysis/constants.py  shared physical constants
```

All companion scripts are **standalone** — no external dependencies, no imports
from other project files. Any script runs directly from the repo root:

```bash
python analysis/demos/<script>.py
```
