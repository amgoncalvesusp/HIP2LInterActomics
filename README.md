# HIP2LInterActomics_GUI

HIP2L InterActomics GUI is a desktop interface for building and exploring protein-ligand and protein-protein intermolecular interaction datasets with LUNA. The application turns docking poses, hydrated complexes, or molecular dynamics trajectory frames into interaction summaries, fingerprint matrices, interpretable heatmaps, PyMOL sessions, and PDF reports.

The software is designed as a practical hub for virtual screening analysis, trajectory inspection, pharmacophoric feature extraction, and data generation for future machine learning models. It helps compare interaction patterns across ligands or frames, filter binding modes, calculate LUNA interaction fingerprints, estimate feature importance with scikit-learn models, and export publication-ready visual summaries.

## Main Capabilities

- Prepare docking complexes by separating receptor, ligand, and crystallographic or modeled waters while removing water lone-pair pseudo-atoms.
- Run LUNA analyses for standard and hydrated protein-ligand datasets, including per-ligand receptor files.
- Support molecular dynamics trajectory-style projects where entries are frames or poses ordered by their numeric identifier.
- Generate interaction frequency plots, prevalent-interaction summaries, ligand-residue heatmaps, interaction-type heatmaps, similarity matrices, and fingerprint importance dashboards.
- Export filtered PyMOL sessions and PDF reports with parameters, interpretation notes, and result figures.
- Use external labels for classification or regression workflows, including ordered heatmaps and model-based feature importance.
- Work from Windows directly and from Linux/WSL when the LUNA environment and graphical dependencies are available.

## Typical Use Cases

- Prioritizing virtual screening hits using interpretable protein-ligand interaction fingerprints.
- Comparing hydrated and non-hydrated docking poses to evaluate water-mediated contacts.
- Summarizing molecular dynamics frames as interaction prevalence across a trajectory.
- Extracting pharmacophoric interaction features for a specific target.
- Creating interaction-encoded datasets for later machine learning model development.

## Repository Notes

The repository contains the GUI source code, runtime helpers, tests, and distribution launchers. Generated analysis folders such as `test_hydrated_linux/`, `test_hydrated_windows/`, and trajectory output directories are ignored because they can be large and are reproducible from the input projects.

The public software and repository name is `HIP2LInterActomics_GUI`. The internal Python package remains `luna_gui` for compatibility with existing imports and saved projects.
