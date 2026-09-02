# HIP²LInterActomics

**Graphical and headless workflows for reproducible analysis, comparison, and interpretation of intermolecular interactions using LUNA.**

HIP²LInterActomics is an open-source application for preparing molecular interaction projects, running [LUNA](https://github.com/keiserlab/LUNA), and analyzing the resulting protein–ligand interaction patterns and interaction fingerprints.

Rather than replacing LUNA's molecular perception and interaction engine, HIP²LInterActomics provides a higher-level workflow around it. The software integrates molecular input preparation, project management, interaction-analysis configuration, interaction fingerprints, similarity analysis, clustering, residue-level statistics, PyMOL visualization, and report generation in a reproducible desktop or headless workflow.

The same project structure can be used from the graphical interface or from the command line, allowing analyses created on a workstation to be reproduced or extended on Linux servers and HPC systems.

---

## Main capabilities

HIP²LInterActomics currently provides five complementary layers of functionality:

1. **Molecular input preparation**
2. **LUNA interaction analysis**
3. **Interaction fingerprint generation and comparison**
4. **Post-processing, visualization, and interpretation**
5. **Graphical, headless, and parameter-sweep execution**

---

## 1. Molecular input preparation

HIP²LInterActomics accepts different organization schemes commonly produced by docking and molecular-dynamics workflows.

### Protein inputs

The receptor can be supplied as:

- a single PDB file used for all ligands;
- a directory containing one PDB file per ligand, pose, complex, or trajectory frame.

When a protein directory is used, receptor and ligand files are paired by their base names.

### Ligand inputs

Supported ligand formats include:

- MOL2
- SDF
- MOL
- PDB
- ENT

Ligands may be supplied as:

- a single multi-molecule file;
- individual molecular files in a directory;
- a previously prepared ligand collection.

Homogeneous folders containing MOL2 files or SDF/MOL files can be consolidated into a single ligand input file. MOL2 consolidation can also remove `LP` pseudo-atoms and renumber atoms and bonds before analysis.

### Complex preprocessing

The built-in complex preparation workflow can process directories containing protein–ligand complexes and generate separated protein and ligand inputs for subsequent LUNA analysis.

This is useful for:

- docking poses stored as protein–ligand complexes;
- collections of independently prepared complexes;
- molecular-dynamics snapshots;
- interaction analyses in which each entry has its own receptor coordinates.

Water molecules can be retained in the protein structures when hydrated interaction analysis is required.

### Ligand selection

After loading an input library, detected ligands can be filtered and individually selected before execution. This allows the same source file to be reused for different subsets without modifying the original molecular file.

---

## 2. LUNA interaction analysis

The interaction-detection stage is performed by LUNA.

HIP²LInterActomics exposes both standard and advanced options while keeping the complete project configuration stored with the analysis.

### Hydrogen and protonation settings

The user can control:

- addition of hydrogens before analysis;
- pH used for hydrogen addition;
- inclusion or exclusion of water molecules.

The default configured pH is **7.4**.

### Interaction configuration

Advanced projects can control LUNA's `InteractionCalculator`, including options related to:

- proximal contacts;
- generic atom–atom contacts;
- dependent interactions;
- water-dependent contacts;
- self-interactions.

A complete LUNA interaction configuration file (`.cfg`) may also be supplied.

HIP²LInterActomics additionally supports a **project-level maximum interaction-distance cap**, allowing upper distance limits to be constrained without replacing shorter limits already defined by LUNA.

### Binding-mode filters

Binding-mode rules can be loaded from `.cfg` files and used to restrict the interactions or poses included in an analysis.

A graphical binding-mode editor is included to facilitate creation and modification of these rules.

### Interaction-aware PyMOL export

PyMOL sessions (`.pse`) can be generated for structural inspection of the analyzed complexes.

Session generation can optionally be restricted to selected interaction classes.

---

## 3. Interaction fingerprints

HIP²LInterActomics supports the three hashed interaction fingerprints implemented by LUNA:

- **EIFP** — Extended Interaction FingerPrint
- **FIFP** — Functional Interaction FingerPrint
- **HIFP** — Hybrid Interaction FingerPrint

The three fingerprints can be calculated individually or in the same project.

Fingerprint parameters exposed by HIP²LInterActomics include:

- number of levels;
- radius step;
- fingerprint length;
- binary or count representation.

When all three fingerprint types are requested, separate outputs are generated for HIFP, EIFP, and FIFP.

### Tanimoto similarity

Interaction fingerprints can be compared using Tanimoto similarity.

HIP²LInterActomics can generate:

- pairwise similarity tables;
- square similarity matrices;
- similarity heatmaps.

Similarity matrices can either be produced during the LUNA workflow or calculated later from an existing fingerprint file.

This makes fingerprint generation and downstream similarity analysis independent steps.

---

## 4. Results and interaction analysis

Completed projects can be opened directly in the **Results** workspace without rerunning LUNA.

The results interface is loaded on demand and provides several complementary analyses.

### Interaction fingerprints

Fingerprint CSV files can be inspected directly inside the application.

When several fingerprint types were generated, HIFP, EIFP, and FIFP can be selected independently.

### Interaction statistics

HIP²LInterActomics summarizes the intermolecular contacts detected across the analyzed entries.

Statistics can be inspected:

- for the complete project;
- for individual ligands or entries;
- by interaction type.

Interaction families can also be temporarily hidden from plots without changing the underlying project data.

### Residue-level analysis

The software builds residue-by-entry interaction matrices that connect the detected interaction types to protein residues.

These data can be visualized as interaction-specific heatmaps, allowing identification of residues that participate recurrently in ligand recognition.

### Combined interaction heatmaps

Interaction types, protein residues, and molecular entries can be combined in larger heatmap representations for comparative analysis of binding patterns across a ligand set.

### Docking-pose and trajectory analysis

A project can be explicitly designated as a:

- docking-pose analysis; or
- molecular-dynamics trajectory/frame analysis.

When this mode is enabled, entries are interpreted as poses or frames rather than unrelated ligands.

The results can then be summarized using entry/frame occurrence and percentage-based interaction statistics, providing a way to evaluate the persistence of contacts across a structural ensemble.

HIP²LInterActomics does **not** perform the molecular-dynamics simulation itself. It analyzes previously generated structural frames or complexes.

### Hierarchical clustering

Ligands or entries can be clustered from their interaction-fingerprint similarity matrix.

The current implementation supports hierarchical linkage methods including:

- average;
- complete;
- single.

The software generates a dendrogram and cluster assignments, which can be exported as CSV.

---

## 5. Fingerprint interpretation and prioritization

One goal of HIP²LInterActomics is to preserve the interpretability provided by LUNA interaction fingerprints instead of treating them only as anonymous numerical vectors.

The results workflow can inspect fingerprint features together with their:

- interaction provenance;
- protein or ligand origin;
- shell/level information;
- occurrence across molecular entries;
- associated residues and interaction types.

Features can be classified according to their structural context, including features associated with protein–ligand noncovalent interactions, ligand- or protein-derived environments, intramolecular interactions, collisions, and features whose assignment is considered unreliable.

### Feature-importance analysis

Optional experimental labels can be supplied through CSV or TSV files.

Labels may represent either:

- a continuous endpoint (**regression**); or
- discrete classes (**classification**).

HIP²LInterActomics can then estimate the relative importance of interaction-fingerprint features and relate prioritized features back to their structural and interaction context.

The implemented workflow uses reproducible random seeds and tree-based scikit-learn models when available, with internal fallback strategies if the machine-learning dependency cannot be used.

These analyses are intended as **exploratory and interpretive tools**. Model-derived feature importance should not be interpreted as causal evidence of molecular activity.

---

## 6. PyMOL integration

Generated `.pse` files can be opened directly from HIP²LInterActomics.

The Results workspace also supports **dynamic binding-mode filtering** of existing projects.

A new binding-mode configuration can be applied after the original analysis to generate a separate set of filtered PyMOL sessions. The original sessions are preserved.

This makes it possible to create multiple structural views of the same project without repeating the complete interaction calculation.

---

## 7. Figures and reports

Plots generated in the Results workspace can be exported for further use.

Currently supported direct figure formats are:

- PNG
- SVG
- PDF

HIP²LInterActomics can also generate project-level:

- **HTML reports**
- **PDF reports**

Reports combine project configuration, interaction summaries, residue statistics, cluster information, and available result figures.

The PDF workflow uses an isolated rendering process so that report-generation failures do not terminate the graphical application.

Scientific interpretation paragraphs included in automatically generated reports are intended to guide result inspection and should not replace expert evaluation of the molecular system.

---

## Graphical workflow

The desktop application is organized into six main sections:

### 1. Setup

Detects or prepares the scientific LUNA environment.

### 2. Project

Defines:

- working directory;
- receptor input;
- ligand input;
- complex preprocessing;
- ligand selection;
- water handling;
- new-project or project-fork mode;
- docking-pose/trajectory mode.

### 3. Analyses

Configures:

- hydrogen handling and pH;
- interaction settings;
- interaction fingerprints;
- Tanimoto similarity;
- PyMOL sessions;
- binding-mode filters;
- optional supervised fingerprint labels.

### 4. Run

Validates the project and executes the scientific workflow.

The interface displays:

- the effective command;
- execution progress;
- a live LUNA log;
- cancellation controls.

### 5. Results

Loads completed projects and provides fingerprint, statistical, residue-level, clustering, PyMOL, and interpretability analyses.

### 6. History

Reopens recently used projects.

---

## Execution modes

HIP²LInterActomics provides three entry points.

| Mode | Command | Purpose |
|---|---|---|
| Graphical desktop | `hip2linteractomics` | Interactive project creation, execution, visualization, and reporting |
| Headless single project | `hipplinteractomics-terminal` | Run or post-process one project from JSON and/or command-line arguments |
| Headless multiple run | `hipplinteractomics-multiple-run` | Generate and execute systematic combinations of fingerprint parameters |

---

## Architecture

HIP²LInterActomics separates interface/orchestration from the chemistry stack.

In the packaged desktop workflow, two environments are used:

### Interface environment

Responsible for:

- PyQt6 user interface;
- project management;
- numerical post-processing;
- plotting;
- result visualization;
- report generation.

Typical components include:

- PyQt6
- Matplotlib
- NumPy
- SciPy
- scikit-learn

### Computation environment (`luna-env`)

Responsible for molecular perception and chemistry-related calculations.

It contains:

- LUNA
- RDKit
- Open Babel
- PyMOL
- Biopython
- supporting scientific dependencies

The interface does not directly import the chemistry stack during normal execution. Scientific calculations are started as subprocesses.

The environments exchange project state and results through the project's working directory.

```text
Interface environment
        |
        | project configuration / entries / results
        v
+-------------------------------+
|       Working directory       |
|                               |
|  project JSON                 |
|  entries                      |
|  results                      |
|  logs                         |
|  generated analysis artifacts |
+-------------------------------+
        ^
        | subprocess execution
        |
Computation environment
        |
        +-- LUNA
        +-- RDKit
        +-- Open Babel
        +-- PyMOL
