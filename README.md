# BioDock
BioDock AI is an automated molecular docking pipeline that retrieves suitable protein structures and prepares both protein receptors and small-molecule ligands for docking.

## Features
* Searches **UniProt** for protein matches based on protein name, gene, and organism.
* Retrieves corresponding structures from the **RCSB Protein Data Bank (PDB)**.
* Compares available structures and selects a suitable structure based on experimental method and resolution.
* Prepares receptor structures for docking.
* Converts ligand structures from SDF to PDBQT using **Meeko**.
* Designed to handle multiple protein matches and structures.
* Analyses predicted molecular interactions using **PLIP**
* Generates a concise biochemical interpretation of the predicted binding interactions.


## Pipeline
```text
Protein input
     ↓
UniProt search
     ↓
Protein profiles
     ↓
RCSB PDB search
     ↓
Candidate structures
     ↓
Best structure selection
     ↓
Receptor preparation
     ↓
Docking-ready receptor + ligand
```

## Requirements

* RDKit
* Requests
* Meeko
* PDBFixer
* OpenMM
* PLIP
* RCSB API

Install the Python dependencies with:

```bash
pip install requests meeko pdbfixer openmm plip rcsb-api
```
## External Tools
BioDock also requires the following external molecular modelling tools:

P2Rank — binding pocket prediction
AutoDock Vina — molecular docking
Open Babel — molecular structure conversion

## Usage

Run the main program:

```bash
python BioDock
```

Enter:

```text
Protein: amylase
Gene:
Scientific Organism name: homo sapiens
```

BioDock will search for matching UniProt profiles, identify associated PDB structures, and select a suitable structure for preparation.

## Status
BioDock is currently under development, with ongoing work on robust receptor preparation, non-standard residue handling, and docking integration.

## Limitations
BioDock is currently under development.

Current limitations include:
* Docking predictions depend heavily on the quality of the available protein structure.
* The highest-ranked P2Rank pocket may not always correspond to the biologically relevant binding site.
* Protein preparation may fail for structures containing unusual residues or complex bonding patterns. Residues with missing atoms and unmatched templates are ignored.
* AutoDock Vina binding affinities are computational estimates and should not be interpreted as experimentally measured binding energies.
* Predicted interactions require experimental or higher-level computational validation.