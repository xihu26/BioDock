# BioDock
BioDock AI is an automated molecular docking pipeline that retrieves suitable protein structures and prepares both protein receptors and small-molecule ligands for docking.

## Features
* Searches **UniProt** for protein matches based on protein name, gene, and organism.
* Retrieves corresponding structures from the **RCSB Protein Data Bank (PDB)**.
* Compares available structures and selects a suitable structure based on experimental method and resolution.
* Prepares receptor structures for docking.
* Converts ligand structures from SDF to PDBQT using **Meeko**.
* Designed to handle multiple protein matches and structures.

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

* Python 3
* RDKit
* Requests
* Meeko
* scipy
* gemmi
* pdbfixer

Install the Python dependencies with:

```bash
pip install rdkit requests meeko scipy gemmi pdbfixer
```

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

## Project Structure

```text
BioDock AI/
├── BioDock
├── Preparation
├── proteins/
├── ligands/
└── README.md
```

## Status
BioDock is currently under development, with ongoing work on robust receptor preparation, non-standard residue handling, and docking integration.
