import subprocess
from rdkit import Chem
import meeko
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def main():
    convert_ligand("cinnamaldehyde") #hashtagged for now
    convert_receptor("1MFU") #protein == protein_id in this file

#Prepares ligand and proteins into pdbqt format for docking.
def convert_ligand(ligand):
    supplier_L = Chem.SDMolSupplier(f"ligands/{ligand}.sdf")
    if supplier_L[0] is None:
        raise ValueError("Could not read ligand.")
    
    molecule_L = Chem.AddHs(supplier_L[0])
    prep = meeko.MoleculePreparation()
    mol_L, valid, error_msg = meeko.PDBQTWriterLegacy.write_string(prep.prepare(molecule_L)[0])

    with open(f"ligands/{ligand}.pdbqt", "w") as file:
        file.write(mol_L)

    if valid == False:
        raise ValueError("Preparation error")

def convert_receptor(protein):
    fixer = PDBFixer(filename=f"proteins/{protein}.pdb")
    fixer.findMissingResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.addMissingHydrogens(7.0) #specified pH to determine which hydrogens are missing and adds them.

    with open(f"proteins/{protein}_fixed.pdb", "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)

    templates = meeko.ResidueChemTemplates.create_from_defaults()
    meeko_path = Path(meeko.__file__).parent
    templates.add_json_file(meeko_path / "data" / "NAKB_templates.json")

#Meeko's receptor-preparation API
    subprocess.run([
    "mk_prepare_receptor.py", "--read_pdb",
    f"proteins/{protein}_fixed.pdb", "-o", f"proteins/{protein}", "-p",
    "--allow_bad_res"], check=True) #ignore or remove residues with missing atoms or unmatched templates

if __name__ == "__main__":
    main()
