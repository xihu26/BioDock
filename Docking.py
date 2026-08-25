import csv
import subprocess
from plip.structure.preparation import PDBComplex

def main():
    protein = "8CB1"
    ligand = "cinnamaldehyde"
    print(find(protein))
    dock(protein, ligand)
    docked_file = convert_docked(ligand, protein)
    print(analyze_interactions(docked_file))

def find(protein):
#Ligand Finding on receptor
    pdb_file = f"proteins/{protein}.pdb"
    subprocess.run(["./p2rank-2.5.1/distro/prank",
        "predict", "-f", pdb_file], check=True)

    with open(f"/Users/lunahuang/BioDock AI/p2rank-2.5.1/distro/test_output/predict_{protein}/{protein}.pdb_predictions.csv", "r") as file:
        reader = csv.DictReader(file, skipinitialspace=True)

        for row in reader:
            if row["rank"].strip() == "1":
                x = float(row["center_x"])
                y = float(row["center_y"])
                z = float(row["center_z"])
                center = (x, y, z)

    with open("proteins/config.txt", "w") as config:
        config.write(f"center_x = {center[0]}\n")
        config.write(f"center_y = {center[1]}\n")
        config.write(f"center_z = {center[2]}\n")
        config.write("size_x = 20\n")
        config.write("size_y = 20\n")
        config.write("size_z = 20\n")

    return f"\nPredicted binding pocket center:\n{center}\n"


def dock(protein, ligand):

    subprocess.run([
    "vina",
    "--receptor", f"proteins/{protein}.pdbqt",
    "--ligand", f"ligands/{ligand}.pdbqt",
    "--config", "proteins/config.txt",
     "--out", f"results/{ligand}_{protein}.pdbqt"
], check=True)

#Analysis
def convert_docked(ligand, protein):
    subprocess.run([
        "obabel",
        f"results/{ligand}_{protein}.pdbqt",
        "-O",
        f"results/{ligand}_{protein}.pdb"
    ], check=True)
    return f"results/{ligand}_{protein}.pdb"

def analyze_interactions(docked_file):
    molecule = PDBComplex()
    molecule.load_pdb(docked_file)
    molecule.analyze()

    for ligand_id, interactions in molecule.interaction_sets.items():
        return (
    f"\nLigand: {ligand_id}\n"
    f"Hydrogen bonds: {len(interactions.hbonds_ldon) + len(interactions.hbonds_pdon)}\n"
    f"Hydrophobic contacts: {len(interactions.all_hydrophobic_contacts)}\n"
    f"π-stacking: {len(interactions.pistacking)}\n"
    f"π-cation interactions: {len(interactions.pication_laro) + len(interactions.pication_paro)}\n"
    f"Salt bridges: {len(interactions.saltbridge_lneg) + len(interactions.saltbridge_pneg)}\n"
    f"Halogen bonds: {len(interactions.halogen_bonds)}\n"
    f"Water bridges: {len(interactions.water_bridges)}\n"
    f"Metal complexes: {len(interactions.metal_complexes)}"
)


if __name__ == "__main__":
    main()