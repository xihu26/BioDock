import csv
import subprocess
from plip.structure.preparation import PDBComplex

def main():
    protein = "3SH4"
    ligand = "donepezil"
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
    ligand_pose = f"results/{ligand}_{protein}_pose.pdb"
    combined_file = f"results/{ligand}_{protein}.pdb"

    #Converting the best docking pose into PDB format using Open Babel
    subprocess.run([
        "obabel",
        f"results/{ligand}_{protein}.pdbqt",
        "-O",
        ligand_pose, #1 = best pose
        "-f", "1",
        "-l", "1" 
    ], check=True)

    # Combine protein and ligand
    with open(f"proteins/{protein}_fixed.pdb", "r") as protein_file:
        protein_lines = protein_file.readlines()

    with open(ligand_pose, "r") as ligand_file:
        ligand_lines = ligand_file.readlines()

    with open(combined_file, "w") as output:
        #protein writing
        for line in protein_lines:
            if line.startswith(("ATOM", "HETATM", "TER")):
                output.write(line)
        output.write("TER\n") #End of molecular chain

        #Write docked ligand
        for line in ligand_lines:
            if line.startswith(("ATOM", "HETATM")):
                output.write(line)
        output.write("END\n")
    return combined_file

def analyze_interactions(docked_file):
    molecule = PDBComplex()
    molecule.load_pdb(docked_file)
    molecule.analyze()

    for ligand_id, interactions in molecule.interaction_sets.items():
        analysis = {
    "Ligand": ligand_id,
    "Hydrogen bonds": int(len(interactions.hbonds_ldon) + len(interactions.hbonds_pdon)),
    "Hydrophobic contacts": int(len(interactions.all_hydrophobic_contacts)),
    "π-stacking": int(len(interactions.pistacking)),
    "π-cation interactions": int(len(interactions.pication_laro) + len(interactions.pication_paro)),
    "Salt bridges": int(len(interactions.saltbridge_lneg) + len(interactions.saltbridge_pneg)),
    "Halogen bonds": int(len(interactions.halogen_bonds)),
    "Water bridges": int(len(interactions.water_bridges)),
    "Metal complexes": int(len(interactions.metal_complexes))
        }

#Verbal interpretation
    ligand_name = analysis["Ligand"]
    hbonds = analysis["Hydrogen bonds"]
    hydrophobic = analysis["Hydrophobic contacts"]
    pistacking = analysis["π-stacking"]
    pication = analysis["π-cation interactions"]
    salt_bridges = analysis["Salt bridges"] #For future refinement purposes

        # Overall biochemical interpretation
    result = (f"\nThe predicted ligand forms {hbonds} Hydrogen bonds and {hydrophobic} hydrophobic contacts with the protein. ")

    if hydrophobic > hbonds and hydrophobic > 0:
        result += (
            "Biochemically, this suggests that the ligand is primarily "
            "stabilised within a hydrophobic region of the binding pocket. "
        )

        if hbonds > 0:
            result += (
                "The hydrogen bonds may help orient the ligand and increase "
                "binding specificity by providing favourable polar interactions."
            )

    elif hbonds > hydrophobic and hbonds > 0:
        result += (
            "Biochemically, the predicted binding appears to rely mainly on "
            "specific polar interactions. Hydrogen bonds can help stabilise "
            "the ligand in a particular orientation within the binding pocket."
        )

    elif hbonds > 0 and hydrophobic > 0:
        result += (
            "The predicted binding appears to involve a combination of "
            "hydrophobic and polar interactions, which may collectively "
            "stabilise the ligand within the binding pocket."
        )

    else:
        result += (
            "Few major non-covalent interactions were detected, suggesting "
            "that the predicted binding mode may be weakly stabilised or "
            "depend on interaction types not captured in this summary."
        )

    return result


if __name__ == "__main__":
    main()