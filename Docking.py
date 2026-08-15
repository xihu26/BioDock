import csv
import subprocess

def main():
    protein = "1C8Q"
    ligand = "cinnamaldehyde"
    print(find(protein))
    dock(protein, ligand)

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
    "--config", "proteins/config.txt"
], check=True)

if __name__ == "__main__":
    main()