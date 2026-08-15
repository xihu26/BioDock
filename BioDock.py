from Retrieval import search_uniprot, search_pdb, best_pdb, download_best, search_ligand
from Preparation import convert_ligand, convert_receptor
from Docking import find, dock

def main():
#Retrieval
  #Protein input
    protein = input("Protein: ")
    gene = input("Gene: ")
    organism = input("Scientific Organism name: ")

#Protein profiling
    profiles = search_uniprot(protein, organism, gene)
    all_profiles = [] #Stores the all protein profiles found with search
    best_profiles = []#Stores the PDB IDs of the best proteins from each profile

    for profile in profiles:
        pdb_ids = search_pdb(profile)
        print("\nACCESSION:", profile["accession"])
        print("\nPDB IDS:", pdb_ids)
        if pdb_ids:
            all_profiles.append({
            "profile": profile,
            "pdb_ids": pdb_ids
        })
    #Finding the best out of the best protein
            structures = best_pdb(pdb_ids)
            best_profiles.append(structures[0].pdb_id)

    if not best_profiles:
        print("No valid PDB structures found for the given protein profiles.")
        return
    
    structures2 = best_pdb(best_profiles)

    for structure in structures2:
        print(f"Trying {structure.pdb_id}\nresolution: {structure.resolution} Å)")

        for structure in structures2:
            try:
                download_best(structure)
                convert_receptor(structure.pdb_id)
                best_protein = structure
                break
            except ValueError:
                continue
        else:
            print("No suitable PDB structure could be prepared.")
            return
    print("Selected protein:", best_protein)

#Ligand search & Preparation
    ligand = input("\nLigand: ")
    search_ligand(ligand)
    convert_ligand(ligand)

#Docking
    print(find(best_protein.pdb_id))
    dock(best_protein.pdb_id, ligand)


if __name__ == "__main__":
    main()