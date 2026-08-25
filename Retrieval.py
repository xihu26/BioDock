import requests
from dataclasses import dataclass
from rcsbapi.data import DataQuery

from Preparation import convert_receptor

@dataclass
class PDBStructure:
    pdb_id: str
    method: str
    resolution: float

def main():
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
        print("ACCESSION:", profile["accession"])
        print("PDB IDS:", pdb_ids)
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
#Ligand search
    ligand = input("Ligand: ")
    search_ligand(ligand)


def search_uniprot(protein, organism, gene):
    query = f"{protein} AND organism_name:\"{organism}\""
    if gene:
        query += f" AND gene:{gene}"
    response = requests.get("https://rest.uniprot.org/uniprotkb/search", params={"query": query, "format": "json", "size": "10"})
    profiles = []
    results = response.json()["results"]

    for result in results:
        profile = {
        "accession": result["primaryAccession"],
        "protein_name": result["proteinDescription"]["recommendedName"]["fullName"]["value"] if "recommendedName" in result["proteinDescription"] else result["proteinDescription"]["submittedName"][0]["fullName"]["value"] if "submittedName" in result["proteinDescription"] else "Unknown",
        "gene": result["genes"][0]["geneName"]["value"] if "genes" in result  else "Unknown",
        "organism": result["organism"]["scientificName"],
        "entry_name": result["uniProtkbId"]
        }
        profiles.append(profile)

    return profiles

def search_pdb(profile):
    pdb_ids = []
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": profile["accession"]
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                        "operator": "exact_match",
                        "value": "UniProt"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        }
    }

    response = requests.post("https://search.rcsb.org/rcsbsearch/v2/query",json=query)
    if response.status_code == 204:
        return []
    else:
        pdb_ids = [pdb["identifier"] for pdb in response.json()["result_set"]]
        return pdb_ids


def best_pdb(pdb_ids):
    query = DataQuery(input_type="entries", input_ids=pdb_ids, return_data_list=["exptl.method", "rcsb_entry_info.resolution_combined", "struct.title"])
    results = query.exec()
    structures = []

    for entry in results["data"]["entries"]:
        resolution = entry["rcsb_entry_info"]["resolution_combined"]
        if resolution is None:
            continue
        structure = PDBStructure(
        pdb_id = entry["rcsb_id"],
        method = entry["exptl"][0]["method"],
        resolution = resolution[0]
        )
        structures.append(structure)

    structures.sort(key=lambda structure: structure.resolution)
    return structures

def download_best(best_protein):
    response = requests.get(f"https://files.rcsb.org/download/{best_protein.pdb_id}.pdb")
    with open(f"proteins/{best_protein.pdb_id}.pdb", "wb") as file:
        file.write(response.content)
        return best_protein


def search_ligand(ligand):
    cid = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{ligand}/cids/JSON").json()["IdentifierList"]["CID"][0]
    response = requests.get(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/SDF?record_type=3d")
    with open(f"ligands/{ligand}.sdf", "wb") as file:
        file.write(response.content)


if __name__ == "__main__":
    main()


    