import pytest
import BioDock
from Retrieval import PDBStructure

#Retrieval
def test_search_uniprot():
    protein1 = "AMPK"
    gene1 = "PRKAA1"
    organism1 = "Homo sapiens"

    protein2 = "amylase"
    gene2 = ""
    organism2 = "Homo sapiens"
    assert BioDock.search_uniprot(protein1, organism1, gene1) == [{'accession': 'Q13131', 'protein_name': "5'-AMP-activated protein kinase catalytic subunit alpha-1", 'gene': 'PRKAA1', 'organism': 'Homo sapiens', 'entry_name': 'AAPK1_HUMAN'}, {'accession': 'Q96E92', 'protein_name': "5'-AMP-activated protein kinase catalytic subunit alpha-1", 'gene': 'PRKAA1', 'organism': 'Homo sapiens', 'entry_name': 'Q96E92_HUMAN'}]
    assert BioDock.search_uniprot(protein2, organism2, gene2) == [{'accession': 'P0DTE8', 'protein_name': 'Alpha-amylase 1C', 'gene': 'AMY1C', 'organism': 'Homo sapiens', 'entry_name': 'AMY1C_HUMAN'}, {'accession': 'P0DTE7', 'protein_name': 'Alpha-amylase 1B', 'gene': 'AMY1B', 'organism': 'Homo sapiens', 'entry_name': 'AMY1B_HUMAN'}, {'accession': 'P0DUB6', 'protein_name': 'Alpha-amylase 1A', 'gene': 'AMY1A', 'organism': 'Homo sapiens', 'entry_name': 'AMY1A_HUMAN'}, {'accession': 'P04746', 'protein_name': 'Pancreatic alpha-amylase', 'gene': 'AMY2A', 'organism': 'Homo sapiens', 'entry_name': 'AMYP_HUMAN'}, {'accession': 'P19961', 'protein_name': 'Alpha-amylase 2B', 'gene': 'AMY2B', 'organism': 'Homo sapiens', 'entry_name': 'AMY2B_HUMAN'}, {'accession': 'Q04446', 'protein_name': '1,4-alpha-glucan-branching enzyme', 'gene': 'GBE1', 'organism': 'Homo sapiens', 'entry_name': 'GLGB_HUMAN'}, {'accession': 'Q7RTS3', 'protein_name': 'Pancreas transcription factor 1 subunit alpha', 'gene': 'PTF1A', 'organism': 'Homo sapiens', 'entry_name': 'PTF1A_HUMAN'}, {'accession': 'Q9ULW5', 'protein_name': 'Ras-related protein Rab-26', 'gene': 'RAB26', 'organism': 'Homo sapiens', 'entry_name': 'RAB26_HUMAN'}, {'accession': 'P35573', 'protein_name': 'Glycogen debranching enzyme', 'gene':'AGL', 'organism': 'Homo sapiens', 'entry_name': 'GDE_HUMAN'}, {'accession': 'Q07837', 'protein_name': 'Amino acid transporter heavy chain SLC3A1', 'gene': 'SLC3A1', 'organism': 'Homo sapiens', 'entry_name': 'SLC31_HUMAN'}]

    with pytest.raises(KeyError):
        BioDock.search_uniprot("", "", "")

def test_search_pdb():
    profiles = [{'accession': 'Q13131', 'protein_name': "5'-AMP-activated protein kinase catalytic subunit alpha-1", 'gene': 'PRKAA1', 'organism': 'Homo sapiens', 'entry_name': 'AAPK1_HUMAN'}, {'accession': 'Q96E92', 'protein_name': "5'-AMP-activated protein kinase catalytic subunit alpha-1", 'gene': 'PRKAA1', 'organism': 'Homo sapiens', 'entry_name': 'Q96E92_HUMAN'}]
    assert BioDock.search_pdb(profiles[0]) == ['4RED', '4RER', '4REW', '5EZV', '6C9F', '6C9G', '6C9H', '6C9J', '7JHG', '7JHH', '7JIJ', '7M74']
    assert BioDock.search_pdb(profiles[1]) == []

def test_best_pdb():
    structures = BioDock.best_pdb([
        '4RED', '4RER', '4REW', '5EZV','6C9F', '6C9G', '6C9H', '6C9J','7JHG', '7JHH', '7JIJ', '7M74'])
    assert structures[0] == PDBStructure(pdb_id='6C9H', method='X-RAY DIFFRACTION',resolution=2.65)
    with pytest.raises(IndexError):
        BioDock.best_pdb([])