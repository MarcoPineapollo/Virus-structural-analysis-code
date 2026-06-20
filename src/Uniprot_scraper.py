# -*- coding: utf-8 -*-
"""
Created on Wed May 28 11:06:08 2025

@author: yinyinye
"""

import requests

def Uniprot_scraper(uniprot_id, fields=None):
    '''
    list of fiels:
    UniProtID: Entry ID of the protein in the UniProtKB database
    ProteinName: Protein name saved in the UniProtKB database
    Organism: The source of organism associated with the protein based on the UniProtKB database
    Lineage: Taxonomic lineage; such as Viruses > Monodnaviria (single-stranded DNA viruses) > Sangervirae > Phixviricota > Malgrandaviricetes > Petitvirales > Microviridae (isometric ssDNA phages) > Bullavirinae > Sinsheimervirus > Escherichia phage phiX174 (Bacteriophage phi-X174)
    GeneName: The name of gene
    TaxonID: Taxonomy ID
    UniProtEntryName
    Sequence: Amino acid sequences of the protein
    PubMedID: Any publication associated with the protein
    '''
    
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
    
    #url = "https://rest.uniprot.org/uniprotkb/A0A0P0BSS5"
    params = {"format": "json"}

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()  # Raise an exception for bad status codes
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return None

    data = response.json()
    result = {}
    
    default_fields = ["ProteinName", "Organism", "Gene"]
    if fields is None:
        fields = default_fields
    
    result["UniProtID"]=data.get("primaryAccession",{})
    
    if "ProteinName" in fields:
        protein_name_rec = data.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
        protein_name_sub = None
        if "submissionNames" in data.get("proteinDescription", {}):
            if len(data["proteinDescription"]["submissionNames"]) > 0:
                protein_name_sub = data["proteinDescription"]["submissionNames"][0].get("fullName", {}).get("value")
        result["ProteinName"] = protein_name_rec if protein_name_rec else protein_name_sub

    if "Organism" in fields:
        result["Organism"] = data.get("organism", {}).get("scientificName")

    if "Lineage" in fields:
        lineage = data.get("organism", {}).get("lineage")
        result["Lineage"] = ">".join(lineage)

    if "GeneName" in fields:
        result["GeneName"] = data.get("genes", [{}])[0].get("geneName", {}).get("value", "None")
        
    if "TaxonID" in fields:
        result["TaxonID"]= data.get("organism", {}).get("taxonId")
        
    if "UniProtEntryName" in fields:
        result["UniProtEntryName"]= data.get("uniProtkbId", {})
    
    if "Sequence" in fields:
        result["Sequence"]=data.get("sequence", {}).get("value")
        
    if "PubMedID" in fields:
        ids = [ref.get("citation", {}).get("id") for ref in data.get("references", []) if ref.get("citation", {}).get("id") is not None]
        result["PubMedID"] = ";".join(map(str, ids))
            
    
    return result