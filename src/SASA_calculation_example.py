# -*- coding: utf-8 -*-
"""
Structural Analysis Pipeline
Created: Apr 7 2026
Author: Chonglin Zhu; Yinyin Ye
"""

import os
import pandas as pd
import numpy as np

import gzip
from pathlib import Path
from pymol import cmd

from Bio.PDB.MMCIFParser import MMCIFParser

from PDB_scraper import PDB_scraper
from Uniprot_scraper import Uniprot_scraper
from virSASA import check_res
from virSASA import virSASA
from virSASA import resSASA_extract
from virSASA import map_chain
from virSASA import download_structure

# =================================================
#Step 1: Extract information from a given PDB ID
# =================================================
PDBID = "2BPA"
field_PDB = ["Auth_Chain", "UniProtID", "PDBMoleculeName", "Organism", "PDBSequence"]
df_PDB = PDB_scraper(PDBID, field_PDB) # use PDB_scraper to extract structural information from the PDB database

# based on the UniProtID of the protein, extract information from UniProt using Uniprot_scraper; 
# if the protein in the PDB structure does not have associated UniprotID, skip this step
uniprotids = df_PDB["UniProtID"].replace(["None", None], np.nan).dropna()
fields_uniprot = ["ProteinName","Lineage", "GeneName", "UniProtEntryName"]
Uniprot_results = []
for uniprotid in uniprotids:
    Uniprot_result = Uniprot_scraper(uniprotid, fields_uniprot)
    if Uniprot_result is not None:
        Uniprot_results.append(Uniprot_result)
df_uniprot = pd.DataFrame(Uniprot_results)

df_PDB_uniprot = pd.merge(df_PDB, df_uniprot, on="UniProtID", how="inner") #combine the PDB and Uniprot information for the PDB structure
df_PDB_uniprot["UniprotID_Protein"] = df_PDB_uniprot["UniProtID"].astype(str) + "_" + df_PDB_uniprot["ProteinName"].astype(str)
# We recommend manually reviewing gene names for all viruses, particularly those undergo polyprotein processing (e.g., picornaviridae)
# We noted that for many picornaviridae proteins, VP1, VP2, VP3, VP4 proteins will be all called as polyprotein

df_PDB_uniprot.to_csv("PDB_UniProt_Metadata.csv", index=False) #save results as csv

# ==================================================================================
#Step 2: Download the assembled virus capsid strucutre (cif.gz) from PDB database
# ==================================================================================
cif_structure = download_structure(PDBID)


# ================================
#Step 3: PDB structure cleanup
# ================================
cif_structure = f'{PDBID}.cif.gz'

parser = MMCIFParser(QUIET=True)

file_path = os.path.join(Path.home(), "Downloads", "PDB_Structures", cif_structure)

with gzip.open(file_path, "rt") as handle:
    structure = parser.get_structure(PDBID, handle)

#remove HOH in the structure
'''
cmd.delete("all")
cmd.load(cif_structure)
cmd.remove("resn HOH")
cmd.save(f"{PDBID}_HOH_clean.cif",f"{PDBID}")
'''

#remove all non amino acid residues in the structure
residues_to_remove = check_res(structure) #return the list of resn to remove
cmd.delete("all")
cmd.load(cif_structure)
for resn in residues_to_remove:
    cmd.remove(f"resn {resn}")
cmd.save(f"{PDBID}_all_clean.cif",f"{PDBID}")


# ================================
#Step 4: SASA calculation
# ================================
cif_structure_clean = f'{PDBID}_all_clean.cif'
# Parse the all clean CIF file
structure_clean = parser.get_structure("structure_clean", cif_structure_clean)

# Calculate SASA at residue (R) level. default parameter: n_point = 100, probe_radius=1.4
structure_sasa = virSASA(structure_clean, level= 'R', n_points=300)

# Combine SASA for all copies of the same protein
sasatable = resSASA_extract(structure_sasa)
sasatable_map = map_chain(sasatable, df_PDB_uniprot)
sasatable_final = sasatable_map.groupby(['chainid','resnum','resname'])['SASA'].sum().reset_index()

# Export the results as a csv file
sasatable_final.to_csv(f'{PDBID}_all_clean_300_SASA.csv', index=False) #modify the output file name if needed



