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
PDBID = "8AT5"
field_PDB = ["Auth_Chain", "UniProtID", "PDBMoleculeName", "Organism", "PDBSequence"]
df_PDB = PDB_scraper(PDBID, field_PDB) # use PDB_scraper to extract structural information from the PDB database

uniprotids = df_PDB["UniProtID"].replace(["None", None], np.nan).fillna(value=np.nan)
fields_uniprot = ["ProteinName","Lineage", "GeneName", "UniProtEntryName"]

Uniprot_results = []
for uniprotid in uniprotids:
    # Check if the ID is valid (handles pure Python None and pandas NaN)
    if pd.notna(uniprotid):
        # based on the UniProtID of the protein, extract information from UniProt using Uniprot_scraper
        Uniprot_result = Uniprot_scraper(uniprotid, fields_uniprot)
        # Safely check if the scraper actually returned data
        if Uniprot_result is not None:
            Uniprot_results.append(Uniprot_result)
        else:
            # If the API request failed, append the None dictionary structure
            Uniprot_results.append({"UniProtID": np.nan, **{field: None for field in fields_uniprot}})
    else:
        # If the input ID itself was missing, append the None dictionary structure
        Uniprot_results.append({"UniProtID": np.nan, **{field: None for field in fields_uniprot}})
   
df_uniprot = pd.DataFrame(Uniprot_results)

df_PDB_uniprot = pd.concat([df_PDB,df_uniprot.drop(columns=['UniProtID'])], axis =1)

valid_proteins = df_PDB_uniprot["ProteinName"].dropna()
if valid_proteins.nunique() == len(df_PDB_uniprot):
    # No duplicates found -> Use ProteinName
    df_PDB_uniprot["UniprotID_Protein"] = (
        df_PDB_uniprot["UniProtID"].astype(str) + "_" + df_PDB_uniprot["ProteinName"].astype(str)
    )
else:
    # Duplicates identified -> Fall back to PDBMoleculeNames
    df_PDB_uniprot["UniprotID_Protein"] = (
        df_PDB_uniprot["UniProtID"].fillna("NA").astype(str) + "_" + df_PDB_uniprot["PDBMoleculeName"].astype(str)
    )
# We recommend manually reviewing information extracted in cols 'ProteinName' and 'PDBMoleculeName' for all viruses
# For many Picornaviridae proteins, VP1, VP2, VP3, VP4 proteins are all called as polyprotein in uniprot database
# if that is the case, PDBMoleculeName (if they are unique) will be used to track different protein

'''
df_PDB_uniprot.to_csv("PDB_UniProt_Metadata.csv", index=False) #save results as csv
'''

# ==================================================================================
#Step 2: Download the assembled virus capsid structure (cif.gz) from PDB database and read the compressed structure 
# ==================================================================================
s = download_structure(PDBID) 
#{PDBID}.cif.gz file is located in C:\Downloads\PDB_Structures

# ================================
#Step 3: PDB structure cleanup
# ================================
file_path = os.path.join(Path.home(), "Downloads", "PDB_Structures", f'{PDBID}.cif.gz')

parser = MMCIFParser(QUIET=True)
with gzip.open(file_path, "rt") as handle:
    structure = parser.get_structure(PDBID, handle) #unzip and read the structure

'''
# remove HOH in the structure
cmd.delete("all")
cmd.load(file_path, f"{PDBID}.cif.gz")
cmd.remove("resn HOH")
cmd.save(f"{PDBID}_HOH_clean.cif",f"{PDBID}")
'''

# remove all non amino acid residues in the structure
residues_to_remove = check_res(structure) #return the list of resn to remove
cmd.delete("all")
cmd.load(file_path, f"{PDBID}.cif.gz")
for resn in residues_to_remove:
    cmd.remove(f"resn {resn}")
cmd.save(f"{PDBID}_all_clean.cif",f"{PDBID}")

'''
# remove chains in the structure
chain_ids_list = [chain.get_id() for chain in structure[0]]
cmd.remove("chain O or chain O-*")
cmd.remove("chain P or chain P-*")
cmd.remove("chain Q or chain Q-*")
cmd.save(f"{PDBID}_all_clean_nospike.cif",f"{PDBID}")
'''

# ================================
#Step 4: SASA calculation
# ================================
cif_structure_clean = f'{PDBID}_all_clean.cif'
# Parse the all clean CIF file
structure_clean = parser.get_structure("structure_clean", cif_structure_clean)

# Calculate SASA at residue (R) level. default parameter: n_point = 100, probe_radius=1.4
structure_sasa = virSASA(structure_clean, level= 'R', n_points=300)

# Combine SASA for all copies of the same chain ID
sasatable = resSASA_extract(structure_sasa)

# replace 'chain ID (Auth_chain)' with 'UnirotID_Protein' in the df_PDB_uniprot
sasatable_map = map_chain(sasatable, df_PDB_uniprot)

# Combine SASA from different copies of the same protein
sasatable_final = sasatable_map.groupby(['chainid','resnum','resname'])['SASA'].sum().reset_index()

# Export the results as a csv file
sasatable_final.to_csv(f'{PDBID}_all_clean_300_SASA.csv', index=False) # Save SASA per residue


# ====================================================================================
#Step 5: summarize total and max SASA of 20 common amino acids for each protein in {PDBID} structure
# ====================================================================================
#PDBID = "2BPA"
df = pd.read_csv(f'{PDBID}_all_clean_300_SASA.csv')

# Map standard abbreviations to the 3-letter codes used in PDB files
residues = [
    'MET', 'CYS', 'HIS', 'TRP', 'LYS',
    'TYR', 'ARG', 'ALA', 'ASN', 'ASP',
    'GLN', 'GLU', 'GLY', 'ILE', 'LEU',
    'PHE', 'PRO', 'SER', 'THR', 'VAL'
]

# Establish a baseline dataframe with all unique chain entries
unique_chains = df['chainid'].unique()
summary_df = pd.DataFrame({'chainid': unique_chains})

# Sequentially process each residue from left to right
for res_name in residues:
    # Filter rows matching the 3-letter code
    res_df = df[df['resname'] == res_name]
    
    # 2. FIXED: Custom logic exclusively for MET
    if res_name == 'MET':
        # Total Met uses all data (including resnum == 1)
        total_sasa = res_df.groupby('chainid')['SASA'].sum().reset_index(name='total')
        
        # Max Met filters out resnum == 1 before finding the maximum
        filtered_for_max = res_df[res_df['resnum'] != 1]
        max_sasa = filtered_for_max.groupby('chainid')['SASA'].max().reset_index(name='maximum')
        
        # Combine the custom total and max vectors together
        agg_df = pd.merge(total_sasa, max_sasa, on='chainid', how='outer')
    else:
        # Standard logic for all other 19 residues (Total and Max include all positions)
        agg_df = res_df.groupby('chainid').agg(
            total=('SASA', 'sum'),
            maximum=('SASA', 'max')
        ).reset_index()
        
    # Clean up column naming/capitalization (e.g., 'MET' -> 'Met')
    display_name = res_name.capitalize()
    agg_df = agg_df.rename(columns={
        'total': f'Total {display_name}',
        'maximum': f'Max {display_name}'
    })
    
    # Merge column profiles back into the primary summary table
    summary_df = pd.merge(summary_df, agg_df, on='chainid', how='left')

# Clean up structural omissions with zero values and save 
summary_df = summary_df.fillna(0.00)
summary_df.to_csv(f'{PDBID}_all_residues_sasa_summary.csv', index=False)
