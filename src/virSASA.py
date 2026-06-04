# -*- coding: utf-8 -*-
"""
Created on Sat May 31 22:07:53 2025

@author: yinyinye
"""
from Bio.PDB.SASA import ShrakeRupley
import pandas as pd
import requests
import os
from pathlib import Path



def download_structure(pdb_id, save_path=None):
    """
    For a given PDBID, download assembled capsid .cif.gz file from ftp
    Default pathway to save the downloaded file: Downloads\PDB_Structures
    """
    if not isinstance(pdb_id, str):
        raise TypeError("pdb_id must be a single string (e.g., '2BPA')")
        
    print("Step 1: Starting download process...")
    
    # If no save_path is provided, default to the user's Downloads folder
    if save_path is None:
        save_path = os.path.join(Path.home(), "Downloads", "PDB_Structures")
        print(f"  No save path provided. Defaulting to: {save_path}")
    
    if not os.path.exists(save_path): os.makedirs(save_path)
    
    pdb_id = pdb_id.strip().upper()
    middle_two = pdb_id[1:3].lower()
    url = f"https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/divided/{middle_two}/{pdb_id.lower()}-assembly1.cif.gz"
    target_file = os.path.join(save_path, f"{pdb_id}.cif.gz")
        
    if os.path.exists(target_file):
        print(f"  Skipping {pdb_id}, already exists.")
        return
        
    try:
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(target_file, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            print(f"  Successfully downloaded {pdb_id}.cif.gz")
        else:
            print(f"  Failed to download {pdb_id}. HTTP Status: {response.status_code}")
   
    except Exception as e:
        print(f"  Error downloading {pdb_id}: {e}")


def check_res(Structure):
    """
    Extracts unique non-standard residue names from a parsed structure.
    Structure:
        parsed cif structure using Bio.PDB.MMCIFParser.MMCIFParser
    Returns:
        List of unique non-standard residue names.
    """
    
    # List of 20 standard amino acids
    standard_aas = ['ALA', 'ILE', 'LYS',  'TRP', 'ASN', 'LEU', 'VAL', 'HIS', 'MET', 'PHE', 'ARG', 'SER', 'GLU', 'GLN', 'CYS', 'ASP', 'GLY', 'TYR', 'PRO', 'THR']

    # Extract unique non-standard residue names
    unique_residues = set()
    for chain in Structure.get_chains():
        for residue in chain:
            if residue.get_resname() not in standard_aas:
                unique_residues.add(residue.get_resname())

    # Return the list of unique non-standard residue names
    return (list(unique_residues)) if unique_residues else "None"


def virSASA(structure, level, n_points=100, probe_radius=1.4):
    
    '''
    structure_clean is the parsed Structure object
    '''
    if level not in ["A", "R"]:
        raise ValueError("Invalid level. Level must be 'A' (atom) or 'R' (residue).")
    sr = ShrakeRupley(n_points=n_points, probe_radius=probe_radius)
    sr.compute(structure[0], level=level)
    return structure

def resSASA_extract(structure, resname=None, chainid=None):
    model = structure[0]
    complete_list = []
    for chain in model:
        for res in chain:
            complete_list.append((chain.get_id().split("_")[0],res.get_id()[1],res.get_resname(),round(res.sasa,2)))

    df = pd.DataFrame(complete_list, columns = ['chainid','resnum','resname','SASA'])
    
    #sum sasa by chain and residue
    df_group = df.groupby(['chainid','resnum','resname'])['SASA'].sum().reset_index()

    if resname is not None and resname not in ['ALA', 'ILE', 'LYS',  'TRP', 'ASN', 'LEU', 'VAL', 'HIS', 'MET', 'PHE', 'ARG', 'SER', 'GLU', 'GLN', 'CYS', 'ASP', 'GLY', 'TYR', 'PRO', 'THR']:
        raise ValueError("Invalid residue names")

    if resname and chainid:
        if chainid not in df_group['chainid'].unique():
            print("Invalid chain id, please use PDB_scraper to confirm the chains in the structure")
        result = df_group[(df_group['resname'] == resname) & (df_group['chainid'] == chainid)]
    elif resname:
        result = df_group[df_group['resname'] == resname]
    elif chainid:
        if chainid not in df_group['chainid'].unique():
            print("Invalid chain id, please use PDB_scraper to confirm the chains in the structure")
        result = df_group[df_group['chainid'] == chainid]
    else:
        result = df_group
        
    return result

def map_chain(sasatable, df_PDB_uniprot):
    mapping = {}
    for index, row in df_PDB_uniprot.iterrows():
        uniprot_id_gene = row['UniprotID_Protein']
        auth_chains = [chain.strip() for chain in row['Auth_Chain'].split(',')]
        for chain in auth_chains:
            mapping[chain] = uniprot_id_gene
    sasatable['chainid'] = sasatable['chainid'].map(mapping)
    return sasatable
    
    