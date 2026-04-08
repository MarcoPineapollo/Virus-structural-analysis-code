# -*- coding: utf-8 -*-
"""
Structural Analysis Pipeline
Created: Apr 7 2026
Author: Chonglin
"""


import os
import requests
import re
import pandas as pd
import warnings
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Polypeptide import is_aa

# --- PyMOL Initialization ---
try:
    import pymol
    pymol.pymol_argv = ['pymol', '-qc'] # -q for quiet, -c for command line only (no GUI)
    pymol.finish_launching()
    from pymol import cmd
except ImportError:
    print("Warning: PyMOL module not found. Step 2 will be skipped.")

# ==========================================
# --- CONFIGURATION AREA ---
# ==========================================

# 1. Input your PDB IDs here
PDB_LIST = ['2MS2','5MUU'] # Input PDB IDs in '2MS2', '5MUU' 

# 2. Set your save folder
SAVE_PATH = r"Structural-Bioinfo-Tools\data" #folder location

# Automatically generate output filenames based on the FIRST PDB ID in your list
PRIMARY_TARGET = PDB_LIST[0]
CLEAN_LOG_EXCEL = os.path.join(SAVE_PATH, f"removed_components.xlsx")
METADATA_CSV = os.path.join(SAVE_PATH, f"PDB_Metadata.csv")

# ==========================================


# --- Step 1: Download (FTP Strategy) ---
def step1_download_structures(pdb_ids, save_path):
    print("Step 1: Starting download process...")
    if not os.path.exists(save_path): os.makedirs(save_path)
    
    for pdb_id in pdb_ids:
        pdb_id = pdb_id.strip().upper()
        middle_two = pdb_id[1:3].lower()
        url = f"https://files.wwpdb.org/pub/pdb/data/assemblies/mmCIF/divided/{middle_two}/{pdb_id.lower()}-assembly1.cif.gz"
        target_file = os.path.join(save_path, f"{pdb_id}.cif.gz")
        
        if os.path.exists(target_file):
            print(f"  Skipping {pdb_id}, already exists.")
            continue
        
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

# --- Step 2: Clean and Track (Handles .gz) ---
def step2_clean_structures(pdb_ids, input_dir, output_excel):
    print("Step 2: Starting cleaning process...")
    standard_aas = ['ALA','ILE','LYS','TRP','ASN','LEU','VAL','HIS','MET','PHE','ARG','SER','GLU','GLN','CYS','ASP','GLY','TYR','PRO','THR']
    aas_selection = '+'.join(standard_aas)
    log_data = []

    for pdb_id in pdb_ids:
        pdb_id = pdb_id.strip().upper()
        cif_path = os.path.join(input_dir, f"{pdb_id}.cif.gz")
        if not os.path.exists(cif_path): continue

        cmd.reinitialize()
        cmd.load(cif_path, pdb_id)
        
        removed_sel = f"{pdb_id} and not (polymer.protein and resn {aas_selection})"
        removed_residues = set()
        try:
            model = cmd.get_model(removed_sel)
            for atom in model.atom:
                removed_residues.add(atom.resn)
        except: pass 
        
        cmd.remove(f"{pdb_id} and not polymer.protein")
        cmd.remove(f"{pdb_id} and polymer.protein and not resn {aas_selection}")
        
        clean_filename = f"{pdb_id}_all_clean.cif"
        cmd.save(os.path.join(input_dir, clean_filename), pdb_id)
        
        log_data.append({
            "PDB_ID": pdb_id, 
            "Removed_Components": ", ".join(sorted(list(removed_residues))) if removed_residues else "None"
        })
        cmd.delete("all")

    if log_data:
        pd.DataFrame(log_data).to_excel(output_excel, index=False)
        print(f"  Clean log saved to {output_excel}")

# --- Step 3: Scrape UniProt Information ---
def step3_fetch_metadata(pdb_id_list, output_csv):
    print("Step 3: Scraping metadata...")
    warnings.simplefilter(action='ignore')

    url = "https://data.rcsb.org/graphql"
    query = """
    query($PDBID: String!) {
      entry(entry_id: $PDBID) {
        polymer_entities {
          rcsb_polymer_entity_container_identifiers { auth_asym_ids }
          rcsb_polymer_entity { pdbx_description }
          uniprots {
            rcsb_uniprot_container_identifiers { uniprot_id }
            rcsb_uniprot_protein { name { value } }
          }
          entity_poly { pdbx_seq_one_letter_code_can }
        }
      }
    }
    """
    
    all_rows = []
    for pdb_id in pdb_id_list:
        pdb_id = pdb_id.strip().upper()
        try:
            response = requests.post(url, json={'query': query, 'variables': {"PDBID": pdb_id}})
            data = response.json().get('data', {}).get('entry', {})
            if not data: continue

            for entity in data.get('polymer_entities', []):
                auth_chains = entity.get('rcsb_polymer_entity_container_identifiers', {}).get('auth_asym_ids', [])
                if not auth_chains: continue

                mol_desc = entity.get('rcsb_polymer_entity', {}).get('pdbx_description', "Unknown")
                uniprots = entity.get('uniprots', [])
                
                u_id = uniprots[0].get('rcsb_uniprot_container_identifiers', {}).get('uniprot_id', "NA") if uniprots else "NA"
                u_name = uniprots[0].get('rcsb_uniprot_protein', {}).get('name', {}).get('value', "NA") if uniprots else "NA"
                sequence = entity.get('entity_poly', {}).get('pdbx_seq_one_letter_code_can', "NA")

                # Label Logic
                mol_desc_upper = mol_desc.upper()
                if re.search(r'\b([VC]P\d+)\b', mol_desc, re.IGNORECASE):
                    suffix = re.search(r'\b([VC]P\d+)\b', mol_desc, re.IGNORECASE).group(1).upper()
                elif "COAT PROTEIN" in mol_desc_upper: suffix = "Coat_Protein"
                elif "CAPSID PROTEIN" in mol_desc_upper: suffix = "Capsid_Protein"
                else: suffix = "_".join(re.sub(r'[^\w\s]', '', mol_desc).split()[:2])
                
                initial_label = f"{u_id}_{suffix}"

                all_rows.append({
                    "PDBID": pdb_id,
                    "Auth_Chain": "'" + "','".join(auth_chains) + "'",
                    "First_Chain": auth_chains[0],
                    "UniProtID": u_id,
                    "ProteinName": u_name,
                    "Sequence": sequence,
                    "Molecule_Name": mol_desc,
                    "UniProtID_list": initial_label
                })
        except Exception as e:
            print(f"  Error fetching metadata for {pdb_id}: {e}")

    df = pd.DataFrame(all_rows)
    if not df.empty:
        def disambiguate(group):
            counts = group['UniProtID_list'].value_counts()
            group['UniProtID_list'] = [
                f"{row['UniProtID_list']}_{row['First_Chain']}" if counts[row['UniProtID_list']] > 1 else row['UniProtID_list'] 
                for _, row in group.iterrows()
            ]
            return group
        
        df = df.groupby('PDBID', group_keys=False).apply(disambiguate)
        if 'First_Chain' in df.columns: df = df.drop(columns=['First_Chain'])
        df.to_csv(output_csv, index=False, encoding='utf-8-sig')
        print(f"  Metadata saved to {output_csv}")

# --- Step 4: Calculate SASA ---
def step4_calculate_sasa(pdb_ids, metadata_csv, struct_dir):
    print("Step 4: Starting SASA calculation...")
    if not os.path.exists(metadata_csv):
        print(f"  Error: Metadata file not found at {metadata_csv}")
        return
        
    df_meta = pd.read_csv(metadata_csv)
    mapping = {}
    for _, row in df_meta.iterrows():
        pid = str(row['PDBID']).strip().upper()
        chains = [c.strip().strip("'").strip('"') for c in str(row['Auth_Chain']).split(',')]
        if pid not in mapping: mapping[pid] = {}
        for c in chains: mapping[pid][c] = row['UniProtID_list']

    parser = MMCIFParser(QUIET=True)
    for sid in pdb_ids:
        sid = sid.strip().upper()
        cif_file = os.path.join(struct_dir, f"{sid}_all_clean.cif")
        if not os.path.exists(cif_file) or sid not in mapping: continue
        
        try:
            struct = parser.get_structure(sid, cif_file)
            sr = ShrakeRupley(n_points=300, probe_radius=1.40)
            sr.compute(struct[0], level="R")
            
            residue_data = []
            for chain in struct[0]:
                orig_chain_id = chain.get_id()
                base = orig_chain_id.split('-')[0]
                stripped = base.rstrip('0123456789')
                
                # Matching Logic
                mapped_id = mapping[sid].get(orig_chain_id) or mapping[sid].get(base) or mapping[sid].get(stripped)
                if not mapped_id and len(base) > 0: mapped_id = mapping[sid].get(base[0])
                
                if not mapped_id: continue

                for res in chain:
                    if is_aa(res, standard=True):
                        sasa_val = round(getattr(res, 'sasa', getattr(res, 'SASA', 0.0)), 2)
                        residue_data.append([mapped_id, res.id[1], res.get_resname(), sasa_val])

            if residue_data:
                df = pd.DataFrame(residue_data, columns=['UniProtid', 'AAid', 'AA', 'SASA'])
                grouped_df = df.groupby(['UniProtid', 'AAid', 'AA'], as_index=False)['SASA'].sum()
                out_csv = os.path.join(struct_dir, f"{sid}_all_clean_300.csv")
                grouped_df.to_csv(out_csv, index=False)
                print(f"  Saved SASA data to {out_csv}")
                
        except Exception as e:
            print(f"  Error calculating SASA for {sid}: {e}")

# ==========================================
# --- EXECUTION SWITCHBOARD ---
# ==========================================
if __name__ == "__main__":
    print(f"\n--- Initializing Pipeline for {PDB_LIST} ---")
    
    #ADD OR REMOVE A '#' TO TURN STEPS ON OR OFF
    
    step1_download_structures(PDB_LIST, SAVE_PATH)
    step2_clean_structures(PDB_LIST, SAVE_PATH, CLEAN_LOG_EXCEL)
    step3_fetch_metadata(PDB_LIST, METADATA_CSV)
    #step4_calculate_sasa(PDB_LIST, METADATA_CSV, SAVE_PATH)
    
    print("\n--- Pipeline Complete ---")
