# -*- coding: utf-8 -*-
"""
Created on Wed Apr  8 15:53:57 2026

@author: dmudo
"""
#!/usr/bin/env python3
"""
Script to extract UniProt information for PDB entries from Table S1.

Requirements:
    pip install pandas openpyxl requests

Usage:
    python uniprot_scraper.py
"""

import pandas as pd
import requests
import time
import os
from pathlib import Path

# ==========================================
# --- CONFIGURATION AREA ---
# ==========================================

# 1. Set the folder location where your input file is stored
# Example: r"C:\Users\Name\Desktop\Data" or "./data"
INPUT_FOLDER = r"D:\Structural Analysis\SASA calculation\20260217\Structural-Bioinfo-Tools\data"

# 2. Set the name of your input Excel file
INPUT_FILENAME = "PDB_list.xlsx"

# 3. Set the output folder (We use the same folder by default)
OUTPUT_FOLDER = INPUT_FOLDER

# 4. Set the name of your output Excel file
OUTPUT_FILENAME = "PDB_UniProt_Info.xlsx"

# ==========================================


def get_uniprot_ids_from_pdb(pdb_id):
    """Get UniProt IDs associated with a PDB entry using RCSB PDB API."""
    uniprot_ids = []
    
    try:
        url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            polymer_entities = data.get('rcsb_entry_container_identifiers', {}).get('polymer_entity_ids', [])
            
            for entity_id in polymer_entities:
                entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
                entity_response = requests.get(entity_url, timeout=30)
                
                if entity_response.status_code == 200:
                    entity_data = entity_response.json()
                    uniprot_refs = entity_data.get('rcsb_polymer_entity_container_identifiers', {}).get('uniprot_ids', [])
                    uniprot_ids.extend(uniprot_refs)
        
        # Alternative: Try EBI PDBe API
        if not uniprot_ids:
            pdbe_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
            pdbe_response = requests.get(pdbe_url, timeout=30)
            
            if pdbe_response.status_code == 200:
                pdbe_data = pdbe_response.json()
                if pdb_id.lower() in pdbe_data:
                    uniprot_data = pdbe_data[pdb_id.lower()].get('UniProt', {})
                    uniprot_ids = list(uniprot_data.keys())
    
    except Exception as e:
        print(f"  Error fetching UniProt IDs for {pdb_id}: {e}")
    
    return list(set(uniprot_ids))

def get_uniprot_info(uniprot_id):
    """Get organism, strains, and taxonomic lineage from UniProt."""
    info = {
        'organism': '',
        'strains': '',
        'taxonomic_lineage': ''
    }
    
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            
            # Get organism info
            organism_data = data.get('organism', {})
            info['organism'] = organism_data.get('scientificName', '')
            
            # Get taxonomic lineage
            lineage = organism_data.get('lineage', [])
            if lineage:
                info['taxonomic_lineage'] = ' > '.join(lineage)
            
            # ===== Get Strains =====
            strains_list = []
            
            # Method 1: Check 'strains' field directly in the root
            if 'strains' in data:
                strains_list.extend(data['strains'])
            
            # Method 2: Check organism section for strains
            if 'strains' in organism_data:
                strains_list.extend(organism_data['strains'])
            
            # Method 3: Check organismHosts for strains
            for host in data.get('organismHosts', []):
                if 'strains' in host:
                    strains_list.extend(host['strains'])
            
            # Method 4: Check comments section
            for comment in data.get('comments', []):
                comment_type = comment.get('commentType', '')
                
                # STRAIN comments
                if comment_type == 'STRAIN':
                    texts = comment.get('texts', [])
                    for t in texts:
                        val = t.get('value', '') if isinstance(t, dict) else str(t)
                        if val:
                            strains_list.append(val)
                
                # Also check BIOPHYSICOCHEMICAL PROPERTIES or other sections
                if 'strains' in comment:
                    strains_list.extend(comment['strains'])
            
            # Method 5: Check features for strain info
            for feature in data.get('features', []):
                if 'strain' in feature.get('description', '').lower():
                    strains_list.append(feature.get('description', ''))
            
            # Method 6: Check references for strain info
            for ref in data.get('references', []):
                source = ref.get('source', {})
                if 'strain' in source:
                    strain_val = source['strain']
                    if isinstance(strain_val, list):
                        strains_list.extend(strain_val)
                    else:
                        strains_list.append(strain_val)
            
            # Combine all strains found
            strains_list = [s for s in strains_list if s]  # Remove empty
            strains_list = list(dict.fromkeys(strains_list))  # Remove duplicates, keep order
            info['strains'] = '; '.join(strains_list)
    
    except Exception as e:
        print(f"  Error fetching info for {uniprot_id}: {e}")
    
    return info

def main():
    print(f"\n--- Initializing UniProt Scraper ---")
    
    # Construct full file paths
    input_file = os.path.join(INPUT_FOLDER, INPUT_FILENAME)
    output_file = os.path.join(OUTPUT_FOLDER, OUTPUT_FILENAME)
    
    # Create output folder if it doesn't exist
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
        print(f"Created output folder at: {OUTPUT_FOLDER}")
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Could not find the input file at '{input_file}'")
        print("Please check your CONFIGURATION AREA settings.")
        return
    
    print(f"Reading {input_file}...")
    df = pd.read_excel(input_file, sheet_name='Sheet1')
    
    if 'Entry ID' not in df.columns:
        print("Error: The input Excel file must contain a column named 'Entry ID'.")
        return

    pdb_ids = df['Entry ID'].unique().tolist()
    print(f"Found {len(pdb_ids)} unique PDB entries")
    
    results = []
    
    for i, pdb_id in enumerate(pdb_ids):
        print(f"Processing {i+1}/{len(pdb_ids)}: {pdb_id}")
        
        uniprot_ids = get_uniprot_ids_from_pdb(pdb_id)
        
        if not uniprot_ids:
            results.append({
                'Entry ID': pdb_id,
                'UniProt ID': '',
                'Organism': '',
                'Strains': '',
                'Taxonomic lineage': ''
            })
            print(f"  No UniProt IDs found")
        else:
            print(f"  Found {len(uniprot_ids)} UniProt IDs: {uniprot_ids}")
            
            for uniprot_id in uniprot_ids:
                info = get_uniprot_info(uniprot_id)
                
                results.append({
                    'Entry ID': pdb_id,
                    'UniProt ID': uniprot_id,
                    'Organism': info['organism'],
                    'Strains': info['strains'],
                    'Taxonomic lineage': info['taxonomic_lineage']
                })
                
                org_display = info['organism'][:30] if info['organism'] else 'N/A'
                strain_display = info['strains'][:20] if info['strains'] else 'N/A'
                print(f"    {uniprot_id}: {org_display}... | Strain: {strain_display}")
        
        time.sleep(0.3)
    
    output_df = pd.DataFrame(results)
    output_df.to_excel(output_file, index=False, sheet_name='Sheet1')
    
    print(f"\n=== COMPLETE ===")
    print(f"Total PDB entries: {len(pdb_ids)}")
    print(f"Total rows in output: {len(results)}")
    print(f"Entries with strain info: {(output_df['Strains'] != '').sum()}")
    print(f"Saved to: {output_file}")

if __name__ == '__main__':
    main()
