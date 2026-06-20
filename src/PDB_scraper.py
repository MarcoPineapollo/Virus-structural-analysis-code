# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:55:22 2025

@author: yinyinye
"""

import requests
import pandas as pd

def PDB_scraper(PDBID, fields=None):
    """
    list of fiels:
    PDBID
    Auth_Chain: the name of chain defined by the author
    UniProtID: the Uniprot ID associated with the protein in the structure
    MutationCount
    PDBMoleculeName: the molecular name saved with the PDB structure (input from the author; may not match with the protein name in Uniprot database)
    Organism: input by the author
    PDBSequence: protein sequence
    PubMedID: publication of the structure
    """
    
    # Define the GraphQL query
    query = """
    query($PDBID: String!) {
      entry(entry_id:$PDBID) {
        pubmed {
          rcsb_pubmed_container_identifiers {
           pubmed_id
          }
        }
        polymer_entities {
          rcsb_polymer_entity_container_identifiers {
           entry_id
           asym_ids
           auth_asym_ids
           uniprot_ids
          }
          uniprots {
            rcsb_uniprot_container_identifiers {
             uniprot_id
            }
            rcsb_uniprot_protein {
              name {
                value
              }
              source_organism{
                scientific_name
              }
            }
          }
          rcsb_polymer_entity {
            pdbx_description
          }
          entity_poly {
            rcsb_entity_polymer_type
            pdbx_seq_one_letter_code_can
            rcsb_sample_sequence_length
            rcsb_mutation_count
          }
        }
      }
    }
    """

    # Define the URL and variables
    url = "https://data.rcsb.org/graphql"
    variables = {"PDBID": PDBID}  
    
    #variables = {"PDBID": '6YFH'}
    
    # Make the request
    try:
        response = requests.post(url, json={"query": query, "variables": variables})
        response.raise_for_status()  # Raise an exception for bad status codes
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}")
        return None

    # Parse the response
    result = response.json()

    # Check if the response contains the expected data
    if 'data' not in result or 'entry' not in result['data']:
        print("Invalid response format")
        return None
    
    default_fields = ["Auth_Chain", "PDBMoleculeName","UniProtID"]
    if fields is None:
        fields = default_fields
    # Extract relevant data
    data = {}
    data["PDBID"] = [entity['rcsb_polymer_entity_container_identifiers']['entry_id'] 
                           for entity in result['data']['entry']['polymer_entities']]
    for field in fields:
        if field == "Auth_Chain":
            auth_chains = [','.join(entity['rcsb_polymer_entity_container_identifiers']['auth_asym_ids']) 
                           for entity in result['data']['entry']['polymer_entities']]
            data[field] = auth_chains
        elif field == "UniProtID":
            uniprot_ids = []
            for entity in result['data']['entry']['polymer_entities']:
                if entity['uniprots']:
                    uniprot_ids.append(entity['uniprots'][0]['rcsb_uniprot_container_identifiers']['uniprot_id'])
                else:
                    uniprot_ids.append(None)
            data[field] = uniprot_ids
        elif field == "MutationCount":
            data[field] = [entity['entity_poly']['rcsb_mutation_count'] 
                           for entity in result['data']['entry']['polymer_entities']]
        elif field == "PDBMoleculeName":
            data[field] = [entity['rcsb_polymer_entity']['pdbx_description'] 
                           for entity in result['data']['entry']['polymer_entities']]
        elif field == "Organism":
            source_organisms = []
            for entity in result['data']['entry']['polymer_entities']:
                if entity['uniprots']:
                    source_organisms.append(entity['uniprots'][0]['rcsb_uniprot_protein']['source_organism']['scientific_name'])
                else:
                    source_organisms.append(None)
            data[field] = source_organisms
        elif field == "PDBSequence":
            data[field] = [entity['entity_poly']['pdbx_seq_one_letter_code_can'] 
                           for entity in result['data']['entry']['polymer_entities']]
        elif field == "PubMedID":
            pubmed_data = result.get('data', {}).get('entry', {}).get('pubmed')
            if pubmed_data is not None:
                data[field] = pubmed_data.get('rcsb_pubmed_container_identifiers', {}).get('pubmed_id')
            else:
                data[field] = None        
        else:
            print(f"Invalid field: {field}")
            return None

    # Create data frame
    df = pd.DataFrame(data)

    return df
