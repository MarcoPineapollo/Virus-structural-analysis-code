# Virus Structural Analysis and UniProt Metadata Scraper

This repository provides a standardized workflow for structural biology data processing. It includes tools for protein assembly retrieval, structural cleaning, solvent accessibility calculation, and UniProt metadata extraction.

## Core Functionalities

* **UniProt Scraper**: Maps PDB entries to UniProt accessions and retrieves organism name, strain information, and taxonomic lineage.
* **Structural Pipeline**: 
    * Downloads `.cif.gz` biological assemblies from the wwPDB.
    * Performs structural cleaning via PyMOL to remove non-protein components.
    * Calculates residue-level SASA using the Shrake-Rupley algorithm.
* **Batch Processing**: Automated handling of multiple entries via Excel-based input.

## Project Structure

* **src/**: Contains Python source code (`protein_pipeline.py` and `uniprot_scraper.py`).
* **data/**: Directory for input datasets and example files.
* **results/**: Directory for output storage and reference datasets.
* **requirements.txt**: List of Python dependencies.

## Installation

### Prerequisites
* Python 3.8 or higher.
* **PyMOL**: Required for the structural cleaning module. Installation via Conda is recommended:
    `conda install -c schrodinger pymol`

### Dependencies
Install the required Python libraries using pip:
`pip install -r requirements.txt`

## Usage

### 1. Metadata Extraction
This script prepares the necessary mapping between PDB entries and UniProt metadata.
* **Input**: An Excel file with an `Entry ID` column located in `data/`.
* **Configuration**: Set `INPUT_FOLDER` and `INPUT_FILENAME` in `src/uniprot_scraper.py`.
* **Execution**: Run `python src/uniprot_scraper.py`.
* **Outcome**: Generates `PDB_UniProt_Info.xlsx` containing organism, strain, and taxonomic data.

### 2. Structural Analysis Pipeline
The pipeline executes the following four steps in sequence. Before running, configure the `PDB_LIST` and `SAVE_PATH` in `src/protein_pipeline.py`.

* **Step 1: Download Structures**
    * Automatically retrieves `.cif.gz` biological assemblies from the wwPDB FTP server.
* **Step 2: Structural Cleaning (Requires PyMOL)**
    * Loads the structure and removes all non-protein components (water, ligands, ions).
    * Filters for standard amino acids to ensure a clean model for analysis.
* **Step 3: Fetch RCSB Metadata**
    * Connects to the RCSB GraphQL API to retrieve polymer entity descriptions and sequence info.
    * Generates a `PDB_Metadata.csv` for chain-to-UniProt mapping.
* **Step 4: SASA Calculation**
    * Uses the Shrake-Rupley algorithm (BioPython) with a 1.4Å probe radius.
    * Outputs residue-level SASA values to individual CSV files (e.g., `XXXX_all_clean_300.csv`).

## Validation
Users can verify the installation by processing the example entry **2MS2** provided in the `data/` folder. Ensure the outputs align with the reference files located in the `results/` directory.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
**Chonglin Zhu**, czhu24@buffalo.edu or chonglinzhu1998@gmail.com, Department of Civil, Structural and Environmental Engineering, University at Buffalo.
**Yinyin Ye**, yinyinye@buffalo.edu, Department of Civil, Structural and Environmental Engineering, University at Buffalo.
