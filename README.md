# SASA Analysis for Icosahedral Virus Structures and PDB/UniProt Metadata Scraper

MarcoPineapollo. (2026). MarcoPineapollo/Virus-structural-analysis-code: Initial release (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.19475950

This repository provides a standardized workflow for structural biology data processing. It includes tools for protein assembly retrieval, structural cleaning, solvent accessibility calculation, and UniProt metadata extraction.

## Core Functionalities
* Scrape information from PDB and Uniprot databases for a given PDB ID or Uniprot ID, respectively.
* Download `.cif.gz` biological assembled viral structures from the PDB.
* Perform structural cleanup via PyMOL command functions to remove HOH and any non-protein components.
* Calculate the residue-level SASA using the Shrake-Rupley algorithm.

## Project Structure
* `src/`: Contains Python source code (Uniprot_scraper.py, PDB_scraper.py, virSASA.py, and SASA_calculation_example.py).
* `example_data/`: Directory for example input structure datasets (.cif.gz).
* `example_results/`: Directory for example output SASA results.

## Installation
**Prerequisites**
* Python 3.8 or higher
* BioPython 1.85 or higher: Can install via Conda, following the instructions at https://biopython.org/wiki/Packages
* PyMOL: Required for the structure cleanup step. Installation via Conda is recommended, following the instructions at https://pymol.org/conda/

## Usage
Please follow the step-by-step instructions detailed in `SASA_calculation_example.py` to calculate SASA per residue using the assembled structure PDBID 2BPA. All required functions for the calculation are defined in `virSASA.py`, `Uniprot_scraper.py`, and `PDB_scraper.py`

## Validation
Users can verify the SASA calculation by processing the example structure 2BPA. The raw cif structure and cif structure after cleanup were shared in the `example_data` folder. The results of SASA (n_points = 300, probe radius = 1.4) and metadata were provided in the `example_results` folder.

## License
This project is licensed under the MIT License - see the LICENSE file for details.

## Contact
Chonglin Zhu, czhu24@buffalo.edu or chonglinzhu1998@gmail.com, Department of Civil, Structural and Environmental Engineering, University at Buffalo. 
Yinyin Ye, yinyinye@buffalo.edu, Department of Civil, Structural and Environmental Engineering, University at Buffalo.
