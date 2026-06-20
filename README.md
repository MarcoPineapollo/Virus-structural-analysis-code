# SASA Analysis for Icosahedral Virus Structures and PDB/UniProt Metadata Scraper

## Core Functionalities
* Retrieve metadata from the Protein Data Bank (PDB) and UniProt databases using specific PDB ID or UniProt ID.
* Download gzipped CIF (`.cif.gz`) files from the RCSB PDB using specific PDB ID.
* Execute structural cleanup via PyMOL command-line functions to remove water molecules and all non-protein components.
* Quantify residue-level Solvent Accessible Surface Area (SASA) utilizing the Shrake-Rupley algorithm.
* Summarize SASA results by each protein for total and maximum SASA values.

## Project Structure
* `src/`: Contains Python source code
  * `Uniprot_scraper.py`
  * `PDB_scraper.py`
  * `virSASA.py`
  * `SASA_calculation_example.py`
* `example_data/`: Directory for example input structure datasets
  * `2BPA.cif.gz`
  * `2BPA_all_clean.cif`
  * `6Q5U.cif.gz`
  * `6Q5U_all_clean.cif`
* `example_results/`: Directory for example output SASA results
  * `2BPA_all_clean_300_SASA.csv`
  * `2BPA_all_residues_sasa_summary.csv`
  * `6Q5U_all_clean_300_SASA.csv`
  * `6Q5U_all_residues_sasa_summary.csv`

## Installation
* **BioPython**: Install via Conda, following the instructions at https://biopython.org/wiki/Packages
* **PyMOL**: Install via Conda, following the instructions at https://pymol.org/conda/

## Prerequisites
* Python 3.11.14
* BioPython 1.85
* PyMOL 3.1.7.2

## Usage
Please follow the step-by-step instructions detailed in `SASA_calculation_example.py` to calculate SASA using the example PDB IDs. All required functions for the calculation are defined in `virSASA.py`, `Uniprot_scraper.py`, and `PDB_scraper.py`.

## Validation
Users can verify the SASA calculation by processing the example PDB IDs. The raw cif structures (`{PDBID}.cif.gz`) and the cleaned-up cif structures (`{PDBID}.all_clean.cif`) were shared in the `example_data` folder. The resulting SASA values per residue (`{PDBID}_all_clean_300_SASA.csv`) and summarized SASA values per protein (`{PDBID}_all_residues_sasa_summary.csv`) were provided in the `sample_results` folder.

## License
This project is licensed under the MIT License - see the `LICENSE` file for details.

## Contact
* **Chonglin Zhu**, czhu24@buffalo.edu or chonglinzhu1998@gmail.com, Department of Civil, Structural and Environmental Engineering, University at Buffalo.
* **Yinyin Ye**, yinyinye@buffalo.edu, Department of Civil, Structural and Environmental Engineering, University at Buffalo.
