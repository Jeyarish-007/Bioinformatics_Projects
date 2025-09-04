# Documentation: PubChem CID Extractor & SDF Downloader

---

## 📑 Table of Contents
1. [Project Overview](#project-overview)  
2. [Aim](#aim)  
3. [Objectives](#objectives)  
4. [Methodology](#methodology)  
5. [Significance](#significance)  
6. [Results](#results)  
7. [Conclusion](#conclusion)  
8. [Usage](#usage)  
9. [Requirements](#requirements)  
10. [License](#license)  
11. [Citation](#citation)  
12. [Contact](#contact)  

---

## Project Overview
The **PubChem CID Extractor & SDF Downloader** is a Python-based workflow for automating compound data enrichment from the **PubChem PUG REST API**.  
It extracts PubChem Compound IDs (CIDs) from phytochemical or metabolite names, downloads 3D SDF structure files, and enriches datasets with hyperlinks for direct access.

This pipeline is designed for researchers in **cheminformatics, drug discovery, natural product research, and molecular docking studies**.

---

## Aim
To develop a **scalable and automated pipeline** for retrieving **PubChem CIDs** and downloading **3D molecular structures (SDF format)** for compound datasets.

---

## Objectives
- Automate the mapping of compound names to PubChem CIDs.  
- Provide options for bulk downloading **3D SDF structure files**.  
- Enrich datasets with **descriptive filenames and PubChem hyperlinks**.  
- Ensure reproducibility and **easy integration** with computational chemistry workflows.  

---

## Methodology
1. **Input Preparation**:  
   - Accepts input datasets in `.xlsx` or `.csv` format containing compound names.  
   
2. **CID Extraction**:  
   - Uses PubChem **PUG REST API** to map compound names to PubChem IDs.  
   - Results are saved as enriched CSV/Excel files with both prefixed (`CID:XXXX`) and numeric IDs.  

3. **Structure Download**:  
   - Downloads **3D SDF files** for valid CIDs.  
   - Handles retries, skips duplicates, and logs failures.  
   
4. **Data Enrichment**:  
   - Filenames include **CompoundName_IMPPATID_CID.sdf** for traceability.  
   - Excel/CSV outputs include **clickable download links** for 3D structures.  

5. **Output**:  
   - Structured storage of SDF files, enriched datasets, and log files.  

---

## Significance
- Enables **large-scale compound library preparation** for docking, QSAR, and MD simulations.  
- Saves time by automating tedious manual lookups.  
- Provides **traceable and shareable datasets** with reproducible compound annotations.  
- Ensures FAIR (Findable, Accessible, Interoperable, Reusable) compliance for compound datasets.  

---

## Results
- Successfully retrieves **PubChem CIDs** for a wide range of phytochemicals.  
- Downloads **3D conformers in SDF format** for valid CIDs.  
- Generates **clickable datasets** with direct links to PubChem.  
- Provides logs for failed or missing compounds for easy re-checking.  

---

## Conclusion
The pipeline streamlines the process of **linking compound datasets with PubChem**, downloading 3D structures, and generating publication-ready enriched datasets.  
It is flexible, modular, and can be integrated into **drug discovery, bioinformatics, and cheminformatics pipelines**.  

---

## Usage
```bash
# Extract CIDs from compound names
python src/extract_cid.py --input data/example_input.xlsx --output data/imppat_enriched.csv

# Download SDF files
python src/download_sdf.py --input data/imppat_enriched.csv --output outputs/sdf_files/

# Generate Excel hyperlinks only
python src/hyperlink_generator.py --input data/example_input.xlsx --output data/with_links.xlsx
```

## Requirements

Python 3.8+

Required libraries:
- pandas
- requests
- openpyxl
- Install dependencies:
```bash
pip install -r requirements.txt
```

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this pipeline in your research, please cite:

> Jeyarish V *PubChem CID Extractor & SDF Downloader*. GitHub Repository, 2025.  
> https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/pubchem-cid-extractor

---

## Contact

GitHub: [@Jeyarish-007](https://github.com/Jeyarish-007)
For issues, suggestions, or collaboration, please open an issue or contact directly.

---
