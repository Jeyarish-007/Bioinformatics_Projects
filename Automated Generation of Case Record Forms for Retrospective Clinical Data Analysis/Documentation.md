# Automated Generation of Case Record Forms for Retrospective Clinical Data Analysis

This project presents an open-source Python workflow to automate the conversion of retrospective clinical Excel records into standardized Case Report Forms (CRFs) for research documentation and analysis.

## Table of Contents

- [Project Overview](#project-overview)
- [Aim](#aim)
- [Objectives](#objectives)
- [Methodology](#methodology)
- [Significance](#significance)
- [Results](#results)
- [Conclusion](#conclusion)
- [Rationale for Data Selection](#rationale-for-data-selection)
- [Usage](#usage)
- [Requirements](#requirements)
- [License](#license)

## Project Overview

This repository provides:
- Sample scripts
- CRF templates
- Example (anonymized) data  
to automate the process of transcribing retrospective patient data from a spreadsheet into structured Word Case Record Forms (CRFs), each corresponding to a unique patient record.

## Aim

To standardize and automate the process of transcribing retrospective patient data from Excel into structured Case Report Forms (CRFs) for improved documentation, analysis, and reproducibility in clinical research.

## Objectives

1. **Data Consolidation:** Gather and compile retrospective clinical data from 350+ patients across relevant demographic and clinical variables.
2. **Process Automation:** Automate population of Word CRF templates from structured Excel data using Python.
3. **Documentation Standardization:** Ensure each patient record is uniformly documented in an individual CRF.
4. **Usability & Reproducibility:** Provide an open-source solution usable for similar data transcription tasks.

## Methodology

1. **Data Collection:**  
   - Retrospective extraction of clinical and demographic data for 350+ patients.

2. **Data Entry:**  
   - Manual verification and input of patient data into Excel; each row = one patient.

3. **CRF Template Creation:**  
   - A CRF table was designed in Microsoft Word.

4. **Python Automation:**  
   - The provided Python script:
     - Reads the Excel file
     - Generates a Word document where each page is a formatted CRF, filled for one patient
     - Applies the required formatting and field structure

5. **Quality Assurance:**  
   - Output checked for accuracy and consistency.

## Significance

- **Efficiency:** Reduces manual effort and human error.
- **Reproducibility:** Easy adaptation for new datasets and research.
- **Standardization:** Uniform formatting ideal for medical/legal compliance and secondary analysis.
- **Open Science:** Fully open and generalizable for other researchers.

## Results

- 350+ CRFs auto-generated, appropriately formatted and separated in a single Word document.
- Dates are formatted dd-mm-yyyy (no time).
- The process saved significant time versus manual data transfer.

## Conclusion

This project demonstrates a robust, reproducible, and user-friendly method for transforming retrospective clinical datasets into formatted Case Record Forms using Python, enabling more efficient research and analysis.

## Rationale for Data Selection

- **Real-world data is essential**: Actual institutional patient data enables authentic insight into clinical outcomes.
- **Compatibility with research practices**: CRFs mirror standard research reporting needs.
- **Demonstrates practical utility**: The approach solves a common bottleneck in retrospective studies.

## Usage

1. **Prepare your Excel file**   
   - Each row = one patient.  
   - Columns: Case No, Age, Gender, Diagnosis, Relapse, Date of sample collection.

2. **Adjust file paths** in the script as needed.

3. **Run the script** (see [`Script.py`](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Automated%20Generation%20of%20Case%20Record%20Forms%20for%20Retrospective%20Clinical%20Data%20Analysis/Script.py)):

   ```bash
   pip install pandas openpyxl python-docx
   ```

   Then run the script in your Python environment:

   ```python
   # Example usage
   Script.py
   ```
4. **Result:**  
   - Word document (`.docx`), each page a filled CRF for one patient.

## Requirements

- Python 3.7+
- `pandas`
- `openpyxl`
- `python-docx`
- Microsoft Word (to view the `.docx` output)

## License

This repository is made available under the MIT License—see [LICENSE](LICENSE) for details.

## Citation

If you use this workflow or code for your research or clinical project, please credit or link back to this repository.

## Contact

*For questions or contributions, please reach out via [GitHub Issues](https://github.com/YOUR-USERNAME/YOUR-REPO/issues) or submit a pull request.*

## Sample Output

**Sample Case Report Form:**

| Case No | [value] |
|---|---|
| Age | [value] |
| Gender | [value] |
| Diagnosis | [value] |
| Relapse | [value] |
| Date of sample collection | [dd-mm-yyyy] |
