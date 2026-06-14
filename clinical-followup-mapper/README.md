# Clinical Follow-Up Mapper

A Streamlit-based application for harmonizing longitudinal patient follow-up records into configurable follow-up intervals while preserving Excel formatting and generating quality-control reports.

---

## Overview

Clinical studies frequently involve repeated patient follow-up visits recorded over extended periods. Organizing these follow-up dates into predefined time intervals is often performed manually, making the process time-consuming and prone to inconsistencies.

Clinical Follow-Up Mapper automates this workflow by:

* Detecting baseline and follow-up dates from Excel files.
* Calculating follow-up durations relative to baseline.
* Grouping follow-up visits into configurable month intervals.
* Preserving original spreadsheet formatting.
* Generating quality-control reports for invalid or inconsistent dates.
* Producing a clean Excel output suitable for downstream analysis.

---

## Features

* Automatic baseline and follow-up date detection
* Robust date validation and error handling
* Dynamic follow-up interval generation
* Excel formatting preservation
* Quality-control reporting
* User-friendly Streamlit web interface
* Downloadable processed Excel output
* No patient records are discarded due to date issues

---

## Repository Structure

```text
Clinical_FollowUp_Mapper/
│
├── app.py
├── README.md
├── Documentation.md
├── requirements.txt
├── LICENSE
│
└── Examples/
    ├── sample_input_upload.xlsx
    └── sample_output_upload.xlsx
```

---

## Installation

Clone the repository:

```bash
git clone https://github.com/Jeyarish-007/Bioinformatics_Projects.git
```

Navigate to the project folder:

```bash
cd Clinical_FollowUp_Mapper
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Usage

Run the application:

```bash
streamlit run app.py
```

Open the provided local URL in your browser.

### Workflow

1. Upload the Excel file.
2. Review generated follow-up intervals.
3. Inspect date validation issues.
4. Download the processed Excel file.

---

## Requirements

* Python 3.8+
* Streamlit
* Pandas
* OpenPyXL
* python-dateutil

---

## Applications

* Clinical oncology studies
* Prospective cohort studies
* Registry management
* Longitudinal biomarker studies
* Biobanking projects
* Follow-up data harmonization

---

## Disclaimer

This software is intended for research and data-management purposes only.

Users are responsible for verifying outputs before clinical, regulatory, or patient-care use.

---

## Citation

If you use this project in your work, please cite:

Jeyarish. Clinical Follow-Up Mapper. GitHub Repository.

---

## License

This project is licensed under the MIT License.
