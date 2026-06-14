# Clinical Follow-Up Mapper

---

## 1. Aim

The Clinical Follow-Up Mapper is designed to automate the organization of longitudinal patient follow-up records into predefined follow-up intervals. The application enables efficient management of clinical datasets while preserving original spreadsheet formatting and generating quality-control reports.

---

## 2. Objective

* Automatically identify baseline and follow-up dates from Excel spreadsheets.
* Calculate follow-up duration relative to baseline dates.
* Categorize follow-up visits into configurable month intervals.
* Detect and report invalid or inconsistent dates.
* Preserve original spreadsheet formatting and structure.
* Generate downloadable Excel reports suitable for clinical research workflows.

---

## 3. Methodology

### 3.1 Input Processing

The application accepts Microsoft Excel files (`.xlsx`, `.xls`) containing patient records.

Required fields include:

* Lab Id
* Date (Baseline)
* Case Number
* Patient Name
* Age
* Sex
* Nature of Sample

Additional follow-up date columns may be present.

---

### 3.2 Date Validation

Each date undergoes validation using Pandas datetime parsing.

Validation checks include:

* Missing dates
* Invalid date formats
* Impossible calendar values
* Out-of-range years
* Follow-up dates occurring before baseline dates

Invalid entries are recorded in a quality-control report.

---

### 3.3 Follow-Up Duration Calculation

For each patient:

[
\text{Follow-Up Duration} =
\text{Date}_{FollowUp}
----------------------

\text{Date}_{Baseline}
]

The duration is converted into months using:

[
\text{Months} =
(\text{Years} \times 12)
+
\text{Months}
+
\text{Partial Month Adjustment}
]

where partial months are rounded upward using a ceiling approach.

---

### 3.4 Follow-Up Interval Mapping

Follow-up visits are grouped into predefined month intervals.

Example:

| Interval   | Column Name          |
| ---------- | -------------------- |
| 1–2 Months | 1-2_months_followups |
| 3–4 Months | 3-4_months_followups |
| 5–6 Months | 5-6_months_followups |
| 7–8 Months | 7-8_months_followups |

Intervals are generated dynamically according to the maximum follow-up duration observed in the dataset.

---

### 3.5 Output Generation

For each patient:

* Baseline information is retained.
* Follow-up dates are assigned to interval-specific columns.
* Original formatting is preserved.
* Quality-control issues are reported separately.

Output is generated as a downloadable Excel workbook.

---

## 4. Workflow

```text
Input Excel File
        │
        ▼
Header Detection
        │
        ▼
Date Validation
        │
        ▼
Follow-Up Calculation
        │
        ▼
Month Interval Assignment
        │
        ▼
Quality Control Report
        │
        ▼
Excel Output Generation
```

---

## 5. Theory

### 5.1 Longitudinal Follow-Up Analysis

Longitudinal clinical studies collect observations from the same patient over multiple time points.

To facilitate analysis, follow-up visits are commonly grouped into predefined intervals.

Examples include:

* 0–3 months
* 3–6 months
* 6–12 months
* Annual follow-up visits

This harmonization improves consistency across datasets and simplifies downstream statistical analysis.

---

### 5.2 Quality-Control Principles

The application follows several quality-control principles:

* No patient records are discarded.
* Invalid dates are flagged rather than removed.
* Data integrity is maintained.
* Traceability is preserved through issue reporting.

---

## 6. Features

* Streamlit-based web interface
* Automated date validation
* Dynamic interval generation
* Preservation of Excel formatting
* Quality-control reporting
* Downloadable processed workbook
* Configurable follow-up intervals
* Robust error handling

---

## 7. Applications

* Clinical oncology studies
* Cohort studies
* Registry management
* Longitudinal biomarker studies
* Translational research
* Follow-up data harmonization
* Clinical trial data organization

---

## 8. Results

The application produces:

1. A harmonized Excel dataset.
2. Dynamically generated follow-up interval columns.
3. A quality-control report identifying problematic records.
4. Preserved spreadsheet formatting for user familiarity.

---

## 9. Conclusion

Clinical Follow-Up Mapper simplifies the organization and quality control of longitudinal patient follow-up datasets. By automating follow-up interval assignment while preserving spreadsheet formatting and reporting date inconsistencies, the application reduces manual effort and improves reproducibility in clinical research workflows.

---

## 10. Summary

| Aspect                | Description                                   |
| --------------------- | --------------------------------------------- |
| Application Type      | Streamlit Web Application                     |
| Programming Language  | Python                                        |
| Input                 | Excel Files (.xlsx, .xls)                     |
| Output                | Harmonized Excel Workbook                     |
| Follow-Up Calculation | Relative to Baseline Date                     |
| Error Handling        | Invalid Date Detection                        |
| Formatting            | Preserved                                     |
| Applications          | Clinical Research, Registries, Cohort Studies |

---

## 11. Requirements

### Software

* Python 3.8+
* Streamlit
* Pandas
* OpenPyXL
* python-dateutil

### Installation

```bash
pip install -r requirements.txt
```

---

## 12. References

* McKinney, W. (2010). Data Structures for Statistical Computing in Python. Proceedings of the 9th Python in Science Conference.
* OpenPyXL Documentation. https://openpyxl.readthedocs.io
* Streamlit Documentation. https://docs.streamlit.io
* Python DateUtil Documentation. https://dateutil.readthedocs.io

---

## 13. Final Notes

* The application is intended for research and data-management workflows.
* Users should review quality-control reports before downstream analysis.
* The software does not replace formal clinical data review procedures.
* Validation should be performed before regulatory or clinical use.
