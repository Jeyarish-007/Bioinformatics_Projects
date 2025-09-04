# Bioinformatics_Projects
Various projects of bioinformatics

## Projects:

[Project 1: FASTA Conversion tool](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/FASTA_conversion.ipynb)

[Project 2: Simple SNP Finder](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Simple_SNP_finder.ipynb)

[Project 3: Date Palm SSR-SNP Detection](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Date_Palm_SSR-SNP_Detection/SSP_Detection)

[Project 4: Mutation Visualization in Aligned DNA Sequences](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Mutation_Visualization_in_Aligned_DNA_Seq)

[Project 5: Smith-Waterman Algorithm Implementation](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Smith-Waterman_Algorithm)

[Project 6: VCF File Processing for Variant Extraction and Filtering](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/VCF%20File%20Processing%20for%20Variant%20Extraction%20and%20Filtering)

[Project 7: DNA Sequence Analysis](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/DNA_Seq_Analysis.ipynb)

[Project 8: RNA Transcription](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/RNA_transcription.ipynb)

[Project 9: GC Calculator](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/GC_Calculator.ipynb)

[Project 10: Clustering Single-Cell RNA Sequencing Data Using GMM](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Clustering_Sc-RNA_Sequencing_Data_Using_GMM)

[Project 11: Codon Usage Analysis](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Codon_Usage_Analysis_of_single_DNA_Sequence.ipynb)

[Project 12: DNA Translation Tool](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Simple_DNA_translation.ipynb)

[Project 13: Simple Primer Design Tool](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Simple_Primer_Designing.ipynb)

[Project 14: Protein Structure Visualizer](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/protein_str_visualization.ipynb)

[Project 15: Flow Cytometry Data Analysis Using Non-Negative Matrix Factorization (NMF)](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Flow%20Cytometry%20Data%20Analysis%20Using%20Non-Negative%20Matrix%20Factorization%20(NMF))

[Project 16: Gene Analysis Using Upset Plot: Techniques, Metastasis Frequency and Sample Size](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Data_Visualization_Techniques/Gene_Analysis_Using_Upset_Plot)

[Project 17: Automated Generation of Case Record Forms for Retrospective Clinical Data Analysis](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Automated%20Generation%20of%20Case%20Record%20Forms%20for%20Retrospective%20Clinical%20Data%20Analysis)

[Project 18: Whole Exome Sequencing (WES) Variant Calling Pipeline](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Whole%20Exome%20Sequencing%20(WES)%20Variant%20Calling%20Pipeline)

[Project 19: PubChem CID Extractor & SDF Downloader](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/pubchem-cid-extractor)

# 1. FASTA Conversion Tool  

[Project 1: FASTA Conversion tool](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/FASTA_conversion.ipynb)

A utility for converting raw DNA sequences to the standardized FASTA format.  
  
## Overview  
  
This tool provides a simple interface to convert DNA sequence data into the widely-used FASTA format, which is required by many bioinformatics analysis tools and databases.  
  
## Features  
  
- Converts raw DNA sequences to FASTA format  
- Allows custom sequence ID and description  
- Creates properly formatted FASTA files for downstream analysis  
- User-friendly interface via Jupyter Notebook  
  
## Usage  
  
1. Run the Jupyter notebook `FASTA_conversion.ipynb`  
2. Enter your DNA sequence when prompted  
3. Provide a sequence ID and description  
4. Specify an output filename  
5. The tool will generate a properly formatted FASTA file  
  
## Implementation Details  
  
The tool uses Biopython's SeqIO module to handle the conversion:  
1. Creates a Seq object from the input DNA sequence  
2. Defines a SeqRecord with the sequence, ID, and description  
3. Writes the SeqRecord to a file in FASTA format  
  
## Dependencies  
  
- Python 3.x  
- Biopython  
- Jupyter Notebook

# 2. Simple SNP Finder  

[Project 2: Simple SNP Finder](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Simple_SNP_finder.ipynb)

A tool for identifying Single Nucleotide Polymorphisms (SNPs) between two DNA sequences.  
  
## Overview  
  
This project provides a straightforward method to detect SNPs, which are single base pair variations between DNA sequences. SNPs are important genetic markers used in various applications including population genetics, disease association studies, and evolutionary biology.  
  
## Features  
  
- Compares two DNA sequences to identify single nucleotide differences  
- Displays the position and nucleotide changes in a tabular format  
- Works with FASTA format input files  
- Simple and efficient implementation  
  
## Usage  
  
1. Run the Jupyter notebook `Simple_SNP_finder.ipynb`  
2. Specify the paths to your two FASTA sequence files  
3. The tool will identify all SNPs and display them in a table showing:  
   - Position of the SNP  
   - Nucleotide in sequence 1  
   - Nucleotide in sequence 2  
  
## Implementation Details  
  
The tool implements a straightforward algorithm:  
1. Reads two DNA sequences from FASTA files using Biopython  
2. Compares each position between the sequences  
3. Records positions where nucleotides differ  
4. Outputs the results in a readable tabular format  
  
## Dependencies  
  
- Python 3.x  
- Biopython  
- Jupyter Notebook

# 3. Date Palm SSR-SNP Detection  

[Project 3: Date Palm SSR-SNP Detection](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Date_Palm_SSR-SNP_Detection/SSP_Detection)
  
A project for detecting Simple Sequence Repeats (SSRs) and Single Nucleotide Polymorphisms (SNPs) in aquaporin genes of Phoenix dactylifera (Date Palm).  
  
## Overview  
  
This project focuses on the in silico detection of genetic markers (SSRs and SNPs) in the aquaporin gene family of date palm. These markers are valuable for genetic studies, breeding programs, and understanding genetic diversity in date palm varieties.  
  
## Objectives  
  
- Retrieve DNA sequences of aquaporin genes from the NCBI database  
- Detect and analyze SSRs in the retrieved sequences using the MISA tool  
- Detect and analyze SNPs in the retrieved sequences  
- Document and share findings  
  
## Methodology  
  
### Sequence Retrieval  
- DNA sequences of aquaporin genes were retrieved from the NCBI Gene database  
- Sequences were downloaded in FASTA format using Biopython's Entrez module  
- 43 aquaporin gene sequences were collected from various chromosomes  
  
### SSR Detection  
- The MISA (Microsatellite Identification Tool) was used to detect SSRs  
- Parameters used:  
  - Minimum repeat number: Mononucleotide (10), Dinucleotide (6), Trinucleotide (5), Tetranucleotide (5), Pentanucleotide (5), Hexanucleotide (5)  
  - Maximum distance between SSRs for compound SSRs: 100 bases  
  
## Results  
  
- A total of 43 aquaporin gene sequences were analyzed  
- 120 SSRs were detected across the sequences  
- The most common SSR types were:  
  - Mono-nucleotide repeats: 50  
  - Di-nucleotide repeats: 40  
  - Tri-nucleotide repeats: 20  
  
## Files  
  
- `Aquaporins_Date_palm_(Input_file).fasta`: Input FASTA file containing the aquaporin gene sequences  
- `aquaporin_sequences.fasta.misa`: Output file from MISA containing SSR regions  
- `aquaporin_sequences.fasta.statistics`: Summary statistics file from MISA  
- `aquaporin_fasta_seq_code.ipynb`: Script to download sequences from NCBI  
  
## Future Work  
  
- Detection of SNPs in the aquaporin gene sequences using tools like Clustal Omega or VarScan

# 4. Mutation Visualization in Aligned DNA Sequences  

[Project 4: Mutation Visualization in Aligned DNA Sequences](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Mutation_Visualization_in_Aligned_DNA_Seq)
  
A tool for identifying and visualizing mutations between DNA sequences.  
  
## Overview  
  
This project provides a Python-based solution for detecting mutations between two DNA sequences and visualizing the differences using both a bar chart and sequence alignment representation.  
  
## Objectives  
  
- Identify point mutations between two given DNA sequences  
- Visualize the mutation locations using a bar chart  
- Highlight mutations in a sequence alignment format  
  
## Methodology  
  
### Data Input  
- Two DNA sequences are defined using `Bio.Seq` from Biopython  
  
### Mutation Identification  
- The tool iterates through both sequences  
- Compares each nucleotide at corresponding positions  
- Stores the index and differing nucleotides in a list  
  
### Visualization  
1. **Bar Chart**:  
   - Uses Matplotlib to plot mutation positions on the x-axis  
   - Labels each mutation with its original and mutated nucleotide  
  
2. **Sequence Alignment**:  
   - Aligns both sequences for comparison  
   - Highlights mutations in red  
   - Displays a `^` symbol under mismatched bases  
  
## Significance  
  
- Helps in understanding genetic variations between two DNA sequences  
- Useful for bioinformatics applications, including mutation analysis in genetic research  
- Provides a clear and visual representation of mutation positions  
  
## Dependencies  
  
- Python 3.x  
- Biopython  
- Matplotlib  
- Jupyter Notebook

# 5. Smith-Waterman Algorithm Implementation  

[Project 5: Smith-Waterman Algorithm Implementation](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Smith-Waterman_Algorithm)

A Python implementation of the Smith-Waterman algorithm for local sequence alignment.  
  
## Overview  
  
The Smith-Waterman algorithm is designed to find the optimal local alignment between two biological sequences (DNA, RNA, or proteins). Unlike global alignment (Needleman-Wunsch), it identifies highly similar subsequences even if the overall sequences are divergent.  
  
## Objectives  
  
- Compute the best local alignment between two sequences  
- Identify conserved regions (e.g., protein domains or gene segments)  
- Provide a scoring system for matches, mismatches, and gaps  
- Serve as a foundation for understanding more complex alignment tools  
  
## Algorithm Details  
  
### Steps  
1. **Initialization**:  
   - Create a scoring matrix of size (m+1) × (n+1) (where m and n are sequence lengths)  
   - Initialize the first row and column to 0 (local alignment allows alignment to start anywhere)  
  
2. **Matrix Filling**:  
   - For each cell, calculate the maximum score from:  
     - 0 (to ensure local alignment)  
     - Diagonal cell + match/mismatch score  
     - Upper cell + gap penalty  
     - Left cell + gap penalty  
  
3. **Traceback**:  
   - Start from the highest-scoring cell  
   - Move diagonally (match/mismatch), up (gap in sequence 2), or left (gap in sequence 1) until reaching 0  
  
### Scoring System  
- Match: +2 (default)  
- Mismatch: -1 (default)  
- Gap penalty: -1 (default)  
  
## Implementation  
  
The algorithm is implemented in Python with the following features:  
- Dynamic programming approach  
- Customizable scoring parameters  
- Traceback for alignment visualization  
- Example usage with sample sequences  
  
## Files  
- `Smith-waterman_algorithm.ipynb`: Jupyter notebook with the implementation and example  
  
## Usage  
  
```python  
# Example usage  [header-7](#header-7)
seq1 = "ACCTA"  
seq2 = "AGCTA"  
      
aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2)  
      
print(f"Sequence 1: {aligned_seq1}")  
print(f"Sequence 2: {aligned_seq2}")  
print(f"Alignment Score: {score}")
```

# 6. VCF File Processing for Variant Extraction and Filtering  

[Project 6: VCF File Processing for Variant Extraction and Filtering](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/VCF%20File%20Processing%20for%20Variant%20Extraction%20and%20Filtering)
  
A tool for manipulating Variant Call Format (VCF) files to extract and filter genetic variants.  
  
## Overview  
  
This project provides utilities for processing VCF files, which are standard formats for storing gene sequence variations. The tool allows for extraction of specific information and filtering of variants based on criteria like depth (DP).  
  
## Features  
  
- Extracts specific information from compressed VCF files (.vcf.gz)  
- Filters variants based on depth (DP) and other criteria  
- Converts data to tab-separated values (TSV) format  
- Handles headers and variant lines separately  
  
## Usage  
  
1. Run the Jupyter notebook `Manipulating_VCF_file.ipynb`  
2. Specify the input VCF.gz file path  
3. Choose extraction or filtering options  
4. The tool will process the file and output results in TSV format  
  
## Implementation Details  
  
The tool uses Python's gzip module to handle compressed VCF files and implements:  
1. Parsing of VCF headers and variant lines  
2. Extraction of CHROM, POS, REF, and ALT fields  
3. Filtering variants based on the DP (depth) field  
4. Output of processed data in TSV format  
  
## Applications  
  
- Quality control of variant calls  
- Preparation of variant data for downstream analysis  
- Extraction of specific variants of interest  
- Conversion between file formats for compatibility with other tools  
  
## Dependencies  
  
- Python 3.x  
- gzip module  
- Jupyter Notebook

# 7. DNA Sequence Analysis 

[Project 7: DNA Sequence Analysis](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/DNA_Seq_Analysis.ipynb)
  
A tool for analyzing nucleotide composition and patterns in DNA sequences.  
  
## Overview  
  
This project provides utilities for analyzing DNA sequences, including nucleotide frequency calculation, GC content analysis, and pattern identification. These analyses are fundamental in genomics research and can provide insights into the functional and evolutionary characteristics of DNA.  
  
## Features  
  
- Calculates nucleotide frequencies (A, T, G, C)  
- Determines GC content percentage  
- Identifies sequence patterns and motifs  
- Provides statistical summaries of DNA composition  
  
## Usage  
  
1. Run the Jupyter notebook `DNA_Seq_Analysis.ipynb`  
2. Input your DNA sequence  
3. The tool will perform various analyses and display the results  
  
## Implementation Details  
  
The tool implements several analysis methods:  
1. Counting occurrences of each nucleotide  
2. Calculating GC content as a percentage  
3. Searching for specific sequence patterns  
4. Generating statistical summaries of sequence composition  
  
## Applications  
  
- Quality assessment of DNA sequences  
- Comparative genomics  
- Primer design for PCR  
- Identification of functional regions in genomes  
  
## Dependencies  
  
- Python 3.x  
- Biopython (optional)  
- Jupyter Notebook

# 8. RNA Transcription  

[Project 8: RNA Transcription](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/RNA_transcription.ipynb)
  
A tool for converting DNA sequences to their corresponding RNA sequences.  
  
## Overview  
  
This project provides a simple utility for transcribing DNA sequences to RNA sequences by replacing thymine (T) with uracil (U) and maintaining the other nucleotides. This process mimics the biological transcription of DNA to RNA.  
  
## Features  
  
- Converts DNA sequences to RNA sequences  
- Handles both standard and non-standard nucleotide codes  
- Provides options for handling ambiguous bases  
- Simple and intuitive interface  
  
## Usage  
  
1. Run the Jupyter notebook `RNA_transcription.ipynb`  
2. Input your DNA sequence  
3. The tool will transcribe the sequence and display the resulting RNA  
  
## Implementation Details  
  
The tool implements the transcription process by:  
1. Reading the input DNA sequence  
2. Replacing all instances of 'T' with 'U'  
3. Maintaining all other nucleotides (A, G, C)  
4. Returning the transcribed RNA sequence  
  
## Applications  
  
- Educational demonstrations of transcription  
- Preparation of sequences for RNA analysis  
- Part of larger genomic analysis pipelines  
- Conversion of DNA database entries to RNA format  
  
## Dependencies  
  
- Python 3.x  
- Biopython (optional)  
- Jupyter Notebook

# 9. GC Calculator

[Project 9: GC Calculator](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/GC_Calculator.ipynb)
  
A simple tool to calculate the GC content in DNA sequences.  
  
## Overview  
  
GC Calculator is a Python-based tool that calculates the percentage of guanine (G) and cytosine (C) nucleotides in a given DNA sequence. The GC content is an important metric in genomics and can provide insights into various biological properties of DNA.  
  
## Features  
  
- Calculates the exact percentage of GC content in any DNA sequence  
- Simple and user-friendly interface  
- Fast computation even for long sequences  
- Provides both raw count and percentage values  
  
## Usage  
  
1. Run the Jupyter notebook `GC_Calculator.ipynb`  
2. Enter your DNA sequence when prompted  
3. The tool will calculate and display the GC content percentage  
  
## Implementation Details  
  
The tool uses a straightforward algorithm:  
1. Count the number of G and C nucleotides in the sequence  
2. Divide by the total sequence length  
3. Multiply by 100 to get the percentage  
  
## Applications  
  
- Genome characterization  
- Primer design for PCR  
- Prediction of DNA stability and melting temperature  
- Taxonomic classification of organisms  
  
## Dependencies  
  
- Python 3.x  
- Jupyter Notebook

# 10. Clustering Single-Cell RNA Sequencing Data Using GMM

[Project 10: Clustering Single-Cell RNA Sequencing Data Using GMM](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Clustering_Sc-RNA_Sequencing_Data_Using_GMM)

## Overview 

Implementation of Gaussian Mixture Models for clustering scRNA-seq data.

## Objective

Cluster and annotate scRNA-seq data to identify cell types based on gene expression.

## Features

- Data preprocessing (normalization, log transformation)
- Highly variable gene selection
- Dimensionality reduction (PCA)
- GMM clustering
- UMAP visualization
- Marker gene identification
- Cell type annotation

## Implementation Workflow

1. Load PBMC dataset
2. Preprocess with normalization
3. Select highly variable genes
4. Reduce dimensions with PCA
5. Apply GMM clustering
6. Visualize with UMAP
7. Identify marker genes
8. Annotate cell types

## Applications

- Single-cell RNA sequencing analysis
- Cell type identification
- Biomarker discovery
- Cellular heterogeneity studies

## Dependencies

- Python 3, NumPy, ScanPy
- scikit-learn, UMAP, Matplotlib

# 11. Codon Usage Analysis

[Project 11: Codon Usage Analysis](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Codon_Usage_Analysis_of_single_DNA_Sequence.ipynb)

## Overview

Analyze codon frequency in DNA sequences.

##Features

- Codon frequency calculation
- Sequence validation
- Usage percentage calculation

##Implementation

- Validates sequence length
- Counts codons using Counter
- Calculates percentages

## Applications

- Codon bias analysis
- Gene expression optimization
- Synthetic gene design

## Dependencies

- Python 3

# 12. Simple DNA Translation Tool

[Project 12: DNA Translation Tool](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Simple_DNA_translation.ipynb)

## Overview

Translate DNA to protein.

## Features

- Complete genetic code
- Stop codon handling

## Applications

- Protein prediction
- Gene annotation

## Dependencies 

- Python 3

# 13. Primer Design Tool

[Project 13: Simple Primer Design Tool](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/Simple_Primer_Designing.ipynb)

## Overview

Design PCR primers.

## Features

- Melting temperature calculation
- GC content analysis
- Primer filtering

## Applications

- PCR design
- Diagnostic assays

## Dependencies

- Python 3

# 14. Protein Structure Visualizer

[Project 14: Protein Structure Visualizer](https://github.com/Jeyarish-007/Bioinformatics_Projects/blob/main/protein_str_visualization.ipynb)

## Overview

Visualize PDB files.

## Features

- Interactive 3D visualization
- Multiple display styles

## Applications

- Drug design
- Structural biology

# 15.  Flow Cytometry Data Analysis Using Non-Negative Matrix Factorization (NMF)

[Project 15: Flow Cytometry Data Analysis Using Non-Negative Matrix Factorization (NMF)](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Flow%20Cytometry%20Data%20Analysis%20Using%20Non-Negative%20Matrix%20Factorization%20(NMF))

## Overview

This project implements an unsupervised machine learning approach using Non-Negative Matrix Factorization (NMF) to identify biologically meaningful subpopulations in high-dimensional flow cytometry data. The workflow includes preprocessing, dimensionality reduction via PCA, clustering with NMF, and visualization of marker contributions and cell populations.

## Features

- FCS File Support : Load and preprocess raw flow cytometry data using FlowCal.
- Data Normalization : Apply min-max scaling for improved NMF performance.
- PCA-Based Component Selection : Determine optimal number of components by retaining 95% variance.
- NMF Clustering : Decompose data into interpretable basis (W) and coefficient (H) matrices.

## Visualization Tools :

- Heatmap of marker contributions per component
- t-SNE projection of clustered cells
- Biological Interpretation : Identify dominant markers and link them to known cell types or artifacts.
  
## Usage

Clone the repository:
```bash
git clone https://github.com/yourusername/flow-cytometry-nmf.git 
cd flow-cytometry-nmf
```
Install dependencies:
```bash
pip install -r requirements.txt
```

Run the Jupyter notebook:
```bash
jupyter notebook NMF_Flow_Cytometry_Analysis.ipynb
```

Follow the steps inside the notebook:
Load FCS file
Preprocess and normalize data
Perform PCA for component selection
Apply NMF clustering
Visualize results with heatmap and t-SNE
Interpret findings

## Implementation Details

### Data Loading & Preprocessing
- Use FlowCal.io.FCSData() to load .fcs files
- Convert to pandas DataFrame and apply min-max normalization
### Dimensionality Reduction
- Apply PCA to determine optimal number of NMF components (retain 95% variance)
### NMF Clustering
- Use sklearn.decomposition.NMF to decompose data into:
- W: Cell-cluster association matrix
- H: Marker contribution matrix
### Visualization
- Plot heatmap of H matrix to show marker contributions
- Use t-SNE for 2D visualization of clusters
### Interpretation
- Analyze dominant markers per component
- Relate clusters to biological cell types or experimental artifacts

## Software & Libraries
- Python 3.10+
- FlowCal – for reading FCS files
- scikit-learn – for PCA, NMF, and t-SNE
- pandas, numpy – for data manipulation
- matplotlib, seaborn – for visualization
- jupyter – for interactive analysis
- Recommended Hardware
- At least 8GB RAM (for datasets with ~500k cells)
- 16GB+ RAM recommended for faster processing

## Dependencies

- Python 3, Biopython
- NGLView, py3Dmol
- Jupyter Notebook

# 16) Gene Analysis Using Upset Plot: Techniques, Metastasis Frequency and Sample Size

[Project 16: Gene Analysis Using Upset Plot: Techniques, Metastasis Frequency and Sample Size](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Data_Visualization_Techniques/Gene_Analysis_Using_Upset_Plot)

To visualize and analyze relationships between gene sequencing techniques (NGS, WGS, WES) and their association with metastasis frequency in cancer-related genes using UpSet plots.

## Overview
This project provides a comprehensive analysis of:
- Distribution of sequencing techniques across 20 cancer-related genes
- Correlation between technique combinations and metastasis frequency
- Sample size distribution across different technique combinations
- Visualization of complex set relationships using UpSet plots

## Features
- **Interactive Visualization**: Combines intersection matrix, bar charts, and strip plots
- **Multi-dimensional Analysis**: Simultaneously visualizes technique combinations, metastasis frequency, and sample sizes
- **Customizable Output**: Adjustable plot parameters and styling
- **Supplementary Statistics**: Automated technique counting and gene-technique mapping

## Usage
```python
# Basic usage:
from upsetplot import UpSet
upset = UpSet(data, show_counts=True)
upset.add_catplot(value='Metric', kind='strip')
upset.plot()
```

## Dependencies

### Core Packages
| Package     | Minimum Version | Purpose                  |
|-------------|-----------------|--------------------------|
| pandas      | ≥1.0            | Data manipulation        |
| matplotlib  | ≥3.0            | Base plotting framework  |
| upsetplot   | ≥0.6            | UpSet visualization      |

### Optional Packages
| Package    | Recommended For                 |
|------------|----------------------------------|
| seaborn    | Enhanced plot styling           |
| jupyter    | Interactive notebook exploration|
| numpy      | Numerical operations            |
| scipy      | Advanced statistical functions  |

### Installation
```bash
# Install core packages
pip install pandas matplotlib upsetplot

# Install optional packages
pip install seaborn jupyter numpy scipy
```

# 17) Automated Generation of Case Record Forms for Retrospective Clinical Data Analysis

[Project 17: Automated Generation of Case Record Forms for Retrospective Clinical Data Analysis](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Automated%20Generation%20of%20Case%20Record%20Forms%20for%20Retrospective%20Clinical%20Data%20Analysis)


To streamline retrospective clinical research documentation by automating the generation of standardized Case Report Forms (CRFs) from patient data stored in Excel, using an open-source Python workflow.

## Overview

This project automates the process of transcribing retrospective patient data into formatted Case Report Forms (CRFs) for clinical analysis. By leveraging Python, `pandas`, and `python-docx`, users can efficiently convert a spreadsheet containing multiple patient records into a single, professional Word document, with each CRF neatly formatted, titled, and ready for analysis or archiving. The approach reduces manual effort, ensures uniformity, and can be easily adapted for similar projects in clinical or research settings.

## Features

- **Automated Document Generation:**  
  Converts an entire Excel sheet of patient records into a Word document with one CRF per page.
- **Customizable CRF Table:**  
  Each page includes a standardized, bordered table with defined clinical fields.
- **Centered and Styled Title:**  
  Every CRF page starts with a centered "Case Report Form" header.
- **Consistent Spacing:**  
  Five blank lines are inserted between the title and the CRF table, improving readability.
- **Date Formatting:**  
  Dates are shown in `dd-mm-yyyy` format, with no time component.
- **Error Handling:**  
  Handles missing or non-standard data gracefully.
- **Easy Adaptation:**  
  The script can be tailored for different fields or templates as needed.

## Usage

1. **Prepare your data:**  
   - Enter each patient record as a row in Excel with columns: Case No, Age, Gender, Diagnosis, Relapse, Date of sample collection.

2. **Clone or download this repository.**

3. **Install the required packages:**
    ```sh
    pip install pandas openpyxl python-docx
    ```

4. **Edit file paths** in the script to point to your Excel data and specify the desired output Word file location.

5. **Run the script:**
    ```sh
    Scripts.py
    ```

6. **Open the generated `.docx` file** to view a compiled set of CRFs—each patient’s data formatted on its own page.

## Implementation Details

- The script reads patient data from an Excel file using `pandas`.
- It iterates through each row (representing a patient), creating a new page in the Word document starting with a centered "Case Report Form" heading.
- Five blank lines are inserted for visual spacing.
- A 6x2 table is generated for each patient, populated with the appropriate fields.
- Table borders are added using low-level XML manipulation for a professional look.
- The `"Date of sample collection"` field is formatted as `dd-mm-yyyy`, omitting any time data.
- Each record starts on a new page (`add_page_break`).
- The final Word document contains one CRF per patient, compiled into a single file.

## Dependencies

- **Python 3.7+**
- [pandas](https://pandas.pydata.org/)  
- [openpyxl](https://openpyxl.readthedocs.io/en/stable/) (for Excel reading support)
- [python-docx](https://python-docx.readthedocs.io/en/latest/) (for Word document generation)

To install dependencies:

```sh
pip install pandas openpyxl python-docx
```

# 18) Whole Exome Sequencing (WES) Variant Calling Pipeline

[Project 18: Whole Exome Sequencing (WES) Variant Calling Pipeline](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/Whole%20Exome%20Sequencing%20(WES)%20Variant%20Calling%20Pipeline)

To provide a reproducible, modular, and step-by-step computational workflow for processing Whole Exome Sequencing (WES) data — from raw FASTQ files to biologically meaningful annotated variant results.

## Overview
This pipeline processes raw paired-end WES data through essential bioinformatics stages including quality control, trimming, alignment, duplicate removal, recalibration, variant calling, filtering, and annotation.  
Each step is implemented as an **individual script** to allow modular execution, tracking, and troubleshooting.

**Workflow Summary:**
1. Quality control of raw reads  
2. Adapter and quality trimming  
3. Alignment to reference genome  
4. Sorting alignments  
5. Duplicate marking  
6. BAM indexing  
7. Base quality score recalibration (BQSR)  
8. Variant calling  
9. Variant filtering  
10. Variant annotation  

## Features
- Modular 10-script design for flexibility
- Compatible with HPC (LSF job scheduling)
- Works with human WES datasets
- Uses well-established tools (FastQC, BWA, Picard, GATK, ANNOVAR/VEP)
- Step-by-step documentation with matching notes
- Easily extendable to other variant callers or annotation tools

## Usage
1. Ensure all required dependencies are installed and available in your environment.
2. Place your input FASTQ files in the designated input directory.
3. Modify the variables in each script to match your:
   - Input FASTQ file names
   - Reference genome path
   - Output directory path
4. Run each script in order:
5. Review outputs after each step before proceeding.

## Implementation Details

- **Data Input:** Paired-end FASTQ files from Illumina sequencing
- **Reference Genome:** hg38 or hg19 (FASTA + index files)
- **Aligner:** BWA-MEM for mapping reads  
- **Processing:** Picard tools for sorting and marking duplicates
- **Recalibration:** GATK BaseRecalibrator and ApplyBQSR  
- **Variant Calling:** GATK HaplotypeCaller in gVCF mode 
- **Filtering:** Hard filters based on variant quality metrics 
- **Annotation:** VEP or ANNOVAR for functional annotation   

## Dependencies

- **FastQC** – Quality control of raw reads  
- **fastp** / **Trimmomatic** – Adapter and quality trimming  
- **BWA** – Read alignment to reference genome  
- **SAMtools** – BAM indexing and basic operations  
- **Picard** – Sorting and marking duplicates  
- **GATK** – BQSR, variant calling, filtering
- **ANNOVAR** / **Ensembl VEP** – Variant annotation
- **Perl/Python** – Required for annotation tools

# 19) PubChem CID Extractor & SDF Downloader

[Project 19: PubChem CID Extractor & SDF Downloader](https://github.com/Jeyarish-007/Bioinformatics_Projects/tree/main/pubchem-cid-extractor)

To provide an automated, scalable, and user-friendly pipeline for extracting PubChem Compound IDs (CIDs) from compound names and downloading corresponding 3D SDF structure files for cheminformatics and bioinformatics research.

## Overview
The project integrates PubChem PUG REST API with Python scripts to enrich compound datasets.
It allows researchers to:
- Map compound names to PubChem CIDs.
- Download 3D molecular structures.
- Generate enriched datasets with hyperlinks and descriptive filenames.
This workflow is particularly suited for natural product research, molecular docking, QSAR modeling, and drug discovery pipelines.  

## Introduction

Chemical compound annotation and structure retrieval are critical in **computational chemistry** and ***bioinformatics**.
Manual retrieval of PubChem CIDs and 3D structures is time-consuming and prone to error.

This pipeline automates the process by:
- Querying PubChem with compound names.
- Extracting corresponding CIDs.
- Downloading 3D conformer files in SDF format.
- Creating enriched datasets with clickable PubChem links.
It ensures data integrity, reproducibility, and ease of integration with downstream analyses.

## Features
- Automated CID extraction from compound names.
- Bulk 3D SDF download with retry and logging for failed cases.
- Descriptive filenames (`CompoundName_IMPPATID_CID.sdf`) for traceability.
- Excel-friendly hyperlinks for direct PubChem access.
- Modular scripts for different use cases:
  - CID extraction
  - SDF downloads
  - Link-only datasets
- Handles `.xlsx` and `.csv` input formats.

## Usage
```bash
# Extract CIDs from compound names
python src/extract_cid.py --input data/example_input.xlsx --output data/imppat_enriched.csv

# Download SDF files from CIDs
python src/download_sdf.py --input data/imppat_enriched.csv --output outputs/sdf_files/

# Quick download from clean CIDs only
python src/quick_download.py --input data/example_input.xlsx --output outputs/sdf_files/

# Enriched dataset with filenames and links
python src/enriched_links.py --input data/example_input.xlsx --output data/with_links.csv

# Generate Excel with hyperlinks only (no downloads)
python src/hyperlink_generator.py --input data/example_input.xlsx --output data/links.xlsx
```

## Implementation Details

- **Language**: `Python 3.8+`

- **Workflow**:
  - Input dataset loaded with Pandas.
  - API calls made to PubChem PUG REST API via requests.
  - Results parsed, validated, and added to DataFrame.
  - Optional: SDF files downloaded in batch with retry mechanism.
  - Outputs saved as CSV/Excel and SDF structure files.
- **Error Handling**:
  - Invalid CIDs skipped.
  - Failed downloads logged in outputs/logs/failed_downloads.txt.
  - Duplicate downloads avoided by checking existing files.

## Dependencies

Core **`Python`** libraries required:
- **pandas** → for data handling (Excel/CSV).
- **requests** → for PubChem API calls.
- **openpyxl** → for Excel file support.
- **os, re, time** → standard library modules.
