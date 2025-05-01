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

## Dependencies

- Python 3, Biopython
- NGLView, py3Dmol
- Jupyter Notebook
