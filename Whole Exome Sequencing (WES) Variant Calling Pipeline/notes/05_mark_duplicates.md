# Step 5 – MarkDuplicates

**Purpose:**  
Identify and mark/remove duplicate reads originating from PCR amplification.

**Inputs:**  
- Sorted BAM file from Step 4

**Outputs:**  
- Duplicate-marked BAM file
- Duplicate metrics text file
- BAM index file

**Why it's important:**  
Duplicate reads can bias variant frequency calculations and lead to false positives.

**Next Step:**  
If following RNA-seq style workflow (as per provided PDF), split reads with **SplitNCigarReads** (Step 6).  
If doing strict WES, you may skip Step 6 (I have done Step 6 for my WES work study).
