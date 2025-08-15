# Step 4 – samtools sort

**Purpose:**  
Sort aligned reads by genomic coordinates.

**Inputs:**  
- BAM file from Step 3

**Outputs:**  
- Sorted BAM file

**Why it's important:**  
Many downstream tools (including GATK) require coordinate-sorted BAM files to process data efficiently.

**Next Step:**  
Mark PCR duplicates using **MarkDuplicates** (Step 5).
