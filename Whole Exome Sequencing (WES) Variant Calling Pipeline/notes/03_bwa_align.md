# Step 3 – BWA Alignment

**Purpose:**  
Align trimmed reads to the reference genome using the BWA-MEM algorithm.

**Inputs:**  
- Paired clean reads from Step 2
- Reference genome FASTA file

**Outputs:**  
- SAM file (raw alignment)
- BAM file (binary alignment)

**Why it's important:**  
Alignment maps sequencing reads to their genomic coordinates, enabling variant calling.

**Next Step:**  
Sort BAM file by genomic coordinates with **samtools sort** (Step 4).
