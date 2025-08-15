# Step 1 – FASTQC

**Purpose:**  
Assess the quality of raw sequencing reads to identify potential issues before downstream processing.

**Inputs:**  
- Raw paired-end FASTQ files (`_R1_` and `_R2_`)

**Outputs:**  
- HTML report  
- ZIP summary report

**Key Quality Metrics Checked:**  
- Per-base sequence quality
- GC content
- Sequence length distribution
- Adapter contamination

**Why it's important:**  
Detecting quality issues early ensures that trimming and filtering steps are targeted, improving downstream alignment and variant calling.

**Next Step:**  
Run **Trimmomatic** (Step 2) to remove adapters and low-quality regions if needed.
