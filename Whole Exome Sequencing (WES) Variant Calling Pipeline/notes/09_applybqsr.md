# Step 9 – ApplyBQSR

**Purpose:**  
Apply base quality score recalibration to BAM files using the model from Step 8.

**Inputs:**  
- BAM file from Step 7
- Recalibration table from Step 8
- Reference genome FASTA

**Outputs:**  
- Recalibrated BAM file

**Why it's important:**  
Improves variant calling accuracy by correcting systematic sequencing errors.

**Next Step:**  
Perform variant calling with **HaplotypeCaller** (Step 10).
