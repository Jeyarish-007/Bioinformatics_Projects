# Step 8 – BaseRecalibrator

**Purpose:**  
Analyze patterns of covariation in the sequence dataset and produce a recalibration table.

**Inputs:**  
- Read-group annotated BAM file from Step 7
- Known sites of variation (dbSNP, Mills indels, etc.)
- Reference genome FASTA

**Outputs:**  
- Recalibration table file (`*.table`)

**Why it's important:**  
Sequencing machines may introduce systematic errors in base quality scores.  
Recalibration models these errors and corrects them in the next step.

**Next Step:**  
Apply the recalibration with **ApplyBQSR** (Step 9).
