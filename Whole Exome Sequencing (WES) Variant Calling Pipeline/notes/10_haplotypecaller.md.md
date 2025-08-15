# Step 10 – HaplotypeCaller

**Purpose:**  
Call SNPs and indels using GATK’s HaplotypeCaller algorithm.

**Inputs:**  
- Recalibrated BAM file from Step 9
- Reference genome FASTA
- Known variant sites (optional)

**Outputs:**  
- VCF file containing variant calls
- VCF index file (`.tbi`)

**Why it's important:**  
This is the final step that produces the list of genetic variants for downstream analysis.

**Next Step:**  
Perform variant annotation (e.g., with VEP, ANNOVAR) and filtering as needed.
