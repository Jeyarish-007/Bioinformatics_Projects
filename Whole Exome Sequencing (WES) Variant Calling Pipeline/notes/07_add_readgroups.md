# Step 7 – AddOrReplaceReadGroups

**Purpose:**  
Add or modify read group information in BAM files.

**Inputs:**  
- Split BAM file from Step 6

**Outputs:**  
- BAM file with correct read group metadata
- BAM index file

**Why it's important:**  
Read groups identify the origin of reads (sample, library, platform) and are required for GATK processing.

**Next Step:**  
Perform base recalibration with **BaseRecalibrator** (Step 8).
