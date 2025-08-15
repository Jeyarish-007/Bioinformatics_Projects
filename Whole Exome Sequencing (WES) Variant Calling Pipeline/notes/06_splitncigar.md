# Step 6 – SplitNCigarReads

**Purpose:**  
Split reads containing Ns in their CIGAR string into separate exon segments.

**Inputs:**  
- Duplicate-marked BAM file from Step 5

**Outputs:**  
- Split BAM file

**Why it's important:**  
This step is normally for RNA-seq data, but included here as per provided commands.  
It prepares the file for base quality score recalibration when dealing with spliced alignments.

**Next Step:**  
Add read group information with **AddOrReplaceReadGroups** (Step 7).
