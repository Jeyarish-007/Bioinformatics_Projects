# Step 2 – Trimmomatic

**Purpose:**  
Remove adapter sequences, trim low-quality bases, and discard very short reads.

**Inputs:**  
- Raw paired-end FASTQ files from Step 1

**Outputs:**  
- Paired clean reads (`*_paired.fastq.gz`)
- Unpaired reads (`*_unpaired.fastq.gz`)

**Main Parameters Used:**  
- **ILLUMINACLIP**: Removes adapters from both ends
- **LEADING**: Trims low-quality bases from the start
- **TRAILING**: Trims low-quality bases from the end
- **SLIDINGWINDOW**: Cuts reads when average quality in a sliding window falls below threshold
- **HEADCROP**: Removes first N bases from each read
- **MINLEN**: Discards reads shorter than threshold

**Why it's important:**  
Trimming improves alignment accuracy and prevents false variant calls due to poor-quality bases.

**Next Step:**  
Align cleaned reads to the reference genome using **BWA** (Step 3).
