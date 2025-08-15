#!/bin/bash
# Step 2: Trimmomatic – Adapter and quality trimming

# ==============================
SAMPLE="SAMPLE_ID"
R1="path/to/${SAMPLE}_R1.fastq.gz"
R2="path/to/${SAMPLE}_R2.fastq.gz"
OUTDIR="path/to/output/trimmomatic"
LOGDIR="path/to/logs"
THREADS=6
TRIMMOMATIC_JAR="path/to/trimmomatic.jar"
ADAPTERS="path/to/TruSeq3-PE-2.fa"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J trimmomatic \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/trimmomatic.out \
     -e ${LOGDIR}/trimmomatic.err \
"module load java && \
java -jar $TRIMMOMATIC_JAR PE -threads $THREADS -phred33 \
  $R1 $R2 \
  $OUTDIR/${SAMPLE}_R1_paired.fastq.gz $OUTDIR/${SAMPLE}_R1_unpaired.fastq.gz \
  $OUTDIR/${SAMPLE}_R2_paired.fastq.gz $OUTDIR/${SAMPLE}_R2_unpaired.fastq.gz \
  ILLUMINACLIP:${ADAPTERS}:2:30:10 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:20:20 HEADCROP:15 MINLEN:100"
