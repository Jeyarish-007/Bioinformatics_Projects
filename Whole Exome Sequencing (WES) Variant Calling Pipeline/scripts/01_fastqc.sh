#!/bin/bash
# Step 1: FASTQC – Quality check of raw reads

# ==============================
# User-defined variables
# ==============================
SAMPLE="SAMPLE_ID"
R1="path/to/${SAMPLE}_R1.fastq.gz"
R2="path/to/${SAMPLE}_R2.fastq.gz"
OUTDIR="path/to/output/fastqc"
LOGDIR="path/to/logs"
THREADS=4
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J fastqc \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 20000000 \
     -o ${LOGDIR}/fastqc.out \
     -e ${LOGDIR}/fastqc.err \
"module load java && \
fastqc -t $THREADS $R1 $R2 -o $OUTDIR"
