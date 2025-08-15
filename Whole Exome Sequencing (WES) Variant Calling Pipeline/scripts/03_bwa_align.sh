#!/bin/bash
# Step 3: BWA – Align reads to reference genome

# ==============================
SAMPLE="SAMPLE_ID"
R1="path/to/${SAMPLE}_R1_paired.fastq.gz"
R2="path/to/${SAMPLE}_R2_paired.fastq.gz"
REF="path/to/reference.fa"
OUTDIR="path/to/output/bwa"
LOGDIR="path/to/logs"
THREADS=6
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J bwa_align \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 50000000 \
     -o ${LOGDIR}/bwa_align.out \
     -e ${LOGDIR}/bwa_align.err \
"module load bwa samtools && \
bwa mem -t $THREADS $REF $R1 $R2 > $OUTDIR/${SAMPLE}.sam && \
samtools view -bS $OUTDIR/${SAMPLE}.sam > $OUTDIR/${SAMPLE}.bam"
