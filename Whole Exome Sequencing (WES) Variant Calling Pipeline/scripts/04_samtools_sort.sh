#!/bin/bash
# Step 4: Sort BAM file

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}.bam"
OUTDIR="path/to/output/sorted_bam"
LOGDIR="path/to/logs"
THREADS=4
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J sort_bam \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/sort_bam.out \
     -e ${LOGDIR}/sort_bam.err \
"module load samtools && \
samtools sort $IN_BAM -o $OUTDIR/${SAMPLE}_sorted.bam"
