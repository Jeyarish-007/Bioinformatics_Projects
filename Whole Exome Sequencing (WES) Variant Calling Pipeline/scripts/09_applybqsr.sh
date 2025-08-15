#!/bin/bash
# Step 9: ApplyBQSR – Apply recalibration

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}_RG.bam"
REF="path/to/reference.fa"
BQSR_TABLE="path/to/${SAMPLE}_recal_data.table"
OUTDIR="path/to/output/recalibrated_bam"
LOGDIR="path/to/logs"
THREADS=6
GATK_JAR="path/to/gatk.jar"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J applybqsr \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/applybqsr.out \
     -e ${LOGDIR}/applybqsr.err \
"java -Xmx25G -jar $GATK_JAR ApplyBQSR \
  -R $REF \
  -I $IN_BAM \
  --bqsr-recal-file $BQSR_TABLE \
  -O $OUTDIR/${SAMPLE}_recalibrated.bam"
