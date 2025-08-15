#!/bin/bash
# Step 5: MarkDuplicates – Remove PCR duplicates

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}_sorted.bam"
OUTDIR="path/to/output/markdup"
LOGDIR="path/to/logs"
THREADS=4
GATK_JAR="path/to/gatk.jar"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J markdup \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/markdup.out \
     -e ${LOGDIR}/markdup.err \
"java -Xmx25G -jar $GATK_JAR MarkDuplicates \
  --INPUT $IN_BAM \
  --OUTPUT $OUTDIR/${SAMPLE}_dedup.bam \
  --METRICS_FILE $OUTDIR/${SAMPLE}_dedup_metrics.txt \
  --REMOVE_DUPLICATES true \
  --CREATE_INDEX true"
