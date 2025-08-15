#!/bin/bash
# Step 8: BaseRecalibrator – Build recalibration model

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}_RG.bam"
REF="path/to/reference.fa"
DBSNP="path/to/dbsnp.vcf.gz"
OUTDIR="path/to/output/bqsr"
LOGDIR="path/to/logs"
THREADS=6
GATK_JAR="path/to/gatk.jar"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J baserecal \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/baserecal.out \
     -e ${LOGDIR}/baserecal.err \
"java -Xmx25G -jar $GATK_JAR BaseRecalibrator \
  -R $REF \
  -I $IN_BAM \
  --known-sites $DBSNP \
  -O $OUTDIR/${SAMPLE}_recal_data.table"
