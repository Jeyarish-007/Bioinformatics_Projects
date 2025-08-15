#!/bin/bash
# Step 10: HaplotypeCaller – Variant calling

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}_recalibrated.bam"
REF="path/to/reference.fa"
DBSNP="path/to/dbsnp.vcf.gz"
OUTDIR="path/to/output/haplotypecaller"
LOGDIR="path/to/logs"
THREADS=6
GATK_JAR="path/to/gatk.jar"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J haplotypecaller \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/haplotypecaller.out \
     -e ${LOGDIR}/haplotypecaller.err \
"java -Xmx25G -jar $GATK_JAR HaplotypeCaller \
  -R $REF \
  -I $IN_BAM \
  -O $OUTDIR/${SAMPLE}_variants.vcf \
  --dbsnp $DBSNP \
  --dont-use-soft-clipped-bases \
  --stand-call-conf 20.0 \
  --create-output-variant-index true \
  --annotation Coverage \
  --annotation-group StandardHCAnnotation"
