#!/bin/bash
# Step 6: SplitNCigarReads – Split reads with Ns in CIGAR string

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}_dedup.bam"
REF="path/to/reference.fa"
OUTDIR="path/to/output/splitncigar"
LOGDIR="path/to/logs"
THREADS=4
GATK_JAR="path/to/gatk.jar"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J splitncigar \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 30000000 \
     -o ${LOGDIR}/splitncigar.out \
     -e ${LOGDIR}/splitncigar.err \
"java -Xmx25G -jar $GATK_JAR SplitNCigarReads \
  -R $REF \
  -I $IN_BAM \
  -O $OUTDIR/${SAMPLE}_SplitNCigar.bam"
