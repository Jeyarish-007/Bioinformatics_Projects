#!/bin/bash
# Step 7: AddOrReplaceReadGroups – Add read group info

# ==============================
SAMPLE="SAMPLE_ID"
IN_BAM="path/to/${SAMPLE}_SplitNCigar.bam"
OUTDIR="path/to/output/addrg"
LOGDIR="path/to/logs"
THREADS=2
PICARD_JAR="path/to/picard.jar"
# ==============================

mkdir -p "$OUTDIR" "$LOGDIR"

bsub -J addrg \
     -n $THREADS \
     -R "rusage[mem=10000]" \
     -M 20000000 \
     -o ${LOGDIR}/addrg.out \
     -e ${LOGDIR}/addrg.err \
"java -Xmx15G -jar $PICARD_JAR AddOrReplaceReadGroups \
  I=$IN_BAM \
  O=$OUTDIR/${SAMPLE}_RG.bam \
  RGID=group1 \
  RGLB=lib1 \
  RGPL=illumina \
  RGPU=unit1 \
  RGSM=$SAMPLE \
  CREATE_INDEX=true"
