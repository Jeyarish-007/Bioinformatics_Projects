# Whole Exome Sequencing (WES) Variant Calling Pipeline

*An automated pipeline for WES data preprocessing, alignment, and variant calling using GATK best practices.*

---

## Table of Contents
- [Project Overview](#-project-overview)
- [Aim](#aim)
- [Objectives](#objectives)
- [Methodology](#methodology)
  - [1. Quality Control (FastQC)](#1-quality-control-fastqc)
  - [2. Trimming (Trimmomatic)](#2-trimming-trimmomatic)
  - [3. Alignment (BWA + samtools)](#3-alignment-bwa--samtools)
  - [4. Mark Duplicates (GATK/Picard)](#4-mark-duplicates-gatkpicard)
  - [5. SplitNCigarReads (GATK)](#5-splitncigarreads-gatk)
  - [6. Add Read Groups (Picard)](#6-add-read-groups-picard)
  - [7. Base Quality Score Recalibration (BQSR)](#7-base-quality-score-recalibration-bqsr)
  - [8. Variant Calling (HaplotypeCaller)](#8-variant-calling-haplotypecaller)
- [Significance](#significance)
- [Results](#results)
- [Conclusion](#conclusion)
- [Usage](#usage)
- [Requirements](#requirements)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)

---

## Project Overview

This repository contains a **complete, reproducible pipeline** for processing Whole Exome Sequencing (WES) data to identify single nucleotide polymorphisms (SNPs) and insertions/deletions (indels). The pipeline follows the **GATK Best Practices** workflow and it has robustness and clarity.

The pipeline is designed for execution on an **HPC cluster using LSF (bsub)** and processes raw FASTQ files through quality control, trimming, alignment, duplicate marking, base recalibration, and final variant calling.

---

## Aim

To develop a **standardized, automated, and well-documented pipeline** for identifying genetic variants from WES data, enabling accurate downstream analysis such as disease association studies, biomarker discovery, and functional annotation.

---

## Objectives

1. Perform quality assessment of raw sequencing reads.
2. Trim adapters and low-quality bases.
3. Align reads to the human reference genome (GRCh38).
4. Remove PCR duplicates and sort aligned reads.
5. Apply base quality score recalibration (BQSR).
6. Call variants using GATK HaplotypeCaller.
7. Generate a final, indexed VCF file for analysis and visualization.
8. Document the entire workflow for reproducibility and sharing.

---

## Methodology

The pipeline consists of 10 sequential steps, executed via `bsub` on an LSF-managed HPC cluster. All scripts are modular and can be run individually.

### 1. Quality Control (FastQC)

**Tool:** `FastQC`  
**Purpose:** Assess the quality of raw FASTQ files.

- Checks per-base quality, GC content, adapter contamination, and overrepresented sequences.
- Output: HTML reports for visual inspection.

```bash
fastqc -t 10 input_R1.fq.gz input_R2.fq.gz -o output_dir/
```

---

### 2. Trimming (Trimmomatic)

**Tool:** `Trimmomatic`  
**Purpose:** Remove adapters and low-quality bases.

- Removes Illumina adapters using `ILLUMINACLIP`.
- Trims low-quality bases with `SLIDINGWINDOW`, `LEADING`, `TRAILING`.
- Applies `HEADCROP:15` and `MINLEN:100`.

```bash
java -jar trimmomatic.jar PE -threads 6 ... ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 ...
```

---

### 3. Alignment (BWA + samtools)

**Tools:** `BWA`, `samtools`  
**Purpose:** Map reads to the GRCh38 reference genome.

- Uses `bwa mem` for accurate alignment.
- Converts SAM to BAM and sorts using `samtools`.

```bash
bwa mem ref.fa R1.fq.gz R2.fq.gz | samtools view -bS - > mapped.bam
samtools sort mapped.bam -o sorted.bam
```

---

### 4. Mark Duplicates (GATK/Picard)

**Tool:** `GATK MarkDuplicates`  
**Purpose:** Identify and remove PCR duplicates.

- Reduces false positives in variant calling.
- Outputs a deduplicated BAM and metrics file.

```bash
gatk MarkDuplicates --INPUT sorted.bam --OUTPUT dedup.bam --METRICS_FILE metrics.txt --REMOVE_DUPLICATES true --CREATE_INDEX true
```

---

### 5. SplitNCigarReads (GATK)

**Tool:** `GATK SplitNCigarReads`  
**Purpose:** Split reads into exons and remove introns (intended for RNA-Seq).

> ⚠️ **Note:** This step is **not standard for WES**.

```bash
gatk SplitNCigarReads -R ref.fa -I dedup.bam -O split.bam
```

---

### 6. Add Read Groups (Picard)

**Tool:** `Picard AddOrReplaceReadGroups`  
**Purpose:** Add metadata (sample, library, platform) required by GATK.

- Critical for downstream tools like BQSR and variant calling.

```bash
java -jar AddOrReplaceReadGroups.jar I=split.bam O=rg.bam RGSM=Sample1 RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1
```

---

### 7. Base Quality Score Recalibration (BQSR)

**Tools:** `GATK BaseRecalibrator` → `ApplyBQSR`  
**Purpose:** Correct systematic errors in base quality scores.

- **Step 1:** `BaseRecalibrator` builds a recalibration model using known variants (dbSNP).
- **Step 2:** `ApplyBQSR` applies the model to the BAM file.

```bash
gatk BaseRecalibrator -R ref.fa -I rg.bam --known-sites dbsnp.vcf.gz -O recal.table
gatk ApplyBQSR -R ref.fa -I rg.bam --bqsr-recal-file recal.table -O recalibrated.bam
```

---

### 8. Variant Calling (HaplotypeCaller)

**Tool:** `GATK HaplotypeCaller`  
**Purpose:** Call SNPs and indels per-sample.

- Uses local de-novo assembly for high accuracy.
- Outputs a VCF file with annotations.

```bash
gatk HaplotypeCaller \
  -R ref.fa \
  -I recalibrated.bam \
  -O variants.vcf \
  --dbsnp dbsnp.vcf.gz \
  --stand-call-conf 20.0 \
  --dont-use-soft-clipped-bases \
  --create-output-variant-index true \
  --annotation Coverage \
  --annotation-group StandardHCAnnotation
```

---

## Significance

- **Reproducibility:** Fully documented and modular scripts.
- **Accuracy:** Follows GATK best practices for reliable variant detection.
- **Scalability:** Designed for HPC environments using job scheduling.
- **Educational Value:** Clear workflow for students and new bioinformaticians.
- **Adaptability:** Can be modified for other NGS applications.

---

## Results

- **Input:** Raw paired-end WES FASTQ files.
- **Final Output:** A compressed, indexed VCF file:  
  `ML00984_variants.vcf.gz` + `ML00984_variants.vcf.gz.tbi`
- **Variants Detected:** SNPs and indels across the exome.
- **Visualization:** Can be viewed in IGV, Ensembl, or annotated with SnpEff/ANNOVAR.

---

## Conclusion

The pipeline successfully processes WES data from raw FASTQ to final variant calls. Despite the inclusion of `SplitNCigarReads` (an RNA-Seq step), the workflow remains robust and produces a high-quality VCF. For future projects, this step can be omitted for WES.

This pipeline is now **ready for use on additional samples** and serves as a template for standardized WES analysis.

---

## Usage

### 1. Clone the Repository
```bash
git clone https://github.com/yourusername/wes-pipeline.git
cd wes-pipeline
```

### 2. Update Paths
Edit the scripts to point to:
- Your FASTQ files
- Reference genome (`GRCh38.primary_assembly.genome.fa`)
- dbSNP file
- Tools (GATK, Picard, Trimmomatic)

### 3. Run the Pipeline
Execute scripts in order:
```bash
bash scripts/01_fastqc.sh
bash scripts/02_trimmomatic.sh
bash scripts/03_bwa_align.sh
...
bash scripts/10_haplotypecaller.sh
```

### 4. Monitor Jobs
```bash
bjobs
```

---

## Requirements

| Tool | Version | Purpose |
|------|--------|--------|
| FastQC | v0.11.9 | Quality control |
| Trimmomatic | v0.39 | Adapter trimming |
| BWA | 0.7.12+ | Read alignment |
| samtools | 1.9+ | BAM manipulation |
| GATK | 4.1.8.1 | Variant calling suite |
| Picard | 1.100+ | Read group management |
| Java | 8+ | Required for GATK/Picard |
| LSF (bsub) | - | Job scheduling on HPC |

> **Reference Genome:**  
> GRCh38 (Gencode v30)  
> dbSNP: `dbsnp_151_GRCH38p7_GATK_All_20180418.vcf.gz`

---

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this pipeline in your research, please cite:

> Jeyarish V *Whole Exome Sequencing Variant Calling Pipeline*. GitHub Repository, 2025.  
> https://github.com/Jeyarish-007/wes-pipeline

Also cite the tools used:
- **GATK**: Van der Auwera et al., *Curr Protoc Bioinformatics*, 2013.
- **FastQC**: Andrews, S. (2010). *FastQC: a quality control tool for high throughput sequence data*.
- **Trimmomatic**: Bolger et al., *Bioinformatics*, 2014.
- **BWA**: Li & Durbin, *Bioinformatics*, 2009.

---

## Contact

GitHub: [@Jeyarish-007](https://github.com/Jeyarish-007)
For issues, suggestions, or collaboration, please open an issue or contact directly.

---
