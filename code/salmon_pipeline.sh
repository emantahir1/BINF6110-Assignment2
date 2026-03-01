#!/bin/bash
# ============================================================
# BINF6110 Assignment 2
# Salmon Quantification Pipeline
# Dataset: Mardanov et al. 2020 - S. cerevisiae strain I-329
# SRA BioProject: PRJNA592304
# Run on: Ubuntu 24 (UTM virtual machine on MacOS)
# Tools: SRA toolkit v3.2.1, Salmon v1.10.2, FastQC, MultiQC
# ============================================================

# ============================================================
# SECTION 1: INSTALL REQUIRED TOOLS
# ============================================================
sudo apt update
sudo apt install sra-toolkit -y
sudo apt install salmon -y
sudo apt install fastqc -y
pip install multiqc --break-system-packages

# ============================================================
# SECTION 2: DOWNLOAD RAW READS FROM SRA
# ============================================================
# 9 samples across 3 biofilm stages (3 replicates each):
# Early biofilm (day 38):   SRR10551663, SRR10551664, SRR10551665
# Thin biofilm (day 83):    SRR10551660, SRR10551661, SRR10551662
# Mature biofilm (day 109): SRR10551657, SRR10551658, SRR10551659
#
# fasterq-dump downloads reads from SRA
# --split-files separates paired-end reads (data confirmed single-end)

mkdir -p ~/assignment2/raw_reads
cd ~/assignment2/raw_reads

fasterq-dump SRR10551665 --outdir . --split-files &
fasterq-dump SRR10551664 --outdir . --split-files &
fasterq-dump SRR10551663 --outdir . --split-files &
fasterq-dump SRR10551662 --outdir . --split-files &
fasterq-dump SRR10551661 --outdir . --split-files &
fasterq-dump SRR10551660 --outdir . --split-files &
fasterq-dump SRR10551659 --outdir . --split-files &
fasterq-dump SRR10551658 --outdir . --split-files &
fasterq-dump SRR10551657 --outdir . --split-files

# Verify all downloads completed
ls -lh ~/assignment2/raw_reads

# ============================================================
# SECTION 3: FASTQC QUALITY CONTROL
# ============================================================
# FastQC assesses raw read quality metrics for each sample including:
# per-base quality scores, sequence length distribution,
# GC content, adapter content, and duplication levels

mkdir -p ~/assignment2/fastqc_results

for sample in SRR10551657 SRR10551658 SRR10551659 \
              SRR10551660 SRR10551661 SRR10551662 \
              SRR10551663 SRR10551664 SRR10551665; do
    echo "Running FastQC on: $sample"
    fastqc ~/assignment2/raw_reads/${sample}_1.fastq \
           -o ~/assignment2/fastqc_results/
done

# ============================================================
# SECTION 4: MULTIQC - AGGREGATE FASTQC RESULTS
# ============================================================
# MultiQC combines all 9 FastQC reports into a single
# interactive HTML report for easy comparison across samples

multiqc ~/assignment2/fastqc_results/ \
        -o ~/assignment2/multiqc_results/

echo "FastQC and MultiQC complete!"

# ============================================================
# SECTION 5: DOWNLOAD REFERENCE TRANSCRIPTOME AND GENOME
# ============================================================
# Reference: S. cerevisiae R64 genome (GCF_000146045.2) from NCBI RefSeq
# Transcriptome (rna.fna.gz): used as the main Salmon index target
# Genome (genomic.fna.gz): used as a decoy to prevent spurious mapping
# of intronic/intergenic reads to the transcriptome (Srivastava et al. 2020)

mkdir -p ~/assignment2/reference
cd ~/assignment2/reference

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# ============================================================
# SECTION 6: BUILD DECOY-AWARE SALMON INDEX
# ============================================================
# Decoy-aware indexing reduces spurious mapping of reads
# originating from non-transcribed genomic regions
# Transcriptome must come FIRST in the concatenated gentrome file

cd ~/assignment2/reference

# Step 1: Extract chromosome names from genome for decoy list
grep "^>" <(gunzip -c GCF_000146045.2_R64_genomic.fna.gz) \
    | cut -d " " -f 1 \
    | sed 's/>//' > decoys.txt

# Verify decoys file
cat decoys.txt

# Step 2: Concatenate transcriptome and genome into gentrome
cat GCF_000146045.2_R64_rna.fna.gz GCF_000146045.2_R64_genomic.fna.gz > gentrome.fna.gz

# Step 3: Build Salmon index
# -t: input gentrome file
# -d: decoy sequence list
# -i: output index directory
# -p: number of threads
salmon index \
    -t gentrome.fna.gz \
    -d decoys.txt \
    -i salmon_index \
    -p 4

# ============================================================
# SECTION 7: QUANTIFY GENE EXPRESSION WITH SALMON v1.10.2
# ============================================================
# -i: path to index
# -l A: automatic library type detection
#        (detected as "U" = unstranded single-end for all samples)
# -r: single-end reads input
# --validateMappings: validates quasi-mappings to improve accuracy
# -p: number of threads
# -o: output directory for each sample

cd ~/assignment2

for sample in SRR10551657 SRR10551658 SRR10551659 \
              SRR10551660 SRR10551661 SRR10551662 \
              SRR10551663 SRR10551664 SRR10551665; do
    echo "Quantifying sample: $sample"
    salmon quant \
        -i reference/salmon_index \
        -l A \
        -r raw_reads/${sample}_1.fastq \
        --validateMappings \
        -p 4 \
        -o quant/${sample}
done

echo "Quantification complete!"

# ============================================================
# SECTION 8: QUALITY CONTROL - CHECK SALMON MAPPING RATES
# ============================================================
# Mapping rates extracted from Salmon output metadata files
# All samples showed acceptable mapping rates (73.7% - 92.1%)

grep "percent_mapped" ~/assignment2/quant/*/aux_info/meta_info.json

# Results:
# SRR10551657 (Mature): 82.95%
# SRR10551658 (Mature): 81.87%
# SRR10551659 (Mature): 83.26%
# SRR10551660 (Thin):   73.75%
# SRR10551661 (Thin):   76.49%
# SRR10551662 (Thin):   92.06%
# SRR10551663 (Early):  85.45%
# SRR10551664 (Early):  84.57%
# SRR10551665 (Early):  85.52%

# ============================================================
# SECTION 9: TRANSFER RESULTS TO MAC FOR R ANALYSIS
# ============================================================
# Quant folders compressed and transferred to Mac for
# downstream analysis in RStudio (deseq2_analysis.R)

cd ~/assignment2
zip -r quant.zip quant/
