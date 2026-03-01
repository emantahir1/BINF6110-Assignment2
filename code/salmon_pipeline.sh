#!/bin/bash
# ============================================================
# BINF6110 Assignment 2
# Salmon Quantification Pipeline
# Dataset: Mardanov et al. 2020 - S. cerevisiae strain I-329
# SRA BioProject: PRJNA592304
# Run on: Ubuntu 24 (UTM virtual machine on MacOS)
# Tools: SRA toolkit v3.2.1, Salmon v1.10.2
# ============================================================

# SECTION 1: DOWNLOAD RAW READS FROM SRA
# 9 samples across 3 biofilm stages (3 replicates each)
# Early biofilm (day 38):   SRR10551663, SRR10551664, SRR10551665
# Thin biofilm (day 83):    SRR10551660, SRR10551661, SRR10551662
# Mature biofilm (day 109): SRR10551657, SRR10551658, SRR10551659
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
# Data confirmed as single-end (only _1.fastq files contained reads)

# SECTION 2: DOWNLOAD REFERENCE TRANSCRIPTOME AND GENOME
# Reference: S. cerevisiae R64 (GCF_000146045.2) from NCBI RefSeq
# Transcriptome used as Salmon index target
# Genome used as decoy to prevent spurious mapping (Srivastava et al. 2020)
mkdir -p ~/assignment2/reference
cd ~/assignment2/reference

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_rna.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna.gz

# SECTION 3: BUILD DECOY-AWARE SALMON INDEX
# Decoy-aware indexing reduces spurious mapping of intronic reads
# Transcriptome must come first in gentrome file
grep "^>" <(gunzip -c GCF_000146045.2_R64_genomic.fna.gz) \
    | cut -d " " -f 1 \
    | sed 's/>//' > decoys.txt

cat GCF_000146045.2_R64_rna.fna.gz GCF_000146045.2_R64_genomic.fna.gz > gentrome.fna.gz

salmon index \
    -t gentrome.fna.gz \
    -d decoys.txt \
    -i salmon_index \
    -p 4

# SECTION 4: QUANTIFY GENE EXPRESSION WITH SALMON v1.10.2
# -l A: automatic library type detection (detected as unstranded single-end)
# --validateMappings: improves accuracy by validating quasi-mappings
# -p 4: 4 threads
cd ~/assignment2

for sample in SRR10551657 SRR10551658 SRR10551659 \
              SRR10551660 SRR10551661 SRR10551662 \
              SRR10551663 SRR10551664 SRR10551665; do
    echo "Processing sample: $sample"
    salmon quant \
        -i reference/salmon_index \
        -l A \
        -r raw_reads/${sample}_1.fastq \
        --validateMappings \
        -p 4 \
        -o quant/${sample}
done

echo "Quantification complete!"

# SECTION 5: QUALITY CONTROL - CHECK MAPPING RATES
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

# Quant folders transferred to Mac via zip for R analysis
# R analysis continues in deseq2_analysis.R
