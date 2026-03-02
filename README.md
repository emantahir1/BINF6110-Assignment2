# BINF6110 Assignment 2: Differential Expression Analysis of Yeast Biofilm Development

**Course:** BINF*6110 - Applied Bioinformatics  
**Dataset:** Mardanov et al. 2020 - *Saccharomyces cerevisiae* strain I-329  
**BioProject:** [PRJNA592304](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA592304)

---

## Table of Contents
1. [Introduction](#introduction)
2. [Methods](#methods)
3. [Results](#results)
4. [Discussion](#discussion)
5. [References](#references)

---
## Introduction

*Saccharomyces cerevisiae* is a single-celled fungal microorganism and one of the most well-studied eukaryotes in biology (Parapouli et al., 2020). It is widely used in the production of fermented foods and beverages including wine, beer, and bread, and serves as a model organism for studying eukaryotic cell biology (Mbuyane et al., 2021). Beyond its role in standard fermentation, certain strains of *S. cerevisiae* are capable of forming a surface biofilm, where cells stick to each other and the surface (Speranza et al., 2020). This is also known as velum, or flor, on the surface of aging wine (David-Vaizant & Alexandre, 2018). This biofilm is formed through a structured, multi-stage developmental process that plays a central role in the production of biologically aged wines such as fino sherry, contributing characteristic aldehyde and acetal aromas to the final product (David-Vaizant & Alexandre, 2018).

As environmental conditions change throughout fermentation, the yeast must transition from active fermentation to a biofilm lifestyle under increasingly nutrient-poor conditions (Mardanov et al., 2020). This transition requires transcriptional reprogramming, thousands of genes must be switched on or off as the yeast adapts its metabolism, cell wall structure, and stress responses (Mardanov et al., 2020). Understanding which genes drive this transition has practical implications for controlling wine quality and for the broader biology of yeast stress adaptation.

A central molecular driver of *S. cerevisiae* biofilm formation is *FLO11*, which encodes a GPI-anchored cell wall glycoprotein that mediates cell-cell and cell-surface adhesion (Bou Zeidan et al., 2013). Its expression is tightly regulated by nutrient-sensing signalling pathways, making it a key link between environmental conditions and biofilm structure (Bou Zeidan et al., 2013). However, biofilm development involves more than adhesion, it requires coordinated changes in metabolism, cell wall composition, and stress tolerance that are not yet fully characterized at the gene level (Irianto et al., 2025).

Mardanov et al. (2020) addressed this by profiling gene expression in *S. cerevisiae* strain I-329, an industrial flor yeast, across three stages of velum development: early biofilm (day 38), thin biofilm (day 83), and mature biofilm (day 109). Their published metadata reveals a clear environmental gradient across these stages: glucose drops from 0.2 g/L in early biofilm to less than 0.1 g/L in mature biofilm, ethanol decreases from 12.4% to 9.6% v/v, and aldehyde concentrations increase from 382.8 to 668.8 mg/L. This provides a well-characterized chemical backdrop against which transcriptional changes can be interpreted (Mardanov et al., 2020). Despite its importance to the wine industry, the gene-level transcriptional changes driving each stage of velum maturation remain incompletely characterized.

In this study, we reanalyze the Mardanov et al. (2020) dataset to identify the key genes and biological pathways driving each stage of velum development in *S. cerevisiae* strain I-329. Using RNA-seq data from nine samples across three developmental stages, we quantify gene expression, identify differentially expressed genes, and perform functional annotation to characterize the transcriptional changes underlying the transition from early to mature biofilm.

---
## Methods

### Data Acquisition

RNA-seq data for *S. cerevisiae* strain I-329 were obtained from the NCBI Sequence Read Archive under BioProject PRJNA592304 (Mardanov et al., 2020). Nine samples were downloaded using fasterq-dump v3.2.1 from the NCBI SRA toolkit, corresponding to three biological replicates at each of three velum developmental stages: early biofilm (SRR10551663–SRR10551665, day 38), thin biofilm (SRR10551660–SRR10551662, day 83), and mature biofilm (SRR10551657–SRR10551659, day 109). Reads were confirmed to be single-end upon download. The *S. cerevisiae* R64 reference transcriptome (GCF_000146045.2_R64_rna.fna.gz) and genome (GCF_000146045.2_R64_genomic.fna.gz) were downloaded from NCBI RefSeq for use in quantification and decoy index construction.

### Quality Control

Raw read quality was assessed for all nine samples using FastQC v0.12.1. FastQC was chosen because it provides a comprehensive and standardized set of quality metrics including per-base quality scores, sequence length distribution, GC content, adapter content, and duplication levels, allowing identification of systematic issues before quantification (Andrews, 2010). Results from all nine samples were aggregated into a single interactive HTML report using MultiQC v1.33, which allows direct visual comparison of quality metrics across all samples simultaneously to identify outliers (Ewels et al., 2016). Mapping rates from Salmon quantification (described below) were additionally extracted from the `aux_info/meta_info.json` output file produced for each sample as a secondary quality metric confirming successful quantification. Adapter trimming was not performed prior to quantification as it has been shown to provide minimal benefit when using pseudoaligners on modern sequencing data, and can, in some cases, introduce biases by artificially truncating reads (Liao & Shi, 2020). GC bias correction was not applied as the dataset consists of standard RNA-seq reads with no known GC content issues reported in the original study.

### Quantification

Reads were quantified using Salmon v1.10.2, a quasi-mapping-based pseudoaligner (Patro et al., 2017). Salmon was chosen over alignment-based tools such as STAR or HISAT2 because the *S. cerevisiae* genome is comprehensively annotated, and the goal of this study is gene-level quantification of known transcripts rather than novel isoform or splice site discovery, which are the primary use cases for full alignment (Patro et al., 2017). Salmon is substantially faster than alignment-based tools while producing equally accurate expression estimates for well-annotated genomes (Patro et al., 2017).

A decoy-aware index was constructed using the full *S. cerevisiae* R64 genome as a decoy sequence. Without a decoy, reads originating from genomic regions that share sequence similarity with transcripts can be incorrectly assigned, inflating expression estimates. Including the genome as a decoy forces Salmon to consider these alternative mapping locations and reduces spurious quantification (Srivastava et al., 2020). Genome chromosome names were extracted to serve as the decoy list, and the transcriptome and genome were concatenated into a single gentrome file with the transcriptome placed first as required by Salmon. Salmon v1.10.2 was then run for each sample with the following parameters:

- `-l A` : automatic library type detection, which determines strandedness and read orientation directly from the data. All samples were confirmed as unstranded single-end (library type U)
- `--validateMappings` : enables selective alignment, which validates quasi-mappings against the full index to improve mapping accuracy by rescuing reads that would otherwise be incorrectly mapped or discarded
- `-p 4` : parallelizes computation across four threads to reduce runtime

Full bash commands are provided in `code/salmon_pipeline.sh`.

### Differential Expression Analysis

All downstream analysis was performed in R v4.5.1. R was chosen over Python-based alternatives because the Bioconductor ecosystem provides the most mature and widely validated suite of tools for RNA-seq analysis, including DESeq2, tximport, and clusterProfiler, all of which have extensive documentation, active maintenance, and peer-reviewed publications supporting their use (Muzellec et al., 2023).

Salmon output was imported using tximport v1.36.1. tximport was used rather than count-based tools such as featureCounts or HTSeq because it is specifically designed to work with pseudoalignment output from tools like Salmon (Soneson et al., 2016). Unlike featureCounts and HTSeq, which require a genome-aligned BAM file, tximport works directly with transcript-level abundance estimates and aggregates them to gene-level counts while propagating quantification uncertainty into the downstream statistical model. This approach has been shown to improve the sensitivity and accuracy of differential expression results compared to simple count summarization (Soneson et al., 2016).

A transcript-to-gene mapping was constructed from the NCBI RefSeq GTF annotation for *S. cerevisiae* R64 using GenomicFeatures v1.60.0. GenomicFeatures was used rather than a pre-built tx2gene table because it constructs the mapping directly from the same annotation file used to build the Salmon index, ensuring complete consistency between the quantification and the gene-level summarization (Lawrence et al., 2013). Transcript version numbers were removed from IDs to ensure compatibility with Salmon output.

Differential expression analysis was performed using DESeq2 v1.48.2 (Love et al., 2014). DESeq2 was chosen over EdgeR because it is specifically designed to handle experiments with small numbers of replicates. With only three samples per group, estimating gene-level variability from the data alone is unreliable. DESeq2 addresses this by pooling information across all genes to produce more stable variability estimates for each individual gene, which reduces false positives. It models count data using a negative binomial distribution, which accounts for the fact that variance in RNA-seq data tends to exceed the mean. A DESeq2 object was created using `DESeqDataSetFromTximport()` with biofilm stage as the experimental factor, with levels ordered Early → Thin → Mature. Three pairwise comparisons were performed: Early vs. Thin, Early vs. Mature, and Thin vs. Mature. Pairwise comparisons were chosen over a likelihood ratio test (LRT) time-course model because the three stages represent biologically distinct conditions rather than evenly spaced time points, and pairwise results allow direct interpretation of which genes change between specific stages.

Genes with Benjamini-Hochberg adjusted p-value < 0.05 were considered statistically significant, controlling the false discovery rate at 5% so that on average no more than 5% of genes called significant are expected to be false positives. An additional |log2 fold change| > 1 threshold was applied to restrict results to genes showing at least a 2-fold change in expression, as smaller fold changes are unlikely to be biologically meaningful in the context of a major developmental transition such as biofilm maturation.

### Visualization

Variance-stabilizing transformation (VST) was applied to normalized count data prior to all visualization using DESeq2's `vst()` function. VST was chosen over simple log2 transformation because it stabilizes variance across the full range of mean expression values, preventing highly expressed genes from dominating distance-based analyses. VST was preferred over DESeq2's alternative rlog transformation because rlog is substantially slower, and for a nine-sample dataset, both transformations produce comparable results, making VST the more practical choice (Love et al., 2014).

**QC Visualization**

Principal component analysis (PCA) was performed on VST-transformed data to assess overall sample relationships and confirm that biofilm stage was the primary source of transcriptional variation. Individual sample labels were added using ggrepel v0.9.1 to allow identification of individual replicates. A sample-to-sample distance matrix was computed from VST-transformed data using Euclidean distance to provide a complementary view of replicate consistency alongside PCA. DESeq2 dispersion estimates were plotted using `plotDispEsts()` to confirm appropriate model fit. A well-fitted model would show gene-level dispersion estimates shrinking towards the fitted trend curve with few outliers.

**Differential Expression Visualization**

A volcano plot was generated for the Early vs. Mature comparison, as it represented the most biologically divergent stages with the largest number of DE genes. The top 15 most significant genes were labelled using ggrepel v0.9.1 to identify key genes directly on the plot. An MA plot was generated for the same comparison using DESeq2's `plotMA()` to visualize the relationship between fold change and mean expression level, confirming that differential expression was not biased towards highly or lowly expressed genes. A heatmap of the top 30 most significant DE genes was produced using pheatmap v1.0.12. Yeast genes are identified in the reference genome by systematic ORF identifiers such as YIR019C, which are not biologically informative on their own. To make results interpretable, ORF IDs were converted to standard gene names (e.g. *FLO11*, *TDH1*) using org.Sc.sgd.db v3.21.0, an R annotation package maintained by the Saccharomyces Genome Database (SGD) that maps every yeast ORF to its curated gene name and functional annotation. This database was used over generic annotation tools because it is the community-standard resource for *S. cerevisiae* gene annotation. Normalized expression counts for *FLO11* (YIR019C) were plotted across all three stages using DESeq2's `plotCounts()` to directly illustrate the progressive upregulation of the key biofilm adhesion gene. All plots were produced using ggplot2 v3.5.1 within R v4.5.1.

### Functional Annotation

Functional annotation was performed using Gene Ontology over-representation analysis (ORA) via clusterProfiler v4.16.0 (Wu et al., 2021). ORA was chosen over gene set enrichment analysis (GSEA) because the goal here is to identify which biological processes are statistically overrepresented in discrete sets of significant DE genes, rather than ranking all genes by a continuous enrichment score. ORA applies a hypergeometric test to determine whether each GO term appears more frequently among DE genes than expected by chance, given the background of all expressed genes in this dataset.

The Biological Process (BP) ontology was used rather than Molecular Function or Cellular Component because BP terms describe the biological roles and pathways genes participate in, which is the most relevant level of annotation for interpreting transcriptional changes accompanying a major developmental transition. Genes from the Early vs. Mature comparison were split into upregulated and downregulated sets and analyzed separately to preserve directionality and distinguish processes that are activated from those that are repressed during biofilm maturation. A background of all expressed genes was used rather than all annotated genes in the genome, as genes not detected in this experiment cannot be considered part of the testable universe. Multiple testing correction was applied using the Benjamini-Hochberg method with padj < 0.05. A secondary qvalueCutoff of 0.2 was applied as an additional filter, which controls the proportion of false positives among significant results using the q-value framework independently of the adjusted p-value, providing greater confidence in reported terms (Wu et al., 2021). Full R code is provided in `code/deseq2_analysis.R`.

---

## Results


---

## Discussion


---

## References
