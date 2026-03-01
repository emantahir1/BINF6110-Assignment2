# ============================================================
# Course: BINF6110 Assignment 2
# Project: Differential Expression Analysis of Yeast Biofilm Development
# Author: Eman Tahir
# Dataset: Mardanov et al. 2020 - S. cerevisiae strain I-329
# SRA BioProject: PRJNA592304
# ============================================================

# ============================================================
# SECTION 1: LOAD LIBRARIES
# ============================================================
library(DESeq2)          # Differential expression analysis (Love et al. 2014)
library(tximport)        # Import Salmon quantification output (Soneson et al. 2015)
library(GenomicFeatures) # Build tx2gene mapping from GTF annotation
library(tidyverse)       # Data manipulation and ggplot2 for visualization
library(pheatmap)        # Heatmap visualization
library(clusterProfiler) # GO enrichment analysis (Wu et al. 2021)
library(org.Sc.sgd.db)   # S. cerevisiae gene annotation database
library(enrichplot)      # Enrichment result visualization
library(AnnotationDbi)   # Database interface for annotation
library(ggrepel)         # Non-overlapping labels for ggplot2

# ============================================================
# SECTION 2: SET UP FILE PATHS AND METADATA
# ============================================================

setwd("~/assignment2")

# Point to Salmon quantification output files for all 9 samples
# Samples are ordered by stage: Mature, Thin, Early
files <- file.path("quant", c(
  "SRR10551657", "SRR10551658", "SRR10551659",  # Mature biofilm (day 109)
  "SRR10551660", "SRR10551661", "SRR10551662",  # Thin biofilm (day 83)
  "SRR10551663", "SRR10551664", "SRR10551665"   # Early biofilm (day 38)
), "quant.sf")

# Verify all files exist before proceeding
file.exists(files)

# Create metadata table linking samples to their biofilm stage
# Levels set in developmental order: Early -> Thin -> Mature
sampleTable <- data.frame(
  sample = c("SRR10551657", "SRR10551658", "SRR10551659",
             "SRR10551660", "SRR10551661", "SRR10551662",
             "SRR10551663", "SRR10551664", "SRR10551665"),
  stage = factor(c("Mature", "Mature", "Mature",
                   "Thin", "Thin", "Thin",
                   "Early", "Early", "Early"),
                 levels = c("Early", "Thin", "Mature"))
)

names(files) <- sampleTable$sample
sampleTable

# ============================================================
# SECTION 3: BUILD TRANSCRIPT-TO-GENE MAPPING (tx2gene)
# ============================================================
# We need to map transcript-level Salmon counts to gene-level counts
# Using the NCBI RefSeq annotation for S. cerevisiae R64 genome

txdb <- makeTxDbFromGFF("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gtf.gz")

# Extract transcript to gene ID mapping
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Remove version numbers from transcript IDs (e.g. NM_001178148.1 -> NM_001178148)
# This is necessary to match the IDs output by Salmon
tx2gene$TXNAME <- sub("\\..*", "", tx2gene$TXNAME)
head(tx2gene)

# ============================================================
# SECTION 4: IMPORT SALMON QUANTIFICATION WITH TXIMPORT
# ============================================================
# tximport imports transcript-level estimates and summarizes to gene level
# This approach is recommended over importing raw counts (Soneson et al. 2015)

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

# Check dimensions: should be genes x samples
dim(txi$counts)  # Expected: ~6000 genes x 9 samples
head(txi$counts)

# ============================================================
# SECTION 5: QUALITY CONTROL - SALMON MAPPING RATES
# ============================================================
# Mapping rates were extracted from Salmon aux_info/meta_info.json files
# All samples showed acceptable mapping rates (73.7% - 92.1%)
# Trimming was not performed prior to quantification, consistent with
# current best practices for RNA-seq

mapping_rates <- data.frame(
  Sample = sampleTable$sample,
  Stage = sampleTable$stage,
  MappingRate = c(82.95, 81.87, 83.26,  # Mature
                  73.75, 76.49, 92.06,  # Thin
                  85.45, 84.57, 85.52)  # Early
)
mapping_rates

# ============================================================
# SECTION 6: DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# ============================================================
# DESeq2 was chosen for its robust handling of small sample sizes (n=3)
# and its well-validated negative binomial model for count data

# Create DESeq2 object with stage as the experimental factor
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~stage)

# Run DESeq2 - estimates size factors, dispersions, and fits the model
dds <- DESeq(dds)

# Check available comparisons
resultsNames(dds)

# Extract pairwise results for all three stage comparisons
# Pairwise comparisons chosen over LRT time-course model for
# clearer interpretation of stage-specific gene expression changes

# Early vs Thin biofilm
res_early_thin <- results(dds, contrast = c("stage", "Thin", "Early"))

# Early vs Mature biofilm (most biologically informative comparison)
res_early_mature <- results(dds, contrast = c("stage", "Mature", "Early"))

# Thin vs Mature biofilm
res_thin_mature <- results(dds, contrast = c("stage", "Mature", "Thin"))

# Summarize significant DE genes (padj < 0.05) for each comparison
summary(res_early_thin, alpha = 0.05)
summary(res_early_mature, alpha = 0.05)
summary(res_thin_mature, alpha = 0.05)

# ============================================================
# SECTION 7: VISUALIZATION
# ============================================================

# Variance stabilizing transformation for visualization
# VST stabilizes variance across the range of mean values
vsd <- vst(dds)

# ----------------------------------------------------------
# FIGURE 1: PCA Plot with sample labels
# ----------------------------------------------------------
# PCA reduces high-dimensional gene expression data to 2D
# allowing visualization of overall sample relationships

pca_data <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = stage, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Yeast Biofilm Samples") +
  theme_bw()
pca_plot

ggsave("~/BINF6110-Assignment2/figures/pca.png",
       pca_plot, width = 8, height = 6, dpi = 300)

# ----------------------------------------------------------
# FIGURE 2: Volcano Plot with labeled top DE genes
# ----------------------------------------------------------
# Volcano plot shows log2 fold change vs statistical significance
# Cutoffs: padj < 0.05 and |log2FC| > 1 (2-fold change)
# Top 15 most significant genes labeled using ggrepel

res_df <- as.data.frame(res_early_mature)
res_df$gene <- rownames(res_df)
res_df$significant <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
  ifelse(res_df$log2FoldChange > 0, "Up in Mature", "Down in Mature"),
  "Not Sig")
res_df <- na.omit(res_df)

# Get top 15 most significant genes for labeling
top_labeled <- res_df %>%
  arrange(padj) %>%
  head(15)

gene_labels <- AnnotationDbi::select(org.Sc.sgd.db,
                                     keys = top_labeled$gene,
                                     columns = c("ORF", "GENENAME"),
                                     keytype = "ORF")

top_labeled <- merge(top_labeled, gene_labels, by.x = "gene", by.y = "ORF")
top_labeled$label <- ifelse(is.na(top_labeled$GENENAME),
                            top_labeled$gene,
                            paste0(top_labeled$GENENAME, " (", top_labeled$gene, ")"))

volcano_plot <- ggplot(res_df, aes(x = log2FoldChange,
                                   y = -log10(padj),
                                   color = significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Up in Mature" = "red",
                                "Down in Mature" = "blue",
                                "Not Sig" = "gray")) +
  geom_label_repel(data = top_labeled,
                   aes(label = label),
                   size = 3,
                   max.overlaps = 20,
                   box.padding = 0.5,
                   show.legend = FALSE) +
  labs(x = "Log2 Fold Change",
       y = "-Log10 adjusted p-value",
       title = "Volcano Plot: Early vs Mature Biofilm",
       color = "Direction") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_bw()
volcano_plot

ggsave("~/BINF6110-Assignment2/figures/volcano.png",
       volcano_plot, width = 10, height = 8, dpi = 300)

# ----------------------------------------------------------
# FIGURE 3: Heatmap - Top 30 DE genes
# ----------------------------------------------------------
# Heatmap shows expression patterns of most significant genes
# across all 9 samples, scaled by row (z-score)

top_genes <- head(order(res_early_mature$padj, na.last = TRUE), 30)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))

# Convert ORF IDs to readable gene names using SGD database
gene_names_map <- AnnotationDbi::select(org.Sc.sgd.db,
                                        keys = rownames(mat),
                                        columns = c("ORF", "GENENAME"),
                                        keytype = "ORF")
matched <- gene_names_map$GENENAME[match(rownames(mat), gene_names_map$ORF)]
rownames(mat_scaled) <- ifelse(is.na(matched),
                               rownames(mat),
                               paste0(matched, " (", rownames(mat), ")"))

# Stage annotation for column colors
annotation_col <- data.frame(Stage = sampleTable$stage)
rownames(annotation_col) <- sampleTable$sample

png("~/BINF6110-Assignment2/figures/heatmap.png", width = 800, height = 900)
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 30 DE Genes: Early vs Mature Biofilm")
dev.off()

# ----------------------------------------------------------
# FIGURE 4: Sample-to-Sample Distance Heatmap
# ----------------------------------------------------------
# Shows pairwise similarity between all 9 samples
# Complements PCA by showing replicate consistency within stages

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(sampleTable$stage, c(1,2,3,1,2,3,1,2,3), sep="_")
colnames(sampleDistMatrix) <- rownames(sampleDistMatrix)

png("~/BINF6110-Assignment2/figures/sample_distance.png", width = 800, height = 700)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         color = colorRampPalette(c("darkblue", "white"))(100),
         main = "Sample-to-Sample Distance Matrix")
dev.off()

# ----------------------------------------------------------
# FIGURE 5: DESeq2 Dispersion Plot
# ----------------------------------------------------------
# Shows per-gene dispersion estimates and fitted trend
# Confirms DESeq2 model fit is appropriate for this dataset

png("~/BINF6110-Assignment2/figures/dispersion.png", width = 800, height = 600)
plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
dev.off()

# ----------------------------------------------------------
# FIGURE 6: MA Plot - Early vs Mature
# ----------------------------------------------------------
# MA plot shows log2 fold change vs mean expression
# Complements volcano plot by showing relationship with expression level

png("~/BINF6110-Assignment2/figures/ma_plot.png", width = 800, height = 600)
plotMA(res_early_mature,
       main = "MA Plot: Early vs Mature Biofilm",
       ylim = c(-10, 10),
       alpha = 0.05)
dev.off()

# ============================================================
# SECTION 8: FUNCTIONAL ANNOTATION - GO ENRICHMENT (ORA)
# ============================================================
# Over-representation analysis (ORA) using Fisher's exact test
# Tests whether GO biological process terms are enriched in our
# DE gene sets compared to all expressed genes as background
# Genes split by direction to distinguish upregulated vs downregulated pathways

# Upregulated in Mature biofilm
upregulated <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(gene)

# Downregulated in Mature biofilm (upregulated in Early)
downregulated <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(gene)

# Background: all expressed genes
all_genes <- res_df %>%
  pull(gene)

# GO ORA for upregulated genes in Mature biofilm
ego_up <- enrichGO(gene = upregulated,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   keyType = "ORF",
                   ont = "BP",
                   pAdjustMethod = "BH",  # Benjamini-Hochberg correction
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)

# GO ORA for downregulated genes in Mature (high in Early)
ego_down <- enrichGO(gene = downregulated,
                     universe = all_genes,
                     OrgDb = org.Sc.sgd.db,
                     keyType = "ORF",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = FALSE)

# ----------------------------------------------------------
# FIGURE 7: GO Dotplot - Upregulated in Mature
# ----------------------------------------------------------
png("~/BINF6110-Assignment2/figures/go_upregulated.png",
    width = 800, height = 600)
dotplot(ego_up, showCategory = 15,
        title = "GO BP: Upregulated in Mature Biofilm")
dev.off()

# ----------------------------------------------------------
# FIGURE 8: GO Dotplot - Downregulated in Mature (Up in Early)
# ----------------------------------------------------------
png("~/BINF6110-Assignment2/figures/go_downregulated.png",
    width = 800, height = 600)
dotplot(ego_down, showCategory = 15,
        title = "GO BP: Downregulated in Mature Biofilm (Upregulated in Early)")
dev.off()

# ============================================================
# SECTION 9: FLO11 EXPRESSION PLOT
# ============================================================
# FLO11 is the key biofilm adhesion gene - plotting its expression
# across all three stages directly illustrates the central finding

flo11_counts <- plotCounts(dds, gene = "YIR019C",
                           intgroup = "stage",
                           returnData = TRUE)

flo11_plot <- ggplot(flo11_counts, aes(x = stage, y = count, fill = stage)) +
  geom_boxplot(alpha = 0.7) +
  geom_point(size = 3) +
  scale_y_log10() +
  labs(title = "FLO11 (YIR019C) Expression Across Biofilm Stages",
       x = "Biofilm Stage",
       y = "Normalized Count (log10)") +
  theme_bw() +
  theme(legend.position = "none")
flo11_plot

ggsave("~/BINF6110-Assignment2/figures/flo11_expression.png",
       flo11_plot, width = 6, height = 5, dpi = 300)

# ============================================================
# SECTION 10: TOP DE GENES TABLE
# ============================================================
# Table of top 20 most significant DE genes from Early vs Mature
# sorted by adjusted p-value

top_genes_table <- as.data.frame(res_early_mature) %>%
  rownames_to_column("ORF") %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(20) %>%
  select(ORF, log2FoldChange, padj)

# Add readable gene names from SGD database
gene_names <- AnnotationDbi::select(org.Sc.sgd.db,
                                    keys = top_genes_table$ORF,
                                    columns = c("ORF", "GENENAME"),
                                    keytype = "ORF")

top_genes_table <- merge(top_genes_table, gene_names, by = "ORF") %>%
  arrange(padj) %>%
  mutate(log2FoldChange = round(log2FoldChange, 3),
         padj = signif(padj, 3))

top_genes_table

# Save table as CSV for inclusion in results
write.csv(top_genes_table,
          "~/BINF6110-Assignment2/results/top_genes.csv",
          row.names = FALSE)

# ============================================================
# SESSION INFO - for reproducibility
# ============================================================
sessionInfo()