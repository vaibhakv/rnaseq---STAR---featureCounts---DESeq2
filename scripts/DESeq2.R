library(tidyverse)
library(Rsubread)
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(ggrepel)

#featureCounts
df <- featureCounts(files = c("data/output/SRR13978640.aln/SRR13978640Aligned.sortedByCoord.out.bam", 
                    "data/output/SRR13978641.aln/SRR13978641Aligned.sortedByCoord.out.bam",  
                    "data/output/SRR13978642.aln/SRR13978642Aligned.sortedByCoord.out.bam",
                    "data/output/SRR13978643.aln/SRR13978643Aligned.sortedByCoord.out.bam",
                    "data/output/SRR13978644.aln/SRR13978644Aligned.sortedByCoord.out.bam",
                    "data/output/SRR13978645.aln/SRR13978645Aligned.sortedByCoord.out.bam"),
                        annot.ext = "genome/Saccharomyces_cerevisiae.R64-1-1.109.gtf", 
                        isGTFAnnotationFile = TRUE,
                        isPairedEnd = TRUE)

raw_counts <- df$counts
colnames(raw_counts) <- paste0("SRR1397864",0:5) 

#metadata
metadata <- tibble(sample = paste0("SRR1397864",0:5), condition = rep(c("SPRC", "control"), each=3))

#making a deseq object from matrix
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ condition )
assay(dds)  

keep <- rowSums(counts(dds)) > 10
dds <-  dds[keep,]

#regularized log transformation
rlog_dds <- rlog(dds)
sampleDist <- dist(t(assay(rlog_dds)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- dds$condition

#heatmap of sample distances
pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         cols = colors)

#Principle component analysis of log transformed counts
plotPCA(rlog_dds)

#normalizing counts 
sizeFactors <- estimateSizeFactors(dds)
normalized_counts <- counts(sizeFactors, normalized=TRUE)

#Variant stabilizing transformation
vst <-  varianceStabilizingTransformation(dds)

plotPCA(vst)

#running DESeq2
res_dds <- DESeq(dds)
res <- results(res_dds)
mcols(res)
summary(res)

#filtering significant results
res_sig <- subset(res, padj < 0.05 & baseMean > 100 & abs(log2FoldChange) > 1)
head(res_sig[order(res_sig$log2FoldChange),])
head(res_sig[order(res_sig$log2FoldChange, decreasing = TRUE),])

#histogram of p-values
hist(res$pvalue)

#MA plot 
res_shrink <- lfcShrink(res_dds, coef = "condition_SPRC_vs_control", type = "apeglm")
plotMA(res_shrink, ylim = c(-3, 3))

#Heatmap of differentially expressed genes
sig_matrix <- normalized_counts[rownames(res_sig),]
sig_matrix.z <- t(apply(sig_matrix, 1, scale))
colnames(sig_matrix.z) <- colnames(sig_matrix)

heat_colors <- colorRampPalette(brewer.pal(6, "PuOr"))(100)

pheatmap(mat = sig_matrix.z, 
         cluster_rows = T, 
         cluster_cols = T, 
         labels_col = colnames(sig_matrix.z), 
         labels_row = rownames(sig_matrix.z),
         color = heat_colors,
         height = 40)

gene_ids_sig <- rownames(res_sig)
gene_ids <- rownames(res)
#bioMart

listEnsembl()
ensembl <- useEnsembl(biomart = "genes")

listDatasets(ensembl) %>% View()
bio.mart <- useMart("ensembl", "scerevisiae_gene_ensembl")

listAttributes(bio.mart) %>% View()
listFilters(bio.mart) %>% View()

bm_attr <- c("ensembl_gene_id", "external_gene_name", 
           "chromosome_name", "start_position", "end_position",
           "description")

bm_res <- getBM(attributes = bm_attr, 
      filters = c("ensembl_gene_id"), 
      values = gene_ids, mart = bio.mart)

bm_res_sig <- getBM(attributes = bm_attr, 
                filters = c("ensembl_gene_id"), 
                values = gene_ids_sig, mart = bio.mart)

bm_res$external_gene_name


res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

res_tb <- res_tb %>% 
  mutate(gene_labels = bm_res$external_gene_name) 

res_tb <- res_tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1) 

ggplot(res_tb, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = gene_labels), alpha = 0.5) +
  geom_vline(xintercept = c(1, -1), linetype = "longdash", alpha = 0.3) +
  geom_hline(yintercept = 2, linetype = "longdash", alpha = 0.3) +
  scale_color_manual(values = c("#333366", "#CC9900")) +
  ggtitle("SPRC vs control - differentially expressed genes") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.2), hjust = 0.5),
        axis.title = element_text(size = rel(1.1)))




