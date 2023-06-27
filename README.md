RNA-Seq – Differential Gene Expression Analysis
================
vaibhakv
2023-06-27

## Experimental setup

- Saccharomyces cerevisiae was used as a model to study the mechanism of
  endogenous H2S that promoted the growth rate of yeast.
- The growth of fungi is controlled by several factors, one of which is
  signaling molecules, such as hydrogen sulfide (H2S), which was
  traditionally regarded as a toxic gas without physiological function.
  However, recent studies have revealed that H2S is produced
  enzymatically and endogenously in several species, where it serves as
  a gaseous signaling molecule performing a variety of critical
  biological functions.
- However, the influence of this endogenous H2S on the biological
  activities occurring within the pathogenic fungi, such as
  transcriptomic and phenotypic alternations, has not been elucidated so
  far. Therefore, the present study was aimed to decipher this concern
  by utilizing S-propargyl-cysteine (SPRC) as a novel and stable donor
  of H2S and Saccharomyces cerevisiae as a fungal model.
- In this analysis, we’ll find out the differentially expressed genes
  between two conditions with three replicates each -
  - 3 samples in the presence of SPRC -
    - SRR13978640
    - SRR13978641
    - SRR13978642
  - 3 samples in the absence of SPRC (control)
    - SRR13978643
    - SRR13978644
    - SRR13978645

``` r
library(tidyverse)
library(Rsubread)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(biomaRt)
library(ggrepel)
```

### featureCounts

Use featureCounts to generate the count table

``` r
df <- featureCounts(files = c("data/output/SRR13978640.aln/SRR13978640Aligned.sortedByCoord.out.bam", 
                    "data/output/SRR13978641.aln/SRR13978641Aligned.sortedByCoord.out.bam",  
                    "data/output/SRR13978642.aln/SRR13978642Aligned.sortedByCoord.out.bam",
                    "data/output/SRR13978643.aln/SRR13978643Aligned.sortedByCoord.out.bam",
                    "data/output/SRR13978644.aln/SRR13978644Aligned.sortedByCoord.out.bam",
                    "data/output/SRR13978645.aln/SRR13978645Aligned.sortedByCoord.out.bam"),
                        annot.ext = "genome/Saccharomyces_cerevisiae.R64-1-1.109.gtf", 
                        isGTFAnnotationFile = TRUE,
                        isPairedEnd = TRUE,
                        nthreads = 4)
```

    ## 
    ##         ==========     _____ _    _ ____  _____  ______          _____  
    ##         =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
    ##           =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
    ##             ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
    ##               ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
    ##         ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
    ##        Rsubread 2.14.2
    ## 
    ## //========================== featureCounts setting ===========================\\
    ## ||                                                                            ||
    ## ||             Input files : 6 BAM files                                      ||
    ## ||                                                                            ||
    ## ||                           SRR13978640Aligned.sortedByCoord.out.bam         ||
    ## ||                           SRR13978641Aligned.sortedByCoord.out.bam         ||
    ## ||                           SRR13978642Aligned.sortedByCoord.out.bam         ||
    ## ||                           SRR13978643Aligned.sortedByCoord.out.bam         ||
    ## ||                           SRR13978644Aligned.sortedByCoord.out.bam         ||
    ## ||                           SRR13978645Aligned.sortedByCoord.out.bam         ||
    ## ||                                                                            ||
    ## ||              Paired-end : yes                                              ||
    ## ||        Count read pairs : yes                                              ||
    ## ||              Annotation : Saccharomyces_cerevisiae.R64-1-1.109.gtf (GTF)   ||
    ## ||      Dir for temp files : .                                                ||
    ## ||                 Threads : 4                                                ||
    ## ||                   Level : meta-feature level                               ||
    ## ||      Multimapping reads : counted                                          ||
    ## || Multi-overlapping reads : not counted                                      ||
    ## ||   Min overlapping bases : 1                                                ||
    ## ||                                                                            ||
    ## \\============================================================================//
    ## 
    ## //================================= Running ==================================\\
    ## ||                                                                            ||
    ## || Load annotation file Saccharomyces_cerevisiae.R64-1-1.109.gtf ...          ||
    ## ||    Features : 7507                                                         ||
    ## ||    Meta-features : 7127                                                    ||
    ## ||    Chromosomes/contigs : 17                                                ||
    ## ||                                                                            ||
    ## || Process BAM file SRR13978640Aligned.sortedByCoord.out.bam...               ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3294536                                              ||
    ## ||    Successfully assigned alignments : 2310619 (70.1%)                      ||
    ## ||    Running time : 0.13 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file SRR13978641Aligned.sortedByCoord.out.bam...               ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 2667430                                              ||
    ## ||    Successfully assigned alignments : 2045225 (76.7%)                      ||
    ## ||    Running time : 0.10 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file SRR13978642Aligned.sortedByCoord.out.bam...               ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 5981649                                              ||
    ## ||    Successfully assigned alignments : 2887042 (48.3%)                      ||
    ## ||    Running time : 0.70 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file SRR13978643Aligned.sortedByCoord.out.bam...               ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3319092                                              ||
    ## ||    Successfully assigned alignments : 1975303 (59.5%)                      ||
    ## ||    Running time : 0.15 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file SRR13978644Aligned.sortedByCoord.out.bam...               ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 3719013                                              ||
    ## ||    Successfully assigned alignments : 2830348 (76.1%)                      ||
    ## ||    Running time : 0.15 minutes                                             ||
    ## ||                                                                            ||
    ## || Process BAM file SRR13978645Aligned.sortedByCoord.out.bam...               ||
    ## ||    Paired-end reads are included.                                          ||
    ## ||    Total alignments : 1144778                                              ||
    ## ||    Successfully assigned alignments : 431397 (37.7%)                       ||
    ## ||    Running time : 0.07 minutes                                             ||
    ## ||                                                                            ||
    ## || Write the final count table.                                               ||
    ## || Write the read assignment summary.                                         ||
    ## ||                                                                            ||
    ## \\============================================================================//

Save the count table into a variable and give the column names

``` r
raw_counts <- df$counts
colnames(raw_counts) <- paste0("SRR1397864",0:5) 
```

Create the metadata for samples

``` r
metadata <- tibble(sample = paste0("SRR1397864",0:5), 
                   condition = rep(c("SPRC", "control"), each=3))
```

### DESeqDataSetFromMatrix

Then we make a DESeq2 object

``` r
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ condition )
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

Visualizing the counts

``` r
head(assay(dds)) 
```

    ##         SRR13978640 SRR13978641 SRR13978642 SRR13978643 SRR13978644 SRR13978645
    ## YDL246C           7           7           5           7          12           3
    ## YDL243C          32          34          57          40          44           8
    ## YDR387C         176         154         179         171         202          41
    ## YDL094C           6           6          12           5           8           1
    ## YDR438W          88          92         100          82         148           9
    ## YDR523C          10           9          14           3          10           4

We can filter out the genes with low counts

``` r
keep <- rowSums(counts(dds)) > 10
dds <-  dds[keep,]
```

Estimating size factors and getting normalized counts

``` r
sizeFactors <- estimateSizeFactors(dds)
normalized_counts <- counts(sizeFactors, normalized=TRUE)
```

We transform the data in order to minimize differences between samples
for rows with small counts, and which helps in normalizing the data with
respect to library size.

- We can use 2 types of transformation
  - `vst`
  - `rlog`

#### Variance stabilizing transformation

``` r
vst <-  varianceStabilizingTransformation(dds)
plotPCA(vst)
```

![](README_files/figure-gfm/vst-1.png)<!-- -->

#### Regularized log transformation

``` r
rlog_dds <- rlog(dds)
plotPCA(rlog_dds)
```

![](README_files/figure-gfm/rlog-1.png)<!-- --> As you can see `rlog`
transformation is more robust than `vst`.

### Heatmap of samples

So for plotting the heatmap of sample distances, we transpose the
`rlog_dds` matrix and calculate distances between them

Then we plot the heatmap.

``` r
sampleDist <- dist(t(assay(rlog_dds)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- dds$condition

pheatmap(sampleDistMatrix, 
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         cols = colors)
```

![](README_files/figure-gfm/sample_heatmap-1.png)<!-- -->

### Running DESeq

We finally run `DESEq` on our original deseq object `dds` and not the
transformed ones. Cause the DESeq does the transformation

``` r
res_dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

``` r
res <- results(res_dds)
```

We ran `DESeq` and stored the results into `res`

Here, you can see what the column names of `res` indicate. And also the
`summary` of each column

##### Column names

``` r
mcols(res)
```

    ## DataFrame with 6 rows and 2 columns
    ##                        type            description
    ##                 <character>            <character>
    ## baseMean       intermediate mean of normalized c..
    ## log2FoldChange      results log2 fold change (ML..
    ## lfcSE               results standard error: cond..
    ## stat                results Wald statistic: cond..
    ## pvalue              results Wald test p-value: c..
    ## padj                results   BH adjusted p-values

##### Summary

``` r
summary(res)
```

    ## 
    ## out of 6222 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 86, 1.4%
    ## LFC < 0 (down)     : 67, 1.1%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 1327, 21%
    ## (mean count < 38)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

#### Extracting significant genes

In order to make a **heatmap** we’ll have to set some **thresholds**
first to extract only the top genes. These include both *upregulated*
and *downregulated*

``` r
res_sig <- subset(res, padj < 0.05 & baseMean > 100 & abs(log2FoldChange) > 1)
```

Most *upregulated*

``` r
head(res_sig[order(res_sig$log2FoldChange),])
```

    ## log2 fold change (MLE): condition SPRC vs control 
    ## Wald test p-value: condition SPRC vs control 
    ## DataFrame with 6 rows and 6 columns
    ##          baseMean log2FoldChange     lfcSE      stat      pvalue        padj
    ##         <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
    ## YOR028C   108.834       -1.55976  0.247070  -6.31304 2.73607e-10 5.15118e-08
    ## YMR081C   136.271       -1.33230  0.223228  -5.96832 2.39707e-09 3.55566e-07
    ## YPL111W   251.685       -1.33022  0.251862  -5.28152 1.28115e-07 1.45842e-05
    ## YKL109W   323.179       -1.15495  0.185954  -6.21096 5.26624e-10 9.54751e-08
    ## YIL101C   122.056       -1.14385  0.225322  -5.07649 3.84466e-07 3.84074e-05
    ## YMR011W   229.655       -1.13112  0.221775  -5.10032 3.39071e-07 3.53139e-05

Most *downregulated*

``` r
head(res_sig[order(res_sig$log2FoldChange, decreasing = TRUE),])
```

    ## log2 fold change (MLE): condition SPRC vs control 
    ## Wald test p-value: condition SPRC vs control 
    ## DataFrame with 6 rows and 6 columns
    ##          baseMean log2FoldChange     lfcSE      stat      pvalue        padj
    ##         <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
    ## YER011W   374.369        2.53137  0.205073  12.34374 5.26512e-35 1.78971e-31
    ## YKL001C   229.786        1.96491  0.200473   9.80137 1.11072e-22 7.76710e-20
    ## YJR010W   972.220        1.92640  0.190171  10.12984 4.07297e-24 3.32287e-21
    ## YJR137C   527.502        1.89133  0.183029  10.33355 4.96890e-25 6.08069e-22
    ## YFR030W   633.872        1.82705  0.178149  10.25576 1.11498e-24 1.09156e-21
    ## YLR092W   456.053        1.68802  0.158546  10.64683 1.80402e-26 2.94355e-23

### MA plot

MA plots are commonly used to represent log fold-change versus mean
expression between two treatments They are scatter plots with **base-2
log fold-change** along the *y-axis* and **normalized mean expression**
along the *x-axis*

First we have to **shrink** the log2fold changes before plotting.

##### What is shrinkage?

In short, it looks at the largest fold changes that are not due to low
counts and uses these to inform a prior distribution. So the large fold
changes from genes with lots of statistical information are not shrunk,
while the imprecise fold changes are shrunk.

This allows you to compare all estimated LFC across experiments.
`DESeq2` does this by default. We’re just doing it to make the MA plot

``` r
res_shrink <- lfcShrink(res_dds, coef = "condition_SPRC_vs_control", type = "apeglm")
```

    ## using 'apeglm' for LFC shrinkage. If used in published research, please cite:
    ##     Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
    ##     sequence count data: removing the noise and preserving large differences.
    ##     Bioinformatics. https://doi.org/10.1093/bioinformatics/bty895

``` r
plotMA(res_shrink, ylim = c(-3, 3))
```

![](README_files/figure-gfm/MA_plot-1.png)<!-- -->

``` r
gene_ids_sig <- rownames(res_sig)
gene_ids <- rownames(res)
```

### bioMart

``` r
head(listEnsembl())
```

    ##         biomart                version
    ## 1         genes      Ensembl Genes 109
    ## 2 mouse_strains      Mouse strains 109
    ## 3          snps  Ensembl Variation 109
    ## 4    regulation Ensembl Regulation 109

``` r
ensembl <- useEnsembl(biomart = "genes")
```

Selecting the **dataset**

``` r
head(listDatasets(ensembl))
```

    ##                        dataset                           description
    ## 1 abrachyrhynchus_gene_ensembl Pink-footed goose genes (ASM259213v1)
    ## 2     acalliptera_gene_ensembl      Eastern happy genes (fAstCal1.2)
    ## 3   acarolinensis_gene_ensembl       Green anole genes (AnoCar2.0v2)
    ## 4    acchrysaetos_gene_ensembl       Golden eagle genes (bAquChr1.2)
    ## 5    acitrinellus_gene_ensembl        Midas cichlid genes (Midas_v5)
    ## 6    amelanoleuca_gene_ensembl       Giant panda genes (ASM200744v2)
    ##       version
    ## 1 ASM259213v1
    ## 2  fAstCal1.2
    ## 3 AnoCar2.0v2
    ## 4  bAquChr1.2
    ## 5    Midas_v5
    ## 6 ASM200744v2

``` r
bio.mart <- useMart("ensembl", "scerevisiae_gene_ensembl")
```

Listing the **attributes** and **filters**

``` r
head(listAttributes(bio.mart))
```

    ##                    name              description         page
    ## 1       ensembl_gene_id           Gene stable ID feature_page
    ## 2 ensembl_transcript_id     Transcript stable ID feature_page
    ## 3    ensembl_peptide_id        Protein stable ID feature_page
    ## 4       ensembl_exon_id           Exon stable ID feature_page
    ## 5           description         Gene description feature_page
    ## 6       chromosome_name Chromosome/scaffold name feature_page

``` r
head(listFilters(bio.mart))
```

    ##                 name                            description
    ## 1    chromosome_name               Chromosome/scaffold name
    ## 2              start                                  Start
    ## 3                end                                    End
    ## 4             strand                                 Strand
    ## 5 chromosomal_region e.g. 1:100:10000:-1, 1:100000:200000:1
    ## 6        with_chembl                      With ChEMBL ID(s)

Getting the data from bioMart

``` r
bm_attr <- c("ensembl_gene_id", "external_gene_name", 
           "chromosome_name", "start_position", "end_position",
           "description")

bm_res <- getBM(attributes = bm_attr, 
      filters = c("ensembl_gene_id"), 
      values = gene_ids, mart = bio.mart)

bm_res_sig <- getBM(attributes = bm_attr, 
                filters = c("ensembl_gene_id"), 
                values = gene_ids_sig, mart = bio.mart)

bm_res_sig
```

    ##    ensembl_gene_id external_gene_name chromosome_name start_position
    ## 1          YAL061W               BDH2               I          33448
    ## 2          YBR213W               MET8              II         650368
    ## 3          YCL027W               FUS1             III          71803
    ## 4          YDR453C               TSA2              IV        1365072
    ## 5          YER011W               TIR1               V         175248
    ## 6          YER042W               MXR1               V         234937
    ## 7          YER069W             ARG5,6               V         295410
    ## 8          YER091C               MET6               V         339864
    ## 9          YFR030W              MET10              VI         213312
    ## 10         YGR060W              ERG25             VII         610564
    ## 11         YGR177C               ATF2             VII         848829
    ## 12         YHR007C              ERG11            VIII         120091
    ## 13         YIL011W               TIR3              IX         333727
    ## 14         YIL074C              SER33              IX         221081
    ## 15         YIL101C               XBP1              IX         175307
    ## 16         YJL088W               ARG3               X         268799
    ## 17         YJR010W               MET3               X         456239
    ## 18         YJR137C               MET5               X         678957
    ## 19         YKL001C              MET14              XI         438777
    ## 20       YKL068W-A                                 XI         309206
    ## 21         YKL109W               HAP4              XI         232227
    ## 22         YKR069W               MET1              XI         571612
    ## 23         YLL061W               MMP1             XII          17956
    ## 24         YLL062C               MHT1             XII          16639
    ## 25         YLR056W               ERG3             XII         253861
    ## 26         YLR092W               SUL2             XII         323544
    ## 27         YLR303W              MET17             XII         732542
    ## 28         YML018C                               XIII         234771
    ## 29         YMR011W               HXT2            XIII         288079
    ## 30         YMR015C               ERG5            XIII         300869
    ## 31         YMR081C               ISF1            XIII         430079
    ## 32         YOL140W               ARG8              XV          58759
    ## 33         YOR010C               TIR2              XV         346195
    ## 34         YOR028C               CIN5              XV         383533
    ## 35         YPL111W               CAR1             XVI         339944
    ## 36         YPL274W               SAM3             XVI          22938
    ## 37         YPR167C              MET16             XVI         876847
    ##    end_position
    ## 1         34701
    ## 2        651192
    ## 3         73341
    ## 4       1365662
    ## 5        176012
    ## 6        235491
    ## 7        298001
    ## 8        342167
    ## 9        216419
    ## 10       611493
    ## 11       850436
    ## 12       121683
    ## 13       334536
    ## 14       222490
    ## 15       177250
    ## 16       269815
    ## 17       457774
    ## 18       683285
    ## 19       439385
    ## 20       309442
    ## 21       233891
    ## 22       573393
    ## 23        19707
    ## 24        17613
    ## 25       254958
    ## 26       326225
    ## 27       733876
    ## 28       235952
    ## 29       289704
    ## 30       302485
    ## 31       431095
    ## 32        60030
    ## 33       346950
    ## 34       384420
    ## 35       340945
    ## 36        24701
    ## 37       877632
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    description
    ## 1                                                                                                                                                                                                                                                                                                                                                                Putative medium-chain alcohol dehydrogenase with similarity to BDH1; transcription induced by constitutively active PDR1 and PDR3 [Source:SGD;Acc:S000000057]
    ## 2                                                                                                                                                                                                                                                                                             Bifunctional dehydrogenase and ferrochelatase; involved in the biosynthesis of siroheme, a prosthetic group used by sulfite reductase; required for sulfate assimilation and methionine biosynthesis [Source:SGD;Acc:S000000417]
    ## 3                                                                                                                                                                                                                                                               Membrane protein localized to the shmoo tip; required for cell fusion; expression regulated by mating pheromone; proposed to coordinate signaling, fusion, and polarization events required for fusion; potential Cdc28p substrate [Source:SGD;Acc:S000000532]
    ## 4                                                                                                                Stress inducible cytoplasmic thioredoxin peroxidase; cooperates with Tsa1p in the removal of reactive oxygen, nitrogen and sulfur species using thioredoxin as hydrogen donor; deletion enhances the mutator phenotype of tsa1 mutants; protein abundance increases in response to DNA replication stress; TSA2 has a paralog, TSA1, that arose from the whole genome duplication [Source:SGD;Acc:S000002861]
    ## 5                                                                                                                                                                                                                                                          Cell wall mannoprotein; expression is downregulated at acidic pH and induced by cold shock and anaerobiosis; abundance is increased in cells cultured without shaking; member of the Srp1p/Tip1p family of serine-alanine-rich proteins [Source:SGD;Acc:S000000813]
    ## 6                                                                                                                                                                                                                             Methionine-S-sulfoxide reductase; involved in the response to oxidative stress; protects iron-sulfur clusters from oxidative inactivation along with MXR2; involved in the regulation of lifespan; reduced activity of human homolog implicated in Alzheimer disease [Source:SGD;Acc:S000000844]
    ## 7  Acetylglutamate kinase and N-acetyl-gamma-glutamyl-phosphate reductase; N-acetyl-L-glutamate kinase (NAGK) catalyzes the 2nd and N-acetyl-gamma-glutamyl-phosphate reductase (NAGSA), the 3rd step in arginine biosynthesis; synthesized as a precursor which is processed in the mitochondrion to yield mature NAGK and NAGSA; enzymes form a metabolon complex with Arg2p; NAGK C-terminal domain stabilizes the enzymes, slows catalysis and is involved in feed-back inhibition by arginine [Source:SGD;Acc:S000000871]
    ## 8                                                                                                                                                                                                                                                                                Cobalamin-independent methionine synthase; involved in methionine biosynthesis and regeneration; requires a minimum of two glutamates on the methyltetrahydrofolate substrate, similar to bacterial metE homologs [Source:SGD;Acc:S000000893]
    ## 9                                                                                                                                                                                                                                                                                                                                                                                                           Subunit alpha of assimilatory sulfite reductase; complex converts sulfite into sulfide [Source:SGD;Acc:S000001926]
    ## 10                                                                                                                                                                      C-4 methyl sterol oxidase; catalyzes the first of three steps required to remove two C-4 methyl groups from an intermediate in ergosterol biosynthesis; mutants accumulate the sterol intermediate 4,4-dimethylzymosterol; human MSMO1 functionally complements the growth defect caused by repression of ERG25 expression [Source:SGD;Acc:S000003292]
    ## 11                                                                                                                                                                                                                                                                                                                                  Alcohol acetyltransferase; may play a role in steroid detoxification; forms volatile esters during fermentation, which is important for brewing and winemaking [Source:SGD;Acc:S000003409]
    ## 12                                                                  Lanosterol 14-alpha-demethylase; catalyzes C-14 demethylation of lanosterol to form 4,4''-dimethyl cholesta-8,14,24-triene-3-beta-ol in ergosterol biosynthesis pathway; transcriptionally down-regulated when ergosterol is in excess; member of cytochrome P450 family; associated and coordinately regulated with the P450 reductase Ncp1p; human CYP51A1 functionally complements the lethality of the erg11 null mutation [Source:SGD;Acc:S000001049]
    ## 13                                                                                                                                                                                                                                                            Cell wall mannoprotein; member of Srp1p/Tip1p family of serine-alanine-rich proteins; expressed under anaerobic conditions and required for anaerobic growth; TIR3 has a paralog, TIR2, that arose from the whole genome duplication [Source:SGD;Acc:S000001273]
    ## 14                                                                                                                 3-phosphoglycerate dehydrogenase and alpha-ketoglutarate reductase; 3PG dehydrogenase that catalyzes the first step in serine and glycine biosynthesis; also functions as an alpha-ketoglutarate reductase, converting alpha-ketoglutarate to D-2-hydroxyglutarate (D-2HG); localizes to the cytoplasm; SER33 has a paralog, SER3, that arose from the whole genome duplication [Source:SGD;Acc:S000001336]
    ## 15            Transcriptional repressor; binds promoter sequences of cyclin genes, CYS3, and SMF2; not expressed during log phase of growth, but induced by stress or starvation during mitosis, and late in meiosis; represses 15% of all yeast genes as cells transition to quiescence; important for maintaining G1 arrest and for longevity of quiescent cells; member of Swi4p/Mbp1p family; phosphorylated by Cdc28p; relative distribution to nucleus increases upon DNA replication stress [Source:SGD;Acc:S000001363]
    ## 16                                                                                                                                                                                                                                                                                                                              Ornithine carbamoyltransferase; also known as carbamoylphosphate:L-ornithine carbamoyltransferase; catalyzes the biosynthesis of the arginine precursor citrulline [Source:SGD;Acc:S000003624]
    ## 17                                                                                                                                                                                                                                                              ATP sulfurylase; catalyzes the primary step of intracellular sulfate activation, essential for assimilatory reduction of sulfate to sulfide, involved in methionine metabolism; human homolog PAPSS2 complements yeast null mutant [Source:SGD;Acc:S000003771]
    ## 18                                                                                                                                                                                                                                                                                                                                                                                      Sulfite reductase beta subunit; involved in amino acid biosynthesis, transcription repressed by methionine [Source:SGD;Acc:S000003898]
    ## 19                                                                                                                                                                                                                                                                                                                                             Adenylylsulfate kinase; required for sulfate assimilation and involved in methionine metabolism; human homolog PAPSS2 complements yeast null mutant [Source:SGD;Acc:S000001484]
    ## 20                                                                                                                                                                                                                                                                                                                                                                                                                 Putative protein of unknown function; identified by homology to Ashbya gossypii [Source:SGD;Acc:S000028524]
    ## 21                                                                                                                                                                                                              Transcription factor; subunit of the heme-activated, glucose-repressed Hap2p/3p/4p/5p CCAAT-binding complex, a transcriptional activator and global regulator of respiratory gene expression; provides the principal activation function of the complex; involved in diauxic shift [Source:SGD;Acc:S000001592]
    ## 22                                                                                                                                                                                                                                                                              S-adenosyl-L-methionine uroporphyrinogen III transmethylase; involved in the biosynthesis of siroheme, a prosthetic group used by sulfite reductase; required for sulfate assimilation and methionine biosynthesis [Source:SGD;Acc:S000001777]
    ## 23                                                                                                                                                                                                                                                                                                                             High-affinity S-methylmethionine permease; required for utilization of S-methylmethionine as a sulfur source; has similarity to S-adenosylmethionine permease Sam3p [Source:SGD;Acc:S000003984]
    ## 24                                                                                                                                                                                                                                                                                                           S-methylmethionine-homocysteine methyltransferase; functions along with Sam4p in the conversion of S-adenosylmethionine (AdoMet) to methionine to control the methionine/AdoMet ratio [Source:SGD;Acc:S000003985]
    ## 25                                                                                                                C-5 sterol desaturase; glycoprotein that catalyzes the introduction of a C-5(6) double bond into episterol, a precursor in ergosterol biosynthesis; transcriptionally down-regulated when ergosterol is in excess; mutants are viable, but cannot grow on non-fermentable carbon sources; substrate of HRD ubiquitin ligase; mutation is functionally complemented by human SC5D [Source:SGD;Acc:S000004046]
    ## 26                                                                                                                                                                                                                                                                                                      High affinity sulfate permease; sulfate uptake is mediated by specific sulfate transporters Sul1p and Sul2p, which control the concentration of endogenous activated sulfate intermediates [Source:SGD;Acc:S000004082]
    ## 27                                                                                                                                                                                                                                                                                                                                                                                            O-acetyl homoserine-O-acetyl serine sulfhydrylase; required for Methionine and cysteine biosynthesis [Source:SGD;Acc:S000004294]
    ## 28                                                                                                Protein of unknown function; green fluorescent protein (GFP)-fusion protein localizes to the membrane of the vacuole; physical interaction with Atg27p suggests a possible role in autophagy; YML018C is not an essential gene; relative distribution to the vacuolar membrane decreases upon DNA replication stress; YML018C has a paralog, THI74, that arose from the whole genome duplication [Source:SGD;Acc:S000004480]
    ## 29                                                                                                                                                                                                                                                                                                                                  High-affinity glucose transporter of the major facilitator superfamily; expression is induced by low levels of glucose and repressed by high levels of glucose [Source:SGD;Acc:S000004613]
    ## 30                                                                                                                                                                                                                                                                                        C-22 sterol desaturase; a cytochrome P450 enzyme that catalyzes the formation of the C-22(23) double bond in the sterol side chain in ergosterol biosynthesis; may be a target of azole antifungal drugs [Source:SGD;Acc:S000004617]
    ## 31                                                                                                                                                                                                                               Serine-rich, hydrophilic protein; overexpression suppresses growth defects of hap2, hap3, and hap4 mutants; expression is under glucose control; cotranscribed with NAM7 in a cyp1 mutant; ISF1 has a paralog, MBR1, that arose from the whole genome duplication [Source:SGD;Acc:S000004686]
    ## 32                                                                                                                                                                                                                                                                                                                                                                             Acetylornithine aminotransferase; catalyzes the fourth step in the biosynthesis of the arginine precursor ornithine [Source:SGD;Acc:S000005500]
    ## 33                                                                                                                                                                                                                                                              Putative cell wall mannoprotein; member of the Srp1p/Tip1p family of serine-alanine-rich proteins; transcription is induced by cold shock and anaerobiosis; TIR2 has a paralog, TIR3, that arose from the whole genome duplication [Source:SGD;Acc:S000005536]
    ## 34                                                                                        Basic leucine zipper (bZIP) transcription factor of the yAP-1 family; physically interacts with the Tup1-Cyc8 complex and recruits Tup1p to its targets; mediates pleiotropic drug resistance and salt tolerance; nuclearly localized under oxidative stress and sequestered in the cytoplasm by Lot6p under reducing conditions; CIN5 has a paralog, YAP6, that arose from the whole genome duplication [Source:SGD;Acc:S000005554]
    ## 35                                                                                                                                                                                                                           Arginase, catabolizes arginine to ornithine and urea; expression responds to both induction by arginine and nitrogen catabolite repression; disruption decreases production of carcinogen ethyl carbamate during wine fermentation and also enhances freeze tolerance [Source:SGD;Acc:S000006032]
    ## 36                                                                                                                                                                                                                                                                                                                           High-affinity S-adenosylmethionine permease; required for utilization of S-adenosylmethionine as a sulfur source; has similarity to S-methylmethionine permease Mmp1p [Source:SGD;Acc:S000006195]
    ## 37                                                                                                                                                                                                                                                                     3'-phosphoadenylsulfate reductase; reduces 3'-phosphoadenylyl sulfate to adenosine-3',5'-bisphosphate and free sulfite using reduced thioredoxin as cosubstrate, involved in sulfate assimilation and methionine metabolism [Source:SGD;Acc:S000006371]

#### Most significant Differentially Expressed Genes

``` r
bm_res_sig$external_gene_name
```

    ##  [1] "BDH2"   "MET8"   "FUS1"   "TSA2"   "TIR1"   "MXR1"   "ARG5,6" "MET6"  
    ##  [9] "MET10"  "ERG25"  "ATF2"   "ERG11"  "TIR3"   "SER33"  "XBP1"   "ARG3"  
    ## [17] "MET3"   "MET5"   "MET14"  ""       "HAP4"   "MET1"   "MMP1"   "MHT1"  
    ## [25] "ERG3"   "SUL2"   "MET17"  ""       "HXT2"   "ERG5"   "ISF1"   "ARG8"  
    ## [33] "TIR2"   "CIN5"   "CAR1"   "SAM3"   "MET16"

### Heatmap

Preparation for the Heatmap

``` r
sig_matrix <- normalized_counts[rownames(res_sig),]
sig_matrix.z <- t(apply(sig_matrix, 1, scale))
rownames(sig_matrix.z) <- bm_res_sig$external_gene_name
colnames(sig_matrix.z) <- colnames(sig_matrix)
```

``` r
heat_colors <- colorRampPalette(brewer.pal(6, "PuOr"))(100)
```

``` r
pheatmap(mat = sig_matrix.z, 
         cluster_rows = T, 
         cluster_cols = T, 
         labels_col = colnames(sig_matrix.z), 
         labels_row = rownames(sig_matrix.z),
         color = heat_colors,
         height = 40)
```

![](README_files/figure-gfm/heatmap-1.png)<!-- -->

Preparation for the **Volcano plot**

``` r
res_tb <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```

Getting the **gene name** instead of ensembl ids

``` r
res_tb <- res_tb %>% 
  mutate(gene_labels = bm_res$external_gene_name) 
```

Setting the **threshold**

``` r
res_tb <- res_tb %>% 
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) > 1) 
```

### Volcano plot

``` r
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
```

    ## Warning: Removed 1327 rows containing missing values (`geom_point()`).

    ## Warning: Removed 1327 rows containing missing values (`geom_text_repel()`).

    ## Warning: ggrepel: 3950 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](README_files/figure-gfm/volcano_plot-1.png)<!-- -->
