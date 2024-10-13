# Breast Cancer Data Preprocessing Script
# This script preprocesses breast cancer RNA-Seq data from TCGA, 
# performs gene annotation, and applies normalization techniques 
# to prepare the data for co-expression network analysis.

### 01. Load required libraries
library(BiocManager)
library(SummarizedExperiment)
library(TCGAbiolinks)
library(EDASeq)
library(dplyr)
library(NOISeq)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(ggbiplot)

### 02. Gene annotation using Ensembl (biomaRt package)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
features <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", 
              "hgnc_symbol", "percentage_gene_gc_content", "gene_biotype")
chrs <- c(1:22, "X", "Y")
annot <- getBM(attributes = features, filters = "chromosome_name", values = chrs, mart = ensembl)
colnames(annot) <- c("ensembl_gene_id", "Chr", "Start", "End", "HGNC_symbol", "GC", "Type")
annot$Length <- abs(annot$End - annot$Start)

### 03. Query RNA-Seq data from TCGA (BRCA)
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
samplesDown <- getResults(query, cols = c("cases"))
half <- samplesDown[1:553]
other_half <- samplesDown[554:1231]

queryDown_half <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                           data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", barcode = half)
queryDown_other_half <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                                 data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts", barcode = other_half)
GDCdownload(query = queryDown_half)
GDCdownload(query = queryDown_other_half)
dataPrep1_1 <- GDCprepare(query = queryDown_half)
dataPrep1_2 <- GDCprepare(query = queryDown_other_half)

### 04. Filter Primary Tumor and Normal samples
BC1.T <- dataPrep1_1[, dataPrep1_1$sample_type == "Primary Tumor"]
BC1.N <- dataPrep1_1[, dataPrep1_1$sample_type == "Solid Tissue Normal"]
BC2.T <- dataPrep1_2[, dataPrep1_2$sample_type == "Primary Tumor"]
BC2.N <- dataPrep1_2[, dataPrep1_2$sample_type == "Solid Tissue Normal"]

### 05. Tumor purity estimation (if necessary)
Purity.BRCA.1 <- TCGAtumor_purity(colnames(BC1.T), 0,0,0,0,0)$pure_barcodes
Purity.BRCA.2 <- TCGAtumor_purity(colnames(BC2.T), 0,0,0,0,0)$pure_barcodes

### 06. Exclude samples based on molecular subtypes
diff.1 <- setdiff(Purity.BRCA.1, TCGA_MolecularSubtype(half)$filtered)
diff.2 <- setdiff(Purity.BRCA.2, TCGA_MolecularSubtype(other_half)$filtered)

### 07. Combine tumor and normal samples
rnas <- cbind(assay(BC1.T)[,diff.1], assay(BC2.T)[,diff.2], assay(BC1.N), assay(BC2.N))

### 08. Identify molecular subtypes
mol_subtypes.1 <- TCGA_MolecularSubtype(colnames(BC1.T)[,diff.1])$subtypes$subtype
mol_subtypes.2 <- TCGA_MolecularSubtype(colnames(BC2.T)[,diff.2])$subtypes$subtype
normal <- data.frame(subtype = rep('normal', 113))
mol_subtypes <- rbind(data.frame(subtype = mol_subtypes.1), data.frame(subtype = mol_subtypes.2), normal)

### 09. Prepare factors for differential analysis
factorBC <- data.frame(Group = "PT", Sample = c(colnames(assay(BC1.T)[,diff.1]), colnames(assay(BC2.T)[,diff.2])))
factorsNormalBC <- data.frame(Group = "Normal", Sample = c(colnames(BC1.N), colnames(BC2.N)))
factors <- rbind(factorBC, factorsNormalBC)
factors <- cbind(factors, mol_subtypes)
rownames(factors) <- factors$Sample

### 10. Filter expression data based on quantile
dataFilt <- TCGAanalyze_Filtering(tabDF = rnas, method = "quantile", qnt.cut = 0.25)
ridx <- rowSums(dataFilt == 0) <= round(dim(rnas)[2]/2)
dataFilt <- dataFilt[ridx, ]
ridx <- rowMeans(dataFilt) >= 50
dataFilt <- dataFilt[ridx, ]
rnas <- rnas[rownames(rnas) %in% rownames(dataFilt), ]

### 11. Match gene IDs and prepare data
rownames(rnas) <- gsub("\\..*", "", rownames(rnas))
inter <- intersect(rownames(rnas), annot$ensembl_gene_id)
rnas1 <- rnas[rownames(rnas) %in% inter, ]
annot1 <- annot[annot$ensembl_gene_id %in% inter, ]
annot1 <- annot1[!duplicated(annot1$ensembl_gene_id), ]

### 12. Normalize expression data
ln.data <- withinLaneNormalization(rnas1, annot1$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data, annot1$GC, which = "full")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full")
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = Ready_factors)
mydata2corr1 <- NOISeq::ARSyNseq(noiseqData, norm = "n", logtransf = FALSE)
rnas2 <- exprs(mydata2corr1)

### 13. Quality control (PCA before and after normalization)
before.pca <- prcomp(t(rnas1), center = TRUE, scale. = TRUE)
ggbiplot(before.pca, var.axes = FALSE, ellipse = TRUE, groups = factors$Group)

after.pca <- prcomp(t(rnas2), center = TRUE, scale. = TRUE)
ggbiplot(after.pca, var.axes = FALSE, ellipse = TRUE, groups = factors$Group)








































  
