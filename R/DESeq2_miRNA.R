setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

## differential expression analysis for 14 TCGA cancer types
# condition: normal/Tumor samples >= 10

library(DESeq2)
library(TCGAbiolinks)
library("BiocParallel")
#register(MulticoreParam(4))

## load the pancancer miRNA count matrix
load("Normal_Tumor_panTCGA_isomir.RData")
dim(panisomir)

## remove the samples which are excludde in the study
# 1. load the sample information object
load("sampleinfo.RData")
head(merged_df)
dim(merged_df)
rownames(merged_df) <- 1:nrow(merged_df)

is.na(match(colnames(panisomir), merged_df$sample_id))

## 2. subset the smaple info table with only 14 cancer types
projects <- c("TCGA-ESCA",
              "TCGA-BLCA",
              "TCGA-KICH",
              "TCGA-UCEC",
              "TCGA-KIRP",
              "TCGA-HNSC",
              "TCGA-LUSC",
              "TCGA-STAD",
              "TCGA-LUAD",
              "TCGA-LIHC",
              "TCGA-PRAD",
              "TCGA-THCA",
              "TCGA-KIRC",
              "TCGA-BRCA")
TCGA14cancers <- merged_df[merged_df$project %in% projects,]
TCGA14cancers
dim(TCGA14cancers)

rownames(TCGA14cancers) <- 1:nrow(TCGA14cancers)
rownames(TCGA14cancers) <- TCGA14cancers$sample_id
head(TCGA14cancers)

# make additional condition column, example TCGA-BRCA-TP, TCGA-BRCA-NT
paste(TCGA14cancers$project, TCGA14cancers$tissue_name, sep = "-")
TCGA14cancers$condition <- paste(TCGA14cancers$project, TCGA14cancers$tissue_name, sep = "-")
head(TCGA14cancers)

#3. now take subset of these 14 cancer types from the pancancer isomir count matrix
count14cancers <- panisomir[,colnames(panisomir) %in% TCGA14cancers$sample_id]
dim(count14cancers)

#4. now rearrange the columns of count matrix in same order as sample info table
count14cancers <- count14cancers[,rownames(TCGA14cancers)]
dim(count14cancers)
count14cancers[1:5,1:5]

match(colnames(count14cancers), rownames(TCGA14cancers))
all(rownames(TCGA14cancers) == colnames(count14cancers))

# 5. run DESeq2 analysis
table(is.na(count14cancers))
count14cancers[1:5,1:5]
count14cancers <- as.matrix(count14cancers)

count14cancers[is.na(count14cancers)] <- 0

dds <- DESeqDataSetFromMatrix(countData = count14cancers,
                              colData = TCGA14cancers,
                              design = ~ condition)
dds

dds$condition

## differential expression analysis
dds <- DESeq(dds)
dds
assay(dds)[1:5,1:5]

## save the deseq object created for further analysis
save(dds, file = "DESeq_run.RData")

## save the normalized counts
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd
assay(vsd)[1:5,1:5]
colData(vsd)
save(vsd, file = "normalized_counts_14_Cancer.RData")

## results for TCGA-BRCA-NT vs TCGA-BRCA-TP
res <- results(dds, contrast = c("condition", "TCGA-BRCA-NT", "TCGA-BRCA-TP"))
res

resOrdered <- res[order(res$padj),]
resOrdered

summary(res)
sum(res$padj < 0.05, na.rm=TRUE)








