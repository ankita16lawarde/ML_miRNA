#differential expression mRNAs/and other RNAs 
## differential expression lncRNAs

## differential expression analysis on RNA-Seq data
setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

library(DESeq2)
library(TCGAbiolinks)


## load the sampleinfo collected for miRNA data
# 1. load the sample information object
load("sampleinfo.RData")

### load RNASeq data
#KICH, BLCA, BRCA, ESCA, HNSC, KIRC, KIRP, LIHC, LUAD, LUSC, PRAD, STAD, THCA, UCEC

project <- c("BRCA", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

## loop over each prject id and save its DESeq object and run from DESeq
for (proj in project) {
  
  load(file = paste0("RNASeq_Data/TCGA-", proj, "RNASeq.RData"))
  # make two rangedsummarizedexperiment objects from data: lncRNA and protein-coding genes and remove miRNA genes
  # extract ensembl id info using rowRanges()
  ##########
  protein_coding <- rowRanges(data)[rowRanges(data)$gene_type == "protein_coding",]
  lncRNA <- rowRanges(data)[rowRanges(data)$gene_type == "lncRNA",]
  protein_coding <- subset(data, rownames(data) %in% protein_coding$gene_id)
  lncRNA <- subset(data, rownames(data) %in% lncRNA$gene_id)
  ###########################################################################################
  # run DESeq on protein coding
  dds_proteincoding <- DESeqDataSet(protein_coding, design = ~shortLetterCode)
  ## differential expression analysis
  dds_proteincoding <- DESeq(dds_proteincoding)
  save(dds_proteincoding, file = paste0(proj ,"_dds_proteincoding.RData"))
  ## save the normalized counts
  vsd <- varianceStabilizingTransformation(dds_proteincoding, blind=FALSE)
  save(vsd, file = paste0(proj, "_normalized_counts.RData"))
  ################################################################################
  #run DESeq fro lncrNA data
  dds_lncRNA <- DESeqDataSet(lncRNA, design = ~shortLetterCode)
  ## differential expression analysis
  dds_lncRNA <- DESeq(dds_lncRNA)
  save(dds_lncRNA, file = paste0(proj, "_dds_lncRNA.RData"))
  ## save the normalized counts
  vsd_lncRNA <- varianceStabilizingTransformation(dds_lncRNA, blind=FALSE)
  save(vsd_lncRNA, file = paste0(proj, "_normalized_lncRNAcounts.RData"))
}

vsd_lncRNA
rowData(vsd_lncRNA)


###############################################################################################################

load("RNASeq_Data/TCGA-BLCARNASeq.RData")
data


load("RNASeq_Data/TCGA-UCECRNASeq.RData")
data

colData(data)[1:5,]

dim(colData(data)[colData(data)$shortLetterCode == "TP",])


# make two rangedsummarizedexperiment objects from data: lncRNA and protein-coding genes and remove miRNA genes
# extract ensembl id info using rowRanges()
table(rowRanges(data)$gene_type)
#protein_coding                          
#19962
#lncRNA 
#16901

head(rowRanges(data)[rowRanges(data)$gene_type == "miRNA",])

##########
protein_coding <- rowRanges(data)[rowRanges(data)$gene_type == "protein_coding",]
table(protein_coding$gene_type)


lncRNA <- rowRanges(data)[rowRanges(data)$gene_type == "lncRNA",]
table(lncRNA$gene_type)

protein_coding <- subset(data, rownames(data) %in% protein_coding$gene_id)
protein_coding


lncRNA <- subset(data, rownames(data) %in% lncRNA$gene_id)
lncRNA
###########################################################################################
colData(data)[1:5,]

## run DESeq on protein coding
dds_proteincoding <- DESeqDataSet(protein_coding, design = ~shortLetterCode)

## differential expression analysis
dds_proteincoding <- DESeq(dds_proteincoding)
dds_proteincoding

save(dds_proteincoding, file = "dds_proteincoding_TCGA-BLCA.RData")

## save the normalized counts
vsd <- varianceStabilizingTransformation(dds_proteincoding, blind=FALSE)
vsd
assay(vsd)[1:5,1:5]
colData(vsd)
save(vsd, file = "normalized_counts_TCGA-BLCA.RData")

## results for TCGA-BRCA-NT vs TCGA-BRCA-TP
res <- results(dds_proteincoding, contrast = c("shortLetterCode", "NT", "TP"))
res

resOrdered <- res[order(res$padj),]
resOrdered

summary(res)
sum(res$padj < 0.05, na.rm=TRUE)


#################################################################################
#run DESeq fro lncrNA data
lncRNA
dds_lncRNA <- DESeqDataSet(lncRNA, design = ~shortLetterCode )

## differential expression analysis
dds_lncRNA <- DESeq(dds_lncRNA)

save(dds_lncRNA, file = "dds_lncRNA_TCGA-BLCA.RData")

## save the normalized counts
vsd_lncRNA <- varianceStabilizingTransformation(dds_lncRNA, blind=FALSE)

assay(vsd)[1:5,1:5]
colData(vsd)
save(vsd_lncRNA, file = "normalized_lncRNAcounts_TCGA-BLCA.RData")

## results for TCGA-BRCA-NT vs TCGA-BRCA-TP
res <- results(dds_lncRNA, contrast = c("shortLetterCode", "NT", "TP"))
res

resOrdered <- res[order(res$padj),]
resOrdered

summary(res)
sum(res$padj < 0.05, na.rm=TRUE)


