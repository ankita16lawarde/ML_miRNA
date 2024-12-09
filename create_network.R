## network creation for each cancer type
# make network of miRNA, mRNA and lncRNAs with padjusted <0.05 & abs(lfc) > 1
# calculate the assortativity coeffecient
# save the list of interacting miRNAs for model building

setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

library(DESeq2)
library(TCGAbiolinks)
library("BiocParallel")
library(igraph)

## loop over for each cancer type

project <- c("BLCA" ,"BRCA", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

#1. load miRNA deseq object for a cancertype
load(file = "DESeq_run.RData")
#dds

# 1. load the sample information object
load("sampleinfo.RData")
rownames(merged_df) <- 1:nrow(merged_df)
merged_df[1:5,1:5]

for (proj in project) {
  ## to create network
  proj
  #paste0("RNASeq_Data/TCGA-", proj, "RNASeq.RData")
  #BLCA
  ## extract results from dds for cancer type tumor - normal comparison
  res_mirna <- results(dds, contrast = c("condition", paste0("TCGA-",proj,"-TP"), paste0("TCGA-",proj,"-NT")))
  #res_mirna
  
  ## take subset of miRNAs with cutoff mentioned before
  res_mirna <- subset(res_mirna, abs(res_mirna$log2FoldChange) > 1 & res_mirna$padj < 0.05)
  #res_mirna
  
  ## take subset of count matrix of BLCA from the toal count matrix of 14 cancer
  #merged_df contains sample info
  selcancer <- merged_df[merged_df$project == paste0("TCGA-", proj),]
  #dim(selcancer)
  # load the normalized count matrix of 14 cancers
  # load the normalized count matrix of 14 cancers
  load("normalized_counts_14_Cancer.RData")
  #vsd
  #dim(vsd)
  #colnames(vsd)
  
  # only BLCA samples
  vsd_mirna <- assay(vsd)[,selcancer$sample_id]
  #dim(vsd_mirna)
  
  ## now subset the dysregulated miRNAS from the matrix
  vsd_mirna <- vsd_mirna[rownames(res_mirna),]
  
  ########################################################################################
  # load mRNA and lncRNA deseq run
  
  # for mRNA
  load(paste0(proj,"_dds_proteincoding.RData"))
  # change rownames to gene symbol: gene:name
  rownames(dds_proteincoding) <- rowData(dds_proteincoding)$gene_name
  
  res_proteincoding <- results(dds_proteincoding, contrast = c("shortLetterCode", "TP", "NT"))
  res_proteincoding <- subset(res_proteincoding, abs(res_proteincoding$log2FoldChange) > 1 & res_proteincoding$padj < 0.05)
  
  ## load normalized count matrix
  load(paste0(proj,"_normalized_counts.RData"))
  rownames(vsd) <- rowData(vsd)$gene_name
  vsd_mrna <- assay(vsd)
  vsd_mrna <- vsd_mrna[rownames(res_proteincoding),] 
  #dim(vsd_mrna)
  
  ## for lncRNA
  load(paste0(proj,"_dds_lncRNA.RData"))
  #dds_lncRNA
  #change rownames to gene_name
  rownames(dds_lncRNA) <- rowData(dds_lncRNA)$gene_name
  
  res_lncRNA <- results(dds_lncRNA, contrast = c("shortLetterCode", "TP", "NT"))
  res_lncRNA <- subset(res_lncRNA, abs(res_lncRNA$log2FoldChange) > 1 & res_lncRNA$padj < 0.05)
  #res_lncRNA
  
  # load lncRNA count matrix
  load(paste0(proj,"_normalized_lncRNAcounts.RData"))
  rownames(vsd_lncRNA) <- rowData(vsd_lncRNA)$gene_name
  vsd_lncRNA <- assay(vsd_lncRNA)
  vsd_lncRNA <- vsd_lncRNA[rownames(res_lncRNA),]
  #dim(vsd_lncRNA)
  
  #########################################################################################
  ## to create network only select common samples from miRNA and mRNA, lncRNA count matrix
  
  ## rowbind mrna and lncrna data
  gene_count <- rbind(vsd_mrna, vsd_lncRNA)
  
  #change colnames : patient barocde
  colnames(vsd_mirna) <- substr(colnames(vsd_mirna),1,12)
  
  ##get the mRNA expression matrix for DEGs
  colnames(gene_count) <- substr(colnames(gene_count),1,12)
  
  #take intersect of colnames
  inter_col <- intersect(colnames(gene_count),colnames(vsd_mirna))
  #inter_col
  #take subset from both expression matrics
  vsd_mirna <- vsd_mirna[,inter_col]
  gene_count <- gene_count[,inter_col]
  
  #dim(vsd_mirna)
  #dim(gene_count)
  
  ##transpose it
  vsd_mirna <- t(vsd_mirna)
  gene_count <- t(gene_count)
  
  ## pearson correlation between miRNA and mRNA
  cor_data <- cor(vsd_mirna, gene_count, method = "pearson")
  
  #cor_data[1:5,1:5]
  
  cor_data[abs(cor_data) < 0.5] <- 0
  
  
  #create netwrok object
  net_mim <- graph_from_biadjacency_matrix(cor_data, weighted = "correlation")
  #net_mim
  # Assign types to nodes
  V(net_mim)$type <- ifelse(V(net_mim)$name %in% colnames(vsd_mirna), "miRNA",
                            ifelse(V(net_mim)$name %in% rownames(vsd_mrna), "mRNA", "lncRNA"))
  
  # Convert the types to a factor to avoid NAs or invalid values
  V(net_mim)$type <- as.factor(V(net_mim)$type)
  
  
  ## delete the nodes with zeor degree
  Isolated = which(degree(net_mim)==0)
  net_mim = delete_vertices(net_mim, Isolated) 
  
  save(net_mim, file = paste0(proj,"_network_gene_name.RData"))
  
  
  fc <- cluster_fast_greedy(net_mim)
  
  V(net_mim)$community <- fc$membership
  
  nodes <- data.frame(id = V(net_mim)$name, group = V(net_mim)$community)
  
  nodes <- nodes[order(nodes$id, decreasing = F),]
  
  nodes <- unique(nodes)
  edges <- as_data_frame(net_mim, what="edges")[1:3]
  
  edges$log2FC_miRNA <- res_mirna$log2FoldChange[match(edges$from, rownames(res_mirna))]
  
  lfc_gene <- res_proteincoding$log2FoldChange[match(edges$to, rownames(res_proteincoding))]
  lfc_gene <- ifelse(is.na(lfc_gene), res_lncRNA$log2FoldChange[match(edges$to, rownames(res_lncRNA))], lfc_gene) 
  
  edges$log2FC_Gene <- lfc_gene
  
  
  ## assotativity coeffieceint
  #assort_corr <- assortativity(net_mim, V(net_mim), directed=F)
  assortativity_coefficient <- assortativity_nominal(net_mim, V(net_mim)$type, directed = FALSE)
  degree_assortativity_coefficient <- assortativity_degree(net_mim, directed = FALSE)
  #assort_corr
  
  edges["assort_nominal",] <- assortativity_coefficient
  edges["assort_degree",] <- degree_assortativity_coefficient
  #edges
  write.csv(edges, file = paste0(proj,"_edges_gene_name.csv"))
  
}





#########################################################################################################
## to create network
#1. load miRNA deseq object for a cancertype

load(file = "DESeq_run.RData")
dds

#BLCA
## extract results from dds for cancer type tumor - normal comparison
res_mirna <- results(dds, contrast = c("condition", "TCGA-BLCA-TP", "TCGA-BLCA-NT"))
res_mirna

## take subset of miRNAs with cutoff mentioned before
res_mirna <- subset(res_mirna, abs(res_mirna$log2FoldChange) > 1 & res_mirna$padj < 0.05)
res_mirna

## take subset of count matrix of BLCA from the toal count matrix of 14 cancer
# load the normalized count matrix of 14 cancers

load("normalized_counts_14_Cancer.RData")

# 1. load the sample information object
load("sampleinfo.RData")
selcancer <- merged_df[merged_df$project == "TCGA-BLCA",]
dim(selcancer)

vsd
dim(vsd)
colnames(vsd)

# only BLCA samples
vsd_mirna <- assay(vsd)[,selcancer$sample_id]
dim(vsd_mirna)

## now subset the dysregulated miRNAS from the matrix
vsd_mirna <- vsd_mirna[rownames(res_mirna),]

########################################################################################
# load mRNA and lncRNA deseq run

# for mRNA
load("BLCA_dds_proteincoding.RData")
dds_proteincoding
table(is.na(rowData(dds_proteincoding)$gene_name))
rowData(dds_proteincoding)[1:5,1:5]
rownames(dds_proteincoding) <- rowData(dds_proteincoding)$gene_name



res_proteincoding <- results(dds_proteincoding, contrast = c("shortLetterCode", "TP", "NT"))
res_proteincoding <- subset(res_proteincoding, abs(res_proteincoding$log2FoldChange) > 1 & res_proteincoding$padj < 0.05)

## load normalized count matrix
load("BLCA_normalized_counts.RData")
vsd_mrna <- assay(vsd)
vsd_mrna <- vsd_mrna[rownames(res_proteincoding),] 
dim(vsd_mrna)

## for lncRNA
load("BLCA_dds_lncRNA.RData")
dds_lncRNA
table(is.na(rowData(dds_lncRNA)$gene_name))

rownames(dds_lncRNA) <- rowData(dds_lncRNA)$gene_name
rowData(dds_lncRNA)[1:5,1:5]


res_lncRNA <- results(dds_lncRNA, contrast = c("shortLetterCode", "TP", "NT"))
res_lncRNA <- subset(res_lncRNA, abs(res_lncRNA$log2FoldChange) > 1 & res_lncRNA$padj < 0.05)
res_lncRNA

# load lncRNA count matrix
load("BLCA_normalized_lncRNAcounts.RData")

vsd_lncRNA <- assay(vsd_lncRNA)
vsd_lncRNA <- vsd_lncRNA[rownames(res_lncRNA),]
dim(vsd_lncRNA)

#########################################################################################
## to create network only select common samples from miRNA and mRNA, lncRNA count matrix

## rowbind mrna and lncrna data
gene_count <- rbind(vsd_mrna, vsd_lncRNA)

#change colnames : patient barocde
colnames(vsd_mirna) <- substr(colnames(vsd_mirna),1,12)

##get the mRNA expression matrix for DEGs
colnames(gene_count) <- substr(colnames(gene_count),1,12)

#take intersect of colnames
inter_col <- intersect(colnames(gene_count),colnames(vsd_mirna))
inter_col
#take subset from both expression matrics
vsd_mirna <- vsd_mirna[,inter_col]
gene_count <- gene_count[,inter_col]

dim(vsd_mirna)
dim(gene_count)

##transpose it
vsd_mirna <- t(vsd_mirna)
gene_count <- t(gene_count)

## pearson correlation between miRNA and mRNA
cor_data <- cor(vsd_mirna, gene_count, method = "pearson")

cor_data[1:5,1:5]

cor_data[abs(cor_data) < 0.5] <- 0

library(igraph)
#create netwrok object
net_mim <- graph_from_biadjacency_matrix(cor_data, weighted = "correlation")
net_mim

## delete the nodes with zeor degree
Isolated = which(degree(net_mim)==0)
net_mim = delete_vertices(net_mim, Isolated) 

save(net_mim, file = "BLCA_network.RData")

load("BLCA_network.RData")

fc <- cluster_fast_greedy(net_mim)
fc
V(net_mim)$community <- fc$membership

nodes <- data.frame(id = V(net_mim)$name, group = V(net_mim)$community)
head(nodes)
## add shape of the nodes
#nodes <- cbind(shape = ifelse(grepl("ENSG", nodes$id, ignore.case = T), "circle", "square"), nodes)
#nodes$shape[nodes$biotype == "lncRNA"] <- "triangle"

nodes <- nodes[order(nodes$id, decreasing = F),]

nodes <- unique(nodes)
edges <- as_data_frame(net_mim, what="edges")[1:3]
edges

write.csv(edges, file = "BLCA_edges.csv")

V(net_mim)$type
net_mim

E(net_mim)$correlation

# Print the vertex IDs
print("Vertex IDs:")
print(V(net_mim))

# Print the vertex names
print("Vertex Names:")
print(V(net_mim)$name)

# Print all vertex attributes
print("Vertex Attributes:")
print(vertex_attr(net_mim))

## assotativity coeffieceint
match(V(net_mim), V(net_mim)$name)
assort_corr <- assortativity(net_mim, V(net_mim), directed=F)
assort_corr
assort_corr <- round(assort_corr, digits = 4)

assortativity_nominal(net_mim, V(net_mim)$community)
assortativity_nominal(net_mim, V(net_mim)$type)

assortativity_degree(net_mim, directed = F)

if (length(V(net_mim)$correlation) == vcount(net_mim)) {
  assortativity_correlation <- assortativity(net_mim, V(net_mim)$correlation, directed = FALSE)
  print(paste("Assortativity (correlation):", assortativity_correlation))
} else {
  print("Vertex correlation attribute length does not match the number of vertices.")
}

if (length(E(net_mim)$correlation) == ecount(net_mim)) {
  assortativity_correlation <- assortativity(net_mim, E(net_mim)$correlation, directed = FALSE)
  print(paste("Assortativity (correlation):", assortativity_correlation))
} else {
  print("Edge correlation attribute length does not match the number of edges.")
}


edges["assort",] <- assort_corr
tail(edges)

## find unique miRNAs per cancer
df <- readxl::read_excel("miRNAinteractions.xlsx", sheet = "miRNAs")
library(dplyr)

colnames(df)

miRNAs_subset_df <- select(df, contains("mi"))
dim(miRNAs_subset_df)
colnames(miRNAs_subset_df)

#unique_values <- setdiff(df$miRNAs_KIRP, unlist(select(df, -miRNAs_UCEC)))
#print(unique_values)


find_unique_values_df <- function(df) {
  unique_values_list <- list()
  
  for (col in colnames(df)) {
    other_cols <- select(df, -all_of(col))
    unique_values <- setdiff(df[[col]], unlist(other_cols))
    unique_values_list[[col]] <- unique_values
  }
  
  # Find the maximum length of the lists to pad shorter lists with NAs
  max_length <- max(sapply(unique_values_list, length))
  unique_values_list_padded <- lapply(unique_values_list, function(x) {
    c(x, rep(NA, max_length - length(x)))
  })
  
  # Convert the list to a dataframe
  unique_values_df <- as.data.frame(unique_values_list_padded)
  return(unique_values_df)
}

unique_values_in_each_column_df <- find_unique_values_df(df)
print(unique_values_in_each_column_df)

colnames(unique_values_in_each_column_df)

miRNAs_subset_df <- select(unique_values_in_each_column_df, contains("mi"))
print(miRNAs_subset_df)
dim(miRNAs_subset_df)
write.csv(miRNAs_subset_df, file = "unique_miRNA_per_cancer.csv")

