setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

library(DESeq2)

# Load your data
## load miRNAs counts from all samples
load("Normal_Tumor_panTCGA_isomir.RData")
dim(panisomir)
panisomir[1:5,1:5]

## load the sample info data
load("sampleinfo.RData")
head(merged_df)
# rename the rows of the sample info df
rownames(merged_df) <- 1:nrow(merged_df)

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
rownames(TCGA14cancers) <- 1:nrow(TCGA14cancers) 
dim(TCGA14cancers)
TCGA14cancers[1:5,]

panisomir[1:5,1:5]
cancer14 <- panisomir[, match(TCGA14cancers$sample_id, colnames(panisomir))]
dim(cancer14)
table(is.na(cancer14))
cancer14[is.na(cancer14)] <- 0

match(colnames(cancer14), TCGA14cancers$sample_id)

#paste0(TCGA14cancers$project,"-",TCGA14cancers$tissue_name)
TCGA14cancers$project_id <- paste0(TCGA14cancers$project,"-",TCGA14cancers$tissue_name)
# Split the data by cancer type
cancer_types <- split(TCGA14cancers, TCGA14cancers$project)
#cancer_types

# Create a list to store counts split by cancer type
counts_by_cancer <- lapply(cancer_types, function(info) {
  cancer14[, info$sample_id, drop=FALSE]
})

## calculate RPM
calculate_rpm <- function(counts) {
  total_reads <- colSums(counts)  # Total reads per sample
  rpm <- t(t(counts) / total_reads * 1e6)  # RPM calculation
  return(rpm)
}

counts_by_cancer$`TCGA-BLCA`[1:5,1:5]

# Check for NA values in counts_by_cancer
sapply(counts_by_cancer, function(x) sum(is.na(x)))

calculate_rpm <- function(counts) {
  total_reads <- colSums(counts)  # Total reads per sample
  rpm <- t(t(counts) / total_reads * 1e6)  # RPM calculation
  return(rpm)
}

# Apply RPM calculation to each cancer type
rpm_by_cancer <- lapply(counts_by_cancer, calculate_rpm)
rpm_by_cancer$`TCGA-BLCA`[1:5,1:5]

dim(rpm_by_cancer$`TCGA-BLCA`)


# Combine RPM-normalized counts into a single matrix
combined_rpm <- do.call(cbind, rpm_by_cancer)
dim(combined_rpm)

head(colnames(combined_rpm))
match(colnames(combined_rpm),TCGA14cancers$sample_id)

combined_rpm <- combined_rpm[,TCGA14cancers$sample_id]
combined_rpm[1:5,1:5]


# Log2 transformation of RPM data
log_transformed_counts <- log2(combined_rpm + 1)
log_transformed_counts[1:5,1:5]

library(preprocessCore)

# Apply quantile normalization
# Quantile normalization on log-transformed data
quantile_normalized_counts <- normalize.quantiles(as.matrix(log_transformed_counts))

quantile_normalized_counts[1:5,1:5]
rownames(quantile_normalized_counts) <- rownames(log_transformed_counts)
colnames(quantile_normalized_counts) <- colnames(log_transformed_counts)
quantile_normalized_counts[1:5,1:5]

# Perform PCA
pca <- prcomp(t(quantile_normalized_counts))
head(pca)

# Prepare PCA data for plotting
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Batch = TCGA14cancers$project)

library(ggplot2)
# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Batch)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Quantile-Normalized Data", x = "PC1", y = "PC2")


save(quantile_normalized_counts, file = "quantile_norm_counts14cancers.RData")

## remove batch effect from QN normalized data
load("quantile_norm_counts14cancers.RData")
dim(quantile_normalized_counts)

## batch effect removal
# Define batch: 1 for training data, 2 for external data
batch <- TCGA14cancers$project
batch

library(sva)

# Now, perform ComBat batch correction with the cancer type information as a covariate
batch_correct_QN_Allmi <- ComBat(dat = quantile_normalized_counts, batch = batch)
batch_correct_QN_Allmi[1:5,1:5]

quantile_normalized_counts[1:5,1:5]

#############################################################################################
mi597 <- read.csv("miRNAs_597_14cancer_mat.csv")
mi597[1:5,1:5]
dim(mi597)

colnames(mi597) <- gsub("\\.", "\\-", colnames(mi597))

quantille_norm597 <- batch_correct_QN_Allmi[colnames(mi597)[2:598],] 
dim(quantille_norm597)

quantille_norm597[1:5,1:5]
set.seed(123)
pca <- prcomp(t(quantille_norm597))
head(pca)


# Prepare PCA data for plotting
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Batch = TCGA14cancers$project)

library(ggplot2)
# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Batch)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Quantile-Normalized Data", x = "PC1", y = "PC2")

set.seed(123)
pca2 <- prcomp(mi597[,2:598])
pca_data2 <- data.frame(PC1 = pca2$x[,1], PC2 = pca2$x[,2], Batch = TCGA14cancers$project)

ggplot(pca_data2, aes(x = PC1, y = PC2, color = Batch)) + 
  geom_point() +
  theme_minimal() +
  labs(title = "PCA of Quantile-Normalized Data", x = "PC1", y = "PC2")

save(quantille_norm597, file = "quantile_norm_597.RData")

## remove the kich an dkirp samples form the quantile_norm 597 matrix
#kich
kidney2_sample_mat[1:5,1:5]

library(dplyr)
quantille_norm597[1:5,1:5]
match(mi597reduced$X, colnames(quantille_norm597))
quantille_norm597_nok <- quantille_norm597[, match(mi597reduced$X, colnames(quantille_norm597))] 
dim(quantille_norm597_nok)

quantille_norm597_nok <- t(quantille_norm597_nok)
quantille_norm597_nok <- as.data.frame(quantille_norm597_nok)

quantille_norm597_nok$response <- mi597reduced$response

dim(quantille_norm597_nok)
is.matrix(quantille_norm597_nok)
quantille_norm597_nok[1:5,1:5]

match(rownames(quantille_norm597_nok), mi597reduced$X)

write.csv(quantille_norm597_nok, file = "quant_norm_597_reduced_mat.csv")


## make csv file from quantile norm dataset,
load("quantile_norm_597.RData")
dim(quantille_norm597)

tail(match(colnames(quantille_norm597), TCGA14cancers$sample_id))

## add last column for the sample type info: response
response <- factor(paste0(TCGA14cancers$project, "-", TCGA14cancers$tissue_name))
response



quantille_norm597[1:5,1:5]
quantille_norm597 <- t(quantille_norm597)
quantille_norm597<- as.data.frame(quantille_norm597)
quantille_norm597$response <- response
dim(quantille_norm597)

write.csv(quantille_norm597, file = "QN_TCGA.csv")
