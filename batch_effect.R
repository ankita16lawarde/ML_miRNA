## combine the validations datasets

setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")


## load the CPTAC dataset 
load("CPTAC/CPTAC_combined.RData")
load("CPTAC/responseCPTAC.RData")
head(responseCPTAC)
CPTAC[1:5,1:5]
dim(responseCPTAC)
dim(CPTAC)

table(is.na(match(responseCPTAC$sample_id, rownames(CPTAC))))
table(is.na(match(rownames(CPTAC), responseCPTAC$sample_id)))

match(rownames(CPTAC), responseCPTAC$sample_id)

CPTAC[151,1:5]

responseCPTAC$sample_id <- rownames(CPTAC)

## load miTED datasets
load("miTED_data_portal/miTED_response.RData")
load("miTED_data_portal/miTED_counts.RData")
miTED_counts[1:5,1:4]

match(rownames(miTED_counts), miTED_response$sample_id)


## load the 10 sample dataset 
load("miTED_data_portal/sample10QN_GSE138734.RData")
load("miTED_data_portal/response_GSE138734.RData")
sample_10_QN[1:5,1:5]
head(response)

## esophagus cancer dataset
load("GSE61047_esca_mat.RData")
microarray_normalized[1:5,1:5]
rownames(microarray_normalized)

response_esca <- cbind("response" = c(rep("TCGA-ESCA-TP",4), rep("TCGA-ESCA-NT", 4)), 
                       "sample_id" = rownames(microarray_normalized))
response_esca <- as.data.frame(response_esca)
response_esca  
  
###########
# load th ekiCH and kipr datasets from TCGA
kirp <- load("kirp_validation_mat.RData")
kirp #"kidney3_sample_mat"
kich <- load("kich_validation_mat.RData")
kich #"kidney2_sample_mat"


load("quantile_norm_counts14cancers.RData")
quantile_normalized_counts[1:5,1:5]
kirp <- quantile_normalized_counts[,kidney3_sample_mat$X]
kirp[1:5,1:5]
colnames(kirp)

kich <- quantile_normalized_counts[,kidney2_sample_mat$X]
colnames(kich)

dim(kich)
kidney <- cbind(kirp, kich)
colnames(kidney)

responsekidney <- cbind("response" = c(rep("TCGA-KIRP-TP",5), rep("TCGA-KIRP-NT", 5), 
                                       rep("TCGA-KICH-TP", 5), rep("TCGA-KICH-NT", 5)), 
                        "sample_id" = colnames(kidney))

## take subset of 597 miRNAs
colnames(CPTAC)
kidney <- kidney[colnames(CPTAC),]
dim(kidney)
kidney <- t(kidney)


# Transpose the datasets
#CPTAC_t <- t(CPTAC)
#miTED_t <- t(miTED_counts)
#sample10_t <- t(sample_10_QN)
#ESCA_t <- t(microarray_normalized)


combined_data <- rbind(CPTAC, miTED_counts, sample_10_QN, microarray_normalized)
dim(combined_data)
combined_data <- rbind(combined_data, kidney)
combined_data[1:5,1:5]
dim(combined_data)
## 


# Create batch variable
batch <- c(rep(1, nrow(CPTAC)), rep(2, nrow(miTED_counts)), 
           rep(3, nrow(sample_10_QN)), rep(4, nrow(microarray_normalized)),
           rep(5, nrow(kidney)))
batch

# Transpose the matrix to get miRNAs in rows and samples in columns
combined_data_norm_t <- t(combined_data)
combined_data_norm_t[1:5,1:5]

combined_data_norm_t <- as.data.frame(combined_data_norm_t)

# Load libraries
library(sva)

# Check if all values are numeric
is_numeric <- sapply(combined_data_norm_t, is.numeric)
is_numeric

if (!all(is_numeric)) {
  stop("Non-numeric values detected in the matrix.")
}

# Check the dimensions of the data matrix and batch vector
dim(combined_data_norm_t)
length(batch)

# Check if all data is numeric
if (!is.numeric(as.matrix(combined_data_norm_t))) {
  combined_data_norm_t <- apply(combined_data_norm_t, 2, as.numeric)
}
combined_data_norm_t[1:5,1:5]

# Check for missing values
any_na <- sum(is.na(combined_data_norm_t))
any_na

combined_data[1:5,1:5]
rownames(combined_data_norm_t) <- colnames(combined_data)

responseCPTAC
response
response_esca
responsekidney
head(miTED_response)
head(miTED_response[,1])

mod_val <- c(responseCPTAC[,1], miTED_response[,1], 
             response[,1], response_esca[,1], responsekidney[,1])
mod_val

mod_val <- model.matrix(~ mod_val)


# Re-run ComBat
combined_data_batch_corrected <- ComBat(dat = as.matrix(combined_data_norm_t), 
                                        batch = batch, par.prior = TRUE, mod = NULL)
#combined_data_batch_corrected
combined_data_batch_corrected[1:5,1:5]

colnames(combined_data_batch_corrected) <- colnames(combined_data_norm_t)
combined_data[1:5,1:5]
rownames(combined_data_batch_corrected) <- colnames(combined_data)
dim(combined_data_batch_corrected)

combined_data_batch_corrected <- t(combined_data_batch_corrected)
combined_data_batch_corrected[1:5,1:5]
save(combined_data_batch_corrected, file = "combined_data_batch_corrected.RData")


load("combined_data_batch_corrected.RData")
dim(combined_data_batch_corrected)

combined_data[1:5,1:5]
combined_data_norm_t[1:5,1:5]

rownames(combined_data_norm_t) <- colnames(combined_data)

# Perform PCA before batch correction
pca_before <- prcomp(t(combined_data_norm_t), scale = TRUE)
pca_after <- prcomp(t(combined_data_batch_corrected), scale = TRUE)
head(pca_after$x)


library(ggplot2)
factor(batch)
# Plot PCA before batch correction
ggplot(data = as.data.frame(pca_before$x), aes(PC1, PC2, color = factor(batch))) +
  geom_point(size = 2) +
  ggtitle("PCA Before Batch Correction")

pca_data <- as.data.frame(pca_after$x)  # Assuming 'pca_after' contains PCA results
pca_data$batch <- factor(batch)         # Add batch info to the data

# Plot PCA after batch correction
ggplot(data = as.data.frame(pca_after$x), aes(PC1, PC2, color = factor(batch))) +
  geom_point(size = 2) +
  ggtitle("PCA After Batch Correction")

ggplot(data = pca_data, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 2) +
  geom_text(aes(label = rownames(pca_data)), vjust = -0.5, size = 2.5) +
  ggtitle("PCA After Batch Correction with Sample Labels") +
  theme_minimal()


# Calculate centroid
centroid <- colMeans(pca_data[, c("PC1", "PC2")])

# Calculate Euclidean distance from each point to the centroid
pca_data$distance_from_centroid <- sqrt((pca_data$PC1 - centroid[1])^2 + (pca_data$PC2 - centroid[2])^2)

# Set a threshold for outliers (e.g., 1.5 times the interquartile range)
threshold <- quantile(pca_data$distance_from_centroid, 0.75) + 1.5 * IQR(pca_data$distance_from_centroid)

# Identify samples that are outliers
outliers <- pca_data[pca_data$distance_from_centroid > threshold, ]

outliers[1:5,1:5]

outlier_samples <- rownames(outliers)
print(outlier_samples)  # This will print the samples not clustered together

write.csv(outlier_samples, "outlier_samples_valdata.csv", row.names = FALSE)

#####################
write.csv(combined_data_batch_corrected, file = "combined_data_batch_corrected.csv")

########################################################################################
responseAll <- rbind(responseCPTAC, miTED_response, response, response_esca, responsekidney)
dim(responseAll)

write.csv(responseAll, file = "responseAll.csv")
### batch correct the validation dataset with training data

## read the training dataset csv file

training_data <- read.csv("quant_norm_597_reduced_mat.csv")

training_data[1:5,1:5]
dim(training_data)
colnames(training_data) <- gsub("\\.","\\-", colnames(training_data))
rownames(training_data) <- training_data$X
training_data <- training_data[,-c(1,599)]

dim(combined_data_norm_t)
combined_data_norm_t[1:5,1:5]
dim(training_data)
colnames(combined_data_norm_t)

# Combine the training and external datasets into one matrix
training_ex_data <- cbind(t(training_data), combined_data_norm_t)
dim(training_ex_data)
training_ex_data[1:5,1:5]

c(rep(1, nrow(CPTAC)), rep(2, nrow(miTED_counts)), 
  rep(3, nrow(sample_10_QN)), rep(4, nrow(microarray_normalized)),
  rep(5, nrow(kidney)))

nrow(CPTAC)
nrow(quantille_norm597)

table(is.na(match(rownames(combined_data_batch_corrected), rownames(quantille_norm597))))
rownames(quantille_norm597)[!is.na(match(rownames(quantille_norm597),rownames(combined_data_batch_corrected)))]


quantille_norm597_reduced <- quantille_norm597[is.na(match(rownames(quantille_norm597),rownames(combined_data_batch_corrected))),]
dim(quantille_norm597_reduced)

nrow(CPTAC)
nrow(miTED_counts)
nrow(sample_10_QN)
nrow(microarray_normalized)

# Define batch: 1 for training data, 2 for external data
batch2 <- c(rep(1, nrow(quantille_norm597_reduced)), rep(2, nrow(CPTAC)), rep(3, nrow(miTED_counts)), 
            rep(4, nrow(sample_10_QN)), rep(5, nrow(microarray_normalized)),
            rep(6, 20))


dim(quantille_norm597_reduced)
dim(combined_data_batch_corrected)

traing_val <- rbind(quantille_norm597_reduced[,1:597],combined_data_batch_corrected)
factor(batch2)

traing_val <- t(traing_val)
combat_corrected_data <- ComBat(dat = traing_val, batch = batch2, ref.batch = 1)


# Perform PCA on your data
pca_result <- prcomp(t(combat_corrected_data), scale. = TRUE)

# Extract PCA scores
pca_scores <- pca_result$x

# Check the proportion of variance explained
summary(pca_result)

# Convert PCA scores to a data frame for ggplot2
pca_df <- as.data.frame(pca_scores)
pca_df$Batch <- factor(batch2)  # Add batch information to the data frame

# Load ggplot2 for plotting
library(ggplot2)

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  labs(title = "PCA of Batch-Corrected Data", x = "PC1", y = "PC2") +
  theme_minimal()

dim(quantille_norm597_reduced)
dim(combat_corrected_data)
external_val_data <- combat_corrected_data[,7106:ncol(combat_corrected_data)]
dim(external_val_data)
external_val_data[1:5,1:5]
save(external_val_data, file = "external_val_data.RData")

cancertype <- read.csv("quant_norm_597_reduced_mat.csv")
cancertype <- cancertype[,c(1,599)]
head(cancertype)
cancertype <- cancertype[,c(2,1)]
colnames(cancertype) <- c("response", "sample_id")

responseAll
head(responseAll)

cancer_types_combined <- rbind(cancertype, responseAll)
head(cancer_types_combined)
dim(cancer_types_combined)
cancer_types_combined <- cancer_types_combined$response
cancer_response_combined <- as.factor(cancer_types_combined)
head(cancer_types_combined)

# Create the model matrix, accounting for cancer type
mod <- model.matrix(~ cancer_types_combined)
head(mod)
dim(mod)

dim(t(training_ex_data))
length(batch2)
nrow(mod)
batch2

training_ex_data <- t(training_ex_data)
dim(training_ex_data)
training_ex_data[1:5,1:5]
batch2
sum(is.na(mod))  # Should return 0
str(mod)
batch2 <- as.factor(batch2)
length(batch2) == nrow(training_ex_data)  # Should return TRUE

levels(batch2)
sum(is.na(training_ex_data))  # Should return 0
sum(is.na(batch2))            # Should return 0
sum(is.na(mod))               # Should return 0
training_ex_data <- as.matrix(training_ex_data)
training_ex_data[1:5,1:5]
str(training_ex_data)  # Ensure this is a numeric matrix

# Example of how to reorder batch2 if needed:
batch2 <- batch2[match(rownames(training_ex_data), names(batch2))]
batch2
head(names(batch2))

dim(training_data)
# Get the sample names from each dataset
sample_names <- c(rownames(training_data), rownames(CPTAC), rownames(miTED_counts),
                  rownames(sample_10_QN), rownames(microarray_normalized), rownames(kidney))

# Assign sample names to batch2
names(batch2) <- sample_names

head(rownames(training_ex_data))
head(names(batch2))

head(batch2)
dim(training_ex_data)
# Ensure row names of data and batch names are identical
all(rownames(training_ex_data) == names(batch2))  # Should return TRUE
batch2 <- batch2[match(rownames(training_ex_data), names(batch2))]
table(batch2)
# Ensure it's a numeric matrix
training_ex_data <- as.matrix(training_ex_data)
str(training_ex_data)  # Should show numeric matrix

# Ensure no missing values
sum(is.na(training_ex_data))  # Should return 0

# Remove miRNAs (columns) that have zero variance across all samples
training_ex_data0var <- training_ex_data[, apply(training_ex_data, 2, var) != 0]

levels(batch2)
table(batch2)
nrow(CPTAC)
nrow(miTED_counts)
nrow(sample_10_QN)
nrow(microarray_normalized)
nrow(kidney)
dim(combined_data_norm_t)
combined_data_norm_t[]

# Now, perform ComBat batch correction with the cancer type information as a covariate
combat_corrected_data <- ComBat(dat = training_ex_data0var, batch = batch2, ref.batch = 1)
traceback()

# # Subset of 100 samples and 100 miRNAs
# subset_data <- training_ex_data[1:100,]
# subset_data[1:4,1:5]
# subset_batch <- batch2[1:100]
# subset_batch
# # Example of checking dimensions
# dim(subset_data)      # Should return (number of genes, number of samples)
# length(subset_batch)  # Should be equal to number of samples
# 
# # Test ComBat on a subset
# combat_corrected_data_subset <- ComBat(dat = subset_data, batch = subset_batch)


training_ex_data[1:5,1:5]
training_ex_data <- t(training_ex_data)
combat_corrected_data <- ComBat(dat = training_ex_data, 
                                batch = batch2,
                                mod = mod, 
                                ref.batch = 1 )

# Perform PCA on your data
pca_result <- prcomp(t(combat_corrected_data), scale. = TRUE)

# Extract PCA scores
pca_scores <- pca_result$x

# Check the proportion of variance explained
summary(pca_result)

# Convert PCA scores to a data frame for ggplot2
pca_df <- as.data.frame(pca_scores)
pca_df$Batch <- factor(batch2)  # Add batch information to the data frame

# Load ggplot2 for plotting
library(ggplot2)

# Plot PCA
ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point(size = 3) +
  labs(title = "PCA of Batch-Corrected Data", x = "PC1", y = "PC2") +
  theme_minimal()



#now take subset of validation set from the combat_corrected_data and use it for validation
dim(combat_corrected_data)
table(batch2)
validation_batch_correct <- combat_corrected_data[, c(7106:ncol(combat_corrected_data))]
dim(validation_batch_correct)

validation_batch_correct <- t(validation_batch_correct)
validation_batch_correct[1:5,1:5]

match(responseAll$sample_id, rownames(validation_batch_correct))
write.csv(validation_batch_correct, file = "validation_batch_correct.csv")


# Assuming 'outlier_samples' is a vector of sample IDs that are considered outliers
# and 'validation_batch_correct' is your dataset

validation_batch_correct <- read.csv("validation_batch_correct.csv")
dim(validation_batch_correct)
validation_batch_correct[1:5,1:5]
rownames(validation_batch_correct) <- validation_batch_correct$X
colnames(validation_batch_correct) <- gsub("\\.","\\-", colnames(validation_batch_correct))

validation_batch_correct <- validation_batch_correct[,-1]

outlier_samples <- read.csv("outlier_samples_valdata.csv")
head(outlier_samples)
outlier_samples <- outlier_samples$x

# Remove outlier samples from the dataset
validation_batch_correct_noout <- validation_batch_correct[!(rownames(validation_batch_correct) %in% outlier_samples), ]
dim(validation_batch_correct_noout)
dim(validation_batch_correct)
dim(responseAll)

head(responseAll)
# Assuming 'outlier_samples' is a vector of sample IDs that are considered outliers
# and 'responseAll' is your dataset

# Remove outlier samples from the responseAll dataset based on the 'sample_id' column
responseAll_noout <- responseAll[!(responseAll$sample_id %in% outlier_samples), ]
dim(responseAll_noout)
head(responseAll_noout)
rownames(responseAll_noout) <- 1:nrow(responseAll_noout)

match(responseAll_noout$sample_id, rownames(validation_batch_correct_noout))
table(responseAll_noout$response)
table(responseAll$response)

write.csv(validation_batch_correct_noout, file = "validation_batch_correct_noout.csv")
write.csv(responseAll_noout, file = "responseAll_noout.csv")


res <- read.csv("responseAll.csv")
res
table(res$response)
