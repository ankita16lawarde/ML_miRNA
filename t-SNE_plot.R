setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")
library(umap)
library(dplyr)
library(Rtsne)
library(ggplot2)
library(DESeq2)
library(ggrepel)

load("normalized_counts_14_Cancer.RData")
assay(vsd)[1:5,1:5]

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

## remove ESCA-NT from the TCGA14cancers
TCGA14cancers$samplestypes <- paste0(TCGA14cancers$project, "-", TCGA14cancers$tissue_name)
TCGA14cancers_noesca <- TCGA14cancers[TCGA14cancers$samplestypes != "TCGA-ESCA-NT",]
dim(TCGA14cancers_noesca)

rownames(TCGA14cancers_noesca) <- 1:nrow(TCGA14cancers_noesca)
#####################################################################################
##################################################################

## now load the interacting miRNAs from each cancertypes into a single object
interactinglist <- readxl::read_excel(path = "miRNA_interactions2.xlsx", sheet = "miRNAs")
interactinglist
colnames(interactinglist)

# make a single list from interactinglis
#interactinglist[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)]
!is.na(interactinglist)

interactinglist <- interactinglist[!is.na(interactinglist)]
#interactinglist
interactinglist <- unique(interactinglist)

## take subset from total count matrix
count_597mi <- assay(vsd)[interactinglist,]
dim(count_597mi)

count_597mi <- count_597mi[,TCGA14cancers_noesca$sample_id]

response <- factor(paste0(TCGA14cancers_noesca$project, "-", TCGA14cancers_noesca$tissue_name))
response

count_597mi[1:5,1:5]

count_597mi <- as.data.frame(t(count_597mi))
## add last column for the sample type info: response

count_597mi$response <- response
dim(count_597mi)
count_597mi[1:5,1:5]

# Perform t-SNE dimensionality reduction
set.seed(123)
#count_597mi

interacting_mi <- count_597mi[,1:597]
colnames(interacting_mi)
dim(interacting_mi)

## now make T-SNE for RFE features
## for boruta
#LDA
rfe <- read.csv("google_colab_output/RFE_features_150new.csv")
head(rfe)
ref_mi <- interacting_mi[,rfe$X0]
dim(ref_mi)
lda <- read.csv("google_colab_output/lda_selected_features.csv")
lda_mi <- interacting_mi[,lda$X0]

rf <- read.csv("google_colab_output/rf_selected_features.csv")
rf_mi <- interacting_mi[,rf$X0]

boruta <- read.csv("google_colab_output/boruta_selected_features.csv")
boruta_mi <- interacting_mi[,boruta$X0]

## read validation dataset
#val_exp <- read.csv("validation_batch_correct.csv")
#val_response <- read.csv("responseAll.csv")
#val_exp[1:5,1:5]

#rownames(val_exp) <- val_exp$X
#val_exp <- val_exp[,-1]
set.seed(123)
tsne_model <- Rtsne(as.matrix(interacting_mi), dims = 2, perplexity = 30, verbose = TRUE)
tsne_model

set.seed(123)
tsne_model <- Rtsne(as.matrix(ref_mi), dims = 2, perplexity = 30, verbose = TRUE)
tsne_model

#tsne_df[1:5,]
tsne_df <- tsne_model$Y
colnames(tsne_df)
colnames(tsne_df) <- c("tSNE1", "tSNE2")

tsne_df <- as.data.frame(tsne_df)
tsne_df$response <- count_597mi$response
tsne_df

tsne_df$response_short <- gsub("TCGA-(.*)", "\\1", tsne_df$response)

# Calculate the centroids for each response (cluster)
centroids <- tsne_df %>%
  group_by(response) %>%
  summarize(tSNE1 = mean(tSNE1), tSNE2 = mean(tSNE2))

# Create a new variable to hold modified cancer names (removing 'TCGA')
centroids$response_short <- gsub("TCGA-(.*)", "\\1", centroids$response)  # Keep -TP or -NT

# Define unique labels for each cancer type to avoid duplication
unique_labels <- centroids %>%
  distinct(response_short, .keep_all = TRUE) %>%
  mutate(label_y = tSNE2 + 3)  # Position labels slightly above the centroids

# Define your existing color palettes
dark_colors <- c("#00008B", "green", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928", "#1E90FF", 
                 "#fb9a99", "#fdbf6f", "mediumpurple2", "#B8860B", "darkgreen", "#FA8072")

neon_colors <- c("#00FFFF", "firebrick", "#FF1493", "#FF4500", "#DA70D6", "#FFD700", "#87CEFA",
                 "slateblue", "#FF69B4", "#FF8C00", "#DDA0DD", "#808000", "deepskyblue", "magenta3")

all_colors <- c(dark_colors, neon_colors)
response_levels <- levels(tsne_df$response)
color_palette <- setNames(all_colors, response_levels)

# Create the t-SNE plot with labeled clusters
plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = response)) +
  geom_point(size = 0.8) +
  scale_color_manual(values = color_palette) +
  labs(title = "RFE features", color = "Tissue") +
  theme_minimal() +
  geom_label_repel(data = unique_labels, aes(label = response_short, color = response),  
                   size = 2, fontface = "bold", show.legend = FALSE, 
                   box.padding = 0.3, point.padding = 0.3) +  # Avoid overlap with labels
  theme(plot.title = element_text(size = 14, face = "bold"), # Title size and style
        legend.position = "none") # Legend position

plot

# Save the plot
ggsave(filename = "t-SNE_RFE_with_colored_labels_no_legend.jpg", plot = plot, width = 12, height = 8, dpi = 300)

dev.off()

# Plot t-SNE with custom colors
plot <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = response)) +
  geom_point(size = 0.6) +
  scale_color_manual(values = color_palette) +
  labs(title = "RFE features", color = "Tissue") +
  theme_minimal()

ggsave(filename = "t-SNE_plot_RFE_features.jpg", plot = plot, width = 8, height = 6, dpi = 300)

dev.off()






## tsne plot for each cancer using the impt miRNAs
# Create a directory to save the plots if it doesn't exist
if (!dir.exists("tSNE_plots_val")) {
  dir.create("tSNE_plots_val")
}

# Loop through each unique cancer type and create a separate t-SNE plot
unique_cancers <- unique(tsne_df$Cancer)
unique_cancers
for (cancer in unique_cancers) {
  plot_data <- tsne_df %>% filter(Cancer == cancer)
  
  plot <- ggplot(plot_data, aes(x = tSNE1, y = tSNE2, color = response)) +
    geom_point(size = 2) + # Modify point size here
    scale_color_manual(values = color_palette) +
    labs(title = paste("t-SNE Plot for", cancer)) +
    theme_minimal()
  
  # Save the plot as a PNG file
  ggsave(filename = paste0("tSNE_plots_val/tSNE_plot_", cancer, ".jpg"), plot = plot, width = 8, height = 6)
}


################## 
## UMAP on the same data
# Perform UMAP on the data
umap_model <- umap(as.matrix(rfe), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_model
# Create a data frame for UMAP results
umap_df <- as.data.frame(umap_model$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$response <- count_597mi$response  # Assuming 'response' is the variable you want to color by

# Plot UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = response)) +
  geom_point() +
  scale_color_discrete(name = "Samples") +
  labs(title = "UMAP Plot of Important Features") +
  theme_minimal()

umap_df$response
substr(umap_df$response, 1,9)
umap_df$Cancer <- substr(umap_df$response, 1,9)
umap_df$Sample <- substr(umap_df$response, 11, 12)

sample(umap_df, 5)

head(umap_df, 5)

# Plot t-SNE
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = response)) +
  geom_point(size = 2) +
  scale_color_discrete(name = "Samples") +
  labs(title = "t-SNE Plot of Important Features")

c("#800000","#FF0000","#FA8072","#FF4500","#B8860B","#808000",
  "#9ACD32","#7FFF00","#006400","#90EE90","#20B2AA","#00FFFF",
  "#7FFFD4","#6495ED","#00BFFF","#1E90FF","#00008B","#4B0082",
  "#7B68EE","#9932CC","#FF00FF","#C71585","#FF69B4","#D2691E")

# Define 14 distinct dark colors
# Define 14 distinct dark colors
dark_colors <- c("#00008B", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928", "#1E90FF", 
                 "#b2df8a", "#fb9a99", "#fdbf6f", "#cab2d6", "#B8860B", "#b3de69", "#FA8072")

# Define neon shades for each dark color
neon_colors <- c("#00FFFF", "#00FF00", "#FF1493", "#FF4500", "#DA70D6", "#FFD700", "#87CEFA",
                 "#7FFF00", "#FF69B4", "#FF8C00", "#DDA0DD", "#808000", "#ADFF2F", "#FFC0CB")

# Combine dark and faint colors
# Combine dark and neon colors
all_colors <- c(dark_colors, neon_colors)

all_colors

# Define the response levels
levels(umap_df$response)
response_levels <- levels(umap_df$response)

# Create a named vector for mapping response to colors
color_palette <- setNames(all_colors, response_levels)

# Plot t-SNE with custom colors
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = response)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = color_palette) +
  labs(title = "Important features") +
  theme_minimal()

# Calculate centroids for each cluster (response)
centroids <- umap_df %>%
  group_by(response) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

library(dplyr)
library(ggrepel)

umap_df[1:5,]
# Calculate centroids for each Cancer type (or any preferred method of labeling)
centroids <- umap_df %>%
  group_by(Cancer) %>%
  summarize(UMAP1 = mean(UMAP1), UMAP2 = mean(UMAP2))

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = response)) +
  geom_point(size = 0.5) +
  geom_text_repel(data = centroids, aes(label = Cancer), 
                  color = "black", 
                  size = 2, 
                  fontface = "bold", 
                  box.padding = 0.5,  # Adjust padding around labels
                  point.padding = 0.5, # Adjust padding from points
                  max.overlaps = Inf) +  # Allow more overlaps to be resolved
  scale_color_manual(values = color_palette, guide = "legend") +  # Use legend for `response`
  labs(title = "Important Features") +
  theme_minimal()
