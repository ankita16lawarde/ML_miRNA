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

## select only UCEC TP samples

## subset from the sampleinfo table of all cancers
ucecTP <- TCGA14cancers[TCGA14cancers$project == "TCGA-UCEC" & TCGA14cancers$tissue_name == "TP",]
dim(ucecTP)

## load the RNA-Seq data and extact clinical info
load("RNASeq_Data/TCGA-UCECRNASeq.RData")
data

colData(data)[1:5,]

dim(colData(data)[colData(data)$shortLetterCode == "TP",])
uece_clinical <- colData(data)[colData(data)$shortLetterCode == "TP",]
uece_clinical$patient[1:5]

match(substr(ucecTP$sample_id,1,12), uece_clinical$patient)

ucecTP2 <- ucecTP$sample_id[!is.na(match(substr(ucecTP$sample_id,1,12), uece_clinical$patient))]
ucecTP2
unique(ucecTP2)

unique(uece_clinical$patient)
dim(ucecTP)

dim(uece_clinical)
# change rownames of clinical data
rownames(uece_clinical) <- substr(rownames(uece_clinical),1,12)
unique(rownames(uece_clinical))
rownames(uece_clinical)

uece_clinical_unique <- unique(uece_clinical)
dim(uece_clinical_unique)

uece_clinical_unique <- uece_clinical[!duplicated(uece_clinical$patient),]
rownames(uece_clinical_unique)

## get the count matrix of miRNAs from the all data
ucecTP_count <- assay(vsd)
ucecTP_count <- ucecTP_count[, ucecTP$sample_id]

dim(ucecTP_count)
# change the colnames
colnames(ucecTP_count) <- substr(colnames(ucecTP_count),1,12)

match(colnames(ucecTP_count), uece_clinical_unique$patient)
match(rownames(uece_clinical_unique), colnames(ucecTP_count))

ucecTP_count <- ucecTP_count[,!is.na(match(colnames(ucecTP_count), uece_clinical_unique$patient))]
uece_clinical_unique <- uece_clinical_unique[!is.na(match(colnames(ucecTP_count), uece_clinical_unique$patient)),]
dim(ucecTP_count)
dim(uece_clinical_unique)

## take common samples only 
common_samples <- intersect(colnames(ucecTP_count), rownames(uece_clinical_unique))
common_samples

ucecTP_count <- ucecTP_count[,common_samples]
uece_clinical_unique <- uece_clinical_unique[common_samples,]

# check the match
match(colnames(ucecTP_count), rownames(uece_clinical_unique))

uece_clinical_unique$vital_status
uece_clinical_unique$days_to_last_follow_up
uece_clinical_unique$paper_os_days
uece_clinical_unique$days_to_death

interactinglist <- readxl::read_excel(path = "miRNA_interactions2.xlsx", sheet = "miRNAs")
interactinglist
colnames(interactinglist)

# make a single list from interactinglis
#interactinglist[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)]
!is.na(interactinglist)

interactinglist <- interactinglist[!is.na(interactinglist)]
#interactinglist
interactinglist <- unique(interactinglist)
interactinglist

ucecTP_count <- ucecTP_count[interactinglist,]

###########################
##### loadd all required data
setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

library(DESeq2)
library(survival)
library(survminer)

load("normalized_counts_14_Cancer.RData")
#assay(vsd)[1:5,1:5]
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

interactinglist <- readxl::read_excel(path = "miRNA_interactions2.xlsx", sheet = "miRNAs")
interactinglist
colnames(interactinglist)

# make a single list from interactinglis
#interactinglist[c(1,3,5,7,9,11,13,15,17,19,21,23,25,27)]
!is.na(interactinglist)

interactinglist <- interactinglist[!is.na(interactinglist)]
#interactinglist
interactinglist <- unique(interactinglist)
interactinglist



#######################################
survival_analysis_all_miRNAs <- function(miRNA_subset, clinical_subset, threshold = 0.25) {
    
    results <- data.frame()
    
    # Function to check if a split satisfies the threshold
    split_threshold <- function(x, threshold) {
      t1 <- table(x)  # Calculate the distribution of hypo and hyper
      if (min(t1) / sum(t1) < threshold) {
        return("no")
      } else {
        return("yes")
      }
    }
    
    for (miRNA in rownames(miRNA_subset)) {
      miRNA_expression <- t(miRNA_subset[miRNA, ])
      clinical_subset$miRNA_exp <- miRNA_expression
      
      # Create factors based on different splitting methods
      mean_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= mean(as.numeric(clinical_subset$miRNA_exp)), "High", "Low"))
      
      # Check for valid factor levels
      if (length(unique(mean_split)) < 2) {
        cat("No valid split for", miRNA, "using mean. Skipping.\n")
        next
      }
      
      mean_split <- relevel(mean_split, ref = "Low")
      
      # Repeat for other splits with similar checks
      median_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= median(as.numeric(clinical_subset$miRNA_exp)), "High", "Low"))
      if (length(unique(median_split)) < 2) {
        cat("No valid split for", miRNA, "using median. Skipping.\n")
        next
      }
      median_split <- relevel(median_split, ref = "Low")
      
      q25_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= quantile(as.numeric(clinical_subset$miRNA_exp), 0.25), "High", "Low"))
      if (length(unique(q25_split)) < 2) {
        cat("No valid split for", miRNA, "using q25. Skipping.\n")
        next
      }
      q25_split <- relevel(q25_split, ref = "Low")
      
      q75_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= quantile(as.numeric(clinical_subset$miRNA_exp), 0.75), "High", "Low"))
      if (length(unique(q75_split)) < 2) {
        cat("No valid split for", miRNA, "using q75. Skipping.\n")
        next
      }
      q75_split <- relevel(q75_split, ref = "Low")
      
      # Store the splits in a list
      split_list <- list(median = median_split, mean = mean_split, q25 = q25_split, q75 = q75_split)
      
      # Check each split method against the threshold
      does_suit_threshold <- sapply(split_list, split_threshold, threshold = threshold)
      
      # Select the first valid split that satisfies the threshold
      split_list_yes <- split_list[which(does_suit_threshold == "yes")]
      
      if (length(split_list_yes) > 0) {
        chosen_split <- split_list_yes[[1]]  # Automatically selects the first valid split
        chosen_split_method <- names(split_list_yes)[1]  # Get the name of the chosen split method
        cat("Selected split method for", miRNA, ":", chosen_split_method, "\n")
      } else {
        cat("No valid split method for miRNA:", miRNA, ". Skipping.\n")
        next  # Skip this miRNA if no valid split method is found
      }
      
      # Check and create necessary columns for survival analysis
      clinical_subset$overall_survival <- ifelse(
        clinical_subset$vital_status == "Dead", 
        clinical_subset$days_to_death, 
        clinical_subset$days_to_last_follow_up
      )
      
      # Create binary 'deceased' column
      clinical_subset$deceased <- ifelse(clinical_subset$vital_status == "Dead", 1, 0)
      
      # Convert overall survival from days to months
      clinical_subset$overall_survival_months <- clinical_subset$overall_survival / 30
      
      # Perform Cox regression to get the Hazard Ratio (HR) and Wald p-value
      cox_fit <- coxph(Surv(overall_survival_months, deceased) ~ miRNA_exp, data = clinical_subset)
      cox_summary <- summary(cox_fit)
      HR <- cox_summary$coefficients[1, "exp(coef)"]
      lower_CI <- cox_summary$conf.int[1, 3]  # lower CI
      upper_CI <- cox_summary$conf.int[1, 4]  # upper CI
      Wald_pvalue <- cox_summary$coefficients[1, "Pr(>|z|)"]
      
      # Proportional hazards (PH) test p-value
      PH_test <- cox.zph(cox_fit)
      PH_test_pvalue <- PH_test$table[1, 3]
      
      # Log-rank test p-value
      logrank_test <- survdiff(Surv(overall_survival_months, deceased) ~ chosen_split, data = clinical_subset)
      LR_test_pvalue <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
      
      # Format the HR and CI as requested
      HR_CI <- paste0(HR, " (", round(lower_CI, 3), ";", round(upper_CI, 3), ")")
      
      # Store results
      results <- rbind(results, data.frame(
        miRNA = miRNA,
        HR_CI = HR_CI,
        Wald_Pvalue = round(Wald_pvalue, 6),
        PH_test_Pvalue = round(PH_test_pvalue, 6),
        LR_test_pvalue = round(LR_test_pvalue, 6),
        Split_Method = chosen_split_method  # Add the split method used
      ))
    }
    
    return(results)
  }
  
projects[1]

for (i in projects) {
  ## subset from the sampleinfo table of all cancers
  projectTP <- TCGA14cancers[TCGA14cancers$project == i & TCGA14cancers$tissue_name == "TP",]
  
  ## load the RNA-Seq data and extact clinical info
  #load("RNASeq_Data/TCGA-UCECRNASeq.RData")
  load(file = paste0("RNASeq_Data/", i, "RNASeq.RData"))
  
  clinical <- colData(data)[colData(data)$shortLetterCode == "TP",]
  
  match(substr(projectTP$sample_id,1,12), clinical$patient)
  
  # change rownames of clinical data
  rownames(clinical) <- substr(rownames(clinical),1,12)
  
  clinical_unique <- clinical[!duplicated(clinical$patient),]
  #rownames(clinical_unique)
  
  ## get the count matrix of miRNAs from the all data
  projectTP_count <- assay(vsd)
  projectTP_count <- projectTP_count[, projectTP$sample_id]
  
  # change the colnames
  colnames(projectTP_count) <- substr(colnames(projectTP_count),1,12)
  
  match(colnames(projectTP_count), clinical_unique$patient)
  match(rownames(clinical_unique), colnames(projectTP_count))
  
  #projectTP_count <- projectTP_count[,!is.na(match(colnames(projectTP_count), clinical_unique$patient))]
  
  #clinical_unique <- clinical_unique[!is.na(match(colnames(projectTP_count), clinical_unique$patient)),]
  ## take common samples only 
  common_samples <- intersect(colnames(projectTP_count), rownames(clinical_unique))
  
  projectTP_count <- projectTP_count[,common_samples]
  clinical_unique <- clinical_unique[common_samples,]
  
  dim(projectTP_count)
  dim(clinical_unique)
  
  projectTP_count <- projectTP_count[interactinglist,]
  
  result_table <- survival_analysis_all_miRNAs(miRNA_subset = projectTP_count, clinical_subset = clinical_unique, threshold = 0.25)
  write.csv(result_table, file = paste0(i,"surv_results.csv"))
}

################################################################

i = "TCGA-ESCA"

projectTP <- TCGA14cancers[TCGA14cancers$project == i & TCGA14cancers$tissue_name == "TP",]
dim(projectTP)

## load the RNA-Seq data and extact clinical info
#load("RNASeq_Data/TCGA-UCECRNASeq.RData")
load(file = paste0("RNASeq_Data/", i, "RNASeq.RData"))

clinical <- colData(data)[colData(data)$shortLetterCode == "TP",]

match(substr(projectTP$sample_id,1,12), clinical$patient)

# change rownames of clinical data
rownames(clinical) <- substr(rownames(clinical),1,12)

clinical_unique <- clinical[!duplicated(clinical$patient),]
#rownames(clinical_unique)

## get the count matrix of miRNAs from the all data
projectTP_count <- assay(vsd)
projectTP_count <- projectTP_count[, projectTP$sample_id]

# change the colnames
colnames(projectTP_count) <- substr(colnames(projectTP_count),1,12)

match(colnames(projectTP_count), clinical_unique$patient)
match(rownames(clinical_unique), colnames(projectTP_count))

#projectTP_count <- projectTP_count[,!is.na(match(colnames(projectTP_count), clinical_unique$patient))]

#clinical_unique <- clinical_unique[!is.na(match(colnames(projectTP_count), clinical_unique$patient)),]



## take common samples only 
common_samples <- intersect(colnames(projectTP_count), rownames(clinical_unique))
common_samples

projectTP_count <- projectTP_count[,common_samples]
clinical_unique <- clinical_unique[common_samples,]

dim(projectTP_count)

dim(clinical_unique)

#projectTP_count <- projectTP_count[interactinglist,]
# After calculating common_samples
common_samples <- intersect(colnames(projectTP_count), rownames(clinical_unique))

if (length(common_samples) == 0) {
  stop("No common samples found.")
}

projectTP_count <- projectTP_count[, common_samples]
clinical_unique <- clinical_unique[common_samples,]

if (any(interactinglist > nrow(projectTP_count))) {
  stop("interactinglist contains out-of-bounds indices.")
}

projectTP_count <- projectTP_count[interactinglist,]

result_table <- survival_analysis_all_miRNAs(miRNA_subset = projectTP_count, clinical_subset = clinical_unique, threshold = 0.25)
write.csv(result_table, file = paste0(i,"surv_results.csv"))




#########################################################
colData(data)

# Example usage:
result_table <- survival_analysis_all_miRNAs(miRNA_subset = ucecTP_count, clinical_subset = uece_clinical_unique, threshold = 0.25)

head(result_table)
order(result_table$LR_test_pvalue, decreasing = FALSE)
result_table[order(result_table$LR_test_pvalue, decreasing = FALSE),]

###### KM plot ##############################

library(survminer)

KM_all_miRNAs <- function(miRNA_subset, clinical_subset, miRNA, threshold = 0.25) {
  
  # Function to check if a split satisfies the threshold
  split_threshold <- function(x, threshold) {
    t1 <- table(x)  # Calculate the distribution of hypo and hyper
    if (min(t1) / sum(t1) < threshold) {
      return("no")
    } else {
      return("yes")
    }
  }
  
  # Extract expression data for the selected miRNA
  miRNA_expression <- t(miRNA_subset[miRNA, ])
  clinical_subset$miRNA_exp <- miRNA_expression
  
  # Create factors based on different splitting methods
  mean_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= mean(as.numeric(clinical_subset$miRNA_exp)), "Low", "High"))
  mean_split <- relevel(mean_split, ref = "Low")
  
  median_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= median(as.numeric(clinical_subset$miRNA_exp)), "Low", "High"))
  median_split <- relevel(median_split, ref = "Low")
  
  q25_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= quantile(as.numeric(clinical_subset$miRNA_exp), 0.25), "Low", "High"))
  q25_split <- relevel(q25_split, ref = "Low")
  
  q75_split <- factor(ifelse(as.numeric(clinical_subset$miRNA_exp) <= quantile(as.numeric(clinical_subset$miRNA_exp), 0.75), "Low", "High"))
  q75_split <- relevel(q75_split, ref = "Low")
  
  # Store the splits in a list
  split_list <- list(median = median_split, mean = mean_split, q25 = q25_split, q75 = q75_split)
  
  # Check each split method against the threshold
  does_suit_threshold <- sapply(split_list, split_threshold, threshold = threshold)
  
  # Select the first valid split that satisfies the threshold
  split_list_yes <- split_list[which(does_suit_threshold == "yes")]
  
  if (length(split_list_yes) > 0) {
    chosen_split <- split_list_yes[[1]]  # Automatically selects the first valid split
    chosen_split_method <- names(split_list_yes)[1]  # Get the name of the chosen split method
    cat("Selected split method for", miRNA, ":", chosen_split_method, "\n")
  } else {
    cat("No valid split method for miRNA:", miRNA, ". Skipping.\n")
    return(NULL)  # Return NULL if no valid split method is found
  }
  
  # Check and create necessary columns for survival analysis
  clinical_subset$overall_survival <- ifelse(
    clinical_subset$vital_status == "Dead", 
    clinical_subset$days_to_death, 
    clinical_subset$days_to_last_follow_up  # Ensure this matches the actual column name in your data
  )
  
  # Create binary 'deceased' column
  clinical_subset$deceased <- ifelse(clinical_subset$vital_status == "Dead", 1, 0)
  # Convert overall survival from days to months
  clinical_subset$overall_survival_months <- clinical_subset$overall_survival / 30
  
  # Use the chosen split (not raw miRNA expression) for survival analysis
  clinical_subset$ID_regulation <- chosen_split
  
  # Perform the Kaplan-Meier fit based on miRNA expression groups
  #fit <- survfit(Surv(overall_survival, deceased) ~ ID_regulation, data = clinical_subset)
  # Perform the Kaplan-Meier fit based on miRNA expression groups
  fit <- survfit(Surv(overall_survival_months, deceased) ~ ID_regulation, data = clinical_subset)
  
  # Fit a Cox proportional hazards model to get hazard ratios
  cox_fit <- coxph(Surv(overall_survival_months, deceased) ~ ID_regulation, data = clinical_subset)
  hr <- exp(coef(cox_fit))  # Calculate hazard ratios
  hr_labels <- paste("HR:", round(hr, 2))  # Create labels for HR
  
  
  # Plot the Kaplan-Meier curve using ggsurvplot
  k <- ggsurvplot(
    fit, 
    fun = "pct",               # Plot survival as a percentage
    data = clinical_subset,    # Clinical data with miRNA groups
    pval = TRUE,               # Display p-value from log-rank test
    surv.median.line = "hv",   # Show horizontal and vertical lines for median survival
    title = paste("Kaplan-Meier Plot for", miRNA),  # Title of the plot
    xlab = "Time (Months)",    # Label for the x-axis
    ylab = "Survival Probability (%)",  # Label for the y-axis
    legend.title = "miRNA Expression",  # Title for the legend
    legend.labs = levels(clinical_subset$ID_regulation)  # Labels for the legend
  )
  
  # Add HR to the plot using annotation
  k$plot <- k$plot + 
    annotate("text", x = 20, y = 10, label = hr_labels, size = 4, vjust = 1)
  
  
  # Return the KM plot object
  return(k)
}



KM_all_miRNAs(ucecTP_count, uece_clinical_unique, miRNA = "hsa-miR-449-3p")

#i <- "TCGA-BLCA"
for(i in projects){
  projectTP <- TCGA14cancers[TCGA14cancers$project == i & TCGA14cancers$tissue_name == "TP",]
  
  ## load the RNA-Seq data and extact clinical info
  #load("RNASeq_Data/TCGA-UCECRNASeq.RData")
  load(file = paste0("RNASeq_Data/", i, "RNASeq.RData"))
  
  clinical <- colData(data)[colData(data)$shortLetterCode == "TP",]
  
  match(substr(projectTP$sample_id,1,12), clinical$patient)
  
  # change rownames of clinical data
  rownames(clinical) <- substr(rownames(clinical),1,12)
  
  clinical_unique <- clinical[!duplicated(clinical$patient),]
  #rownames(clinical_unique)
  
  ## get the count matrix of miRNAs from the all data
  projectTP_count <- assay(vsd)
  projectTP_count <- projectTP_count[, projectTP$sample_id]
  
  # change the colnames
  colnames(projectTP_count) <- substr(colnames(projectTP_count),1,12)
  
  match(colnames(projectTP_count), clinical_unique$patient)
  match(rownames(clinical_unique), colnames(projectTP_count))
  
  #projectTP_count <- projectTP_count[,!is.na(match(colnames(projectTP_count), clinical_unique$patient))]
  
  #clinical_unique <- clinical_unique[!is.na(match(colnames(projectTP_count), clinical_unique$patient)),]
  ## take common samples only 
  common_samples <- intersect(colnames(projectTP_count), rownames(clinical_unique))
  
  projectTP_count <- projectTP_count[,common_samples]
  clinical_unique <- clinical_unique[common_samples,]
  
  dim(projectTP_count)
  dim(clinical_unique)
  
  projectTP_count <- projectTP_count[interactinglist,]
  
  #result_table <- survival_analysis_all_miRNAs(miRNA_subset = projectTP_count, clinical_subset = clinical_unique, threshold = 0.25)
  #write.csv(result_table, file = paste0(i,"surv_results.csv"))
  km_plot <- KM_all_miRNAs(projectTP_count, clinical_unique, miRNA = "hsa-miR-204-5p")
  #km_plot$plot
  ggsave(filename = paste0(i,"KM_plot.jpg"), plot = km_plot$plot, width = 8, height = 8, dpi = 300)

}


#########################################

i = "TCGA-LIHC"

mirna_list <- c("hsa-miR-149-5p",
                "hsa-miR-139-5p",
                "hsa-miR-7-5p",
                "hsa-miR-139-3p",
                "hsa-miR-93-3p")

for(mi in mirna_list){
  projectTP <- TCGA14cancers[TCGA14cancers$project == i & TCGA14cancers$tissue_name == "TP",]
  
  ## load the RNA-Seq data and extact clinical info
  #load("RNASeq_Data/TCGA-UCECRNASeq.RData")
  load(file = paste0("RNASeq_Data/", i, "RNASeq.RData"))
  
  clinical <- colData(data)[colData(data)$shortLetterCode == "TP",]
  
  match(substr(projectTP$sample_id,1,12), clinical$patient)
  
  # change rownames of clinical data
  rownames(clinical) <- substr(rownames(clinical),1,12)
  
  clinical_unique <- clinical[!duplicated(clinical$patient),]
  #rownames(clinical_unique)
  
  ## get the count matrix of miRNAs from the all data
  projectTP_count <- assay(vsd)
  projectTP_count <- projectTP_count[, projectTP$sample_id]
  
  # change the colnames
  colnames(projectTP_count) <- substr(colnames(projectTP_count),1,12)
  
  match(colnames(projectTP_count), clinical_unique$patient)
  match(rownames(clinical_unique), colnames(projectTP_count))
  
  #projectTP_count <- projectTP_count[,!is.na(match(colnames(projectTP_count), clinical_unique$patient))]
  
  #clinical_unique <- clinical_unique[!is.na(match(colnames(projectTP_count), clinical_unique$patient)),]
  ## take common samples only 
  common_samples <- intersect(colnames(projectTP_count), rownames(clinical_unique))
  
  projectTP_count <- projectTP_count[,common_samples]
  clinical_unique <- clinical_unique[common_samples,]
  
  dim(projectTP_count)
  dim(clinical_unique)
  
  projectTP_count <- projectTP_count[interactinglist,]
  
  #result_table <- survival_analysis_all_miRNAs(miRNA_subset = projectTP_count, clinical_subset = clinical_unique, threshold = 0.25)
  #write.csv(result_table, file = paste0(i,"surv_results.csv"))
  km_plot <- KM_all_miRNAs(projectTP_count, clinical_unique, miRNA = mi)
  #km_plot$plot
  ggsave(filename = paste0(i, mi, "_HR_KM_plot.jpg"), plot = km_plot$plot, width = 8, height = 8, dpi = 300)
  
}

projectTP <- TCGA14cancers[TCGA14cancers$project == i & TCGA14cancers$tissue_name == "TP",]

## load the RNA-Seq data and extact clinical info
#load("RNASeq_Data/TCGA-UCECRNASeq.RData")
load(file = paste0("RNASeq_Data/", i, "RNASeq.RData"))

clinical <- colData(data)[colData(data)$shortLetterCode == "TP",]

match(substr(projectTP$sample_id,1,12), clinical$patient)

# change rownames of clinical data
rownames(clinical) <- substr(rownames(clinical),1,12)

clinical_unique <- clinical[!duplicated(clinical$patient),]
#rownames(clinical_unique)

## get the count matrix of miRNAs from the all data
projectTP_count <- assay(vsd)
projectTP_count <- projectTP_count[, projectTP$sample_id]

# change the colnames
colnames(projectTP_count) <- substr(colnames(projectTP_count),1,12)

match(colnames(projectTP_count), clinical_unique$patient)
match(rownames(clinical_unique), colnames(projectTP_count))

#projectTP_count <- projectTP_count[,!is.na(match(colnames(projectTP_count), clinical_unique$patient))]

#clinical_unique <- clinical_unique[!is.na(match(colnames(projectTP_count), clinical_unique$patient)),]
## take common samples only 
common_samples <- intersect(colnames(projectTP_count), rownames(clinical_unique))

projectTP_count <- projectTP_count[,common_samples]
clinical_unique <- clinical_unique[common_samples,]

dim(projectTP_count)
dim(clinical_unique)

projectTP_count <- projectTP_count[interactinglist,]

#result_table <- survival_analysis_all_miRNAs(miRNA_subset = projectTP_count, clinical_subset = clinical_unique, threshold = 0.25)
#write.csv(result_table, file = paste0(i,"surv_results.csv"))
km_plot <- KM_all_miRNAs(projectTP_count, clinical_unique, miRNA = "hsa-miR-30a-5p")
#km_plot$plot
ggsave(filename = paste0(i, "hsa-miR-30a-5p_HR_KM_plot.jpg"), plot = km_plot$plot, width = 8, height = 8, dpi = 300)
