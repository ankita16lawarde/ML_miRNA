## venn diagram with five lists

setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

library(VennDiagram)

# Generate 3 sets of 200 words
set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

# Chart
venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Set 1" , "Set 2 " , "Set 3"),
  filename = '#14_venn_diagramm.png',
  output=TRUE
)


mi <- readxl::read_excel("miRNA_interactions2.xlsx", sheet = "duplicate_removed")
mi

interactingmi597 <- mi$`All 597 miRNAs`
interactingmi597 <- interactingmi597[!is.na(interactingmi597)]
RF299 <- mi$`RF features`
RF299 <- RF299[!is.na(RF299)]
RFE150 <- mi$`RFE features`
RFE150 <- RFE150[!is.na(RFE150)]
Boruta531 <- mi$`Boruta features`
Boruta531 <- Boruta531[!is.na(Boruta531)]
length(Boruta531)
LDA338 <- mi$`LDA features`
LDA338 <- LDA338[!is.na(LDA338)]


# Find common miRNAs across all five lists
common_miRNAs <- Reduce(intersect, list(interactingmi597, RF299, RFE150, Boruta531, LDA338))

common_miRNAs
write.csv(common_miRNAs, file = "common_miRNAs_from_feature_sets.csv")

common_RF_RFE_LDA <- Reduce(intersect, list(RF299, RFE150, LDA338))

common_RF_RFE_LDA 
match(common_miRNAs, common_RF_RFE_LDA)

venn.diagram(
  x = list(interactingmi597, RF299, RFE150, Boruta531, LDA338),
  category.names = c("Interacting_597", "RF_298", "RFE_150", "Boruta_530", "LDA_352"),
  filename = 'mifeatures_venn_diagramm.png',
  output=TRUE,
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  imagetype="png" ,
  height = 1000 , 
  width = 1000 , 
  resolution = 300
)

venn.diagram(
  x = list(interactingmi597, RF299, RFE150, Boruta531, LDA338),
  category.names = c("interacting_597", "RF_298", "RFE_150", "Boruta_530", "LDA_352"),
  filename = 'mifeatures_venn_diagramm.png',
  output = TRUE,
  
  # Numbers
  cex = .3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names and colors
  cat.cex = 0.5,
  cat.col = c("red", "blue", "green", "purple", "orange"),  # Set colors for each category
  imagetype = "png",
  height = 1000, 
  width = 1000, 
  resolution = 300,
  
  # Position the category names outside the circles
  cat.pos = c(-30, 30, 180, 90, 0),  # Adjust positions as necessary
  cat.dist = c(0.05, 0.05, 0.05, 0.05, 0.05)  # Adjust distance from circles
)

# Open a JPEG device
#jpeg('mifeatures_venn_diagramm.jpeg', width = 1000, height = 1000, quality = 100)

venn.diagram(
  x = list(interactingmi597, RF299, RFE150, Boruta531, LDA338),
  category.names = c("597_miRNAs", "RF_298", "RFE_150", "Boruta_530", "LDA_352"),
  
  # Set names and softer colors
  fill = c("#FFB3BA", "#FFDFBA", "#BAFFC9", "#BAE1FF", "#FFC3A0"),  # Softer shades for each circle
  cex = 1.2,  # Increase text size for numbers
  fontface = "bold",
  fontfamily = "sans",
  
  # Positioning
  cat.cex = 1.2,  # Increase size of category names
  cat.col = "black",  # Softer shades for category names
  cat.pos = c(-30, 30, 180, 90, 0),  # Adjust positions
  cat.dist = c(0.07, 0.07, 0.07, 0.07, 0.07),  # Adjust distance from circles
  filename = 'mifeatures_venn_diagramm.jpeg',  # Specify the output file name
  resolution = 300,# Set resolution for JPEG
  show.n = FALSE  # Hide the zero numbers
)

# Define the sets
miRNA_597 <- interactingmi597
RF_299 <- RF299
RFE_150 <- RFE150
Boruta_531 <- Boruta531
LDA_338 <- LDA338


compendium <- read.csv("matched_miRNA_compendium_network.csv")
compendium

cmc_mi <- read.csv("matched_cmc.csv")
cmc_mi <- cmc_mi$x


mi_597 <- compendium$All_597_miRs
isoform <- compendium$matched_isoform
EV <- compendium$matched_EV
#exo <- compendium$matched_exosom
#clinical <- compendium$matched_clinical
#nonco <- compendium$matched_noncoRNA_DB

common_comp_cmc <- Reduce(intersect, list(mi_597, isoform,EV, cmc_mi))
common_comp_cmc


venn.diagram(
  x = list(
    "597_miRNAs" = mi_597,
    "EV_miRNAs" = EV,
    "CMC_miRNAs" = cmc_mi,
    "Compendium_miRNAs" = isoform
  ),
  category.names = c("597_miRNAs", "EV_miRNAs", "CMC_miRNAs", "Compendium_miRNAs"),
  
  # Set names and colors
  fill = c("#FFB3BA", "#FFDFBA", "#BAFFC9", "#BAE1FF"),  # Softer shades for each circle
  cex = 1.2,  # Text size for numbers
  fontface = "bold",
  fontfamily = "sans",
  
  # Positioning
  # Positioning
  cat.cex = 1.2,  # Size of category names
  cat.col = "black",  # Category names in black
  #cat.pos = c(0, 160, 270, 45),  # Adjusted positions
  #cat.dist = c(0.08, 0.08, 0.08, 0.08),  # Increased distance from circles
  show.n = TRUE,  # Show numbers
  filename = "compendium__cmc_venn.jpeg",  # Set to NULL to plot directly
  resolution = 300  # Set resolution for JPEG
)

# Find common miRNAs among all lists
common_miRNAs <- Reduce(intersect, list(mi_597, EV, cmc_mi,isoform))

common_miRNAs



venn.diagram(
  x = list(
    "597_miRNAs" = mi_597,
    "noncoRNA-DB" = nonco
  ),
  category.names = c("597_miRNAs", "noncoRNA-DB"),
  
  # Set names and colors
  fill = c("#FFB3BA", "#BAE1FF"),  # Softer shades for each circle
  cex = 1.2,  # Text size for numbers
  fontface = "bold",
  fontfamily = "sans",
  
  # Positioning
  # Positioning
  cat.cex = 1.2,  # Size of category names
  cat.col = "black",  # Category names in black
  #cat.pos = c(0, 160),  # Adjusted positions
  #cat.dist = c(0.08, 0.08),  # Increased distance from circles
  #show.n = TRUE,  # Show numbers
  filename = "nonco_venn.jpeg",  # Set to NULL to plot directly
  resolution = 300  # Set resolution for JPEG
)

## read the CMC miRNA list file
cmc_mi <- read.csv("matched_cmc.csv")
cmc_mi <- cmc_mi$x

interactingmi597

venn.diagram(
  x = list(
    "597_miRNAs" = interactingmi597,
    "CMC_miRNAs" = cmc_mi
  ),
  category.names = c("597_miRNAs", "CMC_miRNAs"),
  
  # Set names and colors
  fill = c("#FFB3BA", "#BAE1FF"),  # Softer shades for each circle
  cex = 1.2,  # Text size for numbers
  fontface = "bold",
  fontfamily = "sans",
  
  # Positioning
  # Positioning
  cat.cex = 1.2,  # Size of category names
  cat.col = "black",  # Category names in black
  #cat.pos = c(0, 160),  # Adjusted positions
  #cat.dist = c(0.08, 0.08),  # Increased distance from circles
  #show.n = TRUE,  # Show numbers
  filename = "cmc_venn.jpeg",  # Set to NULL to plot directly
  resolution = 300  # Set resolution for JPEG
)



# bar plot 
noncoRNA <- readxl::read_excel(path = "noncoRNA_DB.xlsx", sheet = "barplot")
noncoRNA

# Create the data frame
data <- data.frame(
  Count = c(256, 18),  # counts for 'resistant' and 'sensitive'
  Type = c("Resistant", "Sensitive"),
  Evidence_Code = c("Predicted", "Predicted")  # Assuming same evidence code for simplicity
)

# Add the second set of counts
data2 <- data.frame(
  Count = c(186, 157),  # counts for 'resistant' and 'sensitive' in 'validated'
  Type = c("Resistant", "Sensitive"),
  Evidence_Code = c("Validated", "Validated")  # Again assuming same evidence code for simplicity
)

# Combine both data frames
data_combined <- rbind(data, data2)
data_combined


# Create the bar plot
plot <- ggplot(data_combined, aes(x = Type, y = Count, fill = Evidence_Code)) +
  geom_bar(stat = "identity", position = "dodge") +  # Use dodge to separate bars
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5) +  # Add counts on top of the bars
  labs(title = "Count of miRNAs by Sensitivity and Evidence Code",
       x = "Sensitivity Type",
       y = "Count of miRNAs") +
  scale_fill_manual(values = c("Predicted" = "darkorange", "Validated" = "darkblue")) +  # Custom colors
  theme_minimal() #+
  #theme(axis.text.x = element_text(face = "bold"))  # Make x-axis tick labels bold
ggsave("noncoRNA_barplot.jpeg", plot = plot, width = 10, height = 8, dpi = 300, device = "jpeg")


## pie chart

# Create the data_combined data frame
data_combined <- data.frame(
  Count = c(256, 18, 186, 157),
  Type = c("Resistant", "Sensitive", "Resistant", "Sensitive"),
  Evidence_Code = c("Predicted", "Predicted", "Validated", "Validated")
)

# Aggregate the counts by Type and Evidence_Code
library(dplyr)
agg_data <- data_combined %>%
  group_by(Type, Evidence_Code) %>%
  summarise(Total_Count = sum(Count))

# Create a pie chart for the aggregated data
library(ggplot2)

data_combined
agg_data

# Modify the labels in the 'fill' column to use underscores
#agg_data$fill_label <- gsub(".", "_", interaction(agg_data$Type, agg_data$Evidence_Code))
#agg_data$fill_label


# Create the pie chart with counts added and custom legend title
pie_chart <- ggplot(agg_data, aes(x = "", y = Total_Count, fill = interaction(Type, Evidence_Code, sep = "_"))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(title = "", 
       fill = "Drug-target association") +  # Rename legend title
  theme_void() +
  scale_fill_brewer(palette = "Set3") +  # Choose a color palette
  geom_text(aes(label = Total_Count), position = position_stack(vjust = 0.5), color = "black")

pie_chart


#Save the pie chart as a JPEG file with 300 dpi
ggsave("pie_chart.jpeg", plot = pie_chart, dpi = 300, width = 6, height = 6, units = "in")
