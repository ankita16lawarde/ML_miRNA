setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")

library(dplyr)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "BLCA_edges_gene_name")
inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

# View the result
print(miRNA_groups)
dim(miRNA_groups)

write.csv(miRNA_groups, file = "BLCA_miRNA_gene_list.csv")

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "BRCA_edges_gene_name")
inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_BRCA <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)


inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "ESCA_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_ESCA <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "HNSC_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_HNSC <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)


inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "KIRP_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_KIRP <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "KIRC_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_KIRC <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)


inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "LIHC_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_LIHC <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "LUSC_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_LUSC <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "LUAD_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_LUAD <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)


inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "PRAD_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_PRAD <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)


inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "STAD_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_STAD <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "THCA_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_THCA <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)

inter_net <- readxl::read_excel(path = "interacting_miRNA_table.xlsx", sheet = "UCEC_edges_gene_name")
#inter_net

# Create the comma-separated list of the 'to' column for each miRNA
miRNA_groups_UCEC <- inter_net %>%
  group_by(from) %>%
  summarise(target_genes = paste(to, collapse = ", ")) %>%
  rename(miRNA = from)


# Combine all miRNA groups into a list
miRNA_groups_list <- list(miRNA_groups, miRNA_groups_BRCA, miRNA_groups_ESCA, miRNA_groups_HNSC,
                          miRNA_groups_KICH, 
                          miRNA_groups_KIRC,
                          miRNA_groups_KIRP,
                          miRNA_groups_LIHC,
                          miRNA_groups_LUAD,
                          miRNA_groups_LUSC,
                          miRNA_groups_PRAD,
                          miRNA_groups_STAD,
                          miRNA_groups_THCA,
                          miRNA_groups_UCEC)

# Use Reduce to perform full join across all data frames in the list
combined_miRNA_groups <- Reduce(function(x, y) full_join(x, y, by = "miRNA"), miRNA_groups_list)

print(combined_miRNA_groups)

colnames(combined_miRNA_groups)
colnames(combined_miRNA_groups) <- c("miRNA","BLCA", "BRCA", "ESCA", "HNSC",
                                     "KICH", "KIRC", "KIRP", "LIHC", "LUAD",
                                     "LUSC", "PRAD", "STAD", "THCA", "UCEC")

write.csv(combined_miRNA_groups, file = "combined_miRNA_network_list.csv")
