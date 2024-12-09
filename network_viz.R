setwd("C:/Users/ankita.lawarde/OneDrive - Tartu Ülikool/backup/work/R_work")


## network visualizations
library(igraph)
library(visNetwork)
library(DESeq2)


top_features <- readxl::read_excel(path = "network_interaction_tables_top_all.xlsx", sheet = "feature_imp")
top_features
top_features <- top_features$miRNA
top_features

load("UCEC_network_gene_name.RData")
net_mim


project <- c("BLCA" ,"BRCA", "ESCA", "HNSC", "KICH", "KIRC", "KIRP", "LIHC", "LUAD", "LUSC", "PRAD", "STAD", "THCA", "UCEC")

#1. load miRNA deseq object for a cancertype
load(file = "DESeq_run.RData")
#dds

# 1. load the sample information object
load("sampleinfo.RData")
rownames(merged_df) <- 1:nrow(merged_df)
merged_df[1:5,1:5]

#proj = "UCEC"

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
  
  #save(net_mim, file = paste0(proj,"_network_gene_name.RData"))
  
  
  fc <- cluster_fast_greedy(net_mim)
  
  V(net_mim)$community <- fc$membership
  
  
  # Step 1: Filter nodes of type "miRNA" that are in top_features
  miRNA_nodes <- V(net_mim)[V(net_mim)$type == "miRNA" & V(net_mim)$name %in% top_features]
  
  # Check the matched miRNA nodes
  print(miRNA_nodes)
  
  # Step 2: Create a subgraph with the filtered miRNA nodes
  subnetwork <- induced_subgraph(net_mim, miRNA_nodes)
  
  # Check the number of vertices and edges in the subnetwork
  cat("Number of vertices in the miRNA subnetwork:", vcount(subnetwork), "\n")
  cat("Number of edges in the miRNA subnetwork:", ecount(subnetwork), "\n")
  
  # Check the edges for the selected miRNA nodes in the original network
  edges_check <- E(net_mim)[.inc(miRNA_nodes)]
  cat("Number of edges connected to the selected miRNA nodes in the original network:", length(edges_check), "\n")
  
  # If you want to see which edges are connected to these nodes:
  if (length(edges_check) > 0) {
    print(E(net_mim)[.inc(miRNA_nodes)])
  }
  
  
  # Find all neighbors of the selected miRNA nodes
  neighboring_nodes <- unique(unlist(neighborhood(net_mim, order = 1, nodes = miRNA_nodes)))
  neighboring_nodes
  # Convert neighboring_nodes to a vertex sequence if it's not
  expanded_nodes <- unique(c(as.vector(miRNA_nodes), as.vector(neighboring_nodes)))
  
  # Create a new subgraph that includes both selected nodes and their neighbors
  expanded_subnetwork <- induced_subgraph(net_mim, expanded_nodes)
  
  save(expanded_subnetwork, file = paste0(proj,"expanded_subnet.RData"))
  
  # # Check the number of vertices and edges in the expanded subnetwork
  # cat("Number of vertices in the expanded miRNA subnetwork:", vcount(expanded_subnetwork), "\n")
  # cat("Number of edges in the expanded miRNA subnetwork:", ecount(expanded_subnetwork), "\n")
  # 
  # # Get the number of vertices
  # cat("Number of vertices:", vcount(expanded_subnetwork), "\n")
  # 
  # # Get the vertex attributes
  # vertex_attributes <- vertex.attributes(expanded_subnetwork)
  # print(vertex_attributes)
  # 
  # if (vcount(expanded_subnetwork) > 0) {
  #   nodes_df <- data.frame(id = V(expanded_subnetwork)$name,
  #                          label = V(expanded_subnetwork)$name,
  #                          type = V(expanded_subnetwork)$type)  # Add type column
  #   
  #   # Prepare edges data frame
  #   edges_df <- as.data.frame(as_edgelist(expanded_subnetwork))
  #   colnames(edges_df) <- c("from", "to")
  #   
  #   # Assuming edges data frame already has log2FC information
  #   edges_df$log2FC_miRNA <- res_mirna$log2FoldChange[match(edges_df$from, rownames(res_mirna))]
  #   lfc_gene <- res_proteincoding$log2FoldChange[match(edges_df$to, rownames(res_proteincoding))]
  #   lfc_gene <- ifelse(is.na(lfc_gene), res_lncRNA$log2FoldChange[match(edges_df$to, rownames(res_lncRNA))], lfc_gene) 
  #   edges_df$log2FC_Gene <- lfc_gene
  #   
  #   # Create a new column in nodes_df for log2FC
  #   nodes_df$log2FC <- ifelse(nodes_df$type == "miRNA", 
  #                             edges_df$log2FC_miRNA[match(nodes_df$id, edges_df$from)], 
  #                             edges_df$log2FC_Gene[match(nodes_df$id, edges_df$to)])
  #   
  #   # Handle NA values for log2FC
  #   nodes_df$log2FC[is.na(nodes_df$log2FC)] <- 0  # Replace NA with 0 or any other value you deem appropriate
  #   
  #   # Define border colors based on log2FC values
  #   nodes_df$color.border <- ifelse(nodes_df$log2FC > 0, "green", 
  #                                   ifelse(nodes_df$log2FC < 0, "red", "gray"))
  #   
  #   # Define shapes and colors based on types
  #   shape_mapping <- c("miRNA" = "dot", "mRNA" = "box", "lncRNA" = "ellipse")  # Example shapes
  #   color_mapping <- c("miRNA" = "lightblue", "gene" = "lightyellow", "other" = "lightcoral")  # Example colors
  #   
  #   # Add new columns for shapes and colors
  #   nodes_df$shape <- shape_mapping[nodes_df$type]
  #   nodes_df$color.background <- color_mapping[nodes_df$type]
  #   
  #   # Handle cases where type may not be in shape_mapping or color_mapping
  #   nodes_df$shape[is.na(nodes_df$shape)] <- "dot"  # Default shape for undefined types
  #   nodes_df$color.background[is.na(nodes_df$color.background)] <- "gray"  # Default color for undefined types
  #   
  #   #save(nodes_df, file = paste0(proj, "nodes_df.RData"))
  #   #save(edges_df, file = paste0(proj, "edges_df.RData"))
  #   
  #   # Create the interactive undirected network
  #   sub_net <- visNetwork(nodes_df, edges_df) %>%
  #     visOptions(highlightNearest = TRUE, 
  #                nodesIdSelection = TRUE) %>%
  #     visEdges(arrows = 'none', color = "gray") %>%  # Set arrows to 'none'
  #     visNodes(shape = nodes_df$shape, 
  #              color = list(background = nodes_df$color.background, border = nodes_df$color.border),
  #              borderWidth = 3) %>%  # Set shapes and colors
  #     visLayout(randomSeed = 123)
  # } else {
  #   stop("The expanded subnetwork has no vertices.")
  # }
  
  #save(sub_net, file = paste0(proj, "39_top_sub_network.RData"))
  }

sub_net
ggsave("network_plot.jpeg", plot = sub_net, width = 10, height = 10, dpi = 300)


plot(expanded_subnetwork)

# Define a shape mapping for supported vertex shapes
shape_mapping <- c("miRNA" = "circle", "mRNA" = "square")  # Use supported shapes

# Update your nodes_df accordingly
nodes_df$shape <- shape_mapping[nodes_df$type]
nodes_df

# Replace any unsupported shapes with a default shape
nodes_df$shape[is.na(nodes_df$shape)] <- "circle"  # Default shape for undefined types

# Set vertex attributes for coloring and shapes
V(expanded_subnetwork)$log2FC <- nodes_df$log2FC
#V(expanded_subnetwork)$color.border <- nodes_df$color.border
V(expanded_subnetwork)$vertex.shape <- nodes_df$shape
V(expanded_subnetwork)$color.background <- nodes_df$color.border

V(expanded_subnetwork)$color.background <- ifelse(V(expanded_subnetwork)$color.background == "red", "lightgreen" ,"lightcoral")

# Define border colors based on log2FC values
#border_colors <- ifelse(V(expanded_subnetwork)$log2FC > 0, "green", 
 #                       ifelse(V(expanded_subnetwork)$log2FC < 0, "red", "gray"))  # Green for positive, red for negative, gray for zero

# Define the size of vertex labels and shape border
vertex_size <- 10 # Set a larger vertex size for visibility
label_size <- 2  # Increase label size (default is typically 1)

# Set edge width based on a condition (you can customize this as needed)
E(expanded_subnetwork)$width <- 5  # Set a default edge width (adjust as necessary)

# Plotting the network
jpeg("network_plot.jpeg", width = 3000, height = 3000, quality = 100)  # High-resolution output



# Plotting the network
# Create a thicker border by plotting a second set of vertices with larger size
plot(expanded_subnetwork, 
     vertex.size = 5,                       # Set a constant vertex size for the main vertices
     vertex.label = V(expanded_subnetwork)$name,               # Show vertex labels
     vertex.label.cex = 2,                  # Increase label size (adjust the value)
     vertex.color = V(expanded_subnetwork)$color.background,    # Set vertex background color
     #vertex.color = "transparent",
     #vertex.frame.color = border_colors,       # Set border colors based on log2FC
     vertex.label.color = "black",            # Set label color
     vertex.label.family = "sans",            # Font family (optional)
     #edge.arrow.size = 0.5,                   # Set arrow size
     edge.width = E(expanded_subnetwork)$width,                  # Use edge width defined earlier
     main = "",  # Title
     vertex.shape = V(expanded_subnetwork)$vertex.shape,# Set vertex shapes
)


dev.off() 

library(ggplot2)
library(ggnetwork)

# Assuming expanded_subnetwork is your igraph object
gg_network <- ggnetwork(expanded_subnetwork, layout = "fruchterman")
# View the structure of the ggnetwork object
head(gg_network)

# Convert type to a factor if it isn't already
gg_network$type <- as.factor(gg_network$type)


# Create the ggplot network visualization with node labels
gg_plot <- ggplot(data = gg_network, aes(x = x, y = y)) +
  geom_edges(aes(xend = xend, yend = yend, color = "gray"), 
             arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +  # Add edges
  geom_nodes(aes(fill = ifelse(log2FC > 0, "green", 
                               ifelse(log2FC < 0, "red", "gray")),
                 shape = type),  # Use 'type' directly as shape
             size = 5, stroke = 1.5) +  # Set node aesthetics
  geom_text(aes(label = name),  # Add labels for nodes
            vjust = -1,  # Adjust vertical position (above the node)
            hjust = 0.5,  # Center the label horizontally
            size = 3,  # Size of the text
            color = "black") +  # Color of the text
  scale_shape_manual(values = c("miRNA" = 21, "mRNA" = 22, "lncRNA" = 23)) +  # Define shapes
  theme_blank() +  # Use a blank theme
  labs(title = "Network Visualization") +
  theme(legend.position = "none")

# Print the plot
print(gg_plot)


################
#read the nodes and edge table files
load("LUADnodes_df.RData")
load("LUADedges_df.RData")



nodes_df

head(nodes_df)

nodes_df[nodes_df$id == "hsa-miR-30a-3p"|"miR-30c-2-3p",]

# Subset the edges_df for rows where 'hsa-miR-449b-5p' is either in 'from' or 'to' column
subset_edges <- edges_df[edges_df$from == 'hsa-miR-30a-3p' | edges_df$to == 'hsa-miR-30a-3p', ]

# Display the resulting subset
print(subset_edges)

# Step 2: Get the unique node IDs from both 'from' and 'to' columns of the subset
matching_node_ids <- unique(c(subset_edges$from, subset_edges$to))

# Step 3: Subset the nodes_df based on the matching node IDs
subset_nodes <- nodes_df[nodes_df$id %in% matching_node_ids, ]

# Display the resulting subset of nodes_df
print(subset_nodes)


# Create graph
#g <- graph_from_data_frame(d = edges_df, vertices = nodes_df, directed = FALSE)

g <- graph_from_data_frame(d = subset_edges, vertices = subset_nodes, directed = FALSE)

# Assign node border color based on log2FC value (red for negative, green for positive)
V(g)$color.border <- ifelse(V(g)$log2FC < 0, 'green', 'red')
V(g)$color.border

# Assign node colors based on type (miRNA, mRNA, lncRNA)
V(g)$color <- ifelse(V(g)$type == 'miRNA', 'lightyellow', 
                     ifelse(V(g)$type == 'mRNA', 'lightcoral', 'lightblue'))

# Set node size, shape, and label size
V(g)$size <- 10  # Smaller node size
# Assign shape for miRNA, mRNA, and lncRNA
V(g)$vertex.shape <- ifelse(V(g)$type == 'miRNA', 'circle', 
                     ifelse(V(g)$type == 'mRNA', 'square', 'csquare'))  # Triangle for lncRNA

V(g)$label.cex <- 1  # Smaller label size

# Increase edge length and width
E(g)$width <- 2  # Increase the width of the edges

# Increase number of iterations for better spacing of nodes
layout_fr <- layout_with_fr(g, niter = 200)
#E(g)$color

# Save the figure as JPEG with 300 DPI
jpeg("network_plot_LUAD2.jpeg", width = 12, height = 10, units = "in", res = 300)
plot(g, 
     vertex.size = V(g)$size, 
     vertex.label.cex = V(g)$label.cex, 
     vertex.label.font=2,
     vertex.color = V(g)$color, 
     vertex.shape = V(g)$vertex.shape, 
     vertex.label.color = 'black', 
     vertex.frame.color = V(g)$color.border,  # Set border color
     edge.color = "grey", 
     layout = layout_fr, 
     main = "")
dev.off()


