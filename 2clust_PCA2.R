# Load required libraries
library(ggplot2)
library(tidyverse)
library(cluster)
library(factoextra)
library(ggrepel)


# Read the data from the CSV file
data <- read.csv("transposed_genietop5.csv", row.names = 1)

# Check for zero variance columns and remove them
zero_variance_cols <- sapply(data, function(col) var(col) == 0)
data_filtered <- data[, !zero_variance_cols]

# Perform PCA on the filtered data
pca <- prcomp(data_filtered, scale. = TRUE)

# Summary of PCA
summary(pca)

# Plot PCA results (2D plot of first two principal components)
pca_df <- data.frame(pca$x)

# Perform clustering (Hierarchical Clustering)
diss_matrix <- dist(data_filtered)  # Distance matrix
hclust_res <- hclust(diss_matrix, method = "ward.D2")

# Cut tree to create 4 clusters (You can adjust the number of clusters here)
clusters <- cutree(hclust_res, k = 6)

# Add cluster labels to the PCA data
pca_df$Cluster <- as.factor(clusters)
pca_df$Sample <- rownames(data)  # Add sample names for labeling

# Save PCA plot with clusters and labels to PDF
pdf("PCA_Clustered_Labeled_Plot.pdf")
pca_cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +  # Add labels with slight offset
  labs(title = "PCA with Clusters and Labels", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "purple", "yellow3", "cyan4")) # Cluster colors
print(pca_cluster_plot)
dev.off()

# Save horizontal dendrogram plot to PDF
pdf("Horizontal_Dendrogram_Plot.pdf")
hclust_plot <- fviz_dend(hclust_res, 
                         k = 4,  # Adjust number of clusters
                         rect = TRUE, 
                         rect_fill = TRUE, 
                         rect_border = "black", 
                         main = "Hierarchical Clustering Dendrogram", 
                         xlab = "Samples",
                         horiz = TRUE, 
                         cex = 0.8,  # Label size
                         lwd = 1.2,  # Branch line width
                         k_colors = c("darkred", "darkblue", "darkgreen", "cyan")) +  # FIXED here
  theme(
    plot.title = element_text(size = 6, face = "bold")
  )

print(hclust_plot)
dev.off()

# Save cluster assignments to a CSV file
cluster_df <- data.frame(Sample = rownames(data), Cluster = clusters)
write.csv(cluster_df, "Cluster_Assignments.csv", row.names = FALSE)

# Tabular representation of genes in each cluster
cluster_genes <- lapply(unique(clusters), function(cluster_num) {
  rownames(data)[clusters == cluster_num]
})

# Save the genes in each cluster as CSV files
for (i in unique(clusters)) {
  cluster_genes_df <- data.frame(Gene = cluster_genes[[i]])
  write.csv(cluster_genes_df, paste0("Cluster_", i, "_Genes.csv"), row.names = FALSE)
}



# Save PCA plot with clusters and labels to PDF
pdf("PCA_Clustered_Labeled_Plot2.pdf")

pca_cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 4) +  # Increased point size
  geom_text_repel(
    size = 4, 
    box.padding = 0.5,  # More spacing around text
    point.padding = 0.3, # Space from the dots
    segment.size = 0.5,  # Thicker lines for visibility
    segment.color = "black",
    arrow = arrow(length = unit(0.03, "npc")), # Larger arrow size
    max.overlaps = Inf  # Show all labels, no filtering
  ) +  
  labs(title = "PCA with Clusters and Labels", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  ) +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "purple", "yellow3", "cyan4"))  # Cluster colors

print(pca_cluster_plot)
dev.off()


# Calculate the sum of 1s for each gene (prevalence)
gene_prevalence <- rowSums(data_filtered)

# Add prevalence to the PCA dataframe
pca_df$Prevalence <- gene_prevalence
pdf("PCA_Clustered_Labeled_Plot4.pdf")

# Make sure you have 6 clusters
pca_cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample, size = Prevalence)) +
  geom_point() +
  geom_text_repel(
    size = 4, 
    box.padding = 0.5,  # More spacing around text
    point.padding = 0.3, # Space from the dots
    segment.size = 0.5,  # Thicker lines for visibility
    segment.color = "black",
    arrow = arrow(length = unit(0.03, "npc")), # Larger arrow size
    max.overlaps = Inf  # Show all labels, no filtering
  ) +
  labs(title = "PCA with Clusters, Labels, and Prevalence-Adjusted Dot Size", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  ) +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "purple", "yellow3", "cyan4")) +  # 6 colors for 6 clusters
  scale_size_continuous(range = c(2, 10))  # Adjust the range for dot size based on prevalence

# Print the plot
print(pca_cluster_plot)
dev.off()

# Save PCA plot with clusters and labels to PDF
pdf("PCA_Clustered_Labeled_Plot.pdf")
pca_cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample, size = Prevalence)) +
  geom_point() +  # Adjusted size based on prevalence
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +  # Add labels with slight offset
  labs(title = "PCA with Clusters and Labels", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "purple", "yellow3", "cyan4")) +  # 6 colors for 6 clusters
  scale_size_continuous(range = c(2, 10))  # Adjust the range for dot size based on prevalence

# Print the plot
print(pca_cluster_plot)

# Close the PDF device
dev.off()




# Save PCA plot with clusters, labels, and adjusted dot size to PDF
pdf("PCA_Clustered_Labeled_Plot3.pdf")

pca_cluster_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample, size = Prevalence)) +
  geom_point() +
  geom_text_repel(
    size = 4, 
    box.padding = 0.5,  # More spacing around text
    point.padding = 0.3, # Space from the dots
    segment.size = 0.5,  # Thicker lines for visibility
    segment.color = "black",
    arrow = arrow(length = unit(0.03, "npc")), # Larger arrow size
    max.overlaps = Inf  # Show all labels, no filtering
  ) +
  labs(title = "PCA with Clusters, Labels, and Prevalence-Adjusted Dot Size", x = "PC1", y = "PC2") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16)
  ) +
  scale_color_manual(values = c("darkred", "darkblue", "darkgreen", "purple", "yellow3", "cyan4")) +  # Cluster colors
  scale_size_continuous(range = c(2, 10))  # Adjust the range for dot size based on prevalence

print(pca_cluster_plot)
dev.off()

