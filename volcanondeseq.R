# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(grid) # For unit function
library(biomaRt)

# Read the CSV file
data <- read.csv("ptRNA_input_cdk13_vs_clover_Log2FC.csv")
data <- na.omit(data)

# Extract Ensembl IDs without version numbers and chromosome location
data$EnsemblID <- sub("\\..*", "", data$GeneID)
data$ChrLocation <- sub("^[^.]+\\.", "", data$GeneID)

# Define significance thresholds
padj_threshold <- 0.05

# Filter significant genes based on adjusted p-value
significant_genes <- data %>% filter(padj < padj_threshold)

# Get top 20 highest and lowest log2FoldChange among significant genes
top_genes <- significant_genes %>%
  arrange(desc(log2FoldChange)) %>%
  slice(1:20)

bottom_genes <- significant_genes %>%
  arrange(log2FoldChange) %>%
  slice(1:20)

# Combine top and bottom genes for labeling
label_genes <- bind_rows(top_genes, bottom_genes)

# Set up the Ensembl connection
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl IDs to gene symbols for top and bottom genes
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = label_genes$EnsemblID,
                   mart = ensembl)

# Merge gene symbols back into the label_genes data
label_genes <- merge(label_genes, gene_info, by.x = "EnsemblID", by.y = "ensembl_gene_id", all.x = TRUE)

# Create the volcano plot
volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(padj < padj_threshold, 
                                ifelse(log2FoldChange > 0, "positive", "negative"), 
                                "non-significant")), 
             alpha = 0.8, size = 1, shape = 16) + # Solid dots with higher alpha
  scale_color_manual(values = c("positive" = "red", "negative" = "darkgreen", "non-significant" = "grey")) +
  geom_text_repel(data = label_genes, aes(label = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", 
                                                         GeneID, 
                                                         paste(hgnc_symbol, ChrLocation, sep = "_")), 
                                          fontface = ifelse(padj < padj_threshold, "bold", "plain")), 
                  nudge_x = ifelse(label_genes$log2FoldChange > 0, 5, -5), # Nudge right for positive, left for negative
                  direction = "y", # Align labels vertically
                  hjust = ifelse(label_genes$log2FoldChange > 0, 0, 1), # Left-align for positive, right-align for negative
                  segment.size = 0.2, # Size of the connecting line
                  segment.color = "black", # Color of the connecting line
                  arrow = arrow(type = "closed", length = unit(0.1, "cm"), ends = "last"), # Add round arrowhead
                  size = 2.0) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme(legend.position = "none") # Remove legend for clarity

# Save the plot as a PDF
pdf("ptRNA_input_cdk13_vs_clover_5padj_volcano_plot.pdf")
print(volcano_plot)
dev.off()

# Save significant genes to a CSV file
write.csv(significant_genes, "ptRNA_input_cdk13_vs_clover_5padj_significant_genes.csv", row.names = FALSE)