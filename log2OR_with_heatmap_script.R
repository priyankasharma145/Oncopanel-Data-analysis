# Load required libraries
library(dplyr)
library(ggplot2)

# Load your binary matrix with patient IDs as row names
binary_matrix <- read.csv("dfcitop10.csv", row.names = 1)

# Convert the matrix to binary (1 if value > 0), keeping row names intact
binary_matrix[] <- lapply(binary_matrix, function(col) ifelse(col != 0, 1, 0))

# Extract gene names
gene_names <- colnames(binary_matrix)

# Prepare empty results container
results_df <- data.frame(
  Gene1 = character(),
  Gene2 = character(),
  AdjustedLog2OddsRatio = numeric(),
  FisherEstimate = numeric(),
  FisherPValue = numeric(),
  stringsAsFactors = FALSE
)

# Loop over gene pairs
for (i in seq_along(gene_names)) {
  for (j in seq_along(gene_names)) {
    if (i != j) {
      gene1 <- gene_names[i]
      gene2 <- gene_names[j]
      
      # Build contingency table
      tab <- table(binary_matrix[[gene1]], binary_matrix[[gene2]])
      
      # Only proceed if table is 2x2
      if (all(dim(tab) == c(2, 2))) {
        # Apply continuity correction
        tab <- tab + 0.5
        
        fisher_test <- fisher.test(tab)
        
        adjusted_log2_odds_ratio <- log2((tab[1,1] * tab[2,2]) / (tab[1,2] * tab[2,1]))
        
        results_df <- rbind(results_df, data.frame(
          Gene1 = gene1,
          Gene2 = gene2,
          AdjustedLog2OddsRatio = adjusted_log2_odds_ratio,
          FisherEstimate = fisher_test$estimate,
          FisherPValue = fisher_test$p.value
        ))
      }
    }
  }
}

# View or save results
write.csv(results_df, "gene_pair_analysis_results.csv", row.names = FALSE)


# Save results to a file named "adjusted_log2_FE_testcbio.csv"
write.csv(results_df, file = "dfciadjusted_log2_FE_test_new0.5.csv", row.names = FALSE)

# Create a heatmap with normal coloring
heatmap_plot <- ggplot(results_df, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = AdjustedLog2OddsRatio), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-10, 10)) +
  geom_text(data = results_df %>% filter(FisherPValue < 0.05),
            aes(label = sprintf("%.2f", AdjustedLog2OddsRatio)),
            color = "black", size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white")
  ) +
  labs(title = "Adjusted Co-Occurrence Heatmap", fill = "Adjusted Log2 Odds Ratio")

# Print and save the heatmap as an image
ggsave("combinedadjusted_co_occurrence_heatmap_dfci_2_final_new.png", heatmap_plot, width = 10, height = 8, dpi = 300)
# Create a unique identifier for each pair (to remove duplicates in the plot)
results_df <- results_df %>%
  mutate(Pair = paste0(pmin(Gene1, Gene2), "_", pmax(Gene1, Gene2))) %>%
  distinct(Pair, .keep_all = TRUE)  # Keep only unique pairs

# Create a triangular-style heatmap
heatmap_plot <- ggplot(results_df, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = AdjustedLog2OddsRatio), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-10, 10)) +
  geom_text(data = results_df %>% filter(FisherPValue < 0.05),
            aes(label = sprintf("%.2f", AdjustedLog2OddsRatio)),
            color = "black", size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white")
  ) +
  labs(title = "Adjusted Co-Occurrence Heatmap (Unique Pairs)", fill = "Adjusted Log2 Odds Ratio")

# Save the heatmap
ggsave("combined_adjusted_co_occurrence_heatmap_half.png", heatmap_plot, width = 10, height = 8, dpi = 300)



heatmap_plot <- ggplot(results_df, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = AdjustedLog2OddsRatio), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-10, 10)) +
  geom_text(data = results_df %>% filter(FisherPValue < 0.05),
            aes(label = sprintf("%.2f", AdjustedLog2OddsRatio)),
            color = "black", size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0, size = 16, vjust = 0.5),  # Move labels up
    axis.text.y = element_text(size = 16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white"),
    axis.ticks.length = unit(0, "cm"),  # Hide x-axis ticks
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  # Adjust top margin
  ) +
  scale_x_discrete(position = "top") +  # Move x-axis labels to the top
  labs(title = "Adjusted Co-Occurrence Heatmap", fill = "Adjusted Log2 Odds Ratio")

# Save the updated heatmap
ggsave("combined_adjusted_co_occurrence_heatmap_top_labels.png", heatmap_plot, width = 10, height = 8, dpi = 300)

# Create a unique identifier for each pair (to remove duplicates in the plot)
results_df <- results_df %>%
  mutate(Pair = paste0(pmin(Gene1, Gene2), "_", pmax(Gene1, Gene2))) %>%
  distinct(Pair, .keep_all = TRUE)  # Keep only unique pairs

# Create a triangular-style heatmap
heatmap_plot <- ggplot(results_df, aes(x = Gene1, y = Gene2)) +
  geom_tile(aes(fill = AdjustedLog2OddsRatio), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-4, 8)) +
  geom_text(data = results_df %>% filter(FisherPValue < 0.05),
            aes(label = sprintf("%.2f", AdjustedLog2OddsRatio)),
            color = "black", size = 3) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white")
  ) +
  labs(title = "Adjusted Co-Occurrence Heatmap (Unique Pairs)", fill = "Adjusted Log2 Odds Ratio")

# Save the heatmap
ggsave("combined_adjusted_co_occurrence_heatmap_half.png", heatmap_plot, width = 10, height = 8, dpi = 300)
