# Load necessary libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
library(dplyr)
library(tidyr)

# Read your data
# Replace "your_file.csv" with your actual file path
data <- read.csv("comut_data_mucosal_n84.csv", header = TRUE, stringsAsFactors = FALSE)

# Convert the data from long to wide format
wide_data <- data %>%
  # Group by sample and category to ensure proper aggregation
  group_by(sample, category) %>%
  # Combine all values for the same sample-category pair into a single string
  summarize(combined = paste(value, collapse = ", "), .groups = "drop") %>%
  # Spread the data into wide format
  pivot_wider(names_from = category, values_from = combined, values_fill = "no alteration") %>%
  # Rename the "sample" column to "Patient ID"
  rename(`Patient ID` = sample)

# View the transformed data
print(wide_data)

# Optionally, save the transformed data to a CSV file
write.csv(wide_data, "transformed_data.csv", row.names = FALSE)
