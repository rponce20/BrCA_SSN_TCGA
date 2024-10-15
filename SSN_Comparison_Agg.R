# Code to compare link similarity (10K) between aggregated networks versus 
# unique networks for each phenotype-subtype (Normal/Luminal A, Luminal B, 
# Her2, Basal), to determine heterogeneity.
#
# Input: For each phenotype and subtype
# - Aggregated network files (.tsv)
# - Unique network files (.tsv)
#
# Output: 
# - Plot (~facet) for each comparison - percentage
# - Plot (~facet) for Jaccard Index
################################################################################

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Define function reorder_within to reorder elements within facets
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun, ...)
}

# Define function scale_x_reordered to scale x-axis within facets
scale_x_reordered <- function(...) {
  scale_x_discrete(labels = function(x) gsub("___.*", "", x), ...)
}

# Function to read a network from a .tsv file
read_network <- function(file_path) {
  read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
}

# Function to calculate the proportion of common pairs and the Jaccard Index
compare_networks <- function(aggregated, unique) {
  # Use the first two columns of each network regardless of their names
  aggregated_pairs <- with(aggregated, paste(aggregated[[1]], aggregated[[2]], sep = "-"))
  unique_pairs <- with(unique, paste(unique[[1]], unique[[2]], sep = "-"))
  
  # Calculate common pairs
  common_pairs <- intersect(aggregated_pairs, unique_pairs)
  
  # Calculate the proportion of common pairs
  proportion_common <- length(common_pairs) / length(unique_pairs)
  
  # Calculate Jaccard Index
  jaccard_index <- length(common_pairs) / (length(aggregated_pairs) + length(unique_pairs) - length(common_pairs))
  
  return(list(proportion_common = proportion_common, jaccard_index = jaccard_index))
}

# Function to process a complete subtype
process_subtype <- function(subtype_name, aggregated_network_path, unique_networks_dir) {
  # Read aggregated network
  aggregated_network <- read_network(aggregated_network_path)
  
  # List unique network files in the selected directory
  unique_network_files <- list.files(path = unique_networks_dir, pattern = "^TCGA-.*\\.tsv$", full.names = TRUE)
  
  # Read all unique networks
  unique_networks <- lapply(unique_network_files, read_network)
  names(unique_networks) <- basename(unique_network_files)
  
  # Apply comparison to all unique networks
  comparison_results <- lapply(unique_networks, compare_networks, aggregated = aggregated_network)
  comparison_results_df <- data.frame(
    Sample = names(unique_networks),
    Proportion_Common = sapply(comparison_results, `[[`, "proportion_common"),
    Jaccard_Index = sapply(comparison_results, `[[`, "jaccard_index")
  )
  comparison_results_df$Subtype <- subtype_name
  
  return(comparison_results_df)
}

# Select files for each subtype
cat("Select the aggregated network files for each subtype.\n")

subtype_files <- list()
subtypes <- c("Normal", "Luminal A", "Luminal B", "Her2", "Basal")

for (subtype in subtypes) {
  cat(paste("Select the aggregated network file for", subtype, "\n"))
  aggregated_network_path <- file.choose()
  cat(paste("Select a file in the directory of unique networks for", subtype, "\n"))
  unique_networks_dir <- dirname(file.choose())
  
  subtype_files[[subtype]] <- list(aggregated = aggregated_network_path, unique_dir = unique_networks_dir)
}

# Process all subtypes
all_results <- lapply(names(subtype_files), function(subtype) {
  process_subtype(subtype, subtype_files[[subtype]]$aggregated, subtype_files[[subtype]]$unique_dir)
})

# Combine all results into a single dataframe
combined_results_df <- bind_rows(all_results)

# Sort results by descending order within each subtype for proportion of common pairs
combined_results_df <- combined_results_df %>%
  mutate(Subtype = factor(Subtype, levels = c("Normal", "Luminal A", "Luminal B", "Her2", "Basal"))) %>%
  group_by(Subtype) %>%
  arrange(Subtype, desc(Proportion_Common), .by_group = TRUE)

# Sort results by descending order within each subtype for Jaccard Index
combined_results_df <- combined_results_df %>%
  group_by(Subtype) %>%
  arrange(Subtype, desc(Jaccard_Index), .by_group = TRUE)

# Violin plot for Jaccard Index
ggplot(combined_results_df, aes(x = Subtype, y = Jaccard_Index, fill = Subtype)) +
  geom_violin(trim = FALSE) +  # Add density with violin plot
  geom_boxplot(width = 0.1) +  # Overlay boxplot inside violin
  theme_minimal(base_size = 16) +  
  theme(
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(),  # Remove x-axis ticks
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5), 
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_blank()  
  ) +
  labs(title = "", 
       y = "Jaccard Index") +  # Keep only y-axis label
  scale_fill_manual(values = c("Normal" = "#E69F00", "Luminal A" = "#56B4E9", 
                               "Luminal B" = "#009E73", "Her2" = "coral3", 
                               "Basal" = "#CC79A7"))

# Save plot
ggsave("violin_jaccard_index_updated.jpeg", width = 10, height = 10, dpi = 300, device = "jpeg")
