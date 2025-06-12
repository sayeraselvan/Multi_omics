# Load required libraries
library(readxl)      # For reading Excel files
library(dplyr)       # For data manipulation
library(tidyr)       # For data wrangling
library(tm)          # For text mining
library(cluster)     # For clustering
library(ggplot2)     # For visualization
library(factoextra)  # For PCA and clustering visualization

# Step 1: Load the data
file_path <- "~/Downloads/FAD_Interaction_Modes_Detailed_Grouped_by_PDB_and_Chain.xlsx"
data <- read_excel(file_path)

interaction_details <- data$`FAD Interaction Details`
parsed_interactions <- str_split(interaction_details, ", ")

# Extract motif atoms, interacting atoms, and bond lengths
extract_features <- function(interactions) {
  parsed <- str_match(interactions, "([A-Za-z0-9']+) - ([A-Za-z0-9']+) - ([0-9.]+)")
  motif_atoms <- parsed[, 2]          # Column 2: Motif atoms (e.g., O1A)
  interacting_atoms <- parsed[, 3]   # Column 3: Interacting atoms (e.g., N)
  bond_lengths <- as.numeric(parsed[, 4])  # Column 4: Bond lengths (numeric)
  data.frame(motif_atoms, interacting_atoms, bond_lengths)
}

# Parse all interactions and combine them into a single data frame
all_interactions <- do.call(rbind, lapply(parsed_interactions, extract_features))

# Add chain identifier to the parsed interactions
chain_ids <- rep(paste(data$`PDB ID`, data$Chain, sep = "_"), lengths(parsed_interactions))
all_interactions$chain <- chain_ids

# Split the interaction details into individual interactions
interaction_details <- data$`FAD Interaction Details`
parsed_interactions <- str_split(interaction_details, ", ")

# Extract motif atom, interacting atom, and bond length
extract_features <- function(interactions) {
  parsed <- str_match(interactions, "([A-Za-z0-9']+) - ([A-Za-z0-9']+) - ([0-9.]+)")
  motif_atoms <- parsed[, 2]          # Column 2: Motif atoms (e.g., O1A)
  interacting_atoms <- parsed[, 3]   # Column 3: Interacting atoms (e.g., N)
  bond_lengths <- as.numeric(parsed[, 4])  # Column 4: Bond lengths (numeric)
  data.frame(motif_atoms, interacting_atoms, bond_lengths)
}

# Parse all interactions and combine them into a single data frame
all_interactions <- do.call(rbind, lapply(parsed_interactions, extract_features))

# Add chain identifier to the parsed interactions
chain_ids <- rep(paste(data$`PDB ID`, data$Chain, sep = "_"), lengths(parsed_interactions))
all_interactions$chain <- chain_ids

# Step 3: Normalize the bond lengths
all_interactions$normalized_bond_length <- scale(all_interactions$bond_lengths)

aggregated_interactions <- all_interactions %>%
  unite("interaction_detail", motif_atoms, interacting_atoms, bond_lengths, sep = "_", remove = FALSE) %>%
  group_by(chain) %>%
  summarize(
    interaction_list = paste(interaction_detail, collapse = ";") # Combine all interactions into one string
  )

# Step 2: Vectorize the interactions
# Create a sparse matrix with interaction details
interaction_matrix <- all_interactions %>%
  unite("interaction_pair", motif_atoms, interacting_atoms, sep = "_") %>%
  group_by(chain, interaction_pair) %>%
  summarize(mean_bond_length = mean(normalized_bond_length), .groups = "drop") %>%
  pivot_wider(
    names_from = interaction_pair,
    values_from = mean_bond_length,
    values_fill = 0
  )

# Step 3: Perform clustering
# Convert interaction_matrix to numeric matrix
interaction_matrix_numeric <- interaction_matrix %>%
  select(-chain) %>%
  as.matrix()

# Apply hierarchical clustering
dist_matrix <- dist(interaction_matrix_numeric) # Calculate distance matrix
hc <- hclust(dist_matrix, method = "ward.D2")    # Hierarchical clustering

# Plot the dendrogram
plot(hc, labels = interaction_matrix$chain, main = "Clustering of Chains", xlab = "Chain", ylab = "Height")

## PCA
interaction_matrix_scaled <- scale(interaction_matrix_numeric)

pca_result <- prcomp(interaction_matrix_scaled, center = TRUE, scale. = TRUE)
summary(pca_result)  

# Extract PCA scores (principal components)
pca_scores <- as.data.frame(pca_result$x)

# Add chain identifiers to the PCA scores
pca_scores$chain <- interaction_matrix$chain

# Step 4: Visualize the first two principal components
ggplot(pca_scores, aes(x = PC1, y = PC2, label = chain)) +
  geom_point(aes(color = chain), size = 3) +
  geom_text(vjust = 1.5, size = 3) +
  labs(title = "PCA of Interaction Features", x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal()


### label based on the chain
pca_scores$chain_type <- sub(".*_(.)", "\\1", pca_scores$chain)

ggplot(pca_scores, aes(x = PC1, y = PC2, color = chain_type)) +
  geom_point(size = 3) +
  #geom_text(vjust = 1.5, size = 3) +
  labs(title = "PCA of Interaction Features", x = "Principal Component 1", y = "Principal Component 2") +
  #scale_color_manual(values = c("A" = "blue", "B" = "red", "C" = "green", "M" = "purple", "E" = "orange")) +  # Custom colors for chains A, B, C, M, E
  theme_minimal()
# ggplot(pca_scores %>% filter(chain_type == "K"), aes(x = PC1, y = PC2, color = chain_type)) +
#   geom_point(size = 3) +
#   labs(title = "PCA of Interaction Features", x = "Principal Component 1", y = "Principal Component 2") +
#   theme_minimal()




