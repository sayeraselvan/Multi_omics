##### PCA plot:
xfbe <- fread("D:/Rawdata/af/Base100k_base16k_tb20_bc40_ts20_ts40_ts60_ts80_mpileup_flt.xf.txt")
covbe  <- fread("D:/Rawdata/af/Base100k_base16k_tb20_bc40_ts20_ts40_ts60_ts80_mpileup_flt.cov.txt")
afbe <- fread("D:/Rawdata/af/Base100k_base16k_tb20_bc40_ts20_ts40_ts60_ts80_mpileup_flt.af.txt")
            min(covbe_filtered$V11)
quantile(covbe_filtered$V11, probs = 0.01)

covbe_filtered <- covbe %>%
  mutate(V5 = ifelse(V5 < quantile(V5, probs = 0.01) | V5 > quantile(V5, probs = 0.99), NA, V5)) %>%
  mutate(V6 = ifelse(V6 < quantile(V6, probs = 0.01) | V6 > quantile(V6, probs = 0.99), NA, V6)) %>%
  mutate(V7 = ifelse(V7 < quantile(V7, probs = 0.01) | V7 > quantile(V7, probs = 0.99), NA, V7)) %>%
  mutate(V8 = ifelse(V8 < quantile(V8, probs = 0.01) | V8 > quantile(V8, probs = 0.99), NA, V8)) %>%
  mutate(V9 = ifelse(V9 < quantile(V9, probs = 0.01) | V9 > quantile(V9, probs = 0.99), NA, V9)) %>%
  mutate(V10 = ifelse(V10 < quantile(V10, probs = 0.01) | V10 > quantile(V10, probs = 0.99), NA, V10)) %>%
  mutate(V11 = ifelse(V11 < quantile(V11, probs = 0.01) | V11 > quantile(V11, probs = 0.99), NA, V11)) %>%
  mutate(V12 = ifelse(V12 < quantile(V12, probs = 0.01) | V12 > quantile(V12, probs = 0.99), NA, V12))

covbe_filtered <- na.omit(covbe_filtered)

# covbe_filtered <- covbe %>% filter(V5 >= 50 & V5 <= 200)
# covbe_filtered <- covbe_filtered %>% filter(V6 >= 50 & V6 <= 200)
# covbe_filtered <- covbe_filtered %>% filter(V7 >= 50 & V7 <= 200)
# covbe_filtered <- covbe_filtered %>% filter(V8 >= 50 & V8 <= 200)
# covbe_filtered <- covbe_filtered %>% filter(V9 >= 50 & V9 <= 200)
# covbe_filtered <- covbe_filtered %>% filter(V10 >= 50 & V10 <= 200)
# covbe_filtered <- covbe_filtered %>% filter(V11 >= 50 & V11 <= 200) 29297659

keys_cov <- covbe_filtered[, c("V1", "V2")]
afbe_filtered <- merge(keys_cov, afbe, by = c("V1", "V2"))
xfbe_filtered <- merge(keys_cov, xfbe, by = c("V1", "V2"))

head(afbe_filtered)

afbe_filtered_0.05 <- apply(afbe_filtered[, 5:12], 1, function(row) any(row >= 0.05))

# Identify rows where all values are 0 or 1
afbe_all_zero_or_one <- rowSums(afbe_filtered[, 5:12] == 0) == ncol(afbe_filtered[, 5:12]) |
  rowSums(afbe_filtered[, 5:12] == 1) == ncol(afbe_filtered[, 5:12])

# Filter the data frame to include only those rows
filtered_afbe <- afbe_filtered[afbe_all_zero_or_one, ]


afbe_filtered_step1 <- afbe_filtered[afbe_filtered_0.05, ]
afbe_filtered_0.95 <- apply(afbe_filtered_step1[, 5:12], 1, function(row) any(row <= 0.95))
afbe_filtered_final <- afbe_filtered_step1[afbe_filtered_0.95,]
setnames(afbe_filtered_final2, c("chrom", "pos", "ref", "alt", "base100k", "base16k", "tb20", "bc40", "ts20", "ts40", "ts60", "ts80"))
setnames(xfbe_filtered, c("chrom", "pos", "ref", "alt", "base100k", "base16k", "tb20", "bc40", "ts20", "ts40", "ts60", "ts80"))

keys_afbe_filtered_final <- afbe_filtered_final[, c("chrom", "pos")]
xfbe_filtered_final <- merge(keys_afbe_filtered_final, xfbe_filtered, by = c("chrom", "pos"))

## normalise the af
afbe_filtered_final_subset <- afbe_filtered_final2[, 5:12]
afbe_filtered_final_subset[] <- lapply(afbe_filtered_final_subset, as.numeric)
afbe_filtered_final_transposed <- t(afbe_filtered_final_subset)
afbe_normalized_transposed <- scale(afbe_filtered_final_transposed)
pca_result <- prcomp(afbe_normalized_transposed, center = TRUE, scale. = TRUE)
summary(pca_result)
pca_components <- pca_result$x
screeplot(pca_result, main = "Scree Plot")
biplot(pca_result, main = "PCA Biplot")

library(ggplot2)

# Create the PCA plot
pca_df <- data.frame(PC1 = pca_components[, 1], PC2 = pca_components[, 2])

variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
variance_explained

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = PC1), size = 3, alpha = 0.7) +  # Points with color based on PC1
  labs(
    title = paste("PCA: First vs Second Principal Components\n",
                  "PC1 (", round(variance_explained[1] * 100, 2), "%) vs ",
                  "PC2 (", round(variance_explained[2] * 100, 2), "%)"),
    x = paste("PC1 (", round(variance_explained[1] * 100, 2), "%)", sep = ""),
    y = paste("PC2 (", round(variance_explained[2] * 100, 2), "%)", sep = "")
  ) +
  geom_text_repel(aes(label = rownames(pca_df)), size = 3, box.padding = 0.35, point.padding = 0.5) +  # Avoid label collision
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient for points
  theme_minimal() +  # Minimal theme for cleaner look
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Centered title
    axis.title = element_text(size = 12),  # Axis labels size
    axis.text = element_text(size = 10),  # Axis tick label size
    legend.position = "none"  # Remove legend
  )
library("ggrepel")