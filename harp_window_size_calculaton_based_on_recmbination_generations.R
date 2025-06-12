library(dplyr)

L = 27160941 # Chromosome length in bp
R = 0.0000000304
R = 0.00000

G = 75

Q = 18

# Calculate the scaling factor
scale_factor = L / (R * L * G + 1)

# Calculate the quantile
window_size_bp = qexp(Q / 100, rate = 1 / scale_factor)

# Convert to kb
window_size_kb = window_size_bp / 1000
window_size_bp
window_size_kb
/Users/siva/Desktop/Rawdata/cluster_reports/harp/84_samples_maf_filtration/BC/3R_1mbwindow/BC_R4_F40_trimmedG_sorted_rmdup_filtered_sc.bam.3R.freqs
# 1mb window with 100kb sliding window this BC and base100k is based on the ref from signor
BC_3R <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/178samples_fil/RR_33_MKDN240015467-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.3R.freqs",
                    header = FALSE, 
                    sep = " ")
# Remove the last column from BC_3R cuz it is NA column

BC_3R <- BC_3R[, -ncol(BC_3R)]

snp_table_3R <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/178samples_fil/snp_table_3R.npute", 
                           header = TRUE, 
                           sep = ",")
new_colnames <- colnames(snp_table_3R)
#column name change according to the sample name for the freq file
samples_names <- new_colnames[-(1:2)]
new_colnames <- c("window", "start", "end")
final_colnames <- c(new_colnames, samples_names)
colnames(BC_3R) <- final_colnames
colnames(BC_3R)

### lets try to plot

library(ggplot2)
library(reshape2)

# Reshape the data for plotting
BC_3R_long <- melt(BC_3R, id.vars = c("window", "start", "end"), 
                   variable.name = "haplotype", 
                   value.name = "frequency")

# Plot
ggplot(BC_3R_long, aes(x = start, y = frequency, color = haplotype)) +
  geom_line() +
  labs(title = "Haplotype Frequencies Across 3R",
       x = "Start Position",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.position = "right")


### heatmap
library(ggplot2)
library(reshape2)

ggplot(BC_3R_long, aes(x = start, y = haplotype, fill = frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Haplotype Frequencies",
       x = "Start Position",
       y = "Haplotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


#####
BC_3R1 <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/84_samples/3R_1mbwindow/BC_R4_F40_trimmedG_sorted_rmdup_filtered_sc.bam.3R.freqs",
                     header = FALSE, 
              sep = " ")
# Remove the last column from BC_3R cuz it is NA column
BC_3R1 <- BC_3R1[, -ncol(BC_3R1)]

colnames(BC_3R1) <- final_colnames

# Reshape the data for plotting
BC_3R_long1 <- melt(BC_3R1, id.vars = c("window", "start", "end"), 
                   variable.name = "haplotype", 
                   value.name = "frequency")


### heatmap
library(ggplot2)
library(reshape2)

# Create a heatmap
ggplot(BC_3R_long1, aes(x = start, y = haplotype, fill = frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Haplotype Frequencies",
       x = "Start Position",
       y = "Haplotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

###  84 samples maf filtration 
BC_3R2 <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/84_samples_maf_filtration/BC/3R/BC_R4_F40_trimmedG_sorted_rmdup_filtered_sc.bam.3R.freqs",
                    header = FALSE, 
                     sep = " ")
# Remove the last column from BC_3R cuz it is NA column
BC_3R2 <- BC_3R2[, -ncol(BC_3R2)]

colnames(BC_3R2) <- final_colnames

# Reshape the data for plotting
BC_3R_long2 <- melt(BC_3R2, id.vars = c("window", "start", "end"), 
                    variable.name = "haplotype", 
                    value.name = "frequency")

### heatmap
library(ggplot2)
library(reshape2)

# Create a heatmap
ggplot(BC_3R_long2, aes(x = start, y = haplotype, fill = frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Haplotype Frequencies",
       x = "Start Position",
       y = "Haplotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


#### base 100k
base100k <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/84_samples_maf_filtration/base100k/3R_1mbwindow/Base100k_sorted_rmdup_filtered_sc.bam.3R.freqs",
                     header = FALSE, 
                     sep = " ")
# Remove the last column from BC_3R cuz it is NA column
base100k <- base100k[, -ncol(base100k)]

colnames(base100k) <- final_colnames

# Reshape the data for plotting
base100k_long <- melt(base100k, id.vars = c("window", "start", "end"), 
                    variable.name = "haplotype", 
                    value.name = "frequency")

### heatmap
library(ggplot2)
library(reshape2)

# Create a heatmap
ggplot(base100k_long, aes(x = start, y = haplotype, fill = frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Haplotype Frequencies",
       x = "Start Position",
       y = "Haplotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))



###### for 31.10.2024 presentation
snp_table_2L <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/snp_table_2L", 
                           header = TRUE, 
                           sep = ",")

new_colnames2L <- colnames(snp_table_2L)
#column name change according to the sample name for the freq file
samples_names_2L <- new_colnames2L[-(1:2)]
new_colnames <- c("window", "start", "end")
final_colnames_2L <- c(new_colnames, samples_names_2L)


DSIM_12 <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/1mbwindow/filtered_2L_DSIM_12.bam.2L.freqs",
                    header = FALSE, 
                    sep = " ")
RR33 <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/RR_1mbwindow/RR_33_MKDN240015467-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.2L.freqs",
                   header = FALSE, 
                   sep = " ")
RR34 <- read.table("~/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/RR_34_1mbwindow/RR_34_MKDN240015468-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.2L.freqs",
                   header = FALSE, 
                   sep = " ")
RR90 <- read.table("~/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/RR_90_1mbwindow/RR_90_MKDN240015469-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.2L.freqs",
                   header = FALSE, 
                   sep = " ")
Sz225 <- read.table("~/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/Sz225_1mbwindow/Sz225_pass_trimmed_sorted_rmdup_filtered.bam.2L.freqs",
                    header = FALSE, 
                    sep = " ")
RR91 <- read.table("~/Desktop/Rawdata/cluster_reports/harp/31_10_2024/harp/90samples/RR_91_1mbwindow/RR_91_MKDN240015470-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.2L.freqs",
           header = FALSE, 
           sep = " ")

DSIM_12 <- DSIM_12[, -ncol(DSIM_12)]
colnames(DSIM_12) <- final_colnames_2L
RR33 <- RR33[, -ncol(RR33)]
colnames(RR33) <- final_colnames_2L
RR90 <- RR90[, -ncol(RR90)]
colnames(RR90) <- final_colnames_2L
RR91 <- RR91[, -ncol(RR91)]
colnames(RR91) <- final_colnames_2L
Sz225 <- Sz225[, -ncol(Sz225)]
colnames(Sz225) <- final_colnames_2L
RR34 <- RR34[, -ncol(RR34)]
colnames(RR34) <- final_colnames_2L

DSIM_12_long2 <- melt(DSIM_12, id.vars = c("window", "start", "end"), 
                    variable.name = "haplotype", 
                    value.name = "frequency")
RR33_long <- melt(RR33, id.vars = c("window", "start", "end"), 
                  variable.name = "haplotype", 
                  value.name = "frequency")

RR90_long <- melt(RR90, id.vars = c("window", "start", "end"), 
                  variable.name = "haplotype", 
                  value.name = "frequency")

RR91_long <- melt(RR91, id.vars = c("window", "start", "end"), 
                  variable.name = "haplotype", 
                  value.name = "frequency")

Sz225_long <- melt(Sz225, id.vars = c("window", "start", "end"), 
                   variable.name = "haplotype", 
                   value.name = "frequency")

RR34_long <- melt(RR34, id.vars = c("window", "start", "end"), 
                  variable.name = "haplotype", 
                  value.name = "frequency")

ggplot(RR33_long, aes(x = start, y = haplotype, fill = frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Heatmap of Haplotype Frequencies",
       x = "Start Position",
       y = "Haplotype") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))


##### heterozygous and missing data per sample plot:
#inbred lines 90 samples from our reference:
het_data <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/SCRIPTS/inbred_lines_90_haplotypes.het.in.het", header = TRUE)
het_data$O_HET <- het_data$N_SITES - het_data$O.HOM.
het_data$Heterozygosity_Percent <- (het_data$O_HET / het_data$N_SITES) * 100

ggplot(het_data, aes(x = INDV, y = Heterozygosity_Percent)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Heterozygosity Percentage per Sample", x = "Sample", y = "Heterozygosity (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

missing_data <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/SCRIPTS/inbred_lines_90_haplotypes.missing.in.imiss", header = TRUE)
missing_data$Missing_Percent <- missing_data$F_MISS * 100
ggplot(missing_data, aes(x = INDV, y = Missing_Percent)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(title = "Missing Data Percentage per Sample", x = "Sample", y = "Missing Data (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#184 samples from our reference:
het_data_184 <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/SCRIPTS/signor_final_184samples_gatk_all_chr_genotyped.het.in.het", header = TRUE)
het_data_184$O_HET <- het_data_184$N_SITES - het_data_184$O.HOM.
het_data_184$Heterozygosity_Percent <- (het_data_184$O_HET / het_data_184$N_SITES) * 100

ggplot(het_data_184, aes(x = INDV, y = Heterozygosity_Percent)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Heterozygosity Percentage per Sample", x = "Sample", y = "Heterozygosity (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

missing_data_184 <- read.table("/Users/siva/Desktop/Rawdata/cluster_reports/SCRIPTS/signor_final_184samples_gatk_all_chr_genotyped.missing.in.imiss", header = TRUE)
missing_data_184$Missing_Percent <- missing_data_184$F_MISS * 100
ggplot(missing_data_184, aes(x = INDV, y = Missing_Percent)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(title = "Missing Data Percentage per Sample", x = "Sample", y = "Missing Data (%)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### heterozygoisty per variant plot

file_path <- "/Users/siva/Desktop/Rawdata/cluster_reports/biallelic_snps_inbred_lines_90_heterozygosity_stats.txt"
heterozygosity <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(heterozygosity)
#heterozygosity <- heterozygosity[-1, ]
library(ggplot2)

heterozygosity$POS <- as.numeric(heterozygosity$POS)
heterozygosity$nHet <- as.numeric(heterozygosity$nHet)
heterozygosity$nhetpercent <- heterozygosity$nHet /90

chromosomes <- unique(heterozygosity$CHR)
print(chromosomes)
chrom = "X"


chrom_data <- heterozygosity[heterozygosity$CHR == chrom & heterozygosity$nHet != 0, ]

p <- ggplot(chrom_data, aes(x = POS, y = nHet)) +
    geom_point(alpha = 0.6, size = 0.5) +  
    labs(title = paste("Manhattan Plot of Heterozygous Counts - Chromosome", chrom),
         x = "Position on Chromosome",
         y = "Heterozygous Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
  
print(p)

filtered_data <- chrom_data[chrom_data$nHet %in% 0:90, ]
nhet_counts <- table(filtered_data$nHet)
print(nhet_counts)

## count table
chromosomes <- unique(heterozygosity$CHR)  # Get the list of unique chromosomes

filtered_chrom_data_list <- list()

for (chrom in chromosomes) {
  chrom_data <- heterozygosity[heterozygosity$CHR == chrom & heterozygosity$nHet != 0, ]
  filtered_chrom_data_list[[chrom]] <- chrom_data  # Store the filtered data in the list
}
filtered_2L <- nhet_counts_df[nhet_counts_df$Chromosome == "2L", ]
filtered_2R <- nhet_counts_df[nhet_counts_df$Chromosome == "2R", ]
filtered_3L <- nhet_counts_df[nhet_counts_df$Chromosome == "3L", ]
filtered_3R <- nhet_counts_df[nhet_counts_df$Chromosome == "3R", ]
filtered_X <- nhet_counts_df[nhet_counts_df$Chromosome == "X", ]

### shanon index::
# Load data
# Remove the first two columns (chromosome and position columns) 
haplotype_freqs <- base100k[, -c(1,2,3,94)]

# Function to calculate Shannon Diversity Index for a vector of frequencies
shannon_index <- function(freqs) {
  # Normalize frequencies to sum to 1 for each row
  freqs <- freqs / sum(freqs)
  # Remove zero frequencies to avoid log(0)
  freqs <- freqs[freqs > 0]
  # Calculate Shannon Index
  -sum(freqs * log(freqs))
}

# Apply the function to each row to calculate Shannon Index for each window
shannon_indices <- apply(haplotype_freqs, 1, shannon_index)
print(shannon_indices)
# Combine the results with the original data to see window positions with Shannon Index
result <- data.frame(Chromosome = base100k[,1], Position = base100k[,2], Position_end = base100k[,3], Shannon_Index = shannon_indices)

# Print the first few rows
head(result)
resultRR90 <- result
resultRR91 <- result
resultDSIM12 <- result
resultBC_3R <- result
resultbase100k <- result
# Load necessary library
library(ggplot2)

# Assuming `result` contains the calculated Shannon indices and genomic positions
# Preview the result data frame
head(result)

# Plot Shannon Index along the genome
ggplot(resultDSIM12, aes(x = Position, y = Shannon_Index, color = Chromosome)) +
  geom_line() +  # Line plot to connect points
  geom_point(size = 1) +  # Optional: add points for each window
  labs(
    title = "Shannon Diversity Index Along the Genome",
    x = "Genomic Position",
    y = "Shannon Diversity Index"
  ) +
  facet_wrap(~ Chromosome, scales = "free_x") +  # Free x-axis per chromosome, fixed y-axis
  scale_y_continuous(
    limits = c(min(resultDSIM12$Shannon_Index), max(resultDSIM12$Shannon_Index))
  ) +  # Consistent y-axis limits across panels
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


### other diversity index::

# Function to calculate Pi (Nucleotide Diversity)
pi_diversity <- function(freqs) {
  freqs <- freqs / sum(freqs)  # Normalize frequencies
  n <- length(freqs)
  sum(choose(n, 2) * (freqs * (1 - freqs)))  # Pi formula
}

# Function to calculate Tajima's D
# Uses the 'ape' library for Tajima's D
tajimas_d <- function(freqs) {
  # Tajima's D calculation requires a genetic distance matrix or sequence data,
  # so it's typically computed on a distance matrix of sequences.
  # For simplicity, assume that `freqs` represent the genetic distances.
  # Using 'ape' library for Tajima's D calculation (replace with actual data)
  tajima.test(freqs)$D
}
# 
# # Function to calculate Fst (Fixation Index) between two populations
# fst_index <- function(freqs_pop1, freqs_pop2) {
#   # Assuming freqs_pop1 and freqs_pop2 are the haplotype frequencies for two populations
#   p1 <- freqs_pop1 / sum(freqs_pop1)
#   p2 <- freqs_pop2 / sum(freqs_pop2)
#   
#   # Fst formula calculation
#   Hs <- mean(p1 * (1 - p1)) + mean(p2 * (1 - p2))  # Within-population diversity
#   Ht <- (mean(p1) * (1 - mean(p1)) + mean(p2) * (1 - mean(p2)))  # Total population diversity
#   Fst <- (Ht - Hs) / Ht  # Fst formula
#   return(Fst)
# }

#### 
# Load required libraries
library(dplyr)
library(ggplot2)

# Define a function to calculate Shannon Diversity Index
shannon_diversity <- function(freqs) {
  freqs <- freqs[freqs > 0]  # Remove zero frequencies
  -sum(freqs * log(freqs))  # Shannon Index
}

# Define a function to calculate Simpson's Diversity Index
simpsons_diversity <- function(freqs) {
  freqs <- freqs[freqs > 0]  # Remove zero frequencies
  sum(freqs^2)  # Simpson's Index
}

# Define a function to calculate Simpson's Reciprocal Index
simpsons_reciprocal <- function(freqs) {
  freqs <- freqs[freqs > 0]  # Remove zero frequencies
  1 / sum(freqs^2)  # Reciprocal of Simpson's Index
}

# Define a function to calculate Gini-Simpson Index
gini_simpson <- function(freqs) {
  freqs <- freqs[freqs > 0]  # Remove zero frequencies
  1 - sum(freqs^2)  # Gini-Simpson Index (same as Simpson's)
}

# Define a function to calculate Pielou's Evenness Index
pielous_evenness <- function(freqs) {
  freqs <- freqs[freqs > 0]  # Remove zero frequencies
  shannon_index <- -sum(freqs * log(freqs))  # Shannon Index
  evenness <- shannon_index / log(length(freqs))  # Pielou's Evenness
  return(evenness)
}

# Apply all the functions to calculate the indices for each row
haplotype_freqs$Shannon_Index <- apply(haplotype_freqs, 1, function(row) shannon_diversity(row))
haplotype_freqs$Simpsons_Index <- apply(haplotype_freqs, 1, function(row) simpsons_diversity(row))
haplotype_freqs$Simpsons_Reciprocal <- apply(haplotype_freqs, 1, function(row) simpsons_reciprocal(row))
haplotype_freqs$Gini_Simpson <- apply(haplotype_freqs, 1, function(row) gini_simpson(row))
haplotype_freqs$Pielous_Evenness <- apply(haplotype_freqs, 1, function(row) pielous_evenness(row))

head(haplotype_freqs)
haplotype_freqs <- haplotype_freqs %>%
  bind_cols(DSIM_12[, 1:3])  

library(tidyverse)

# Reshape the data to long format
haplotype_freqs_long <- haplotype_freqs %>%
  pivot_longer(
    cols = c(Shannon_Index, Simpsons_Index, Simpsons_Reciprocal, Gini_Simpson, Pielous_Evenness),
    names_to = "Index",
    values_to = "Value"
  )

# Plot
ggplot(haplotype_freqs_long, aes(x = start, y = Value, color = Index, group = Index)) +
  geom_line() +  # Line plot to connect points
  geom_point(size = 1) +  # Optional: add points for each window
  labs(
    title = "Diversity Indices Along the Genome",
    x = "Genomic Position",
    y = "Index Value"
  ) +
  facet_wrap(~ window, scales = "free_x") +  # Free x-axis per chromosome, fixed y-axis
  scale_y_continuous(
    limits = c(min(haplotype_freqs_long$Value), max(haplotype_freqs_long$Value))
  ) +  # Consistent y-axis limits across panels
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### BC_3R and base 100k indexes

resultbase100k$Source <- "Base100k"
resultBC_3R$Source <- "BC_3R"

# Combine the data frames
combined_df <- bind_rows(resultbase100k, resultBC_3R)

# Plot using ggplot2
ggplot(combined_df, aes(x = Position, y = Shannon_Index, color = Source)) +
  geom_line() +  # Line plot to connect points
  geom_point(size = 1) +  # Optional: add points for each window
  labs(
    title = "Shannon Diversity Index Along the Genome",
    x = "Genomic Position",
    y = "Shannon Diversity Index"
  ) +
  facet_wrap(~ Chromosome, scales = "free_x") +  # Free x-axis per chromosome
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### haplotype count distributioN:
sum_values <- apply(BC_3R[, 4:87], 1, sum, na.rm = TRUE)
BC_3R$sum_column <- sum_values
haplotype_count_distribution_function <- function(row_values, threshold) {
  sorted_values <- sort(row_values, decreasing = TRUE)
  cumulative_sum <- 0
  count <- 0
  while (cumulative_sum < threshold && count < length(sorted_values)) {
    cumulative_sum <- cumulative_sum + sorted_values[count + 1]
    count <- count + 1
  }
  return(c(count, cumulative_sum))
}

threshold <- 0.8  # Define your threshold
base100khapcount <- apply(base100k[, 4:87], 1, haplotype_count_distribution_function, threshold)
BC_3Rhapcount <- apply(BC_3R[, 4:87], 1, haplotype_count_distribution_function, threshold)

BC_3Rhapcount_df <- t(data.frame(BC_3Rhapcount))
colnames(BC_3Rhapcount_df) <- c("num_values_added", "cumulative_sum")

head(base100k)

BC_3Rhapcount_df_genomic <- cbind(BC_3R[, 1:3], BC_3Rhapcount_df)

head(base100khapcount_df_genomic)
BC_3Rhapcount <- apply(BC_3R[, 4:87], 1, haplotype_count_distribution_function, threshold)

BC_3Rhapcount_df <- t(data.frame(BC_3Rhapcount))
colnames(BC_3Rhapcount_df) <- c("num_values_added", "cumulative_sum")

head(base100k)

BC_3Rhapcount_df <- cbind(BC_3R[, 1:3], BC_3Rhapcount_df)

head(base100khapcount_df_genomic)
ggplot(BC_3Rhapcount_df_genomic, aes(x = start, y = num_values_added)) +
  geom_line() + 
  geom_point() + 
  labs(title = "hap count distribution",
       x = "genomic pos",
       y = "No of values added") +
  theme_minimal()

ggplot(base100khapcount_df_genomic, aes(x = num_values_added)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Number of Values Added",
       x = "Number of Values Added",
       y = "Frequency") +
  theme_minimal()


# Alternatively, you can use a density plot
# Plot the density of 'num_values_added'
ggplot(base100khapcount_df_genomic, aes(x = num_values_added)) +
  geom_density(fill = "skyblue", alpha = 0.7) +
  labs(title = "Density of Number of Values Added",
       x = "Number of Values Added",
       y = "Density") +
  theme_minimal()

#### 

base100khapcount_df_genomic$source <- "Base100k"
BC_3Rhapcount_df_genomic$source <- "BC_3R"
combined_df <- rbind(base100khapcount_df_genomic, BC_3Rhapcount_df_genomic)
ggplot(combined_df, aes(x = start, y = num_values_added, color = source)) +
  geom_line() + 
  geom_point() +
  labs(title = "Hap count Distribution for Base100k and BC_3R",
       x = "Genomic Position",
       y = "Number of Values Added") +
  scale_color_manual(values = c("Base100k" = "blue", "BC_3R" = "red")) +
  theme_minimal()

##### install poolseq
install.packages(c("foreach", "matrixStats"))
install.packages("/Users/siva/Downloads/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
library(poolSeq)
install.packages("/Users/siva/Downloads/ACER-ACER_1.0.3.tar.gz", repos=NULL, type="source")
library(ACER)
######### now we 
