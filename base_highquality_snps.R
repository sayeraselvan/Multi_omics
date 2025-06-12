library(poolSeq)

Base100k_Base16k_TE_masked_snps_only <- read.sync(file="~/Desktop/Rawdata/Base100k_Base16k_maxdepth8k_sts_minq40.mpileup_indel_TE_masked_snp_only.sync", gen=c(0,0), repl=c(1,2))
Base100k_Base16k_TE_masked_snps <- read.sync(file="~/Desktop/Rawdata/Base100k_Base16k_maxdepth8k_sts_minq40.mpileup_indel_TE_masked_snp.sync", gen=c(0,0), repl=c(1,2))

af_Base100k_Base16k_TE_masked_snps_only <- af(sync = Base100k_Base16k_TE_masked_snps_only, gen=c(0,0), repl=c(1,2))
af_Base100k_Base16k_TE_masked_snps <- af(sync = Base100k_Base16k_TE_masked_snps, gen=c(0,0), repl=c(1,2))

cov_Base100k_Base16k_TE_masked_snps_only <- coverage(sync = Base100k_Base16k_TE_masked_snps_only, gen=c(0,0), repl=c(1,2))
cov_Base100k_Base16k_TE_masked_snps <- coverage(sync = Base100k_Base16k_TE_masked_snps, gen=c(0,0), repl=c(1,2))

file_manage_and_rename <- function(data) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input must be a matrix or data frame.")
  }
  chrom_pos <- rownames(data)
  if (is.null(chrom_pos)) {
    stop("Input does not have row names for chromosomal positions.")
  }
  data_df <- data.frame(chrom_pos = chrom_pos, data)
  data_df <- data_df %>%
    separate(chrom_pos, into = c("CHROM", "POS"), sep = "\\.") %>%
    mutate(POS = as.numeric(POS))
  colnames(data_df) <- c("CHROM", "POS", "Base100k", "Base16k")
  rownames(data_df) <- NULL
  return(data_df)
}

af_Base100k_Base16k_TE_masked_snps_only_df <- file_manage_and_rename(af_Base100k_Base16k_TE_masked_snps_only)
af_Base100k_Base16k_TE_masked_snps_df  <- file_manage_and_rename(af_Base100k_Base16k_TE_masked_snps)

cov_Base100k_Base16k_TE_masked_snps_only_df <- file_manage_and_rename(cov_Base100k_Base16k_TE_masked_snps_only)
cov_Base100k_Base16k_TE_masked_snps_df <- file_manage_and_rename(cov_Base100k_Base16k_TE_masked_snps)


###### FINAL SNP FILTERED, LETS compute the coverage
FINAL_SYNC_BASE_EVOLVED <- read.sync(file="~/Desktop/Rawdata/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync", 
                                                  gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
af_FINAL_SYNC_BASE_EVOLVED <- af(sync = FINAL_SYNC_BASE_EVOLVED, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
cov_FINAL_SYNC_BASE_EVOLVED <- coverage(sync = FINAL_SYNC_BASE_EVOLVED,  gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))

colnames(cov_FINAL_SYNC_BASE_EVOLVED) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")
colnames(af_FINAL_SYNC_BASE_EVOLVED) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")

head(cov_FINAL_SYNC_BASE_EVOLVED)

# coverage_columns <- cov_evol_df_no_zeros[, -c(1, 2)]
cov_FINAL_SYNC_BASE_EVOLVED[cov_FINAL_SYNC_BASE_EVOLVED == 0] <- NA
mean_coverage_per_sample <- colMeans(cov_FINAL_SYNC_BASE_EVOLVED, na.rm = TRUE)
overall_average_coverage <- sum(mean_coverage_per_sample)

print(mean_coverage_per_sample)
print(overall_average_coverage)
0.02*overall_average_coverage

tail_1_percent <- apply(cov_FINAL_SYNC_BASE_EVOLVED, 2, function(column) {
  quantile(column, probs = 0.02, na.rm = TRUE)
})
print(tail_1_percent)

#### after indel and TE masking, we have the snps, 

sync_data <- read.table("~/Desktop/Rawdata/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked.sync", header = FALSE, stringsAsFactors = FALSE)
colnames(sync_data) <- c("CHR", "POS", "REF" ,"Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")
colnames(triallelic) <- c("CHR", "POS", "REF" ,"Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")
colnames(tetraallelic) <- c("CHR", "POS", "REF" ,"Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")

count_non_zero <- function(allele_str) {
  counts <- as.numeric(unlist(strsplit(allele_str, ":"))[1:4])  # A:T:C:G
  sum(counts > 0)
}
biallelic <- sync_data[apply(sync_data[, 4:ncol(sync_data)], 1, function(row) {
  any(sapply(row, count_non_zero) == 2)
}), ]

triallelic <- sync_data[apply(sync_data[, 4:ncol(sync_data)], 1, function(row) {
  any(sapply(row, count_non_zero) == 3)
}), ]

tetraallelic <- sync_data[apply(sync_data[, 4:ncol(sync_data)], 1, function(row) {
  any(sapply(row, count_non_zero) == 4)
}), ]


write.table(biallelic, "~/Desktop/Rawdata/biallelic.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(triallelic, "~/Desktop/Rawdata/fet/triallelic.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(tetraallelic, "~/Desktop/Rawdata/fet/tetraallelic.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

triallelic_chrom_pos <- triallelic %>% select(V1, V2)
tetraallelic_chrom_pos <- tetraallelic %>% select(V1, V2)

excluded_chrom_pos <- bind_rows(triallelic_chrom_pos, tetraallelic_chrom_pos)
filtered_sync_data <- sync_data %>%
  anti_join(excluded_chrom_pos, by = c("V1", "V2"))

write.table(filtered_sync_data, "~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
head(filtered_sync_data)

### coverage distribution:
FINAL_MASKED_SYNC_BASE_EVOLVED <- read.sync(file="~/Desktop/Rawdata/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked.sync", 
                                     gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
af_FINAL_MASKED_SYNC_BASE_EVOLVED <- af(sync = FINAL_MASKED_SYNC_BASE_EVOLVED, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
cov_FINAL_MASKED_SYNC_BASE_EVOLVED <- coverage(sync = FINAL_MASKED_SYNC_BASE_EVOLVED, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))

colnames(cov_FINAL_MASKED_SYNC_BASE_EVOLVED) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")
colnames(af_FINAL_MASKED_SYNC_BASE_EVOLVED) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")

cov_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes <- cov_FINAL_MASKED_SYNC_BASE_EVOLVED
cov_BC_40_no_zeroes <- cov_BC_40
cov_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes[cov_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes == 0] <- NA
cov_BC_40_no_zeroes[cov_BC_40_no_zeroes == 0] <- NA

coverage_2nd_percentiles <- sapply(cov_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, function(x) quantile(x, probs = 0.02, na.rm = TRUE))
coverage_2nd_percentiles


library(data.table)
dt <- as.data.table(cov_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes)
dt <- as.data.table(cov_BC_40_no_zeroes)

coverage_2nd_percentiles <- dt[, lapply(.SD, function(x) quantile(x, probs = 0.01, na.rm = TRUE))]
coverage_2nd_percentiles

na_rows_count <- sum(apply(cov_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, 1, function(row) any(is.na(row))))
na_rows_count

mean_coverage_per_sample <- colMeans(dt, na.rm = TRUE)
overall_average_coverage <- sum(mean_coverage_per_sample)

print(mean_coverage_per_sample)
print(overall_average_coverage)

coverage_values <- coverage_values[!is.na(coverage_values)]

library(ggplot2)

ggplot(data.frame(Coverage = coverage_values), aes(x = Coverage)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Coverage Distribution",
    x = "Coverage",
    y = "Frequency"
  ) +
  theme_minimal()

##### only the biallelic snps for the whole dataset

bi_FINAL_MASKED_SYNC_BASE_EVOLVED <- read.sync(file="~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync", 
                                            gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
af_bi_FINAL_MASKED_SYNC_BASE_EVOLVED <- af(sync = FINAL_MASKED_SYNC_BASE_EVOLVED, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED <- coverage(sync = FINAL_MASKED_SYNC_BASE_EVOLVED, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))

colnames(af_bi_FINAL_MASKED_SYNC_BASE_EVOLVED) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")
colnames(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")

cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes <- cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED
cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes[cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes == 0] <- NA

coverage_2nd_percentiles <- sapply(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, function(x) quantile(x, probs = 0.02, na.rm = TRUE))
coverage_2nd_percentiles


library(data.table)
dt <- as.data.table(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes)

coverage_2nd_percentiles <- dt[, lapply(.SD, function(x) quantile(x, probs = 0.01, na.rm = TRUE))]
coverage_2nd_percentiles

na_rows_count <- sum(apply(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, 1, function(row) any(is.na(row))))
na_rows_count

mean_coverage_per_sample <- colMeans(dt, na.rm = TRUE)
overall_average_coverage <- sum(mean_coverage_per_sample)

print(mean_coverage_per_sample)
print(overall_average_coverage)

coverage_values <- coverage_values[!is.na(coverage_values)]

library(ggplot2)

ggplot(data.frame(Coverage = coverage_values), aes(x = Coverage)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  labs(
    title = "Coverage Distribution",
    x = "Coverage",
    y = "Frequency"
  ) +
  theme_minimal()



##### FET plot::
#### Fisher exact test for BC and SC:  
Base100k_TB20_fisher <- read.table("~/Desktop/Rawdata/fet/FBase100k_TB20_no_na.fet")
Base100k_BC40_fisher <- read.table("~/Desktop/Rawdata/fet/FBase100k_BC40_no_na.fet")

Base16k_TS20_fisher <- read.table("~/Desktop/Rawdata/fet/FBase16k_TS20_no_na.fet")
Base16k_TS40_fisher <- read.table("~/Desktop/Rawdata/fet/FBase16k_TS40_no_na.fet")
Base16k_TS60_fisher <- read.table("~/Desktop/Rawdata/fet/FBase16k_TS60_no_na.fet")
Base16k_TS80_fisher <- read.table("~/Desktop/Rawdata/fet/FBase16k_TS80_no_na.fet")

TB20_BC40_fisher <- read.table("~/Desktop/Rawdata/fet/FTB20_BC40_no_na.fet")
TS20_TS80_fisher <- read.table("~/Desktop/Rawdata/fet/FTS20_TS80_no_na.fet")

Base100k_Base16k_fisher1 <- read.table("~/Desktop/Rawdata/fet/FBase100k_Base16k_mincount1_no_na.fet")
Base100k_Base16k_fisher <- read.table("~/Desktop/Rawdata/fet/FBase100k_Base16k_no_na.fet")

## min count 10
Base100k_TB20_fisher10 <- read.table("~/Desktop/Rawdata/fet/mincount_10_fet/FBase100k_TB20_no_na.fet")
Base100k_BC40_fisher10 <- read.table("~/Desktop/Rawdata/fet/mincount_10_fet/FBase100k_BC40_no_na.fet")

Base16k_TS20_fisher10 <- read.table("~/Desktop/Rawdata/fet/mincount_10_fet/FBase16k_TS20_no_na.fet")
Base16k_TS40_fisher10 <- read.table("~/Desktop/Rawdata/fet/mincount_10_fet/FBase16k_TS40_no_na.fet")
Base16k_TS60_fisher10 <- read.table("~/Desktop/Rawdata/fet/mincount_10_fet/FBase16k_TS60_no_na.fet")
Base16k_TS80_fisher10 <- read.table("~/Desktop/Rawdata/fet/mincount_10_fet/FBase16k_TS80_no_na.fet")


library(ggplot2)

chromosome_of_interest <- "2L"  

start_position <- 1 
end_position <- 20000000
filtered_data <- BCfisher[Base100k_BC40_fisher$V1 == chromosome_of_interest & 
                            BCfisher$V2 >= start_position & 
                            BCfisher$V2 <= end_position & 
                            BCfisher$V3, ]
filtered_data <- Base100k_Base16k_fisher[Base100k_Base16k_fisher$V1 == chromosome_of_interest, ]

#filtered_data <- Base16k_TS80_fisher
filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "BH")

#fdrs<-p.adjust(pvalues, method="BH")

ggplot(filtered_data, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.2, size = 0.1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(filtered_data, aes(x = V2, y = -log10(FDR))) +
      geom_point(alpha = 0.5, size = 0.2, color = "blue") +
      labs(
        title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
        #subtitle = paste("Positions", start_position, "to", end_position),
        x = "Genomic Position",
        y = "-log10(P-value)"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

#### vcf filtered snps and the sync file - gonna perform fisher exact test

Mafvcf <- read.sync(file="~/Desktop/Rawdata/fet/Mafvcf_Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync", 
                                            gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
af_Mafvcf <- af(sync = Mafvcf, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
cov_Mafvcf <- coverage(sync = Mafvcf, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))

colnames(cov_Mafvcf) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")
colnames(af_Mafvcf) <- c("Base100k","Base16k","TB20","BC40","TS20","TS40","TS60","TS80")

cov_Mafvcf_no_zeroes <- cov_Mafvcf

cov_Mafvcf_no_zeroes[cov_Mafvcf_no_zeroes == 0] <- NA
mean_coverage_per_sample <- colMeans(cov_Mafvcf_no_zeroes, na.rm = TRUE)
overall_average_coverage <- sum(mean_coverage_per_sample)

print(mean_coverage_per_sample)
print(overall_average_coverage)
0.02*overall_average_coverage

tail_1_percent <- apply(cov_Mafvcf_no_zeroes, 2, function(column) {
  quantile(column, probs = 0.01, na.rm = TRUE)
})
print(tail_1_percent)


### fisher test results:
chromposvcf <- read.table("~/Desktop/Rawdata/fet/chrom_pos_filteredvcf.txt")
chromposmafvcf <- read.table("~/Desktop/Rawdata/fet/chrom_pos_filteredmafvcf.txt")
filter_by_chrompos <- function(data, chrompos) {
  data$pair <- paste(data$V1, data$V2, sep = "_")
  chrompos$pair <- paste(chrompos$V1, chrompos$V2, sep = "_")
  filtered_data <- data[data$pair %in% chrompos$pair, ]
  filtered_data$pair <- NULL
  chrompos$pair <- NULL
  
  return(filtered_data)
}

datasets <- list(
  Base100k_BC40_fisher = Base100k_BC40_fisher,
  Base100k_TB20_fisher = Base100k_TB20_fisher,
  Base16k_TS20_fisher = Base16k_TS20_fisher,
  Base16k_TS40_fisher = Base16k_TS40_fisher,
  Base16k_TS60_fisher = Base16k_TS60_fisher,
  Base16k_TS80_fisher = Base16k_TS80_fisher
)

chromposvcf_results <- lapply(datasets, filter_by_chrompos, chrompos = chromposvcf)
chromposmafvcf_results <- lapply(datasets, filter_by_chrompos, chrompos = chromposmafvcf)
names(chromposvcf_results) <- paste0("chromposvcf_", names(datasets))
names(chromposmafvcf_results) <- paste0("chromposmafvcf_", names(datasets))
list2env(chromposvcf_results, envir = .GlobalEnv)
list2env(chromposmafvcf_results, envir = .GlobalEnv)

chromosome_of_interest <- "2L"
filtered_data <- chromposmafvcf_Base16k_TS80_fisher[chromposmafvcf_Base16k_TS80_fisher$V1 == chromosome_of_interest, ]

filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "BH")

#fdrs<-p.adjust(pvalues, method="BH")

ggplot(filtered_data, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.2, size = 0.1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(filtered_data, aes(x = V2, y = -log10(FDR))) +
  geom_point(alpha = 0.5, size = 0.2, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

### MAF change
library(ggplot2)
head(af_Mafvcf)
af_Maf_BC <- af_Mafvcf[, c("Base100k", "TB20", "BC40")]
af_Maf_SC <- af_Mafvcf[, c("Base16k", "TS20", "TS40", "TS60", "TS80")]

af_Mafvcf_long <- reshape2::melt(af_Maf_BC_filtered, variable.name = "Group", value.name = "MAF")
af_Mafvcf_long <- af_Mafvcf_long[!is.na(af_Mafvcf_long$MAF), ]
af_Mafvcf_long <- af_Mafvcf_long[!is.na(af_Mafvcf_long$MAF) & af_Mafvcf_long$MAF > 0, ]

ggplot(af_Mafvcf_long, aes(x = MAF, colour = Var2, linetype = "4")) +
  geom_density(size = 1) +
  labs(
    title = "Density Plot of Minor Allele Frequencies",
    x = "Minor Allele Frequency (MAF)",
    y = "Density"
  ) +
  theme_minimal()

head(chromposmafvcf_Base100k_BC40_fisher)
chrompos <- paste(chromposmafvcf_Base100k_BC40_fisher$V1, chromposmafvcf_Base100k_BC40_fisher$V2, sep = ".")

# Extract row names (chromosome and position) from af_Maf_BC
af_Maf_BC_chrompos <- rownames(af_Maf_BC)

# Filter af_Maf_BC where the chrom.pos exists in chromposmafvcf_Base100k_BC40_fisher
af_Maf_BC_filtered <- af_Maf_BC[af_Maf_BC_chrompos %in% chrompos, ]

# View the filtered data
head(af_Maf_BC_filtered)

alleles(sync = Mafvcf)
pol_Mafvcf <- polarization(Mafvcf)


#### Rising alelle:
sync_data_bi <- sync_data %>%
  anti_join(triallelic, by = c("CHR", "POS")) %>%
  anti_join(tetraallelic, by = c("CHR", "POS"))
head(sync_data_bi)

get_major_minor <- function(sync_column) {
  counts <- as.numeric(strsplit(sync_column, ":")[[1]][1:4])  
  alleles <- c("A", "T", "C", "G")                  
  sorted_indices <- order(-counts) 
  major <- alleles[sorted_indices[1]]
  minor <- ifelse(length(sorted_indices) > 1 && counts[sorted_indices[2]] > 0, 
                  alleles[sorted_indices[2]], NA) 
  return(c(major, minor))
}

sync_data_bi <- sync_data_bi %>%
  rowwise() %>%
  mutate(
    Major_Base100k = get_major_minor(Base100k)[1],
    Minor_Base100k = get_major_minor(Base100k)[2],
    Major_Base16k = get_major_minor(Base16k)[1],
    Minor_Base16k = get_major_minor(Base16k)[2]
  )

head(sync_data_bi)
head(test)
set.seed(123) 
n_random_rows <- 1000
random_indices <- sample(1:nrow(sync_data_bi), n_random_rows)
random_rows <- sync_data_bi[random_indices, ]

print(random_rows)
is_polymorphic <- function(base_col) {
  counts <- strsplit(base_col, ":") %>% lapply(as.numeric)
  sapply(counts, function(x) sum(x > 0) >= 2)
}
random_rows_polymorphic <- random_rows %>%
  filter(is_polymorphic(Base100k) & is_polymorphic(Base16k))
random_rows_polymorphic <- random_rows_polymorphic %>%
  rowwise() %>%
  mutate(
    Major_Base100k = get_major_minor(Base100k)[1],
    Minor_Base100k = get_major_minor(Base100k)[2],
    Major_Base16k = get_major_minor(Base16k)[1],
    Minor_Base16k = get_major_minor(Base16k)[2]
  )

random_rows_polymorphic_found <- random_rows_polymorphic %>%
  filter(Major_Base100k == Major_Base16k & Minor_Base100k == Minor_Base16k)

print(filtered_rows)

### major and minor allele frequency count
library(dplyr)
library(tidyr)
library(stringr)

calculate_allele_frequencies <- function(base_column, major_allele_base100k) {
  counts <- str_split(base_column, ":", simplify = TRUE) %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric))
  
  colnames(counts) <- c("A", "T", "C", "G", "N", "del")
  
  # Total counts per row
  total_counts <- rowSums(counts)
  
  # Major allele frequency (from Base100k's major allele)
  major_freq <- ifelse(total_counts == 0, 0, 
                       counts[cbind(seq_len(nrow(counts)), match(major_allele_base100k, colnames(counts)))] / total_counts)
  
  major_freq
}

# Function to calculate minor allele frequencies
calculate_allele_frequencies_minor <- function(base_column, minor_allele_base100k) {
  counts <- str_split(base_column, ":", simplify = TRUE) %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric))
  
  colnames(counts) <- c("A", "T", "C", "G", "N", "del")
  
  # Total counts per row
  total_counts <- rowSums(counts)
  
  # Minor allele frequency (from Base100k's minor allele)
  minor_freq <- ifelse(total_counts == 0, 0, 
                       counts[cbind(seq_len(nrow(counts)), match(minor_allele_base100k, colnames(counts)))] / total_counts)
  
  minor_freq
}

# Columns to calculate allele frequencies for
allele_count_columns <- c("Base100k", "Base16k", "TB20", "BC40", "TS20", "TS40", "TS60", "TS80")

# Process rows to calculate allele frequencies
random_rows_polymorphic_found_AF <- random_rows_polymorphic_found %>%
  rowwise() %>%
  # Calculate major allele frequencies
  mutate(across(
    all_of(allele_count_columns),
    ~ calculate_allele_frequencies(.x, Major_Base100k),
    .names = "{.col}_Major_AF"
  )) %>%
  mutate(across(
    all_of(allele_count_columns),
    ~ calculate_allele_frequencies_minor(.x, Minor_Base100k),
    .names = "{.col}_Minor_AF"
  )) %>%
  ungroup()
head(random_rows_polymorphic_found_AF)

test1 <- random_rows_polymorphic_found_AF %>%
  select(Base16k_Minor_AF, TS80_Minor_AF, TS20_Minor_AF, TS40_Minor_AF, TS60_Minor_AF, Base100k_Minor_AF, BC40_Minor_AF)

test1 <- test1 %>%
  mutate(
    allele_trend = if_else(Base16k_Minor_AF < TS20_Minor_AF, "Minor Rising", "Major Rising")
  )

head(test1)
library(ggplot2)

ggplot(test1, aes(x = Base16k_Minor_AF, y = TS80_Minor_AF)) +
  geom_point(aes(color = allele_trend), size = 1) +
  labs(
    title = "Minor Allele Frequencies: Base16k vs TS80",
    x = "Base16k Minor Allele Frequency",
    y = "TS80 Minor Allele Frequency"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Minor Rising" = "blue", "Major Rising" = "red"))

test1$Ts80minordiff <- test1$TS80_Minor_AF - test1$Base16k_Minor_AF  
test1$Ts60minordiff <- test1$TS60_Minor_AF - test1$Base16k_Minor_AF 
test1$Ts40minordiff <- test1$TS40_Minor_AF - test1$Base16k_Minor_AF 
test1$Ts20minordiff <- test1$TS20_Minor_AF - test1$Base16k_Minor_AF 
test1$TsBC40minordiff <- test1$BC40_Minor_AF - test1$Base100k_Minor_AF 

test1_long_diff <- test1 %>%
  select(TS80_Minor_AF, TS20_Minor_AF) %>%
  mutate(Line = row_number()) %>% # Add a Line identifier
  pivot_longer(
    cols = starts_with("TS"), 
    names_to = "TimePoint",
    values_to = "Minor_Diff"
  )

# Plot the minor allele frequency differences
ggplot(test1_long_diff, aes(x = Minor_Diff, color = TimePoint)) +
  geom_density(alpha = 0.4, size = 1) +
  labs(
    title = "Density of Minor Allele Frequency Differences Over Time Points",
    x = "Minor Allele Frequency Difference",
    y = "Density",
    fill = "Time Point",
    color = "Time Point"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )



hetest1_long <- test1 %>%
  select(Base16k_Minor_AF, TS80_Minor_AF, allele_trend) %>%
  gather(key = "TimePoint", value = "Minor_AF", -allele_trend)


# Plot Minor Allele Frequency (MAF) for Base16k and TS80
ggplot(test1_long, aes(x = TimePoint, y = Minor_AF, color = TimePoint)) +
  geom_point(aes(color = allele_trend), size = 4) +
  geom_line(aes(group = allele_trend), size = 1) +  # Connecting lines between points for trends
  labs(
    title = "Minor Allele Frequency (MAF) Changes Over Time",
    x = "Time Point",
    y = "Minor Allele Frequency (MAF)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Minor Rising" = "blue", "Major Rising" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(test1, aes(x = minordiff, linetype = "4")) +
  geom_density(size = 1) +
  labs(
    title = "Density Plot of Minor Allele Frequencies",
    x = "Minor Allele Frequency (MAF)",
    y = "Density"
  ) +
  theme_minimal()

### adapter chi sq test:
adapter_chiqs_test_syncdata <- sync_data[
  paste(sync_data$CHR, sync_data$POS) %in% paste(chromposmafvcf$V1, chromposmafvcf$V2), 
]

non_matching_rows_adapter_chiqs_test_syncdata <- chromposmafvcf[
  !(paste(chromposmafvcf$V1, chromposmafvcf$V2) %in% paste(sync_data$CHR, sync_data$POS)), 
]

head(non_matching_rows_adapter_chiqs_test_syncdata)
head(adapter_chiqs_test_syncdata)

adapter_chiqs_test_syncdata_Base16k_TS80 <- adapter_chiqs_test_syncdata[c(1,2,3,5,11)]
adapter_chiqs_test_syncdata_Base16k_TS80 <- adapter_chiqs_test_syncdata_Base16k_TS80[
  paste(adapter_chiqs_test_syncdata_Base16k_TS80$CHR, adapter_chiqs_test_syncdata_Base16k_TS80$POS) %in% paste(chromposmafvcf_Base16k_TS80_fisher$V1, chromposmafvcf_Base16k_TS80_fisher$V2), 
]
write.table(adapter_chiqs_test_syncdata_Base16k_TS80, "~/Desktop/Rawdata/fet/adapter_chiqs_test_syncdata_Base16k_TS80.sync", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
adapter_chiqs_test_syncdata_Base16k_TS80_sync <- read.sync(file="~/Desktop/Rawdata/fet/adapter_chiqs_test_syncdata_Base16k_TS80.sync", gen=c(0,80), repl=c(1,1))

af_adapter_chiqs_test_syncdata_Base16k_TS80 <- af(adapter_chiqs_test_syncdata_Base16k_TS80_sync, gen=c(0,80), repl=c(1,1))
cov_adapter_chiqs_test_syncdata_Base16k_TS80 <- coverage(adapter_chiqs_test_syncdata_Base16k_TS80_sync, gen=c(0,80), repl=c(1,1))

head(af_adapter_chiqs_test_syncdata_Base16k_TS80)
adapted.chisq.test(freq = af_adapter_chiqs_test_syncdata_Base16k_TS80, coverage =cov_adapter_chiqs_test_syncdata_Base16k_TS80 , Ne = 47, gen = c(0,80), poolSize = rep(250,2), 
                   mincov = 1, MeanStart = TRUE, IntGen = FALSE, TA = FALSE, RetVal = 0)


### testing: FBase16k_TS80.fet _no_na and the biallelic SNPS after maf from the vcf
FBASE16K_TS80 <- read.table("~/Desktop/Rawdata/fet/biallelic_maf/FBase16k_TS80_no_na.fet")


library(ggplot2)

chromosome_of_interest <- "2L"  

start_position <- 1 
end_position <- 20000000

filtered_data <- FBASE16K_TS80[FBASE16K_TS80$V1 == chromosome_of_interest, ]
filtered_data_sig <- adapter_chiqs_test_syncdata[
  paste(adapter_chiqs_test_syncdata$CHR, adapter_chiqs_test_syncdata$POS) %in% paste(FBASE16K_TS80$V1, FBASE16K_TS80$V2), 
]

#filtered_data <- Base16k_TS80_fisher
filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "BH")

#fdrs<-p.adjust(pvalues, method="BH")

ggplot(filtered_data, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.2, size = 0.1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### polymorphic in all populations:
sync_file <- "~/Desktop/Rawdata/fet/Mafvcf_Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync"
Mafvcf_sync <- read.table(sync_file, header = FALSE, stringsAsFactors = FALSE)

is_biallelic <- function(allele_counts) {
  counts <- as.numeric(unlist(strsplit(allele_counts, ":")))
  non_zero_counts <- counts[counts > 0]
  length(non_zero_counts) <= 2
}

biallelic_sites <- Mafvcf_sync[
  apply(Mafvcf_sync[, 4:ncol(Mafvcf_sync)], 1, function(row) {
    all(sapply(row, is_biallelic))
  }),
]
write.table(biallelic_sites, "~/Desktop/Rawdata/fet/Mafvcf_Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync_onlybialleic_in_all_populations.sync", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

## now we are only looking at the biallelic sites of chromposmaf base16k vs TS80k, we take random 1000 snps and then we calculate the maf and plot the maf trajectory
head(chromposmafvcf_Base16k_TS80_fisher)
dim(chromposmafvcf_Base16k_TS80_fisher)

head(sync_data_bi)

main_data_based_on_base16k_TS80k <- sync_data_bi[
  paste(sync_data_bi$CHR, sync_data_bi$POS) %in% paste(chromposmafvcf_Base16k_TS80_fisher$V1, chromposmafvcf_Base16k_TS80_fisher$V2), 
]

main_data_based_on_base16k_TS80k_2L <- main_data_based_on_base16k_TS80k[main_data_based_on_base16k_TS80k$CHR == "2L", ]
set.seed(123) 
n_random_rows <- 1000
random_indices <- sample(1:nrow(main_data_based_on_base16k_TS80k_2L), n_random_rows)
random_rows <- main_data_based_on_base16k_TS80k_2L[random_indices, ]

print(random_rows)
is_polymorphic <- function(base_col) {
  counts <- strsplit(base_col, ":") %>% lapply(as.numeric)
  sapply(counts, function(x) sum(x > 0) >= 2)
}
random_rows_polymorphic <- random_rows %>%
  filter(is_polymorphic(Base100k) & is_polymorphic(Base16k))
random_rows_polymorphic <- random_rows_polymorphic %>%
  rowwise() %>%
  mutate(
    Major_Base100k = get_major_minor(Base100k)[1],
    Minor_Base100k = get_major_minor(Base100k)[2],
    Major_Base16k = get_major_minor(Base16k)[1],
    Minor_Base16k = get_major_minor(Base16k)[2]
  )

random_rows_polymorphic_found <- random_rows_polymorphic %>%
  filter(Major_Base100k == Major_Base16k & Minor_Base100k == Minor_Base16k)

library(dplyr)
library(tidyr)
library(stringr)

calculate_allele_frequencies <- function(base_column, major_allele_base100k) {
  counts <- str_split(base_column, ":", simplify = TRUE) %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric))
  
  colnames(counts) <- c("A", "T", "C", "G", "N", "del")
  
  # Total counts per row
  total_counts <- rowSums(counts)
  
  # Major allele frequency (from Base100k's major allele)
  major_freq <- ifelse(total_counts == 0, 0, 
                       counts[cbind(seq_len(nrow(counts)), match(major_allele_base100k, colnames(counts)))] / total_counts)
  
  major_freq
}

# Function to calculate minor allele frequencies
calculate_allele_frequencies_minor <- function(base_column, minor_allele_base100k) {
  counts <- str_split(base_column, ":", simplify = TRUE) %>%
    as.data.frame() %>%
    mutate(across(everything(), as.numeric))
  
  colnames(counts) <- c("A", "T", "C", "G", "N", "del")
  
  # Total counts per row
  total_counts <- rowSums(counts)
  
  # Minor allele frequency (from Base100k's minor allele)
  minor_freq <- ifelse(total_counts == 0, 0, 
                       counts[cbind(seq_len(nrow(counts)), match(minor_allele_base100k, colnames(counts)))] / total_counts)
  
  minor_freq
}

# Columns to calculate allele frequencies for
allele_count_columns <- c("Base100k", "Base16k", "TB20", "BC40", "TS20", "TS40", "TS60", "TS80")

# Process rows to calculate allele frequencies
random_rows_polymorphic_found_AF <- random_rows_polymorphic_found %>%
  rowwise() %>%
  # Calculate major allele frequencies
  mutate(across(
    all_of(allele_count_columns),
    ~ calculate_allele_frequencies(.x, Major_Base100k),
    .names = "{.col}_Major_AF"
  )) %>%
  mutate(across(
    all_of(allele_count_columns),
    ~ calculate_allele_frequencies_minor(.x, Minor_Base100k),
    .names = "{.col}_Minor_AF"
  )) %>%
  ungroup()
head(random_rows_polymorphic_found_AF)

test1 <- random_rows_polymorphic_found_AF %>%
  select(CHR, POS, REF, Base100k_Minor_AF, TB20_Minor_AF, BC40_Minor_AF, Base16k_Minor_AF, TS20_Minor_AF, TS40_Minor_AF, TS60_Minor_AF, TS80_Minor_AF)

test1$Variant_ID <- paste("Variant", seq_len(nrow(test1)))
long_data <- test1 %>%
  pivot_longer(
    cols = c(Base100k_Minor_AF, TB20_Minor_AF, BC40_Minor_AF, Base16k_Minor_AF, TS20_Minor_AF, TS40_Minor_AF, TS60_Minor_AF, TS80_Minor_AF),
    names_to = "Generation",   # Create a new column for the generation
    values_to = "Minor_AF"     # Allele frequency values
  ) %>%
  mutate(
    Generation = factor(Generation, 
                        levels = c("Base100k_Minor_AF", "TB20_Minor_AF", "BC40_Minor_AF", "Base16k_Minor_AF", "TS20_Minor_AF", 
                                   "TS40_Minor_AF", "TS60_Minor_AF", "TS80_Minor_AF"),
                        labels = c("Base100k", "TB20", "BC40", "Base16k", "TS20", "TS40", "TS60", "TS80"))
  )

long_data_BC <- long_data %>%
  filter(Generation %in% c("Base100k", "TB20", "BC40"))
long_data_SC <- long_data %>%
  filter(Generation %in% c("Base16k", "TS20", "TS40", "TS60", "TS80"))

ggplot(long_data_SC, aes(x = Generation, y = Minor_AF, group = Variant_ID, color = Variant_ID)) +
  geom_line(alpha = 0.7) +
  geom_point() +
  labs(
    title = "Allele Frequency Trajectory Across Generations",
    x = "Generation",
    y = "Minor Allele Frequency",
    color = "Variant ID"
  ) +
  # scale_x_discrete(labels = c("Base16k", "TS20", "TS40", "TS60", "TS80")) +  # Customizing x-axis labels
  theme_minimal() +
  theme(legend.position = "none")

fixation_threshold <- 1
loss_threshold <- 0

# Classify alleles based on thresholds
test1 <- test1 %>%
  mutate(
    Status = case_when(
      TS80_Minor_AF >= fixation_threshold ~ "Fixing",
      TS80_Minor_AF <= loss_threshold ~ "Lost",
      TRUE ~ "Polymorphic"  # Neither fixing nor lost
    )
  )

allele_proportions <- test1 %>%
  group_by(Status) %>%
  summarise(Count = n()) %>%
  mutate(Proportion = Count / sum(Count))

print(allele_proportions)
dim(test1)

## manhattan plot, marking the fixation allele and the lost allele;
test1_manhattan <- chromposmafvcf_Base16k_TS80_fisher[
  paste(chromposmafvcf_Base16k_TS80_fisher$V1, chromposmafvcf_Base16k_TS80_fisher$V2) %in% paste(test1$CHR, test1$POS), 
]

test1_manhattan$FisherPValue <- as.numeric(sub(".*=", "", test1_manhattan$V6))
test1_manhattan$PValue <- 10^(-as.numeric(sub(".*=", "", test1_manhattan$V6)))
test1_manhattan$FDR <- p.adjust(test1_manhattan$PValue, method = "BH")

ggplot(test1_manhattan, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.2, size = 1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### now we plot the fixation allele, lost alleles and segregating allele
fixation_threshold <- 1
loss_threshold <- 0

test1 <- test1 %>%
  mutate(
    Status = case_when(
      TS80_Minor_AF >= fixation_threshold ~ "Fixing",
      TS80_Minor_AF <= loss_threshold ~ "Lost",
      TRUE ~ "Polymorphic"  # Neither fixing nor lost
    )
  )

allele_summary <- test1 %>%
  group_by(Status) %>%
  summarise(
    Count = n(),
    Proportion = round(Count / nrow(test1) * 100, 2)  # Convert to percentage
  )

allele_summary_text <- paste(
  allele_summary$Status, 
  paste0("(", allele_summary$Count, ", ", allele_summary$Proportion, "%)"), 
  collapse = "; "
)

test1_manhattan <- test1_manhattan %>%
  mutate(
    CHR = V1,
    POS = V2
  ) %>%
  left_join(test1 %>% select(CHR, POS, Status, TS80_Minor_AF), by = c("CHR", "POS"))

head(test1_manhattan)
color_mapping <- c("Fixing" = "red", "Lost" = "blue", "Polymorphic" = "green")
color_mapping <- c("Polymorphic" = "blue")

ggplot(test1_manhattan, aes(x = POS, y = FisherPValue, color = Status)) +
  geom_point(alpha = 0.7, size = 0.1) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "Manhattan Plot Highlighting Minor and Major Rising Alleles",
    subtitle = paste("Allele Status Summary: ", allele_summary_text),
    x = "Genomic Position",
    y = "-log10(P-value)",
    color = "Allele Status"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### polymorphic only
# Filter the dataset to include only polymorphic sites
polymorphic_data <- test1_manhattan %>% filter(Status == "Polymorphic")

ggplot(polymorphic_data, aes(x = POS, y = FisherPValue, color = Status)) +
  geom_point(alpha = 0.7, size = 0.5) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "Manhattan Plot Highlighting Minor and Major Rising Alleles",
    #subtitle = paste("Allele Status Summary: ", allele_summary_text),
    x = "Genomic Position",
    y = "-log10(P-value)",
    color = "Allele Status"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#### now it is minor rising and major rising:

test1 <- test1 %>%
  mutate(
    Status = case_when(
      TS80_Minor_AF > Base16k_Minor_AF ~ "Minor Rising",  # Minor allele frequency is increasing
      TS80_Minor_AF < Base16k_Minor_AF ~ "Major Rising",  # Minor allele frequency is decreasing
      TRUE ~ "No Change"  # Handle edge cases where they are equal
    )
  )

allele_summary <- test1 %>%
  group_by(Status) %>%
  summarise(
    Count = n(),
    Proportion = round(Count / nrow(test1) * 100, 2)  # Convert to percentage
  )

allele_summary_text <- paste(
  allele_summary$Status, 
  paste0("(", allele_summary$Count, ", ", allele_summary$Proportion, "%)"), 
  collapse = "; "
)

test1_manhattan <- test1_manhattan %>%
  mutate(
    CHR = V1,
    POS = V2
  ) %>%
  left_join(test1 %>% select(CHR, POS, Status, TS80_Minor_AF, Base16k_Minor_AF), by = c("CHR", "POS"))
color_mapping <- c("Minor Rising" = "orange", "Major Rising" = "purple", "No Change" = "gray")

ggplot(test1_manhattan, aes(x = POS, y = FisherPValue, color = Status)) +
  geom_point(alpha = 0.7, size = 0.1) +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "Manhattan Plot Highlighting Minor and Major Rising Alleles",
    subtitle = paste("Allele Status Summary: ", allele_summary_text),
    x = "Genomic Position",
    y = "-log10(P-value)",
    color = "Allele Status"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(allele_summary)

##########  lets do the whole 2L chromsome:
FBASE16K_TS80 <- read.table("~/Desktop/Rawdata/fet/biallelic_maf/FBase16k_TS80_no_na.fet")
FBASE16K_TS80_maindata <- sync_data_bi[
  paste(sync_data_bi$CHR, sync_data_bi$POS) %in% paste(FBASE16K_TS80$V1, FBASE16K_TS80$V2), 
]

FBASE16K_TS80_maindata_2L <- FBASE16K_TS80_maindata[FBASE16K_TS80_maindata$CHR == "2L", ]

write.table(FBASE16K_TS80_maindata_2L, "~/Desktop/Rawdata/fet/biallelic_maf/FBASE16K_TS80_maindata_2L.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

FBASE16K_TS80_maindata_2L_sync <- read.sync(file="~/Desktop/Rawdata/fet/biallelic_maf/FBASE16K_TS80_maindata_2L.sync", 
                                            gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))

af_FBASE16K_TS80_maindata_2L_sync <- af(sync = FBASE16K_TS80_maindata_2L_sync, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
cov_FBASE16K_TS80_maindata_2L_sync <- coverage(sync = FBASE16K_TS80_maindata_2L_sync, gen=c(0,10,20,30,40,50,60,70), repl=c(1,1,1,1,1,1,1,1))
cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80 <- cov_FBASE16K_TS80_maindata_2L_sync[,c(2,8)]

set.seed(123) 
n_random_rows <- 1000
random_indices <- sample(1:nrow(main_data_based_on_base16k_TS80k_2L), n_random_rows)
random_rows <- main_data_based_on_base16k_TS80k_2L[random_indices, ]

print(random_rows)
is_polymorphic <- function(base_col) {
  counts <- strsplit(base_col, ":") %>% lapply(as.numeric)
  sapply(counts, function(x) sum(x > 0) >= 2)
}
random_rows_polymorphic <- FBASE16K_TS80_maindata_2L %>%
  filter(is_polymorphic(Base100k) & is_polymorphic(Base16k))
random_rows_polymorphic <- random_rows_polymorphic %>%
  rowwise() %>%
  mutate(
    Major_Base100k = get_major_minor(Base100k)[1],
    Minor_Base100k = get_major_minor(Base100k)[2],
    Major_Base16k = get_major_minor(Base16k)[1],
    Minor_Base16k = get_major_minor(Base16k)[2]
  )

random_rows_polymorphic_found <- random_rows_polymorphic %>%
  filter(Major_Base100k == Major_Base16k & Minor_Base100k == Minor_Base16k)


allele_count_columns <- c("Base100k", "Base16k", "TB20", "BC40", "TS20", "TS40", "TS60", "TS80")

random_rows_polymorphic_found_AF <- random_rows_polymorphic_found %>%
  rowwise() %>%
  mutate(across(
    all_of(allele_count_columns),
    ~ calculate_allele_frequencies(.x, Major_Base100k),
    .names = "{.col}_Major_AF"
  )) %>%
  mutate(across(
    all_of(allele_count_columns),
    ~ calculate_allele_frequencies_minor(.x, Minor_Base100k),
    .names = "{.col}_Minor_AF"
  )) %>%
  ungroup()
head(random_rows_polymorphic_found_AF)

#library(readr)
#write_tsv(random_rows_polymorphic_found_AF, "~/Desktop/Rawdata/fet/biallelic_maf/2L_base16k_ts80_random_rows_polymorphic_found_AF.tsv")


#### we found discrepancy in the dataset, that is there is so many polymorphic sites emerging 
# Define the threshold
threshold <- 20

# Identify rows where any column has values <= threshold
rows_with_values_le_20 <- cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80[
  apply(cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80, 1, function(row) any(row <= threshold, na.rm = TRUE)), 
  , drop = FALSE
]

# Display the rows
cat("Rows with values <= 20:\n")
print(rows_with_values_le_20)

rows_to_keep <- !(rownames(cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80) %in% rownames(rows_with_values_le_20))
cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil <- cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80[rows_to_keep, ]

cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil <- as.data.frame(cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil)

cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil$Ratio <- 
  cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil$F10.R1.cov / 
  cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil$F70.R1.cov

ggplot(cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil, aes(x = Ratio)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(
    title = "Density Plot of Ratio",
    x = "Ratio (F10.R1.cov / F70.R1.cov)",
    y = "Density"
  ) +
  theme_minimal()
cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_0.8_TO_1.2 <- cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil[
  cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil$Ratio >= 0.8 &
    cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_fil$Ratio <= 1.2, ]

## get only these rows

FBASE16K_TS80$combined <- paste(FBASE16K_TS80$V1, FBASE16K_TS80$V2, sep = ".")

row_names_to_extract <- rownames(cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80_0.8_TO_1.2)
filtered_FBASE16K_TS80 <- FBASE16K_TS80[FBASE16K_TS80$combined %in% row_names_to_extract, ]
head(filtered_FBASE16K_TS80)

filtered_FBASE16K_TS80$FisherPValue <- as.numeric(sub(".*=", "", filtered_FBASE16K_TS80$V6))
filtered_FBASE16K_TS80$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_FBASE16K_TS80$V6)))
filtered_FBASE16K_TS80$FDR <- p.adjust(filtered_FBASE16K_TS80$PValue, method = "BH")

#fdrs<-p.adjust(pvalues, method="BH")

ggplot(filtered_FBASE16K_TS80, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.2, size = 0.2, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


####only bialleic no no new risation of new alleless ### this is performed in the bash command line and then it is read

head(random_rows_polymorphic_found_AF)

FBASE16K_TS80_only2 <- read.table("~/Desktop/Rawdata/fet/biallelic_maf/FBase16k_TS80_no_na_filter.fet")
FBASE16K_TS80_only2_maindata <- sync_data_bi[
  paste(sync_data_bi$CHR, sync_data_bi$POS) %in% paste(FBASE16K_TS80_only2$V1, FBASE16K_TS80_only2$V2), 
]

FBASE16K_TS80_only2_maindata_2L <- FBASE16K_TS80_only2_maindata[FBASE16K_TS80_only2_maindata$CHR == "3R", ]

write.table(FBASE16K_TS80_only2_maindata_2L, "~/Desktop/Rawdata/fet/biallelic_maf/FBASE16K_TS80_only2_maindata_3R.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

FBASE16K_TS80_only2_maindata_2L_sync <- read.sync(file="~/Desktop/Rawdata/fet/biallelic_maf/FBASE16K_TS80_only2_maindata_3R.sync", 
                                                  gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))

af_FBASE16K_TS80_maindata_2L_sync <- af(sync = FBASE16K_TS80_only2_maindata_2L_sync, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_FBASE16K_TS80_maindata_2L_sync <- coverage(sync = FBASE16K_TS80_only2_maindata_2L_sync, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))

af_FBASE16K_TS80_maindata_2L_sync_base16k_ts80 <- af_FBASE16K_TS80_maindata_2L_sync[,c(4,8)]
cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80 <- cov_FBASE16K_TS80_maindata_2L_sync[,c(4,8)]

colnames(af_FBASE16K_TS80_maindata_2L_sync) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
colnames(cov_FBASE16K_TS80_maindata_2L_sync) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")


filtered_data <- FBASE16K_TS80_only2[FBASE16K_TS80_only2$V1 == "3R", ]
filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "BH")

#fdrs<-p.adjust(pvalues, method="BH")

ggplot(filtered_data, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.5, size = 0.1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
### adapter fisher test
FBASE16K_TS80_maindata_2L_sync_base16k_ts80_adapterchisq <- adapted.chisq.test(freq = af_FBASE16K_TS80_maindata_2L_sync_base16k_ts80, coverage =cov_FBASE16K_TS80_maindata_2L_sync_base16k_ts80 , Ne = 47, gen = c(0,80), poolSize = rep(250,2), 
                   mincov = 1, MeanStart = TRUE, IntGen = FALSE, TA = FALSE, RetVal = 0)
df <- as.data.frame(FBASE16K_TS80_maindata_2L_sync_base16k_ts80_adapterchisq)
df$chr_pos <- rownames(df)
df <- tidyr::separate(df, chr_pos, into = c("chr", "pos"), sep = "\\.")
df$pos <- as.numeric(df$pos)
df$logp <- -log10(df[, 1])  # Assuming p-values are in the first column

library(ggplot2)

# Basic Manhattan plot
ggplot(df, aes(x = pos, y = logp)) +
  geom_point(color = "blue", size = 0.1, alpha = 0.5) +  # Scatter plot
  labs(
    title = "Manhattan Plot",
    x = "Genomic Position",
    y = expression(-log[10](p))
  ) +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")  # Add significance threshold

##### rising allele:
FBASE16K_TS80_only2_maindata_2L_sync_rising <- read.sync(file="~/Desktop/Rawdata/fet/biallelic_maf/FBASE16K_TS80_only2_maindata_2L.sync", 
                                                  gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2), polarization = "rising")
af_FBASE16K_TS80_maindata_2L_sync_rising <- af(sync = FBASE16K_TS80_only2_maindata_2L_sync_rising, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_FBASE16K_TS80_maindata_2L_sync_rising <- coverage(sync = FBASE16K_TS80_only2_maindata_2L_sync_rising, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))


#### now we compare only the alleles from the fisher exact test with the ones from the same polymorphic in both bases

filtered_no_new_alleles.snps.sync_3.2 <- read.sync(file="~/Desktop/Rawdata/filtered_no_new_alleles.snps.sync", 
                                                   gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
af_filtered_no_new_alleles.snps.sync_3.2 <- af(sync = filtered_no_new_alleles.snps.sync_3.2, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_filtered_no_new_alleles.snps.sync_3.2 <- coverage(sync = filtered_no_new_alleles.snps.sync_3.2, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
colnames(af_filtered_no_new_alleles.snps.sync_3.2) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
colnames(cov_filtered_no_new_alleles.snps.sync_3.2) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")

dt <- as.data.table(cov_filtered_no_new_alleles.snps.sync_3.2, keep.rownames = TRUE)
Base16k_TS80_fisher

chrompos <- paste(Base16k_TS80_fisher$V1, Base16k_TS80_fisher$V2, sep = ".")  

Base16k_TS80_fisher_filtered <- Base16k_TS80_fisher[chrompos %in% dt$rn, ]

filtered_data <- Base16k_TS80_fisher_filtered[Base16k_TS80_fisher_filtered$V1 == "2L", ]
filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "BH")

#fdrs<-p.adjust(pvalues, method="BH")

ggplot(filtered_data, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.1, size = 0.1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    #subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### adapted fisher exact test:
filtered_no_new_alleles.snps.sync_3.2_table <- read.table(file="~/Desktop/Rawdata/filtered_no_new_alleles.snps.sync")
filtered_no_new_alleles.snps.sync_3.2_table_2L <- filtered_no_new_alleles.snps.sync_3.2_table[filtered_no_new_alleles.snps.sync_3.2_table$V1 == "2L", ]

write.table(filtered_no_new_alleles.snps.sync_3.2_table_2L, "~/Desktop/Rawdata/fet/biallelic_maf/filtered_no_new_alleles.snps.sync_3.2_table_2L.sync", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

filtered_no_new_alleles.snps.sync_3.2_table_2L_sync <- read.sync(file="~/Desktop/Rawdata/fet/biallelic_maf/filtered_no_new_alleles.snps.sync_3.2_table_2L.sync", 
                                                  gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2), polarization = "rising")

af_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync <- af(sync = filtered_no_new_alleles.snps.sync_3.2_table_2L_sync, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync <- coverage(sync = filtered_no_new_alleles.snps.sync_3.2_table_2L_sync, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))

af_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync_base16k_ts80 <- af_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync[,c(4,8)]
cov_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync_base16k_ts80 <- cov_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync[,c(4,8)]

colnames(af_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
colnames(cov_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")


filtered_no_new_alleles.snps.sync_3.2_2L_adapterchisqtest <- adapted.chisq.test(freq = af_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync_base16k_ts80, coverage = cov_filtered_no_new_alleles.snps.sync_3.2_table_2L_sync_base16k_ts80 , Ne = 47, gen = c(0,80), poolSize = rep(250,2), 
                                                                               mincov = 1, MeanStart = TRUE, IntGen = FALSE, TA = FALSE, RetVal = 0)

df <- as.data.frame(filtered_no_new_alleles.snps.sync_3.2_2L_adapterchisqtest)
df$chr_pos <- rownames(df)
df <- tidyr::separate(df, chr_pos, into = c("chr", "pos"), sep = "\\.")
df$pos <- as.numeric(df$pos)
df$logp <- -log10(df[, 1])  # Assuming p-values are in the first column

library(ggplot2)
ggplot(df, aes(x = pos, y = logp)) +
  geom_point(color = "blue", size = 0.1, alpha = 0.2) +  # Scatter plot
  labs(
    title = "Manhattan Plot",
    x = "Genomic Position",
    y = expression(-log[10](p))
  ) +
  theme_minimal()
  #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")  # Add significance threshold

#### manhattan plot with all chromosomes:

Base100k_TB20_fisher_mincov57_mincount1 <- read.table("~/Desktop/Rawdata/FBase100k_TB20_mincov17_mincount1_no_na.fet")
filtered_data <- Base100k_TB20_fisher_mincov57_mincount1
filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "BH")

colnames(filtered_data)[1:2] <- c("Chromosome", "Position")

chromosome_levels <- c("2L", "2R", "3L", "3R", "X")
filtered_data <- filtered_data %>%
  filter(Chromosome %in% chromosome_levels)

# Assign factor levels to chromosomes for proper order
filtered_data$Chromosome <- factor(filtered_data$Chromosome, levels = chromosome_levels)

# Add a cumulative position column for clustering chromosomes
chromosome_offsets <- filtered_data %>%
  group_by(Chromosome) %>%
  summarize(max_pos = max(Position, na.rm = TRUE)) %>%
  mutate(offset = cumsum(lag(max_pos, default = 0)))

filtered_data <- filtered_data %>%
  left_join(chromosome_offsets, by = "Chromosome") %>%
  mutate(CumulativePosition = Position + offset)

# Create Manhattan plot
ggplot(filtered_data, aes(x = CumulativePosition, y = FisherPValue, color = Chromosome)) +
  geom_point(alpha = 0.05, size = 0.2) +
  labs(title = "Manhattan Plot", x = "Chromosomes", y = "-log10(pvalue)") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.spacing = unit(0.5, "lines"),   # Adjust spacing between panels
    legend.position = "none"             # Remove legend
  ) +
  scale_color_manual(values = c("red", "blue", "red", "blue", "red")) + # Set colors
  scale_x_continuous(breaks = chromosome_offsets$offset + chromosome_offsets$max_pos / 2,
                     labels = chromosome_levels) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 20))
  


### fisher exact test fro bcf file format: this is done with reference to ref and alt allele:
# af_txt <- read.table("/Users/siva/Desktop/Rawdata/af/base100k_bc40_mpileup.alf.txt")
# cov_text <- read.table("/Users/siva/Desktop/Rawdata/af/base100k_bc40_mpileup.cov.txt")
AF_alt <- read.table("/Users/siva/Desktop/Rawdata/af/base100k_bc40_mpileup.AF_alf.txt")
# 
# cov_txt <- cov_text %>% filter(V5 >= quantile(V5, probs = 0.01) & V5 <= quantile(V5, probs = 0.99))
# cov_txt <- cov_txt %>% mutate(combined = paste(V1, V2, sep = "_"))
# 
# af_txt <- af_txt %>% mutate(combined = paste(V1, V2, sep = "_"))
# af_text <- af_txt %>% filter(combined %in% cov_txt$combined)
# 
# rm(af_txt)
# rm(cov_text)
# rm(AF_alt)
# # Find the rows where either V5 or V6 is 0 in cov_txt
# #rows_to_remove <- cov_txt$V5 == 0 | cov_txt$V6 == 0
# 
# # Remove these rows from both af_text and cov_txt
# #af_text <- af_text[!rows_to_remove, ]
# #cov_txt <- cov_txt[!rows_to_remove, ]
# 
bigcage <- adapted.chisq.test(freq = af_matrix_filtered, coverage = cov_matrix_filtered, Ne = 2000 , gen = c(0,40), poolSize = rep(250,2), 
                    mincov = 1, MeanStart = TRUE, IntGen = FALSE, TA = FALSE, RetVal = 0)
 
# # # multiallelic states, keep only one
# # af_text$combined <- paste(af_text$V1, af_text$V2, sep = "_")
# # af_text_no_duplicates <- af_text[!duplicated(af_text$combined), ]
# # cov_txt$combined <- paste(cov_txt$V1, cov_txt$V2, sep = "_")
# # cov_txt_no_duplicates <- cov_txt[!duplicated(af_text$combined), ]
# # af_text_no_duplicates$combined <- NULL
# # cov_txt_no_duplicates$combined <- NULL
# # head(af_text_no_duplicates)
# # head(cov_txt_no_duplicates)
# 
# # remove multi allelic sites
# af_text$combined <- paste(af_text$V1, af_text$V2, sep = "_")
# duplicates_af <- af_text[duplicated(af_text$combined), ]
# af_text_no_duplicates <- af_text[!af_text$combined %in% duplicates_af$combined, ]
# cov_txt$combined <- paste(cov_txt$V1, cov_txt$V2, sep = "_")
# cov_txt_no_duplicates <- cov_txt[!cov_txt$combined %in% duplicates_af$combined, ]
# 
# af_text_no_duplicates$combined <- NULL
# cov_txt_no_duplicates$combined <- NULL
# head(af_text_no_duplicates)
# head(cov_txt_no_duplicates)
# 
# af_text_no_duplicates$Chr.pos <- paste(af_text_no_duplicates$V1, af_text_no_duplicates$V2, sep = ".")
# af_matrix <- af_text_no_duplicates[, c("Chr.pos", "V5", "V6")]
# rownames(af_matrix) <- af_text_no_duplicates$Chr.pos  
# af_matrix <- as.matrix(af_matrix[, -1]) 
# head(af_matrix)
# 
# cov_txt_no_duplicates$Chr.pos <- paste(cov_txt_no_duplicates$V1, cov_txt_no_duplicates$V2, sep = ".")
# cov_matrix <- cov_txt_no_duplicates[, c("Chr.pos", "V5", "V6")]
# cov_matrix <- as.matrix(cov_matrix[, -1])  
# rownames(cov_matrix) <- cov_txt_no_duplicates$Chr.pos 
# head(cov_matrix)
af_txt <- read.table("/Users/siva/Desktop/Rawdata/af/base100k_bc40_mpileup.alf.txt", header = FALSE)
cov_txt <- read.table("/Users/siva/Desktop/Rawdata/af/base100k_bc40_mpileup.cov.txt", header = FALSE)

colnames(af_txt) <- c("V1", "V2", "V3", "V4", "V5", "V6")
colnames(cov_txt) <- c("V1", "V2", "V3", "V4", "V5", "V6")

af_txt$combined <- paste(af_txt$V1, af_txt$V2, sep = "_")
cov_txt$combined <- paste(cov_txt$V1, cov_txt$V2, sep = "_")
duplicates_af <- af_txt[duplicated(af_txt$combined), ]
af_txt_no_duplicates <- af_txt[!af_txt$combined %in% duplicates_af$combined, ]
cov_txt_no_duplicates <- cov_txt[!cov_txt$combined %in% duplicates_af$combined, ]

af_txt_no_duplicates$combined <- NULL
cov_txt_no_duplicates$combined <- NULL

cov_txt_filtered <- cov_txt_no_duplicates %>% filter(V5 >= quantile(V5, probs = 0.05) & V5 <= quantile(V5, probs = 0.95))
cov_txt_filtered <- cov_txt_filtered %>% filter(V6 >= quantile(V6, probs = 0.05) & V6 <= quantile(V6, probs = 0.95))

cov_txt_filtered <- cov_txt_filtered %>%
  mutate(combined = paste(V1, V2, sep = "_"))

af_txt_no_duplicates <- af_txt_no_duplicates %>%
  mutate(combined = paste(V1, V2, sep = "_"))
af_txt_filtered <- af_txt_no_duplicates %>%
  filter(combined %in% cov_txt_filtered$combined)
af_txt_filtered <- af_txt_filtered %>%
  select(-combined)

cov_txt_filtered <- cov_txt_filtered %>%
  select(-combined)
head(af_txt_filtered)
head(cov_txt_filtered)

af_txt_filtered$Chr.pos <- paste(af_txt_filtered$V1, af_txt_filtered$V2, sep = ".")
af_matrix <- af_txt_filtered[, c("Chr.pos", "V5", "V6")]
rownames(af_matrix) <- af_txt_filtered$Chr.pos
af_matrix <- as.matrix(af_matrix[, -1])
head(af_matrix)

cov_txt_filtered$Chr.pos <- paste(cov_txt_filtered$V1, cov_txt_filtered$V2, sep = ".")
cov_matrix <- cov_txt_filtered[, c("Chr.pos", "V5", "V6")]
cov_matrix <- as.matrix(cov_matrix[, -1])
rownames(cov_matrix) <- cov_txt_filtered$Chr.pos
head(cov_matrix)

### now it is the removal of alleles: based on the 0.05 and 0.95
af_matrix_df <- as.data.frame(af_matrix)
af_matrix_filtered <- af_matrix_df %>%
  filter(!(V5 < 0.05 & V6 < 0.05) & !(V5 > 0.95 & V6 > 0.95))

filtered_row_names <- rownames(af_matrix_filtered)
cov_matrix_filtered <- cov_matrix[rownames(cov_matrix) %in% filtered_row_names, ]
head(af_matrix_filtered)
head(cov_matrix_filtered)


df <- as.data.frame(bigcage)
df$chr_pos <- rownames(df)
df_cleaned <- na.omit(df)
head(df_cleaned)
df_cleaned$chr <- sub("\\..*", "", df_cleaned$chr_pos)  # Extract the part before the dot
df_cleaned$pos <- as.numeric(sub(".*\\.", "", df_cleaned$chr_pos))  # Extract and convert the part after the dot
head(df_cleaned)
df <- df_cleaned
df$logp <- -log10(df[, 1])  # Assuming p-values are in the first column
filtered_data <- df

colnames(filtered_data)[colnames(filtered_data) == "chr"] <- "Chromosome"
colnames(filtered_data)[colnames(filtered_data) == "pos"] <- "Position"

chromosome_levels <- c("2L", "2R", "3L", "3R", "X")
filtered_data <- filtered_data %>%
  filter(Chromosome %in% chromosome_levels)
filtered_data$Chromosome <- factor(filtered_data$Chromosome, levels = chromosome_levels)
chromosome_offsets <- filtered_data %>%
  group_by(Chromosome) %>%
  summarize(max_pos = max(Position, na.rm = TRUE)) %>%
  mutate(offset = cumsum(lag(max_pos, default = 0)))

filtered_data <- filtered_data %>%
  left_join(chromosome_offsets, by = "Chromosome") %>%
  mutate(CumulativePosition = Position + offset)
ggplot(filtered_data, aes(x = CumulativePosition, y = logp, color = Chromosome)) +
  geom_point(alpha = 0.2, size = 0.1) +
  labs(title = "Manhattan Plot 57", x = "Chromosomes", y = "-log10(pvalue)") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), # Remove vertical grid lines
    panel.spacing = unit(0.5, "lines"),   # Adjust spacing between panels
    legend.position = "none"             # Remove legend
  ) +
  scale_color_manual(values = c("red", "blue", "red", "blue", "red")) + # Set colors
  scale_x_continuous(breaks = chromosome_offsets$offset + chromosome_offsets$max_pos / 2,
                     labels = chromosome_levels) +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, by = 20))

## we check SNPs
# chromosome_value <- filtered_data[which(filtered_data$logp > 93), "Chromosome"]
# position_value <- filtered_data[which(filtered_data$logp > 93), "Position"]
# 
# # Create the corresponding chr.pos value
# chr_pos_value <- paste(chromosome_value, position_value, sep = ".")
# af_matrix_filtered[rownames(af_matrix_filtered) %in% chr_pos_value, ]
# cov_matrix_filtered[rownames(cov_matrix_filtered) %in% chr_pos_value, ]


#### pca for the snp: presence, absence of the SNPs
AF_bin <- data.frame(
  base = ifelse(AF_alt$V5 > 0, 1, 0),
  evolved = ifelse(!is.na(AF_alt$V6) & AF_alt$V6 > 0, 1, 0)
)


pca_result <- prcomp(AF_bin_transposed_non_constant, scale. = TRUE)  # Apply PCA
summary(pca_result)

# Create a data frame for plotting with PCA results
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

# Assuming the first half of rows in AF_bin correspond to "Base" and the second half to "Evolved"
n_base <- floor(nrow(AF_bin_transposed_non_constant) / 2)  # Assuming 50% base and 50% evolved
pca_data$population <- rep(c("Base", "Evolved"), each = n_base)

# Plot the PCA results
library(ggplot2)

ggplot(pca_data, aes(x = PC1, y = PC2, color = population)) +
  geom_point(alpha = 0.5) +  # Points with transparency
  labs(title = "PCA of SNP Differentiation Between Base and Evolved Populations",
       x = "Principal Component 1", y = "Principal Component 2") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "red"))


##### PCA plot:
xfbe <- fread("/Users/siva/Desktop/Rawdata/af/Base100k_base16k_tb20_bc40_ts20_ts40_ts60_ts80_mpileup_flt.xf.txt")
covbe  <- fread("/Users/siva/Desktop/Rawdata/af/Base100k_base16k_tb20_bc40_ts20_ts40_ts60_ts80_mpileup_flt.cov.txt")
afbe <- fread("/Users/siva/Desktop/Rawdata/af/Base100k_base16k_tb20_bc40_ts20_ts40_ts60_ts80_mpileup_flt.af.txt")

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

afbe_filtered_0.01 <- apply(afbe_filtered[, 5:12], 1, function(row) any(row >= 0.01))
afbe_filtered_step2 <- afbe_filtered[afbe_filtered_0.01, ]
afbe_filtered_0.99 <- apply(afbe_filtered_step2[, 5:12], 1, function(row) any(row <= 0.99))
afbe_filtered_final2 <- afbe_filtered_step2[afbe_filtered_0.99,]

#
afbe_filtered_final_dif <- afbe_filtered_final
afbe_filtered_final_dif[, 5:12] <- lapply(afbe_filtered_final_dif[, 5:12], function(x) as.numeric(as.character(x)))
afbe_filtered_final_dif$base16k_ts20 <- abs(afbe_filtered_final_dif$base16k - afbe_filtered_final_dif$ts20)
afbe_filtered_final_dif$base16k_ts40 <- abs(afbe_filtered_final_dif$base16k - afbe_filtered_final_dif$ts40)
afbe_filtered_final_dif$base16k_ts60 <- abs(afbe_filtered_final_dif$base16k - afbe_filtered_final_dif$ts60)
afbe_filtered_final_dif$base16k_ts80 <- abs(afbe_filtered_final_dif$base16k - afbe_filtered_final_dif$ts80)
head(afbe_filtered_final_dif
)

