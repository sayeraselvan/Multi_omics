install.packages("~/Downloads/hapsim_0.31.tar.gz")
library(hapsim)
data(ACEdata)
library(poolSeq)
library(ACER)
library("data.table")
library("foreach")
library("doParallel")
library("poolSeq")
library("tidyr")
library("tidyverse")
library("dplyr")
library("ggplot2")
library("ggpubr")

allelefreqs(ACEdata)

data(ACEdata)
x <- allelefreqs(ACEdata)
hist(x$freqs)

divlocus(ACEdata)
haplodata(ACEdata)
freqs <- haplofreqs(ACEdata, 1, 52)
freqs

BC_Sync40 <- read.sync(file="/Users/siva/Desktop/Rawdata/Base100k_Base16k_maxdepth8k_sts_minq40_snp.sync", gen=c(0,20), repl=c(1,1))
BC_Sync20 <- read.sync(file="/Users/siva/Desktop/Rawdata/Base100k_Base16k_maxdepth8k_sts_minq20_snp.sync", gen=c(0,20), repl=c(1,1))

BC_Sync40_snp1 <- read.sync(file="/Users/siva/Desktop/Rawdata/Base100k_Base16k_maxdepth8k_sts_minq40_snp1.sync", gen=c(0,20), repl=c(1,1))
BC_Sync20_snp1 <- read.sync(file="/Users/siva/Desktop/Rawdata/Base100k_Base16k_maxdepth8k_sts_minq20_snp1.sync", gen=c(0,20), repl=c(1,1))

SC_Sync <- read.sync(file="/Users/siva/Desktop/Rawdata/evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sst_minq20_snp.sync", gen=c(0,10,20,30,40,50), repl=c(1,1,1,1,1,1))

af_BC_40 <- af(sync = BC_Sync40, repl = c(1,1), gen = c(0,20))
cov_BC_40 <- coverage(sync = BC_Sync40, repl = c(1,1), gen = c(0,20))
af_BC_20 <- af(sync = BC_Sync20, repl = c(1,1), gen = c(0,20))
cov_BC_20 <- coverage(sync = BC_Sync20, repl = c(1,1), gen = c(0,20))

af_BC_40_snp1 <- af(sync = BC_Sync40_snp1, repl = c(1,1), gen = c(0,20))
cov_BC_40_snp1 <- coverage(sync = BC_Sync40_snp1, repl = c(1,1), gen = c(0,20))
af_BC_20_snp1 <- af(sync = BC_Sync20_snp1, repl = c(1,1), gen = c(0,20))
cov_BC_20_snp1 <- coverage(sync = BC_Sync20_snp1, repl = c(1,1), gen = c(0,20))

af_evol <- af(sync = SC_Sync, gen=c(0,10,20,30,40,50), repl=c(1,1,1,1,1,1))
cov_evol <- coverage(sync = SC_Sync, gen=c(0,10,20,30,40,50), repl=c(1,1,1,1,1,1))

head(cov_BC_gen)
head(af_BC_gen)
head(af_SC_gen)
head(cov_SC_gen)

rows_with_zero <- which(rowSums(cov_BC_40 == 0) > 0)
cov_BC_40_0 <- cov_BC_40[-rows_with_zero, ]
af_BC_40_0 <- af_BC_40[-rows_with_zero, ]

rows_with_zero <- which(rowSums(cov_BC_20 == 0) > 0)
cov_BC_20_0 <- cov_BC_20[-rows_with_zero, ]
af_BC_20_0 <- af_BC_20[-rows_with_zero, ]

rows_with_zero <- which(rowSums(cov_BC_40_snp1 == 0) > 0)
cov_BC_40_snp1_0 <- cov_BC_40_snp1[-rows_with_zero, ]
af_BC_40_snp1_0 <- af_BC_40_snp1[-rows_with_zero, ]

rows_with_zero <- which(rowSums(cov_BC_20_snp1 == 0) > 0)
cov_BC_20_snp1_0 <- cov_BC_20_snp1[-rows_with_zero, ]
af_BC_20_snp1_0 <- af_BC_20_snp1[-rows_with_zero, ]

# rows_with_zero_SC <- which(rowSums(cov_SC_gen == 0) > 0)
# cov_SC_gen <- cov_SC_gen[-rows_with_zero_SC, ]
# af_SC_gen <- af_SC_gen[-rows_with_zero_SC, ]

head(cov_BC_gen)
head(af_BC_gen)

# 
# file_manage <- function(af_gen0) {
#   df_gen0 <- data.frame(chrom_pos = names(af_gen0), af = af_gen0)
#   df_gen0 <- df_gen0 %>% separate(chrom_pos, into = c("CHROM", "POS"), sep = "\\.")
#   rownames(df_gen0) <- NULL
#   return(df_gen0)
# }
# file_manage_cov <- function(af_gen0) {
#   df_gen0 <- data.frame(chrom_pos = names(af_gen0), cov = af_gen0)
#   df_gen0 <- df_gen0 %>% separate(chrom_pos, into = c("CHROM", "POS"), sep = "\\.")
#   rownames(df_gen0) <- NULL
#   return(df_gen0)
# }

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

file_manage_and_rename_SC <- function(data) {
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
  colnames(data_df) <- c("CHROM", "POS", "TB20", "BC40", "TS20", "TS40", "TS60", "TS80")
  rownames(data_df) <- NULL
  return(data_df)
}

# af_BC_40 <- af(sync = BC_Sync40, repl = c(1,1), gen = c(0,20))
# cov_BC_40 <- coverage(sync = BC_Sync40, repl = c(1,1), gen = c(0,20))
# af_BC_20 <- af(sync = BC_Sync20, repl = c(1,1), gen = c(0,20))
# cov_BC_20 <- coverage(sync = BC_Sync20, repl = c(1,1), gen = c(0,20))
# 
# af_BC_40_snp1 <- af(sync = BC_Sync40_snp1, repl = c(1,1), gen = c(0,20))
# cov_BC_40_snp1 <- coverage(sync = BC_Sync40_snp1, repl = c(1,1), gen = c(0,20))
# af_BC_20_snp1 <- af(sync = BC_Sync20_snp1, repl = c(1,1), gen = c(0,20))
# cov_BC_20_snp1 <- coverage(sync = BC_Sync20_snp1, repl = c(1,1), gen = c(0,20))

cov_BC_40_df <- file_manage_and_rename(cov_BC_40)
cov_BC_20_df <- file_manage_and_rename(cov_BC_20)
cov_BC_40_snp1_df <- file_manage_and_rename(cov_BC_40_snp1)
cov_BC_20_snp1_df <- file_manage_and_rename(cov_BC_20_snp1)

cov_evol_df <- file_manage_and_rename_SC(cov_evol)

input_objects <- list(
  "cov_BC_20" = "cov_BC_20_df",
  "cov_BC_40_snp1" = "cov_BC_40_snp1_df",
  "cov_BC_20_snp1" = "cov_BC_20_snp1_df"
)
input_objects <- list(
  "cov_BC_40_0" = "cov_BC_40_0_df",
  "cov_BC_20_0" = "cov_BC_20_0_df",
  "cov_BC_40_snp1_0" = "cov_BC_40_0_snp1_df",
  "cov_BC_20_snp1_0" = "cov_BC_20_0_snp1_df"
)

for (input_name in names(input_objects)) {
  output_name <- input_objects[[input_name]]
  assign(output_name, file_manage_and_rename(get(input_name)))
}

af_BC_40_df <- file_manage_and_rename(af_BC_gen)
af_BC_20_df <- file_manage_and_rename(af_BC_gen)
af_BC_40_snp1
af_BC_20_snp1

library(ggplot2)
ggplot(cov_SC_gen_df, aes(x = average_coverage)) +
  geom_histogram(binwidth = 1, fill = 'skyblue', color = 'skyblue', alpha = 0.7) +
  labs(title = "Coverage Distribution for Small cage ", 
       x = "Coverage Depth", 
       y = "Frequency") +
  theme_minimal()

cov_BC_40_df$average_coverage <- rowMeans(cov_BC_40_df[, c('Base100k', 'Base16k')])
coverage_BC_2_percentile <- quantile(cov_BC_40_df$`Base16k`, 0.02)
coverage_BC_2_percentile

cov_SC_gen_df$average_coverage <- rowMeans(cov_SC_gen_df[, c('0_r1', '20_r1', '40_r1','60_r1','80_r1')])
coverage_SC_2_percentile <- quantile(cov_SC_gen_df$average_coverage, 0.02)
coverage_SC_2_percentile
coverage_SC_2_percentile <- quantile(cov_SC_gen_df$'0_r1', 0.02)

coverage_evol_2_percentile <- quantile(cov_evol_df$`TS20`, 0.02)
coverage_evol_2_percentile

cov_evol_df_no_zeros <- cov_evol_df
cov_evol_df_no_zeros[cov_evol_df_no_zeros == 0] <- NA
coverage_2nd_percentiles <- sapply(cov_evol_df_no_zeros[, -c(1,2)], function(x) quantile(x, probs = 0.01, na.rm = TRUE))
coverage_2nd_percentiles

cov_BC_20_df_no_zeros <- cov_BC_20_df
cov_BC_20_df_no_zeros[cov_BC_20_df_no_zeros == 0] <- NA
coverage_2nd_percentiles <- sapply(cov_BC_20_df_no_zeros[, -c(1,2)], function(x) quantile(x, probs = 0.02, na.rm = TRUE))
coverage_2nd_percentiles

cov_BC_40_df_no_zeros <- cov_BC_40_df
cov_BC_40_df_no_zeros[cov_BC_40_df_no_zeros == 0] <- NA
coverage_2nd_percentiles <- sapply(cov_BC_40_df_no_zeros[, -c(1,2)], function(x) quantile(x, probs = 0.02, na.rm = TRUE))
coverage_2nd_percentiles


####
coverage_columns <- cov_BC_20_df_no_zeros[, -c(1, 2)]
mean_coverage_per_sample <- colMeans(coverage_columns, na.rm = TRUE)
overall_average_coverage <- mean(mean_coverage_per_sample)

coverage_columns <- cov_evol_df_no_zeros[, -c(1, 2)]
mean_coverage_per_sample <- colMeans(coverage_columns, na.rm = TRUE)
overall_average_coverage <- sum(mean_coverage_per_sample)

# Print the average coverage for each sample and the overall average coverage
print(mean_coverage_per_sample)
print(overall_average_coverage)


