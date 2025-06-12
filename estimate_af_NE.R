## Allele frequency and FIsher exact test
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

BC_Sync <- read.sync(file="/Users/siva/Desktop/Rawdata/BC_base100k_gen20_40_rep4_snp.sync", gen=c(0,20,40), repl=c(1,1,1))
SC_Sync <- read.sync(file="/Users/siva/Desktop/Rawdata/SC_base16k_gen20_40_60_80_rep4_snp.sync", gen=c(0,20,40,60,80), repl=c(1,1,1,1,1))

af_BC_gen <- af(sync = BC_Sync, repl = c(1,1,1), gen = c(0,20,40))
cov_BC_gen <- coverage(sync = BC_Sync, repl = c(1,1,1), gen = c(0,20,40))

af_SC_gen <- af(sync = SC_Sync, repl = c(1,1,1,1,1), gen = c(0,20,40,60,80))
cov_SC_gen <- coverage(sync = SC_Sync, repl = c(1,1,1,1,1), gen = c(0,20,40,60,80))

head(cov_BC_gen)
head(af_BC_gen)
head(af_SC_gen)
head(cov_SC_gen)

## to remove the sites with 0 coverage
rows_with_zero <- which(rowSums(cov_BC_gen == 0) > 0)
cov_BC_gen <- cov_BC_gen[-rows_with_zero, ]
af_BC_gen <- af_BC_gen[-rows_with_zero, ]

rows_with_zero_SC <- which(rowSums(cov_SC_gen == 0) > 0)
cov_SC_gen <- cov_SC_gen[-rows_with_zero_SC, ]
af_SC_gen <- af_SC_gen[-rows_with_zero_SC, ]

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
    colnames(data_df) <- c("CHROM", "POS", "0_r1", "20_r1", "40_r1")
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
  colnames(data_df) <- c("CHROM", "POS", "0_r1", "20_r1", "40_r1", "60_r1", "80_r1")
  rownames(data_df) <- NULL
  return(data_df)
}

af_BC_gen_df <- file_manage_and_rename(af_BC_gen)
cov_BC_gen_df <- file_manage_and_rename(cov_BC_gen)

af_SC_gen_df <- file_manage_and_rename_SC(af_SC_gen)
cov_SC_gen_df <- file_manage_and_rename_SC(cov_SC_gen)

### vidya and yiwen script:

Ne.estimator <- function(dp, af, replicates, nb_SNPs, nb_rounds, time,
                         poolSize, nb_replicates, census, chrom) {
  chr_interest <- chrom
  dt.dp <- as.data.table(dp[dp$CHROM == chr_interest,])
  dt.af <- as.data.table(af[af$CHROM == chr_interest,])
  setkey(dt.af, CHROM, POS)
  setkey(dt.dp, CHROM, POS)
  # dt.af[, var:= rowVars(as.matrix(dt.af[,6:13]))]
  # dt.af <- dt.af[var!=0, -14]
  # dt.dp <-dt.dp[dt.af[,.(CHROM,POS)]]
  ne_estimates <- NULL
  for (r in replicates){ 
    print(r)
    # if (times[1]==0) {
    pref.af_i <- paste0(time[1], "_r", r) #time point i
    pref.dp_i <- paste0(time[1], "_r", r) #time point i
    # }
    # else{
    #   pref.af_i <- paste0("xf.","f",times[1], ".r", r) #time point i
    #   pref.dp_i <- paste0("dp.","f",times[1], ".r", r) #time point i
    # }
    pref.af_j <- paste0(time[2], "_r", r) #time point j with j>i
    pref.dp_j <- paste0(time[2], "_r", r) #time point j with j>i
    
    for(j in seq_len(nb_rounds)){
      cov_trial <- sample_n(dt.dp, size = nb_SNPs)
      SNP <- cov_trial[,.(CHROM,POS)] 
      af_trial <- dt.af[POS %in% SNP$POS,]
      pi <- unlist(subset(af_trial, select = pref.af_i))
      pj <- unlist(subset(af_trial, select = pref.af_j))
      covi <- unlist(subset(cov_trial, select = pref.dp_i))
      covj <- unlist(subset(cov_trial, select = pref.dp_j))
      ne <- estimateNe(p0 = pi, pt = pj, cov0 = covi, covt = covj, t = time[2]-time[1], 
                       ploidy = 2, truncAF = 0.05, method = "P.planI", poolSize = poolSize, Ncensus = census)
      ne_estimates <- rbind(ne_estimates, data.frame(replicate = r, start = time[1], end = time[2], trial = j,  ne = ne))
    } 
  }
  #print(ne_estimates)
  ne_estimates <- na.omit(ne_estimates)
  ne <- ne_estimates$ne
  if(length(which(ne<0))>0){
    ne_estimates <- ne_estimates[-which(ne<0)]
  }
  ne_estimates$replicate <- as.factor(ne_estimates$replicate)
  # if(length(which(ne<0))>0){ne_estimates <- ne_estimates[-which(ne<0)]}; 
  # plot <- ggplot(ne_estimates) + geom_boxplot(aes(replicate, ne, group = replicate))
  # print(plot)
  median_ne <- t(sapply(replicates, function(x) {idx <- which(ne_estimates$replicate == x);
  val <- sort(ne_estimates$ne[idx]);
  med <- median(val);
  va <- var(val);
  up <- val[qbinom(1-0.025, length(idx), 0.5)+1]; # CI 95% of a median
  down <- val[qbinom(0.025, length(idx), 0.5)]; # CI 95% of a median
  round(c(down,med,up,va))}))
  NE_result <- data.frame(replicate = replicates, type = rep(chrom, nb_replicates), 
                          start = rep(time[1], nb_replicates), end = rep(time[2], nb_replicates), 
                          median_ne = median_ne[, 2], var_ne = median_ne[, 4],
                          CI = paste0(median_ne[,1], "-",median_ne[,3]))
  NE_result$replicate <- as.factor(NE_result$replicate)
  return(NE_result)
}

# clean_1 is removing only the 0s in the dataframe of af
af_BC_gen_df[, 3:ncol(af_BC_gen_df)] <- lapply(af_BC_gen_df[, 3:ncol(af_BC_gen_df)], as.numeric)
rows_af_BC_zero <- apply(af_BC_gen_df[, 3:ncol(af_BC_gen_df)], 1, function(row) any(row == 0))
af_BC_gen_df_clean_1 <- af_BC_gen_df[!rows_af_BC_zero, ]
cov_BC_gen_df_clean_1 <- cov_BC_gen_df[!rows_af_BC_zero,]
af_BC_gen_df_no_clean_1 <- af_BC_gen_df[rows_af_BC_zero, ]
head(af_BC_gen_df_clean_1)
head(rows_af_BC_zero)

# clean_2 is removing only the 1 in the dataframe
rows_af_BC_one <- apply(af_BC_gen_df_clean_1[, 3:ncol(af_BC_gen_df_clean_1)], 1, function(row) any(row == 1))
af_BC_gen_df_clean_2 <- af_BC_gen_df_clean_1[!rows_af_BC_one, ]
cov_BC_gen_df_clean_2 <- cov_BC_gen_df_clean_1[!rows_af_BC_one,]
af_BC_gen_df_no_clean_2 <- af_BC_gen_df_no_clean_1[rows_af_BC_one, ]
head(af_BC_gen_df_clean_2)
head(rows_af_BC_one)

# clean 3 is removing less than 0.05 in the dataframe
rows_af_BC_0.05 <- apply(af_BC_gen_df_clean_2[, 3:ncol(af_BC_gen_df_clean_2)], 1, function(row) any(row < 0.05))
af_BC_gen_df_clean_3 <- af_BC_gen_df_clean_2[!rows_af_BC_0.05, ]
cov_BC_gen_df_clean_3 <- cov_BC_gen_df_clean_2[!rows_af_BC_0.05,]
af_BC_gen_df_no_clean_3 <- af_BC_gen_df_no_clean_2[rows_af_BC_0.05, ]
head(af_BC_gen_df_clean_3)
head(rows_af_BC_0.05)

# clean 4 is removing more than 0.95 in the dataframe
rows_af_BC_0.95 <- apply(af_BC_gen_df_clean_3[, 3:ncol(af_BC_gen_df_clean_3)], 1, function(row) any(row > 0.95))
af_BC_gen_df_clean_4 <- af_BC_gen_df_clean_3[!rows_af_BC_0.95, ]
cov_BC_gen_df_clean_4 <- cov_BC_gen_df_clean_3[!rows_af_BC_0.95,]
af_BC_gen_df_no_clean_4 <- af_BC_gen_df_no_clean_3[rows_af_BC_0.95, ]
head(af_BC_gen_df_clean_4)
head(rows_af_BC_0.95)

# SC clean
# clean_1 is removing only the 0s in the dataframe of af
af_SC_gen_df[, 3:ncol(af_SC_gen_df)] <- lapply(af_SC_gen_df[, 3:ncol(af_SC_gen_df)], as.numeric)
rows_af_SC_zero <- apply(af_SC_gen_df[, 3:ncol(af_SC_gen_df)], 1, function(row) any(row == 0))
af_SC_gen_df_clean_1 <- af_SC_gen_df[!rows_af_SC_zero, ]
cov_SC_gen_df_clean_1 <- cov_SC_gen_df[!rows_af_SC_zero,]
af_SC_gen_df_no_clean_1 <- af_SC_gen_df[rows_af_SC_zero, ]
head(cov_SC_gen_df_clean_1)
head(rows_af_SC_zero)

# clean_2 is removing only the 1 in the dataframe
rows_af_SC_one <- apply(af_SC_gen_df_clean_1[, 3:ncol(af_SC_gen_df_clean_1)], 1, function(row) any(row == 1))
af_SC_gen_df_clean_2 <- af_SC_gen_df_clean_1[!rows_af_SC_one, ]
cov_SC_gen_df_clean_2 <- cov_SC_gen_df_clean_1[!rows_af_SC_one,]
af_SC_gen_df_no_clean_2 <- af_SC_gen_df_no_clean_1[rows_af_SC_one, ]
head(af_SC_gen_df_clean_2)
head(rows_af_SC_one)

# clean 3 is removing less than 0.05 in the dataframe
rows_af_SC_0.05 <- apply(af_SC_gen_df_clean_2[, 3:ncol(af_SC_gen_df_clean_2)], 1, function(row) any(row < 0.05))
af_SC_gen_df_clean_3 <- af_SC_gen_df_clean_2[!rows_af_SC_0.05, ]
cov_SC_gen_df_clean_3 <- cov_SC_gen_df_clean_2[!rows_af_SC_0.05,]
af_SC_gen_df_no_clean_3 <- af_SC_gen_df_no_clean_2[rows_af_SC_0.05, ]
head(af_SC_gen_df_clean_3)
head(rows_af_SC_0.05)

# clean 4 is removing more than 0.95 in the dataframe
rows_af_SC_0.95 <- apply(af_SC_gen_df_clean_3[, 3:ncol(af_SC_gen_df_clean_3)], 1, function(row) any(row > 0.95))
af_SC_gen_df_clean_4 <- af_SC_gen_df_clean_3[!rows_af_SC_0.95, ]
cov_SC_gen_df_clean_4 <- cov_SC_gen_df_clean_3[!rows_af_SC_0.95,]
af_SC_gen_df_no_clean_4 <- af_SC_gen_df_no_clean_3[rows_af_SC_0.95, ]
head(af_SC_gen_df_clean_4)
head(rows_af_SC_0.95)

# Example Output
head(rows_zero)

#Ne.estimator(DP.SO, AF.SO, replicates = 1:4, nb_SNPs = 1000, nb_rounds = 1000,
#            time = c(1,10), poolSize = c(500,500), nb_replicates = 4, census = 1250,
#             chrom = "2L")

#Ne.estimator(DP.OS, AF.OS, replicates = 1:4, nb_SNPs = 1000, nb_rounds = 1000,
#             time = c(1,10), poolSize = c(500,500), nb_replicates = 4, census = 1250,
#             chrom = "2L")

Ne.estimator(dp = cov_SC_gen_df_clean_4, af = af_SC_gen_df_clean_4, replicates = 1:1, nb_SNPs = 1000, nb_rounds = 100,
             time = c(60,80), poolSize = c(250,250), nb_replicates = 1, census = 800,
             chrom = "2L")
Ne.estimator(dp = cov_BC_gen_df_clean_4, af = af_BC_gen_df_clean_4, replicates = 1:1, nb_SNPs = 1000, nb_rounds = 100,
             time = c(20,40), poolSize = c(250,250), nb_replicates = 1, census = 100000,
             chrom = "X")

af_BC_gen_df_clean4_matrix <- as.matrix(af_BC_gen_df_clean_4[, c("0_r1", "20_r1", "40_r1")])
rownames(af_BC_gen_df_clean4_matrix) <- paste(af_BC_gen_df_clean_4$CHROM, af_BC_gen_df_clean_4$POS, sep = ".")
colnames(af_BC_gen_df_clean4_matrix) <- colnames(af_BC_gen)

cov_BC_gen_df_clean4_matrix <- as.matrix(cov_BC_gen_df_clean_4[, c("0_r1", "20_r1", "40_r1")])
rownames(cov_BC_gen_df_clean4_matrix) <- paste(cov_BC_gen_df_clean_4$CHROM, cov_BC_gen_df_clean_4$POS, sep = ".")
colnames(cov_BC_gen_df_clean4_matrix) <- colnames(cov_BC_gen)

af_SC_gen_df_clean4_matrix <- as.matrix(af_SC_gen_df_clean_4[, c("0_r1", "20_r1", "40_r1", "60_r1", "80_r1")])
rownames(af_SC_gen_df_clean4_matrix) <- paste(af_SC_gen_df_clean_4$CHROM, af_SC_gen_df_clean_4$POS, sep = ".")
colnames(af_SC_gen_df_clean4_matrix) <- colnames(af_SC_gen)

cov_SC_gen_df_clean4_matrix <- as.matrix(cov_SC_gen_df_clean_4[, c("0_r1", "20_r1", "40_r1", "60_r1", "80_r1")])
rownames(cov_SC_gen_df_clean4_matrix) <- paste(cov_SC_gen_df_clean_4$CHROM, cov_SC_gen_df_clean_4$POS, sep = ".")
colnames(cov_SC_gen_df_clean4_matrix) <- colnames(cov_SC_gen)

#p_values <- adapted.cmh.test(freq=afMat, coverage=covMat, Ne=rep(300, 3), gen=c(0,10), repl=1:3, poolSize=rep(1000, ncol(afMat)))
#p_values <- adapted.chisq.test(freq=afMat, coverage=covMat, Ne=300, gen=c(0,10), poolSize=rep(1000, ncol(afMat)))

af_SC_gen_df_clean4_matrix_0_80 <- af_SC_gen_df_clean4_matrix[,c(1,5)]
cov_SC_gen_df_clean4_matrix_0_80 <-  cov_SC_gen_df_clean4_matrix[,c(1,5)]

SC_p_values_0_80 <- adapted.chisq.test(freq = af_SC_gen_df_clean4_matrix_0_80 , coverage = cov_SC_gen_df_clean4_matrix_0_80 , Ne = 47, gen = c(0,80), poolSize=rep(250,2))

SC_p_values_0_80_sim <- adapted.chisq.test(freq = af_SC_gen_df_clean4_matrix_0_80 , coverage = cov_SC_gen_df_clean4_matrix_0_80 , Ne = 300, gen = c(0,80), poolSize=rep(250,2))

#### NE Loop:
# Define the parameters and combinations of timepoints
time_combinations <- list(c(0, 20), c(0, 40), c(0, 60), c(0, 80), 
                          c(20, 40), c(20, 60), c(20, 80), 
                          c(40, 60), c(40, 80))

results <- list()

for (i in seq_along(time_combinations)) {
  current_time <- time_combinations[[i]]
  ne_result <- Ne.estimator(
    dp = cov_SC_gen_df_clean_4,
    af = af_SC_gen_df_clean_4,
    replicates = 1:1,
    nb_SNPs = 1000,
    nb_rounds = 100,
    time = current_time,
    poolSize = c(250, 250),
    nb_replicates = 1,
    census = 800,
    chrom = "3R"
  )
  results[[paste0("time_", current_time[1], "_", current_time[2])]] <- ne_result
}

results

af_BC_gen_df_clean_4
af_SC_gen_df_clean_4

sum(cov_BC_gen_df$`20_r1` == 8)
cov_value <- 
sum(cov_BC_gen_df$`0_r1` == cov_value | cov_BC_gen_df$`20_r1` == cov_value | cov_BC_gen_df$`40_r1` == cov_value)
hi <- 0
for (cov_value in 1:49) {
  #result <- sum(cov_BC_gen_df$`0_r1` == cov_value | cov_BC_gen_df$`20_r1` == cov_value | cov_BC_gen_df$`40_r1` == cov_value)
  result <- sum(cov_BC_gen_df$`0_r1` == cov_value | cov_BC_gen_df$`20_r1` == cov_value)
  
  hi <- hi + result
  print(hi)
}
cov_value <- 8000
sum(cov_BC_gen_df$`0_r1` > cov_value | cov_BC_gen_df$`20_r1` > cov_value | cov_BC_gen_df$`40_r1` > cov_value)
sum(cov_BC_gen_df$`0_r1` > cov_value)

#`0_r1`

library(ggplot2)
ggplot(cov_SC_gen_df, aes(x = average_coverage)) +
  geom_histogram(binwidth = 1, fill = 'skyblue', color = 'skyblue', alpha = 0.7) +
  labs(title = "Coverage Distribution for Small cage ", 
       x = "Coverage Depth", 
       y = "Frequency") +
  theme_minimal()

cov_BC_gen_df$average_coverage <- rowMeans(cov_BC_gen_df[, c('0_r1', '20_r1', '40_r1')])
coverage_BC_2_percentile <- quantile(cov_BC_gen_df$`20_r1`, 0.02)
coverage_BC_2_percentile

cov_SC_gen_df$average_coverage <- rowMeans(cov_SC_gen_df[, c('0_r1', '20_r1', '40_r1','60_r1','80_r1')])
coverage_SC_2_percentile <- quantile(cov_SC_gen_df$'80_r1', 0.01)
coverage_SC_2_percentile
coverage_SC_2_percentile <- quantile(cov_SC_gen_df$'0_r1', 0.02)



#####

bi_FINAL_MASKED_SYNC_BASE_EVOLVED <- read.sync(file="~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync", 
                                               gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
af_bi_FINAL_MASKED_SYNC_BASE_EVOLVED <- af(sync = bi_FINAL_MASKED_SYNC_BASE_EVOLVED, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED <- coverage(sync = bi_FINAL_MASKED_SYNC_BASE_EVOLVED, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))

bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table <- read.table(file="~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic.sync")
                                               
colnames(af_bi_FINAL_MASKED_SYNC_BASE_EVOLVED) <- c("Base100k","TB20","BC40","Base16k", "TS20","TS40","TS60","TS80")
colnames(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED) <- c("Base100k","TB20","BC40","Base16k", "TS20","TS40","TS60","TS80")

cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes <- cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED
cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes[cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes == 0] <- NA

# coverage_2nd_percentiles <- sapply(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, function(x) quantile(x, probs = 0.01, na.rm = TRUE))
# coverage_2nd_percentiles


library(data.table)
dt <- as.data.table(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes)
dt <- as.data.table(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, keep.rownames = TRUE)
setnames(dt, "rn", "RowName")
head(dt)
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


library(data.table)

dt <- as.data.table(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes)

coverage_2nd_percentiles <- dt_no_zeroes[, lapply(.SD, function(x) quantile(x, probs = 0.01, na.rm = TRUE))]
print(coverage_2nd_percentiles)

# Count rows with NA values
na_rows_count <- sum(apply(cov_bi_FINAL_MASKED_SYNC_BASE_EVOLVED_no_zeroes, 1, function(row) any(is.na(row))))
print(na_rows_count)

# Remove rows where any column has a zero value
dt_no_zeroes <- dt[apply(dt, 1, function(row) all(row != 0))]
print(dt_no_zeroes)

mean_coverage_per_sample <- colMeans(dt_no_zeroes[,-1])
overall_average_coverage <- sum(mean_coverage_per_sample)

print(mean_coverage_per_sample)
print(overall_average_coverage)


### now we make the sync file only for the stuff where we focus on the biallelic and sites with no 0 coverage
chrompos <- paste(bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table$V1, bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table$V2, sep = ".")  

bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table_filtered <- bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table[chrompos %in% dt_no_zeroes$RowName, ]

write.table(bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table_filtered, "~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic_no0coverage_14million_snps.sync", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

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




### coverage plots
library(ggplot2)
library(data.table)

# Assuming dt_no_zeroes is your data.table
# Define the percentiles and average for each column
percentiles <- list(
  `Base100k` = c(0.01, 0.02, 0.98, 0.99),
  `TB20` = c(0.01, 0.02, 0.98, 0.99),
  `BC40` = c(0.01, 0.02, 0.98, 0.99),
  `Base16k` = c(0.01, 0.02, 0.98, 0.99),
  `TS20` = c(0.01, 0.02, 0.98, 0.99),
  `TS40` = c(0.01, 0.02, 0.98, 0.99),
  `TS60` = c(0.01, 0.02, 0.98, 0.99),
  `TS80` = c(0.01, 0.02, 0.98, 0.99)
)

# Compute the values of the percentiles and averages
percentile_values <- lapply(names(percentiles), function(col_name) {
  quantiles <- quantile(dt_no_zeroes[[col_name]], probs = percentiles[[col_name]], na.rm = TRUE)
  avg <- mean(dt_no_zeroes[[col_name]], na.rm = TRUE)
  c(quantiles, average = avg)
})

# Convert the result to a data.frame for easier manipulation
percentile_df <- as.data.frame(t(simplify2array(percentile_values)))
colnames(percentile_df) <- c("1st Percentile", "2nd Percentile", "98th Percentile", "99th Percentile", "Average")
percentile_df$Sample <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
# this is for the 14 million snps
# percentile_df <- data.table(
#   `1st Percentile` = c(57, 17, 25, 56, 14, 19, 16, 18),
#   `2nd Percentile` = c(63, 40, 57, 67, 35, 48, 40, 44),
#   `98th Percentile` = c(123, 168, 224, 136, 148, 197, 175, 185),
#   `99th Percentile` = c(128, 175, 234, 141, 155, 206, 182, 193),
#   `Average` = c(95.71544, 126.64230, 166.94786, 105.72177, 111.37802, 145.01925, 131.82221, 140.85954),
#   Sample = c("Base100k", "TB20", "BC40", "Base16k", "TS20", "TS40", "TS60", "TS80")
# )
#base16k_data <- percentile_df[Sample == "TS80"]
#sample_data <- data.table(value = dt_no_zeroes$TS80)  # Replace with actual data

percentile_df <- data.table(
  `1st Percentile` = c(49, 6, 8, 43, 5, 7, 6, 6),
  `2nd Percentile` = c(59, 17, 25, 56, 14, 19, 16, 17),
  `98th Percentile` = c(123, 166, 223, 135, 147, 196, 173, 183),
  `99th Percentile` = c(127, 173, 233, 140, 153, 204, 180, 190),
  `Average` = c(94.93376, 120.59545, 161.62808, 103.42605, 106.24906, 140.26912, 125.82939, 134.57591),
  Sample = c("Base100k", "TB20", "BC40", "Base16k", "TS20", "TS40", "TS60", "TS80")
)

base16k_data <- percentile_df[Sample == "BC40"]

# Print the extracted row for verification
print(base16k_data)

sample_data <- data.table(value = dt$BC40)  # Replace with actual data
# 
# # Create a function to plot the coverage distribution for the 'Base16k' sample
# plot_coverage_distribution <- function(sample_data, base16k_data) {
#   ggplot(sample_data, aes(x = value)) +
#     geom_density(fill = "skyblue", color = "black", alpha = 0.7) +
#     geom_vline(aes(xintercept = base16k_data$`1st Percentile`), color = "red", linetype = "dashed") + # 1st Percentile
#     geom_vline(aes(xintercept = base16k_data$`2nd Percentile`), color = "orange", linetype = "dashed") + # 2nd Percentile
#     geom_vline(aes(xintercept = base16k_data$`98th Percentile`), color = "green", linetype = "dashed") + # 98th Percentile
#     geom_vline(aes(xintercept = base16k_data$`99th Percentile`), color = "purple", linetype = "dashed") + # 99th Percentile
#     geom_vline(aes(xintercept = base16k_data$Average), color = "blue", linetype = "solid", size = 1) + # Average
#     labs(title = "Coverage Distribution for Base16k",
#          x = "Coverage", y = "Frequency") +
#     theme_minimal() +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     theme(legend.position = "none")
# }
# 
# # Generate the plot for 'Base16k'
# p <- plot_coverage_distribution(sample_data, base16k_data)
# 
# # Display the plot
# print(p)

first_percentile <- base16k_data$`1st Percentile`
ninety_ninth_percentile <- base16k_data$`99th Percentile`
filtered_data <- sample_data[value >= first_percentile & value <= ninety_ninth_percentile]

plot(density(filtered_data$value), main = "Coverage Density for Base16k",
     xlab = "Coverage", ylab = "Density", col = "blue")

ggplot(filtered_data, aes(x = value)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Coverage Density for BC40",
       x = "Coverage",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 250)) +  # Set xlim with a range between 10 and 250
  ylim(0, 0.05)  # Adjust y-axis to fit the density curve

### only the biallelic matched from the two base populations out of 14 million snps, we have 4 million snps

biallelic_matched.snps.sync_4mil <- read.sync(file="~/Desktop/Rawdata/biallelic_matched.snps.sync", 
                                               gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
af_biallelic_matched.snps.sync_4mil <- af(sync = biallelic_matched.snps.sync_4mil, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_biallelic_matched.snps.sync_4mil <- coverage(sync = biallelic_matched.snps.sync_4mil, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
colnames(af_biallelic_matched.snps.sync_4mil) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
colnames(cov_biallelic_matched.snps.sync_4mil) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")


## out of this, we need to have same polymorphic state snps in both base and evolved, out of 4 millions snps, we have 3.2 million snps
filtered_no_new_alleles.snps.sync_3.2 <- read.sync(file="~/Desktop/Rawdata/filtered_no_new_alleles.snps.sync", 
                                            gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
af_filtered_no_new_alleles.snps.sync_3.2 <- af(sync = filtered_no_new_alleles.snps.sync_3.2, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_filtered_no_new_alleles.snps.sync_3.2 <- coverage(sync = filtered_no_new_alleles.snps.sync_3.2, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
colnames(af_filtered_no_new_alleles.snps.sync_3.2) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
colnames(cov_filtered_no_new_alleles.snps.sync_3.2) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")

dt <- as.data.table(cov_filtered_no_new_alleles.snps.sync_3.2)

coverage_2nd_percentiles <- dt[, lapply(.SD, function(x) quantile(x, probs = 0.98, na.rm = TRUE))]
print(coverage_2nd_percentiles)

mean_coverage_per_sample <- colMeans(dt)
overall_average_coverage <- sum(mean_coverage_per_sample)
print(mean_coverage_per_sample)
print(overall_average_coverage)


#### Now lets go to STS (allele frequnecy spectrum of te base and evolved)
filtered_no_new_alleles.snps.sync_3.2Pris <- read.sync(file="~/Desktop/Rawdata/filtered_no_new_alleles.snps.sync",
                                                   gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2), polarization = "rising")
af_filtered_no_new_alleles.snps.sync_3.2Pris <- af(sync = filtered_no_new_alleles.snps.sync_3.2Pris, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
cov_filtered_no_new_alleles.snps.sync_3.2Pris <- coverage(sync = filtered_no_new_alleles.snps.sync_3.2Pris, gen=c(0,0,20,40,20,40,60,80), repl=c(1,2,1,1,2,2,2,2))
colnames(af_filtered_no_new_alleles.snps.sync_3.2Pris) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")
colnames(cov_filtered_no_new_alleles.snps.sync_3.2Pris) <- c("Base100k","TB20","BC40","Base16k","TS20","TS40","TS60","TS80")

Base100k <- af_filtered_no_new_alleles.snps.sync_3.2[,4]
af_filtered_no_new_alleles_filtered <- af_filtered_no_new_alleles.snps.sync_3.2[
  Base100k >= 0.01 & Base100k <= 0.99, ]

af_filtered_no_new_alleles_filtered <- af_filtered_no_new_alleles.snps.sync_3.2[
  apply(af_filtered_no_new_alleles.snps.sync_3.2, 1, function(row) all(row >= 0.01 & row <= 0.99)), ]

af_filtered_no_new_alleles_filtered <- af_filtered_no_new_alleles.snps.sync_3.2Pris[
  apply(af_filtered_no_new_alleles.snps.sync_3.2Pris, 1, function(row) all(row != 0 & row != 1)), ]

af_filtered_no_new_alleles_filtered_df <- as.data.table(af_filtered_no_new_alleles_filtered, keep.rownames = TRUE)
af_filtered_no_new_alleles_filtered_df[, Chromosome := gsub("^([A-Za-z0-9]+)\\..*$", "\\1", rn)]
head(af_filtered_no_new_alleles_filtered_df)
unique(af_filtered_no_new_alleles_filtered_df$Chromosome)

chromosome_2L <- af_filtered_no_new_alleles_filtered_df[Chromosome == "2L"]

ggplot(chromosome_2L, aes(x = Base100k)) +
  geom_histogram(bins = 100, fill = "blue", color = "black", alpha = 0.5) +
  labs(title = "Site Frequency Spectrum for Chromosome 2L",
       x = "Frequency",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) + 
  scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, by = 200))  # Adjust the y-axis limits and breaks


af_filtered_selected_chromosomes <- af_filtered_no_new_alleles_filtered_df[
  af_filtered_no_new_alleles_filtered_df$Chromosome %in% c("2L", "2R", "3L", "3R", "X"), ]

# Plot the Site Frequency Spectrum (SFS) for these chromosomes
ggplot(af_filtered_selected_chromosomes, aes(x = BC40)) +
  geom_histogram(bins  = 100, fill = "blue", color = "black", alpha = 0.5) +
  labs(title = "SFS for BC40: chrom 2L, 2R, 3L, 3R, and X",
       x = "Frequency",
       y = "Density") +
  theme_minimal() +
  scale_x_continuous(limits = c(0, 1)) + 
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 600, by = 200)) +
  facet_wrap(~Chromosome, nrow = 5)  
# # Load necessary library
# library(ggplot2)
# 
# # Plot the SFS for Base100k per chromosome
# ggplot(af_filtered_no_new_alleles_filtered, aes(x = Base100k)) +
#   geom_histogram(bins = 1, fill = "blue", color = "black", alpha = 0.5) +
#   facet_wrap(~ Chromosome, scales = "free_y") + # Separate plot for each chromosome
#   labs(title = "Site Frequency Spectrum for Base100k by Chromosome",
#        x = "Frequency",
#        y = "Count") +
#   theme_minimal()
# 
# 
# # View the filtered data
# print(af_filtered_no_new_alleles_filtered)
# 
# # Base100k <- af_filtered_no_new_alleles.snps.sync_3.2Pref[,8]
# 
# Base100k_log <- log10(Base100k + 0.001)
# 
# # Plot the SFS
# ggplot(data.frame(Frequency = Base100k), aes(x = Frequency)) +
#   geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.7) +
#   labs(
#     title = "Site Frequency Spectrum (SFS) for Base100k",
#     x = "Allele Frequency",
#     y = "Count"
#   ) +
#   theme_minimal()
# 
# library(ggplot2)
# 
# ggplot(data.frame(Frequency = Base100k), aes(x = Frequency)) +
#   geom_density(fill = "blue", alpha = 0.5) +
#   labs(title = "Density Plot for Base100k",
#        x = "Frequency",
#        y = "Density") +
#   theme_minimal()
# library(ggplot2)
# 
# ggplot(data.frame(Frequency = Base100k), aes(x = Frequency)) +
#   geom_histogram(binwidth = 0.005, fill = "blue", color = "black", alpha = 0.7) +
#   scale_y_continuous(trans = "log10") +  # Apply log10 transformation to the y-axis
#   labs(title = "Density Plot for Base100k with Log-Scaled Y-Axis",
#        x = "Frequency",
#        y = "Density (Log-Scaled)") +
#   theme_minimal()
# ggplot(data.frame(Frequency = Base100k), aes(x = Frequency)) +
#   geom_histogram(binwidth = 0.005, fill = "blue", color = "black", alpha = 0.7) +
#   scale_y_sqrt() +  # Apply square root transformation to the y-axis
#   labs(title = "Density Plot for Base100k with Square Root-Scaled Y-Axis",
#        x = "Frequency",
#        y = "Density (Sqrt-Scaled)") +
#   theme_minimal()


input_file <- "~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic_no0coverage_14million_snps.sync"
data <- read.table(input_file, header = FALSE, stringsAsFactors = FALSE)

is_biallelic_ATCG <- function(column) {
  # Split the allele counts
  allele_counts <- strsplit(column, ":")
  
  # Count rows where exactly two non-zero values are found in the first 4 positions (A:T:C:G)
  biallelic <- sapply(allele_counts, function(row) {
    counts <- as.numeric(row[1:4])  # Only consider A, T, C, G
    sum(counts > 0) == 2  # Check if exactly 2 alleles are non-zero
  })
  
  return(biallelic)
}


# Apply the function to column 4
biallelic_column4 <- is_biallelic_ATCG(data$V6)

# Count the number of biallelic sites
biallelic_count <- sum(biallelic_column4)

# Output the result
cat("Number of biallelic sites in column 4:", biallelic_count, "\n")

columns_to_process <- 4:ncol(data)
for (col in columns_to_process) {
  biallelic_sites <- is_biallelic_ATCG(data[[col]])
  biallelic_count <- sum(biallelic_sites)
  cat("Number of biallelic sites in column", col, "(A:T:C:G only):", biallelic_count, "\n")
}

###/Users/siva/Desktop/Rawdata/
# Base100k_TB20_fisher_mincov57_mincount1 <- read.table("~/Desktop/Rawdata/FBase100k_TB20_mincov57_mincount1_no_na.fet")
# filtered_data <- Base100k_TB20_fisher_mincov57_mincount1
Base100k_TB20_fisher_mincov57_mincount1 <- read.table("~/Desktop/Rawdata/base16k_ts80_no_na.fet")
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

#### calculate individual af 
biallellic_all <- read.table("~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic_no0coverage_14million_snps.sync")
biallellic_one <- biallellic_all[c(1,2,3,5)]
head(biallellic_one)
write.table(biallellic_one, "~/Desktop/Rawdata/fet/Base16k_14millionssnps.sync", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

base16k_14mil <- read.sync(file = "~/Desktop/Rawdata/fet/Base16k_14millionssnps.sync", gen = 0, rep = 1) 
af_base16k_14mil <- af(base16k_14mil, gen = 0, rep = 1)
cov_base16k_14mil <- coverage(base16k_14mil, gen = 0, rep = 1)
head(af_base16k_14mil)

## exploratory data analysis:
# vcf file to the 14million SNPs:
bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table_filtered <- read.table("~/Desktop/Rawdata/fet/Base100k_base16k_maxdepth8k_sts_filtered_minq20_from_minq40_snp_only.Evolved_R4_TB20_BC40_TS20_TS40_TS60_TS80_sts_minq20.sync_mincount21_indel_TE_masked_biallelic_no0coverage_14million_snps.sync")
dim(bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table_filtered)
chromposvcf <- read.table("~/Desktop/Rawdata/fet/chrom_pos_filteredvcf.txt")
chromposmafvcf <- read.table("~/Desktop/Rawdata/fet/chrom_pos_filteredmafvcf.txt")

filtered_data <- bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table_filtered %>%
  filter(paste(V1, V2) %in% paste(chromposvcf$V1, chromposvcf$V2))
filtered_datamaf <- bi_FINAL_MASKED_SYNC_BASE_EVOLVED_table_filtered %>%
  filter(paste(V1, V2) %in% paste(chromposmafvcf$V1, chromposmafvcf$V2))
dim(filtered_data)
dim(filtered_datamaf)


### whole minq 20 -> sync file: 

base100k <- read.table("~/Desktop/Rawdata/chrompos/base100k_mincount5m.chrompos")
base16k <- read.table("~/Desktop/Rawdata/chrompos/base100k_mincount5m.chrompos")
ts20 <- read.table("~/Desktop/Rawdata/chrompos/ts20_mincount5m.chrompos")
ts40 <- read.table("~/Desktop/Rawdata/chrompos/ts40_mincount5m.chrompos")
ts60 <- read.table("~/Desktop/Rawdata/chrompos/ts60_mincount5m.chrompos")
ts80 <- read.table("~/Desktop/Rawdata/chrompos/ts80_mincount5m.chrompos")
tb20 <- read.table("~/Desktop/Rawdata/chrompos/tb20_mincount5m.chrompos")
bc40 <- read.table("~/Desktop/Rawdata/chrompos/bc40_mincount5m.chrompos")

# Filter for specific chromosomes
ts20_filtered <- subset(ts20, V1 %in% c("2L", "2R", "3L", "3R", "X"))
base100k_filtered <- subset(base100k, V1 %in% c("2L", "2R", "3L", "3R", "X"))
base16k_filtered <- subset(base16k, V1 %in% c("2L", "2R", "3L", "3R", "X"))
ts40_filtered <- subset(ts40, V1 %in% c("2L", "2R", "3L", "3R", "X"))
ts60_filtered <- subset(ts60, V1 %in% c("2L", "2R", "3L", "3R", "X"))
ts80_filtered <- subset(ts80, V1 %in% c("2L", "2R", "3L", "3R", "X"))
tb20_filtered <- subset(tb20, V1 %in% c("2L", "2R", "3L", "3R", "X"))
bc40_filtered <- subset(bc40, V1 %in% c("2L", "2R", "3L", "3R", "X"))
base100k_positions <- paste(base100k_filtered$V1, base100k_filtered$V2, sep = ".")
bc40_positions <- paste(bc40_filtered$V1, bc40_filtered$V2, sep = ".")
tb20_positions <- paste(tb20_filtered$V1, tb20_filtered$V2, sep = ".")
# 
# # Step 5: Find common positions between the datasets
# common_base100k_bc40 <- intersect(base100k_positions, bc40_positions)
# common_base100k_tb20 <- intersect(base100k_positions, tb20_positions)
# common_bc40_tb20 <- intersect(bc40_positions, tb20_positions)
# 
# # Find the total number of unique positions for each pair
# only_base100k <- setdiff(base100k_positions, union(bc40_positions, tb20_positions))
# only_bc40 <- setdiff(bc40_positions, union(base100k_positions, tb20_positions))
# only_tb20 <- setdiff(tb20_positions, union(base100k_positions, bc40_positions))
# 
# # Step 6: Create a data frame for the pie chart
# data <- data.frame(
#   Category = c("Base100k only", "BC40 only", "TB20 only", "Base100k & BC40", "Base100k & TB20", "BC40 & TB20", "All 3"),
#   Count = c(length(only_base100k), length(only_bc40), length(only_tb20),
#             length(common_base100k_bc40), length(common_base100k_tb20), length(common_bc40_tb20),
#             length(intersect(common_base100k_bc40, tb20_positions)))
# )


base100k_positions <- paste(base100k_filtered$V1, base100k_filtered$V2, sep = ".")
bc40_positions <- paste(bc40_filtered$V1, bc40_filtered$V2, sep = ".")
tb20_positions <- paste(tb20_filtered$V1, tb20_filtered$V2, sep = ".")

# Step 5: Find the intersections and unique positions
only_base100k <- setdiff(base100k_positions, union(bc40_positions, tb20_positions))
only_bc40 <- setdiff(bc40_positions, union(base100k_positions, tb20_positions))
only_tb20 <- setdiff(tb20_positions, union(base100k_positions, bc40_positions))

common_base100k_bc40 <- intersect(base100k_positions, bc40_positions)
common_base100k_tb20 <- intersect(base100k_positions, tb20_positions)
common_bc40_tb20 <- intersect(bc40_positions, tb20_positions)

# Step 6: Prepare data for the Venn diagram
venn_data <- list(
  Base100k = base100k_positions,
  BC40 = bc40_positions,
  TB20 = tb20_positions
)

# Step 7: Plot the Venn diagram
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("Base100k", "BC40", "TB20"),
  filename = NULL,  # Do not save to file
  output = TRUE,
  imagetype = "png",  # Type of image
  height = 500, width = 500, resolution = 300,
  col = "transparent", fill = c("#F1C40F", "#E74C3C", "#8E44AD"),
  alpha = 0.5,
  label.col = "black",
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  cat.col = c("#F1C40F", "#E74C3C", "#8E44AD"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  rotation = 1
)

# Display the Venn diagram
grid.draw(venn.plot)


base16k_positions <- paste(base16k_filtered$V1, base16k_filtered$V2, sep = ".")
ts20_positions <- paste(ts20_filtered$V1, ts20_filtered$V2, sep = ".")
ts40_positions <- paste(ts40_filtered$V1, ts40_filtered$V2, sep = ".")
ts60_positions <- paste(ts60_filtered$V1, ts60_filtered$V2, sep = ".")
ts80_positions <- paste(ts80_filtered$V1, ts80_filtered$V2, sep = ".")

# Step 5: Find the intersections and unique positions
only_base16k <- setdiff(base16k_positions, union(union(ts20_positions, ts40_positions), union(ts60_positions, ts80_positions)))
only_ts20 <- setdiff(ts20_positions, union(union(base16k_positions, ts40_positions), union(ts60_positions, ts80_positions)))
only_ts40 <- setdiff(ts40_positions, union(union(base16k_positions, ts20_positions), union(ts60_positions, ts80_positions)))
only_ts60 <- setdiff(ts60_positions, union(union(base16k_positions, ts20_positions), union(ts40_positions, ts80_positions)))
only_ts80 <- setdiff(ts80_positions, union(union(base16k_positions, ts20_positions), union(ts40_positions, ts60_positions)))

common_base16k_ts20 <- intersect(base16k_positions, ts20_positions)
common_base16k_ts40 <- intersect(base16k_positions, ts40_positions)
common_base16k_ts60 <- intersect(base16k_positions, ts60_positions)
common_base16k_ts80 <- intersect(base16k_positions, ts80_positions)

common_ts20_ts40 <- intersect(ts20_positions, ts40_positions)
common_ts20_ts60 <- intersect(ts20_positions, ts60_positions)
common_ts20_ts80 <- intersect(ts20_positions, ts80_positions)

common_ts40_ts60 <- intersect(ts40_positions, ts60_positions)
common_ts40_ts80 <- intersect(ts40_positions, ts80_positions)

common_ts60_ts80 <- intersect(ts60_positions, ts80_positions)

common_all <- Reduce(intersect, list(base16k_positions, ts20_positions, ts40_positions, ts60_positions, ts80_positions))

# Step 6: Prepare data for the Venn diagram
# venn_data <- list(
#   Base16k = base16k_positions,
#   TS20 = ts20_positions,
#   TS40 = ts40_positions,
#   TS60 = ts60_positions,
#   TS80 = ts80_positions
# )
# 
# # Step 7: Plot the Venn diagram
# venn.plot <- venn.diagram(
#   x = venn_data,
#   category.names = c("Base16k", "TS20", "TS40", "TS60", "TS80"),
#   filename = NULL,  # Do not save to file
#   output = TRUE,
#   imagetype = "png",  # Type of image
#   height = 500, width = 500, resolution = 300,
#   col = "transparent", fill = c("#F1C40F", "#E74C3C", "#8E44AD", "#3498DB", "#2ECC71"),
#   alpha = 0.5,
#   label.col = "black",
#   cex = 1.5,
#   fontface = "bold",
#   fontfamily = "sans",
#   cat.col = c("#F1C40F", "#E74C3C", "#8E44AD", "#3498DB", "#2ECC71"),
#   cat.cex = 1.5,
#   cat.fontface = "bold",
#   cat.fontfamily = "sans",
#   rotation = 1
# )
# 
# # Display the Venn diagram
# grid.draw(venn.plot)

# Install and load UpSetR
if (!require(UpSetR)) install.packages("UpSetR")
library(UpSetR)

# Prepare data as a list of positions
list_data <- list(
  Base16k = base16k_positions,
  TS20 = ts20_positions,
  TS40 = ts40_positions,
  TS60 = ts60_positions,
  TS80 = ts80_positions
)

# Create a binary matrix for UpSetR
# Each column represents a dataset, and rows indicate membership in those datasets
intersection_matrix <- fromList(list_data)

# Plot the UpSet plot
upset(intersection_matrix, 
      nsets = 5,  # Number of sets to display
      order.by = "freq",  # Order intersections by frequency
      main.bar.color = "skyblue",  # Color of the main bars
      sets.bar.color = "darkblue",  # Color of the set bars
      text.scale = 1.2)  # Scale text size


list_data <- list(
  #Base100k = base100k_positions,
  #Base16k = base16k_positions,
  Base = base16k_positions,
  BC40 = bc40_positions,
  TB20 = tb20_positions,
  TS20 = ts20_positions,
  TS40 = ts40_positions,
  TS60 = ts60_positions,
  TS80 = ts80_positions
)

# Create a binary matrix for UpSetR
# Each column represents a dataset, and rows indicate membership in those datasets
intersection_matrix <- fromList(list_data)

# Plot the UpSet plot
upset(intersection_matrix, 
      nsets = length(list_data),  # Use all sets
      order.by = "freq",  # Order intersections by frequency
      main.bar.color = "skyblue",  # Color of the main bars
      sets.bar.color = "darkblue",  # Color of the set bars
      text.scale = 1.1,  # Scale text size
      keep.order = TRUE,  # Keep the order of datasets as in the list
      mb.ratio = c(0.7, 0.3))  


row_counts <- data.frame(
  Dataset = c("base16k", "ts20", "ts40", "ts60", "ts80", "base100k", "tb20", "bc40"),
  Rows = c(
    nrow(base16k_filtered), 
    nrow(ts20_filtered), 
    nrow(ts40_filtered), 
    nrow(ts60_filtered), 
    nrow(ts80_filtered), 
    nrow(base100k_filtered), 
    nrow(tb20_filtered), 
    nrow(bc40_filtered)
  )
)



print(row_counts)


# present in evlved but not in base
# filtered_data <- list(ts20_filtered, ts40_filtered, ts60_filtered, ts80_filtered, tb20_filtered, bc40_filtered)
# filtered_positions <- unique(do.call(rbind, filtered_data))

base_filtered_positions <- unique(rbind(base100k_filtered, base16k_filtered))
filtered_data <- list(ts20_filtered, ts40_filtered, ts60_filtered, ts80_filtered, tb20_filtered, bc40_filtered)
intersection_filtered_positions <- Reduce(function(x, y) merge(x, y, by = c("V1", "V2")), filtered_data)
intersection_filtered_positions$V1 <- as.factor(intersection_filtered_positions$V1)

write.table(intersection_filtered_positions, "~/Desktop/Rawdata/chrompos/intersection_of_all_evolvedm.chrompos", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# Output the number of unique positions
cat("Number of positions common to all filtered datasets:", nrow(intersection_filtered_positions), "\n")

### bc cage alone

common_tb20_bc40_positions <- intersect(tb20_positions, bc40_positions)

common_ts20_ts40_filtered <- ts20_filtered[tb20_positions %in% common_tb20_bc40_positions, ]
common_ts20_ts40_filtered <- ts40_filtered[bc40_positions %in% common_tb20_bc40_positions, ]

# Find common positions between ts20, ts40, and base100k
base100k_positions <- paste(base100k_filtered$V1, base100k_filtered$V2, sep = ".")
common_ts20_ts40_base100k_positions <- intersect(common_ts20_ts40_positions, base100k_positions)

# Filter out rows from ts20_filtered, ts40_filtered, and base100k_filtered where the positions match
common_ts20_ts40_base100k_filtered <- ts20_filtered[tb20_positions %in% common_ts20_ts40_base100k_positions, ]
common_ts20_ts40_base100k_filtered <- ts40_filtered[bc40_positions %in% common_ts20_ts40_base100k_positions, ]
common_ts20_ts40_base100k_filtered <- base100k_filtered[base100k_positions %in% common_ts20_ts40_base100k_positions, ]

# Save the results to new files
write.table(common_ts20_ts40_filtered, "common_ts20_ts40_filtered.sync", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(common_ts20_ts40_base100k_filtered, "common_ts20_ts40_base100k_filtered.sync", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Print number of common sites
cat("Number of common sites between ts20_filtered and ts40_filtered:", nrow(common_ts20_ts40_filtered), "\n")
cat("Number of common sites between ts20_filtered, ts40_filtered, and base100k_filtered:", nrow(common_ts20_ts40_base100k_filtered), "\n")
# 
# # Step 7: Plot the pie chart
# ggplot(data, aes(x = "", y = Count, fill = Category)) +
#   geom_bar(stat = "identity", width = 1) +
#   coord_polar("y") +
#   theme_minimal() +
#   labs(title = "Interaction Pie Chart: Variant Positions Between Base100k, BC40, and TB20") +
#   theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
#   scale_fill_manual(values = c("#F1C40F", "#E74C3C", "#8E44AD", "#2ECC71", "#3498DB", "#1ABC9C", "#9B59B6")) +
#   theme(legend.title = element_blank())

