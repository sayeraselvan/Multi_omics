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

BC_Sync <- read.sync(file="/Users/siva/Desktop/Rawdata/BC_base_40gen_R4_snp1_filtered.sync", gen=c(0, 40), repl=c(1, 1))
SC_Sync <- read.sync(file="/Users/siva/Desktop/Rawdata/SC_base_40gen_R4_snp1_filtered.sync", gen=c(0, 40), repl=c(1, 1))

allele_BC <- alleles(sync = BC_Sync) 
allele_SC <- alleles(sync = SC_Sync) 

#colRanks(matrix, ties.method = "average")

# file management
af_BC_gen0 <- af(sync = BC_Sync, repl = 1, gen = 0)
af_BC_gen40 <- af(sync = BC_Sync, repl = 1, gen = 40)

af_SC_gen0 <- af(sync = SC_Sync, repl = 1, gen = 0)
af_SC_gen40 <- af(sync = SC_Sync, repl = 1, gen = 40)

cov_BC_gen0 <- coverage(sync = BC_Sync, repl = 1, gen = 0)
cov_BC_gen40 <- coverage(sync = BC_Sync, repl = 1, gen = 40)

cov_SC_gen0 <- coverage(sync = SC_Sync, repl = 1, gen = 0)
cov_SC_gen40 <- coverage(sync = SC_Sync, repl = 1, gen = 40)

file_manage <- function(af_gen0) {
  df_gen0 <- data.frame(chrom_pos = names(af_gen0), af = af_gen0)
  df_gen0 <- df_gen0 %>% separate(chrom_pos, into = c("chrom", "pos"), sep = "\\.")
  rownames(df_gen0) <- NULL
  return(df_gen0)
}

file_manage_cov <- function(af_gen0) {
  df_gen0 <- data.frame(chrom_pos = names(af_gen0), cov = af_gen0)
  df_gen0 <- df_gen0 %>% separate(chrom_pos, into = c("chrom", "pos"), sep = "\\.")
  rownames(df_gen0) <- NULL
  return(df_gen0)
}

af_BC_gen0_df <- file_manage(af_BC_gen0)
af_BC_gen40_df <- file_manage(af_BC_gen40)
af_SC_gen0_df <- file_manage(af_SC_gen0)
af_SC_gen40_df <- file_manage(af_SC_gen40)
cov_BC_gen0_df <- file_manage_cov(cov_BC_gen0)
cov_BC_gen40_df <- file_manage_cov(cov_BC_gen40)
cov_SC_gen0_df <- file_manage_cov(cov_SC_gen0)
cov_SC_gen40_df  <- file_manage_cov(cov_SC_gen40)



estimateNe(p0=AF_BC_gen0_2L, pt=AF_BC_gen40_2L, cov0=Cov_BC_gen0_2L, covt=Cov_BC_gen40_2L, t=40, Ncensus=100000, poolSize=c(250, 250), truncAF = NA, method ="P.planI")
?estimateNe

checkSNP(p0=AF_BC_gen0_2L, pt=AF_BC_gen40_2L, cov0=Cov_BC_gen0_2L, covt=Cov_BC_gen40_2L, truncAF = NA)

af_BC_gen0_ex <- af_BC_gen0[af_BC_gen0 > 0.01 & af_BC_gen0 < 0.99]
af_BC_gen40_ex <- af_BC_gen40[af_BC_gen40 > 0.01 & af_BC_gen40 < 0.99]

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
  print(ne_estimates)
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


af_BC_gen0_2L <- subset(af_BC_gen0_df, chrom == "2L")
af_BC_gen40_2L <- subset(af_BC_gen40_df, chrom == "2L")
cov_BC_gen0_2L <- subset(cov_BC_gen0_df, chrom == "2L")
cov_BC_gen40_2L <- subset(cov_BC_gen40_df, chrom == "2L")

af_BC_gen0_40_2L <- merge(af_BC_gen0_2L, af_BC_gen40_2L, by = c("chrom", "pos"), suffixes = c("_gen0", "_gen40"))

af_BC_gen0_40_2L$time_0_r1 <- af_BC_gen0_40_2L$af_gen0
af_BC_gen0_40_2L$time_40_r1 <- af_BC_gen0_40_2L$af_gen40

cov_BC_gen0_40_2L <- merge(cov_BC_gen0_2L, cov_BC_gen40_2L, by = c("chrom", "pos"), suffixes = c("_gen0", "_gen40"))

cov_BC_gen0_40_2L$time_0_r1 <- cov_BC_gen0_40_2L$cov_gen0
cov_BC_gen0_40_2L$time_40_r1 <- cov_BC_gen0_40_2L$cov_gen40

colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "chrom"] <- "CHROM"
colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "pos"] <- "POS"

colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "chrom"] <- "CHROM"
colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "pos"] <- "POS"

colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "0_r 1"] <- "0_r1"
colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "40_r 1"] <- "40_r1"

colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "time_0_r1"] <- "0_r1"
colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "time_40_r1"] <- "40_r1"

colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "time_0_r1"] <- "0_r1"
colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "time_40_r1"] <- "40_r1"

cov_BC_gen0_40_2L$CHROM <- as.character(cov_BC_gen0_40_2L$CHROM)
af_BC_gen0_40_2L$CHROM <- as.character(af_BC_gen0_40_2L$CHROM)

#Ne.estimator(DP.SO, AF.SO, replicates = 1:4, nb_SNPs = 1000, nb_rounds = 1000,
#            time = c(1,10), poolSize = c(500,500), nb_replicates = 4, census = 1250,
#             chrom = "2L")

#Ne.estimator(DP.OS, AF.OS, replicates = 1:4, nb_SNPs = 1000, nb_rounds = 1000,
#             time = c(1,10), poolSize = c(500,500), nb_replicates = 4, census = 1250,
#             chrom = "2L")


Ne.estimator(dp = cov_BC_gen0_40_2L_NE_NE, af = af_BC_gen0_40_2L_NE_NE, replicates = 1:1, nb_SNPs = 15000, nb_rounds = 10,
             time = c(0,40), poolSize = c(250,250), nb_replicates = 1, census = 100000,
             chrom = "2L")

#head(cov_BC_gen0_40_2L)

zero_rows <- cov_BC_gen0_40_2L$`40_r1` == 0
zero_rows <- cov_BC_gen0_40_2L$`0_r1` == 0

cov_BC_gen0_40_2L <- cov_BC_gen0_40_2L[!zero_rows, ]
af_BC_gen0_40_2L <- af_BC_gen0_40_2L[!zero_rows, ]

maf_rows <- af_BC_gen0_40_2L$`0_r1` > 0.05 & af_BC_gen0_40_2L$`0_r1` < 0.95
maf_rows_1 <- af_BC_gen0_40_2L_NE$`40_r1` > 0.05 & af_BC_gen0_40_2L_NE$`40_r1` < 0.95

cov_BC_gen0_40_2L_NE <- cov_BC_gen0_40_2L[maf_rows, ]
af_BC_gen0_40_2L_NE <- af_BC_gen0_40_2L[maf_rows, ]

cov_BC_gen0_40_2L_NE_NE <- cov_BC_gen0_40_2L_NE[maf_rows_1, ]
af_BC_gen0_40_2L_NE_NE <- af_BC_gen0_40_2L_NE[maf_rows_1, ]

###### SC::

af_BC_gen0_2L <- subset(af_SC_gen0_df, chrom == "X")
af_BC_gen40_2L <- subset(af_SC_gen40_df, chrom == "X")
cov_BC_gen0_2L <- subset(cov_SC_gen0_df, chrom == "X")
cov_BC_gen40_2L <- subset(cov_SC_gen40_df, chrom == "X")

af_BC_gen0_40_2L <- merge(af_BC_gen0_2L, af_BC_gen40_2L, by = c("chrom", "pos"), suffixes = c("_gen0", "_gen40"))

af_BC_gen0_40_2L$time_0_r1 <- af_BC_gen0_40_2L$af_gen0
af_BC_gen0_40_2L$time_40_r1 <- af_BC_gen0_40_2L$af_gen40

cov_BC_gen0_40_2L <- merge(cov_BC_gen0_2L, cov_BC_gen40_2L, by = c("chrom", "pos"), suffixes = c("_gen0", "_gen40"))

cov_BC_gen0_40_2L$time_0_r1 <- cov_BC_gen0_40_2L$cov_gen0
cov_BC_gen0_40_2L$time_40_r1 <- cov_BC_gen0_40_2L$cov_gen40

colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "chrom"] <- "CHROM"
colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "pos"] <- "POS"

colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "chrom"] <- "CHROM"
colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "pos"] <- "POS"

colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "0_r 1"] <- "0_r1"
colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "40_r 1"] <- "40_r1"

colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "time_0_r1"] <- "0_r1"
colnames(af_BC_gen0_40_2L)[colnames(af_BC_gen0_40_2L) == "time_40_r1"] <- "40_r1"

colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "time_0_r1"] <- "0_r1"
colnames(cov_BC_gen0_40_2L)[colnames(cov_BC_gen0_40_2L) == "time_40_r1"] <- "40_r1"

cov_BC_gen0_40_2L$CHROM <- as.character(cov_BC_gen0_40_2L$CHROM)
af_BC_gen0_40_2L$CHROM <- as.character(af_BC_gen0_40_2L$CHROM)

#Ne.estimator(DP.SO, AF.SO, replicates = 1:4, nb_SNPs = 1000, nb_rounds = 1000,
#            time = c(1,10), poolSize = c(500,500), nb_replicates = 4, census = 1250,
#             chrom = "2L")

#Ne.estimator(DP.OS, AF.OS, replicates = 1:4, nb_SNPs = 1000, nb_rounds = 1000,
#             time = c(1,10), poolSize = c(500,500), nb_replicates = 4, census = 1250,
#             chrom = "2L")


Ne.estimator(dp = cov_BC_gen0_40_2L_NE_NE, af = af_BC_gen0_40_2L_NE_NE, replicates = 1:1, nb_SNPs = 5000, nb_rounds = 10,
             time = c(0,40), poolSize = c(250,250), nb_replicates = 1, census = 800,
             chrom = "X")

#head(cov_BC_gen0_40_2L)

zero_rows <- cov_BC_gen0_40_2L$`40_r1` == 0
cov_BC_gen0_40_2L <- cov_BC_gen0_40_2L[!zero_rows, ]
af_BC_gen0_40_2L <- af_BC_gen0_40_2L[!zero_rows, ]

zero_rows <- cov_BC_gen0_40_2L$`0_r1` == 0

cov_BC_gen0_40_2L <- cov_BC_gen0_40_2L[!zero_rows, ]
af_BC_gen0_40_2L <- af_BC_gen0_40_2L[!zero_rows, ]

maf_rows <- af_BC_gen0_40_2L$`0_r1` > 0.1 & af_BC_gen0_40_2L$`0_r1` < 0.9
maf_rows_1 <- af_BC_gen0_40_2L_NE$`40_r1` > 0.1 & af_BC_gen0_40_2L_NE$`40_r1` < 0.9

cov_BC_gen0_40_2L_NE <- cov_BC_gen0_40_2L[maf_rows, ]
af_BC_gen0_40_2L_NE <- af_BC_gen0_40_2L[maf_rows, ]

cov_BC_gen0_40_2L_NE_NE <- cov_BC_gen0_40_2L_NE[maf_rows_1, ]
af_BC_gen0_40_2L_NE_NE <- af_BC_gen0_40_2L_NE[maf_rows_1, ]

af_traj_X_SC_chrom <- as.character(af_BC_gen0_40_2L_NE_NE$CHROM)
af_traj_X_SC_pos <- as.numeric(af_BC_gen0_40_2L_NE_NE$POS)
?af.traj


af.traj_output <- af.traj(sync = SC_Sync, chr = af_traj_X_SC_chrom, pos = af_traj_X_SC_pos, repl = 1)


df <- data.frame(
  position = rownames(af.traj_output),
  F0 = af.traj_output[, 1],
  F40 = af.traj_output[, 2]
)
rownames(df) = NULL

df$absolute_diff <- df$F40 - df$F0

df1 <- df[df$absolute_diff > 0.4, ]    

library(ggplot2)

# Assuming df is your data
df_long <- reshape2::melt(df1, id.vars = "position", variable.name = "Time", value.name = "Allele_Frequency")

# Plot
ggplot(df_long, aes(x = Time, y = Allele_Frequency, group = position)) +
  geom_line() +
  theme_minimal() +
  labs(x = "Allele Frequency", y = "Time", title = "Allele Frequency Trajectory")


#### Fisher exact test for BC and SC:  
BCfisher <- read.table("~/Desktop/Rawdata/BC_base100k_gen20_40_rep4_withref_no_na.fet")
BCfisher <- read.table("~/Desktop/Rawdata/BC_base100k_gen20_40_rep4_withref_mincount_2_min20_max_1%_no_na.fet")
SCfisher <- read.table("/Users/siva/Desktop/Rawdata/SC_base16k_gen20_40_60_80_rep4_withref.fet")


library(ggplot2)

chromosome_of_interest <- "2L"  
start_position <- 1 
end_position <- 20000000
filtered_data <- BCfisher[BCfisher$V1 == chromosome_of_interest & 
                            BCfisher$V2 >= start_position & 
                            BCfisher$V2 <= end_position & 
                            BCfisher$V3, ]
#filtered_data <- BCfisher[BCfisher$V1 == chromosome_of_interest, ]

filtered_data$FisherPValue <- as.numeric(sub(".*=", "", filtered_data$V6))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$FDR <- p.adjust(filtered_data$PValue, method = "fdr")

#fdrs<-p.adjust(pvalues, method="BH")

# Ensure V6 is numeric

ggplot(filtered_data, aes(x = V2, y = FisherPValue)) +
  geom_point(alpha = 0.2, size = 0.1, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggplot(filtered_data, aes(x = V2, y = -log10(FDR))) +
  geom_point(alpha = 0.2, size = 0, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(sliding_results, aes(x = Start, y = Log10_GeoMeanP)) +
  geom_point(alpha = 0.2, size = 0.5, color = "blue") +
  labs(
    title = paste("Manhattan Plot for Chromosome", chromosome_of_interest),
    subtitle = paste("Positions", start_position, "to", end_position),
    x = "Genomic Position",
    y = "-log10(P-value)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
library("EnvStats")
chr1<- BCfisher[BCfisher$V1 == "2L",]
# Chr1_ave <- aggregate(p~FeatureID, data=chr1, FUN=function(x) c(mean=geoMean(x)))
# Chr1_ave$CHR <- rep(1,nrow(Chr1_ave))
# Chr1_ave$CHR2 <- rep("Chromosome_1",nrow(Chr1_ave))

window_size <- 10000
chr1$PValue <- 10^(-as.numeric(sub(".*=", "", chr1$V6)))
filtered_data$PValue <- 10^(-as.numeric(sub(".*=", "", filtered_data$V6)))
filtered_data$Window <- floor(filtered_data$V2 / window_size)

chr1$Window <- floor(chr1$V2 / window_size)

library(dplyr)
window_summary <- filtered_data %>%
  group_by(V1, Window) %>%
  summarise(
    Start = min(V2),  
    End = max(V2),
    GeometricMeanP = exp(mean(log(PValue)))  # Geometric mean of p-values
  ) %>%
  mutate(Log10_GeoMeanP = -log10(GeometricMeanP))  # Convert back to -log10

# View the summary
print(window_summary)

window_size <- 100
step_size <- 50

sliding_window <- function(data, window_size, step_size) {
  results <- data.frame(
    Start = numeric(),
    End = numeric(),
    Log10_GeoMeanP = numeric()
  )
  for (start in seq(min(data$V2), max(data$V2) - window_size, by = step_size)) {
    end <- start + window_size
    window_data <- data %>% filter(V2 >= start & V2 < end)
    if (nrow(window_data) > 0) {
      geometric_mean_p <- exp(mean(log(window_data$FDR)))
      results <- rbind(
        results,
        data.frame(
          Start = start,
          End = end,
          Log10_GeoMeanP = -log10(geometric_mean_p)
        )
      )
    }
  }
  
  return(results)
}

sliding_results <- sliding_window(filtered_data, window_size, step_size)

print(sliding_results)








######## histogram R insert size::
library(ggplot2)

input_file <- "/Users/siva/Desktop/Rawdata/pdf/base100k.txt"
data <- read.table(input_file, header = TRUE, sep = "\t")
file_base_name <- tools::file_path_sans_ext(basename(input_file))

data$cumulative_fraction <- cumsum(data$All_Reads.fr_count) / sum(data$All_Reads.fr_count)

ggplot(data) +
  # Bar plot for read count
  geom_bar(aes(x = insert_size, y = All_Reads.fr_count), stat = "identity", fill = "red", alpha = 0.5) +
  # Line plot for cumulative fraction
  geom_line(aes(x = insert_size, y = cumulative_fraction * max(All_Reads.fr_count)), 
            color = "black", linetype = "dashed", linewidth = 1) +  # Updated `linewidth`
  # Scale for cumulative fraction (secondary y-axis)
  scale_y_continuous(
    name = "Read Count",
    sec.axis = sec_axis(~ . / max(data$All_Reads.fr_count), 
                        name = "Cumulative Fraction of Reads > Insert Size", 
                        labels = scales::percent),
    limits = c(0, 400000)  # Set y-axis range from 0 to 350,000
  ) +
  # Axis labels and title
  labs(
    title = paste("Insert Size Histogram", file_base_name),
    x = "Insert Size",
    y = "Read Count"
  ) +
  # Set x-axis limits
  xlim(0, 800) +
  # Customize the theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title.y.right = element_text(color = "black")  # Customize secondary axis
  ) 


# Plot 2: Histogram with cumulative fraction
# ggplot(data) +
#   # Bar plot for read count
#   geom_bar(aes(x = insert_size, y = All_Reads.fr_count), stat = "identity", fill = "red", alpha = 0.5) +
#   # Line plot for cumulative fraction
#   geom_line(aes(x = insert_size, y = cumulative_fraction * max(All_Reads.fr_count)), 
#             color = "black", linetype = "dashed", linewidth = 1) +  # Updated `linewidth`
#   # Scale for cumulative fraction (secondary y-axis)
#   scale_y_continuous(
#     name = "Read Count",
#     sec.axis = sec_axis(~ . / max(data$All_Reads.fr_count), 
#                         name = "Cumulative Fraction of Reads > Insert Size", 
#                         labels = scales::percent)
#   ) +
#   # Axis labels and title
#   labs(
#     title = paste("Insert Size Histogram", file_base_name),
#     x = "Insert Size",
#     y = "Read Count"
#   ) +
#   # Set x-axis limits
#   xlim(0, 800) +
#   # Customize the theme
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 14),
#     axis.text = element_text(size = 12),
#     axis.title.y.right = element_text(color = "black")  # Customize secondary axis
#   ) 










