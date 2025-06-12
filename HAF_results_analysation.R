# PLOT ALL HARP RESULTS Base100k, BC, TS, RR crosses in the most efficient way
library(ACER)
library(poolSeq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

snptablecol_2L <- c(
  "2L", "Ref", "Sz5", "Sz6", "Sz8", "Sz9", "Sz11", "Sz15", "Sz22", "Sz24", "Sz27", 
  "Sz30", "Sz33", "Sz39", "Sz40", "Sz45", "Sz47", "Sz52", "Sz53", "Sz54", "Sz56", 
  "Sz62", "Sz68", "Sz75", "Sz82", "Sz83", "Sz87", "Sz90", "Sz91", "Sz101", "Sz102", 
  "Sz106", "Sz110", "Sz120", "Sz124", "Sz125", "Sz129", "Sz130", "Sz137", "Sz138", 
  "Sz141", "Sz142", "Sz143", "Sz146", "Sz149", "Sz153", "Sz160", "Sz166", "Sz169", 
  "Sz175", "Sz182", "Sz185", "Sz186", "Sz191", "Sz193", "Sz196", "Sz200", "Sz201", 
  "Sz202", "Sz205", "Sz208", "Sz209", "Sz213", "Sz214", "Sz217", "Sz218", "Sz220", 
  "Sz225", "Sz232", "Sz236", "Sz238", "Sz241", "Sz244", "Sz249", "Sz254", "Sz257", 
  "Sz259", "Sz262", "Sz271", "Sz275", "Sz281", "Sz282", "Sz287", "Sz288", "Sz289", 
  "DSIM_98", "DSIM_158", "DSIM_159", "DSIM_284", "DSIM_116", "DSIM_286", "DSIM_12"
)
snptablecol_2R <- c(
  "2R", "Ref", "Sz5", "Sz6", "Sz8", "Sz9", "Sz11", "Sz15", "Sz22", "Sz24", "Sz27", 
  "Sz30", "Sz33", "Sz39", "Sz40", "Sz45", "Sz47", "Sz52", "Sz53", "Sz54", "Sz56", 
  "Sz62", "Sz68", "Sz75", "Sz82", "Sz83", "Sz87", "Sz90", "Sz91", "Sz101", "Sz102", 
  "Sz106", "Sz110", "Sz120", "Sz124", "Sz125", "Sz129", "Sz130", "Sz137", "Sz138", 
  "Sz141", "Sz142", "Sz143", "Sz146", "Sz149", "Sz153", "Sz160", "Sz166", "Sz169", 
  "Sz175", "Sz182", "Sz185", "Sz186", "Sz191", "Sz193", "Sz196", "Sz200", "Sz201", 
  "Sz202", "Sz205", "Sz208", "Sz209", "Sz213", "Sz214", "Sz217", "Sz218", "Sz220", 
  "Sz225", "Sz232", "Sz236", "Sz238", "Sz241", "Sz244", "Sz249", "Sz254", "Sz257", 
  "Sz259", "Sz262", "Sz271", "Sz275", "Sz281", "Sz282", "Sz287", "Sz288", "Sz289", 
  "DSIM_98", "DSIM_158", "DSIM_159", "DSIM_284", "DSIM_116", "DSIM_286", "DSIM_12"
)
snptablecol_3L <- c(
  "3L", "Ref", "Sz5", "Sz6", "Sz8", "Sz9", "Sz11", "Sz15", "Sz22", "Sz24", "Sz27", 
  "Sz30", "Sz33", "Sz39", "Sz40", "Sz45", "Sz47", "Sz52", "Sz53", "Sz54", "Sz56", 
  "Sz62", "Sz68", "Sz75", "Sz82", "Sz83", "Sz87", "Sz90", "Sz91", "Sz101", "Sz102", 
  "Sz106", "Sz110", "Sz120", "Sz124", "Sz125", "Sz129", "Sz130", "Sz137", "Sz138", 
  "Sz141", "Sz142", "Sz143", "Sz146", "Sz149", "Sz153", "Sz160", "Sz166", "Sz169", 
  "Sz175", "Sz182", "Sz185", "Sz186", "Sz191", "Sz193", "Sz196", "Sz200", "Sz201", 
  "Sz202", "Sz205", "Sz208", "Sz209", "Sz213", "Sz214", "Sz217", "Sz218", "Sz220", 
  "Sz225", "Sz232", "Sz236", "Sz238", "Sz241", "Sz244", "Sz249", "Sz254", "Sz257", 
  "Sz259", "Sz262", "Sz271", "Sz275", "Sz281", "Sz282", "Sz287", "Sz288", "Sz289", 
  "DSIM_98", "DSIM_158", "DSIM_159", "DSIM_284", "DSIM_116", "DSIM_286", "DSIM_12"
)
snptablecol_3R <- c(
  "3R", "Ref", "Sz5", "Sz6", "Sz8", "Sz9", "Sz11", "Sz15", "Sz22", "Sz24", "Sz27", 
  "Sz30", "Sz33", "Sz39", "Sz40", "Sz45", "Sz47", "Sz52", "Sz53", "Sz54", "Sz56", 
  "Sz62", "Sz68", "Sz75", "Sz82", "Sz83", "Sz87", "Sz90", "Sz91", "Sz101", "Sz102", 
  "Sz106", "Sz110", "Sz120", "Sz124", "Sz125", "Sz129", "Sz130", "Sz137", "Sz138", 
  "Sz141", "Sz142", "Sz143", "Sz146", "Sz149", "Sz153", "Sz160", "Sz166", "Sz169", 
  "Sz175", "Sz182", "Sz185", "Sz186", "Sz191", "Sz193", "Sz196", "Sz200", "Sz201", 
  "Sz202", "Sz205", "Sz208", "Sz209", "Sz213", "Sz214", "Sz217", "Sz218", "Sz220", 
  "Sz225", "Sz232", "Sz236", "Sz238", "Sz241", "Sz244", "Sz249", "Sz254", "Sz257", 
  "Sz259", "Sz262", "Sz271", "Sz275", "Sz281", "Sz282", "Sz287", "Sz288", "Sz289", 
  "DSIM_98", "DSIM_158", "DSIM_159", "DSIM_284", "DSIM_116", "DSIM_286", "DSIM_12"
)
snptablecol_X <- c(
  "X", "Ref", "Sz5", "Sz6", "Sz8", "Sz9", "Sz11", "Sz15", "Sz22", "Sz24", "Sz27", 
  "Sz30", "Sz33", "Sz39", "Sz40", "Sz45", "Sz47", "Sz52", "Sz53", "Sz54", "Sz56", 
  "Sz62", "Sz68", "Sz75", "Sz82", "Sz83", "Sz87", "Sz90", "Sz91", "Sz101", "Sz102", 
  "Sz106", "Sz110", "Sz120", "Sz124", "Sz125", "Sz129", "Sz130", "Sz137", "Sz138", 
  "Sz141", "Sz142", "Sz143", "Sz146", "Sz149", "Sz153", "Sz160", "Sz166", "Sz169", 
  "Sz175", "Sz182", "Sz185", "Sz186", "Sz191", "Sz193", "Sz196", "Sz200", "Sz201", 
  "Sz202", "Sz205", "Sz208", "Sz209", "Sz213", "Sz214", "Sz217", "Sz218", "Sz220", 
  "Sz225", "Sz232", "Sz236", "Sz238", "Sz241", "Sz244", "Sz249", "Sz254", "Sz257", 
  "Sz259", "Sz262", "Sz271", "Sz275", "Sz281", "Sz282", "Sz287", "Sz288", "Sz289", 
  "DSIM_98", "DSIM_158", "DSIM_159", "DSIM_284", "DSIM_116", "DSIM_286", "DSIM_12"
)

snptablecol_freq_2L <- snptablecol_2L[-(1:2)]
snptablecol_freq_2R <- snptablecol_2R[-(1:2)]
snptablecol_freq_3L <- snptablecol_3L[-(1:2)]
snptablecol_freq_3R <- snptablecol_3R[-(1:2)]
snptablecol_freq_X <- snptablecol_X[-(1:2)]

threecolumns <- c("window", "start", "end")
snptablecol_freq_2L <- c(threecolumns, snptablecol_freq_2L)
snptablecol_freq_2R <- c(threecolumns, snptablecol_freq_2R)
snptablecol_freq_3L <- c(threecolumns, snptablecol_freq_3L)
snptablecol_freq_3R <- c(threecolumns, snptablecol_freq_3R)
snptablecol_freq_X <- c(threecolumns, snptablecol_freq_X)

chromosomes <- c("2L", "2R", "3L", "3R", "X")

#hapfreqbasepath <- "/Users/siva/Desktop/Rawdata/harp_results/RR33/400/%s/RR_33_MKDN240015467-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.%s.freqs"
#for (chrom in chromosomes) {
#  hapfreqfilepath <- sprintf(hapfreqbasepath, chrom, chrom)
#  var_name <- paste0("BCR4F40_", chrom) #change the name here and the hapfreqbasepath
#  assign(var_name, read.table(hapfreqfilepath, header = FALSE, sep = " "))
#}

base_paths <- list(
  "Base16k_400kb" = "/Users/siva/Desktop/Rawdata/harp_results/Base16k/400/%s/Base_16k_MKDN240015458-1A_227Y7WLT4_L1_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "Base100_400kb" = "/Users/siva/Desktop/Rawdata/harp_results/Base100k/400/%s/Base_100k_MKDN240015457-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "BCR4F40_115kb" = "/Users/siva/Desktop/Rawdata/harp_results/BCR4F40/115/%s/BC_R4_F40_EKDN220048411-1A_HLTJ3DSX5_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "RR33_400kb" = "/Users/siva/Desktop/Rawdata/harp_results/RR33/400/%s/RR_33_MKDN240015467-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.%s.freqs",
  "RR34_400kb" = "/Users/siva/Desktop/Rawdata/harp_results/RR34/400/%s/RR_34_MKDN240015468-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.%s.freqs",
  "RR90_400kb" = "/Users/siva/Desktop/Rawdata/harp_results/RR90/400/%s/RR_90_MKDN240015469-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.%s.freqs",
  "RR91_400kb" = "/Users/siva/Desktop/Rawdata/harp_results/RR91/400/%s/RR_91_MKDN240015470-1A_227Y7WLT4_L5_trimmedG_sorted_rmdup_filtered_sc.bam.%s.freqs",
  "TBR4F20_180kb" = "/Users/siva/Desktop/Rawdata/harp_results/TBR4F20/180/%s/TB_R4_F20_MKDN240015459-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F20_180kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F20/180/%s/TS_R4_F20_MKDN240015460-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F40_115kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F40/115/%s/TS_R4_F40_EKDN220048410-1A_HLTJ3DSX5_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F60_85kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F60/85/%s/TS_R4_F60_MKDN240015461-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F80_67kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F80/67/%s/TS_R4_F80_MKDN240015462-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs"
)
for (folder_name in names(base_paths)) {
  path_template <- base_paths[[folder_name]]
  for (chrom in chromosomes) {
    file_path <- sprintf(path_template, chrom, chrom)
    var_name <- paste0(folder_name, "_", chrom)
    assign(var_name, read.table(file_path, header = FALSE, sep = " "))
  }
}

# After running this script, you'll have data frames like Base16k_2L, Base100k_3R, we change the colnames
for (folder_name in names(base_paths)) {
  for (chrom in chromosomes) {
    var_name <- paste0(folder_name, "_", chrom) 
    temp_df <- get(var_name)
    temp_df <- temp_df[, -ncol(temp_df)]
    assign(var_name, temp_df)
    snp_col_name <- paste0("snptablecol_freq_", chrom)
    colnames(temp_df) <- get(snp_col_name)  
    assign(var_name, temp_df)
  }
}


for (folder_name in names(base_paths)) {
  for (chrom in chromosomes) {
    var_name <- paste0(folder_name, "_", chrom)  
    if (exists(var_name)) {
      temp_df <- get(var_name)
      melted_df <- melt(temp_df, id.vars = c("window", "start", "end"), 
                        variable.name = "haplotype", 
                        value.name = "frequency")
      
      heatmap_plot <- ggplot(melted_df, aes(x = start, y = haplotype, fill = frequency)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(title = paste("Haplotype frequency ", var_name),
           x = "Genomic position",
           y = "Founder haplotype") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
    
    ggsave(paste0("~/Desktop/harp_results/",var_name,".png"), plot = heatmap_plot, width = 12, height = 10)
    }
  }
}



### 100kb regions:
base_paths_1 <- list(
  "Base16k_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/Base16k/100/%s/Base_16k_MKDN240015458-1A_227Y7WLT4_L1_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "Base100k_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/Base100k/100/%s/Base_100k_MKDN240015457-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "BCR4F40_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/BCR4F40/100/%s/BC_R4_F40_EKDN220048411-1A_HLTJ3DSX5_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TBR4F20_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/TBR4F20/100/%s/TB_R4_F20_MKDN240015459-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F20_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F20/100/%s/TS_R4_F20_MKDN240015460-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F40_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F40/100/%s/TS_R4_F40_EKDN220048410-1A_HLTJ3DSX5_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F60_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F60/100/%s/TS_R4_F60_MKDN240015461-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs",
  "TSR4F80_100kb" = "/Users/siva/Desktop/Rawdata/harp_results/TSR4F80/100/%s/TS_R4_F80_MKDN240015462-1A_227Y7WLT4_L2_trimmedG_sorted_rmdup_filtered_recal_reads_sc.bam.%s.freqs"
)
for (folder_name in names(base_paths_1)) {
  path_template <- base_paths_1[[folder_name]]
  for (chrom in chromosomes) {
    file_path <- sprintf(path_template, chrom, chrom)
    var_name <- paste0(folder_name, "_", chrom)
    assign(var_name, read.table(file_path, header = FALSE, sep = " "))
  }
}


# After running this script, you'll have data frames like Base16k_2L, Base100k_3R, we change the colnames

for (folder_name in names(base_paths_1)) {
  for (chrom in chromosomes) {
    var_name <- paste0(folder_name, "_", chrom) 
    temp_df <- get(var_name)
    temp_df <- temp_df[, -ncol(temp_df)]
    assign(var_name, temp_df)
    snp_col_name <- paste0("snptablecol_freq_", chrom)
    colnames(temp_df) <- get(snp_col_name)  
    assign(var_name, temp_df)
  }
}
for (folder_name in names(base_paths_1)) {
  for (chrom in chromosomes) {
    var_name <- paste0(folder_name, "_", chrom)  
    if (exists(var_name)) {
      temp_df <- get(var_name)
      melted_df <- melt(temp_df, id.vars = c("window", "start", "end"), 
                        variable.name = "haplotype", 
                        value.name = "frequency")
      
      heatmap_plot <- ggplot(melted_df, aes(x = start, y = haplotype, fill = frequency)) +
        geom_tile() +
        scale_fill_gradient(low = "white", high = "blue", limits = c(0, 1)) +
        labs(title = paste("Haplotype frequency ", var_name),
             x = "Genomic position",
             y = "Founder haplotype") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8))
      
      ggsave(paste0("~/Desktop/harp_results/100kb_equal/",var_name,".png"), plot = heatmap_plot, width = 12, height = 10)
    }
  }
}
# base_paths and base_paths_1, shannon index calculation:

for (folder_name in names(base_paths_1)) {
  for (chrom in chromosomes) {
    var_name <- paste0(folder_name, "_", chrom)  
    if (exists(var_name)) {
      temp_df <- get(var_name)
      haplotype_freqs <- temp_df[, -c(1,2,3)]
      shannon_index <- function(freqs) {
        freqs <- freqs / sum(freqs)
        freqs <- freqs[freqs > 0]
        -sum(freqs * log(freqs))
      }
      
      shannon_indices <- apply(haplotype_freqs, 1, shannon_index)
      #print(shannon_indices)
      result <- data.frame(Chromosome = temp_df[,1], Position = temp_df[,2], Position_end = temp_df[,3], Shannon_Index = shannon_indices, Source = var_name)
      var_name1 <- paste0(var_name, "_shannon_index")  
      assign(var_name1, result)
    }
  }
}

### shannon diversity index plots 

combined_df <- bind_rows(Base100k_100kb_X_shannon_index,TBR4F20_100kb_X_shannon_index, BCR4F40_100kb_X_shannon_index)
combined_df <- bind_rows(Base100_400kb_2L_shannon_index,Base100k_100kb_2L_shannon_index,TBR4F20_100kb_2L_shannon_index, BCR4F40_100kb_2L_shannon_index, TBR4F20_180kb_2L_shannon_index, BCR4F40_115kb_2L_shannon_index)

ggplot(combined_df, aes(x = Position, y = Shannon_Index, color = Source)) +
  geom_line() +  
  geom_point(size = 0) +
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  labs(
    title = "",
    x = "Genomic Position",
    y = "Shannon Diversity Index"
  ) +
  facet_wrap(~ Chromosome, scales = "free_x") +  # Free x-axis per chromosome
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 1)
  ) + theme_minimal() + scale_y_continuous(breaks = seq(0, 4.5, by = 0.5), limits = c(0, 4.5))

combined_df <- bind_rows(Base100k_100kb_2L_shannon_index, BCR4F40_100kb_2L_shannon_index, TBR4F20_100kb_2L_shannon_index, Base16k_100kb_2L_shannon_index,TSR4F20_100kb_2L_shannon_index, TSR4F40_100kb_2L_shannon_index, TSR4F60_100kb_2L_shannon_index,TSR4F80_100kb_2L_shannon_index)
combined_df <- bind_rows(BCR4F40_100kb_2L_shannon_index, TSR4F40_100kb_2L_shannon_index)

combined_df <- bind_rows(Base16k_100kb_X_shannon_index,
                         TSR4F20_100kb_X_shannon_index, 
                         TSR4F40_100kb_X_shannon_index, 
                         TSR4F60_100kb_X_shannon_index,
                         TSR4F80_100kb_X_shannon_index)

ggplot(combined_df, aes(x = Position, y = Shannon_Index, color = Source)) +
  geom_line() +  
  geom_point(size = 0) +
  labs(
    title = "",
    x = "Genomic Position",
    y = "Shannon Diversity Index"
  ) +
  facet_wrap(~ Chromosome, scales = "free_x") +  # Free x-axis per chromosome
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(hjust = 1)
  ) + theme_minimal() + scale_y_continuous(breaks = seq(0, 4.5, by = 0.5), limits = c(0, 4.5))

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

threshold <- 0.8  ## very important

for (folder_name in names(base_paths_1)) {
  for (chrom in chromosomes) {
    var_name <- paste0(folder_name, "_", chrom)  
    if (exists(var_name)) {
      temp_df <- get(var_name)
      var_name_values <- apply(temp_df[, 4:93], 1, haplotype_count_distribution_function, threshold)
      var_name_values <- t(data.frame(var_name_values)) 
      colnames(var_name_values) <- c("num_values_added", "cumulative_sum")
      var_name_values_df <- cbind(temp_df[, 1:3], var_name_values, Source = var_name)
      result_var_name <- paste0(var_name, "_hapcountdf")
      assign(result_var_name, var_name_values_df)
    }
  }
}



#### Plots

combined_df <- rbind(Base16k_100kb_X_hapcountdf, 
                     TSR4F20_100kb_X_hapcountdf, 
                     TSR4F40_100kb_X_hapcountdf, 
                     TSR4F60_100kb_X_hapcountdf, 
                     TSR4F80_100kb_X_hapcountdf)
combined_df <- rbind(Base100k_100kb_X_hapcountdf,
                     TBR4F20_100kb_X_hapcountdf, 
                     BCR4F40_100kb_X_hapcountdf)
ggplot(combined_df, aes(x = start, y = num_values_added, color = Source)) +
  geom_line() + 
  #geom_point(size = 0) +
  labs(title = "",
       x = "Genomic Position",
       y = "No. of Haplotypes") +
  # scale_color_manual(values = c("Base100k" = "blue", "BC_3R" = "red")) +
theme_minimal() + scale_y_continuous(breaks = seq(0, 50, by = 5), limits = c(0, 50))

ggplot(combined_df, aes(x = num_values_added, color = Source)) +
  #geom_histogram(binwidth = 1,alpha = 0.7) +
  geom_histogram(binwidth = 1) +
  
  labs(title = "Distribution of Number of Values Added",
       x = "Number of Values Added",
       y = "Frequency") +
  theme_minimal()

ggplot(combined_df, aes(x = num_values_added, color = Source)) +
  geom_density(fill = "white", alpha = 0.7) +
  labs(title = "",
       x = "Number of Haplotypes",
       y = "Density") +
  theme_minimal()

# Fit a linear model with start position and Source as factors
regression_model <- lm(num_values_added ~ Source, data = combined_df)
summary(regression_model)
ggplot(combined_df, aes(x = start, y = num_values_added, color = Source)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
  labs(title = "",
       x = "Genomic Position",
       y = "No. of Haplotypes") +
  theme_minimal() + 
  scale_y_continuous(breaks = seq(0, 60, by = 5), limits = c(0, 60))

### sfs plot for haplotype 
# calculate_sfs <- function(haplotype_freqs) {
#   freq_vector <- unlist(haplotype_freqs)
#   freq_vector <- freq_vector[freq_vector > 0 & !is.na(freq_vector)]
#   freq_vector <- round(freq_vector, 2)
#   sfs <- table(freq_vector)
#   sfs_df <- as.data.frame(sfs)
#   colnames(sfs_df) <- c("Frequency", "Count")
#   return(sfs_df)
# }
# 
# for (folder_name in names(base_paths_1)) {
#   for (chrom in chromosomes) {
#     var_name <- paste0(folder_name, "_", chrom)
#     if (exists(var_name)) {
#       temp_df <- get(var_name)
#       var_name_sfs <- paste0(var_name, "_sfs")
#       haplotype_freqs <- temp_df[, 4:93]
#       sfs <- calculate_sfs(haplotype_freqs)
#       sfs$Source <- var_name
#       assign(var_name_sfs, sfs)
#     }
#   }
# }
# combined_sfs_df <- rbind(Base16k_100kb_3L_sfs, TSR4F20_100kb_3L_sfs, TSR4F40_100kb_3L_sfs, TSR4F60_100kb_3L_sfs, TSR4F80_100kb_3L_sfs)
# ###
