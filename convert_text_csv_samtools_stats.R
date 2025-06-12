
# Load necessary library
library(dplyr)

input_file <- "/Users/siva/Desktop/Rawdata/cluster_reports/concat_bam_stats.txt"  # Change to your file path
output_file <- "/Users/siva/Desktop/Rawdata/signor_stats.csv"  # Change to your desired output path

# Read the text file
lines <- readLines(input_file)

stats_list <- list()
current_file <- NULL

# Loop through each line to extract statistics
for (line in lines) {
  # Check for the file header
  if (grepl("# Stats for", line)) {
    if (!is.null(current_file)) {
      stats_list[[current_file]] <- current_stats
    }
    current_file <- sub("# Stats for ", "", line)
    current_file <- trimws(current_file)
    current_stats <- list()
  } else if (grepl(":", line)) {
    # Extract the statistic name and value, ignoring comments
    parts <- strsplit(line, ":")[[1]]
    stat_name <- trimws(parts[1])
    
    # Remove any comments from the value part
    value_part <- trimws(parts[2])
    value_numeric <- as.numeric(sub("#.*", "", value_part))  # Remove comments and convert to numeric
    
    if (!is.na(value_numeric)) {
      current_stats[[stat_name]] <- value_numeric
    }
  }
}



#### this is for signor paper stats

for (line in lines) {
  # Check for the file header
  if (grepl("Processing BAM file:", line)) {
    if (!is.null(current_file)) {
      stats_list[[current_file]] <- current_stats
    }
    current_file <- sub("Processing BAM file: ", "", line)
    current_file <- trimws(current_file)
    current_stats <- list()
  } else if (grepl(":", line)) {
    # Extract the statistic name and value, ignoring comments
    parts <- strsplit(line, ":")[[1]]
    stat_name <- trimws(parts[1])
    
    # Remove any comments from the value part
    value_part <- trimws(parts[2])
    value_numeric <- as.numeric(sub("#.*", "", value_part))  # Remove comments and convert to numeric
    
    if (!is.na(value_numeric)) {
      current_stats[[stat_name]] <- value_numeric
    }
  }
}
###
# Don't forget to save the last file's statistics
if (!is.null(current_file)) {
  stats_list[[current_file]] <- current_stats
}

# Convert list to data frame
stats_df <- bind_rows(lapply(stats_list, as.data.frame), .id = "File")
stats_df$File <- sub("/home/vetlinux05/Siva/raw_data/2023/sam_files/", "", stats_df$File)
stats_df$File <- sub("/home/vetlinux05/Siva/raw_data/2024/X204SC24073525-Z01-F001/01.RawData/sam_files/", "", stats_df$File)
stats_df$File <- sub(".*?(/)", "\\", stats_df$File)

# Save the data frame to a CSV file
write.csv(stats_df, output_file, row.names = FALSE)

# Print the resulting data frame to console
print(stats_df)
####
