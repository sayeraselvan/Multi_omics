#!/bin/bash

input_dir="/netscratch/dep_coupland/grp_fulgione/siva/outputgemma/FDR/gene_list_per_variable"
output_dir="/netscratch/dep_coupland/grp_fulgione/siva/outputgemma/FDR/GO_plots"

# Make sure the output directory exists
mkdir -p "$output_dir"

for input_file in "$input_dir"/*.txt; do
    # Extract the file name (without extension) for naming the output files
    file_name=$(basename -- "$input_file")
    file_name_wo_front="${file_name%spl_gene_names_}"
    file_name_wo_back="${file_name_wo_front%.*}"

    script="/netscratch/dep_coupland/grp_fulgione/siva/scripts/GOenrichment.R"
    Rscript "$script" "$input_file" "$output_dir/$file_name_wo_back"
done

