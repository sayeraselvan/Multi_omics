#!/bin/bash

if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <csv file> <excel sheet>"
	exit 1
fi

input_file="$1"
excel_sheet="$2"

rscript="/netscratch/dep_coupland/grp_fulgione/siva/scripts/lfmm.R"

Rscript "$rscript" "$input_file" "$excel_sheet"
