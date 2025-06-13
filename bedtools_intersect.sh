#!/bin/bash

input="/netscratch/dep_coupland/grp_fulgione/siva/outputgemma/file_list.txt"

while read line; do
	BIOCLIM=$(echo $line | cut -d',' -f1)
	bedtools intersect -a /netscratch/dep_coupland/grp_fulgione/siva/outputgemma/FDR/AA_chr/${BIOCLIM}.FDR.sig.bed -b /netscratch/dep_coupland/grp_fulgione/siva/Arabis_alpina_mpipz_v5.1_annotation.gff -wb | tee /netscratch/dep_coupland/grp_fulgione/siva/outputgemma/FDR/bedtools_intersect_wb/${BIOCLIM}_bedtools_intersect_FDR_wb_out.txt
done < $input

#!/bin/bash

input1="/netscratch/dep_coupland/grp_fulgione/siva/outputgemma/file_list1.txt"

while read line; do
	BIOCLIM=$(echo $line | cut -d',' -f1)
	bedtools intersect -a /netscratch/dep_coupland/grp_fulgione/siva/outputgemma/FDR/AA_chr/${BIOCLIM}.FDR.sig.bed -b /netscratch/dep_coupland/grp_fulgione/siva/Arabis_alpina_mpipz_v5.1_annotation.gff -wo | tee /netscratch/dep_coupland/grp_fulgione/siva/outputgemma/FDR/bedtools_intersect_wo/${BIOCLIM}_bedtools_intersect_FDR_wo_out.txt
done < $input1
