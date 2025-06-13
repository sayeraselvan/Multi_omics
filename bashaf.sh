#!/bin/bash

#conda init bash

env_name="r_env"

#conda activate $env_name

if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <file_path> <threshold_value>"
	exit 1
fi

file_path="$1"
threshold_value="$2"
filter="/netscratch/dep_coupland/grp_fulgione/siva/scripts/filter.py"
af="/netscratch/dep_coupland/grp_fulgione/siva/scripts/af3.py"
af1="/netscratch/dep_coupland/grp_fulgione/siva/scripts/af1.py"

python "$filter" "$file_path" "$threshold_value"

rm chrom_pos1.txt chrom_pos2.txt chrom_pos3.txt chrom_pos_chr1.txt chrom_pos_chr2.txt chrom_pos_chr3.txt 
rm chrom_pos.bed chrom_pos_chr.bed

input="/netscratch/dep_coupland/grp_fulgione/siva/scripts/awkfiles.txt"
input1="/netscratch/dep_coupland/grp_fulgione/siva/scripts/extraction.txt"
input2="/netscratch/dep_coupland/grp_fulgione/siva/scripts/extraction1.txt"
chrpos="chrom_pos.txt"

while read -r col1 col2 _; do
	echo "Submit ${col2}"
	bcftools view -R $chrpos ${col1} > ${col2}
done < "$input1"
echo "all variant extractions are done"

while read -r col1 _; do
	echo "csv converion ${col1}"
	bcftools query -f '%CHROM,%POS,%REF,%ALT[,\t%GT]\n' ${col1} > ${col1%.txt}.csv
done < "$input2"

#while read -r col1 col2 _; do
#	echo "Submit ${col2}"
#	awk 'NR==FNR{a[$1$2]; next} ($1$2 in a)' $chrpos ${col1} > ${col2}
#done < "$input"

spain="spain"
scad="scad"
alps="alps"
txt="txt"
ind="ind"

mkdir "$spain"
mkdir "$scad"
mkdir "$alps"
mkdir "$ind"

filename="bran_spv.csv"
pathy=$(readlink -f "$filename")
export FILE_PATH="${pathy%bran_spv.csv}"

python $af

mv bran_spv.csv E10_exp2v.csv E17_e15v.csv E1_E6_exp2v.csv E20v.csv E21v.csv E22v.csv E2v.csv E3_exp2v.csv E4_exp2_spv.csv E5v.csv est_E18v.csv leit_spainv.csv ord_e19_spv.csv pajaresv.csv pena_e9v.csv PORTv.csv spain_pyrenesv.csv vega_spav.csv yeclav.csv $spain

mv nor1v.csv nor2v.csv nor4v.csv nor5v.csv nor6v.csv nor7v.csv par_nor6v.csv s1v.csv s2v.csv S3v.csv S5v.csv sve2v.csv sve3v.csv sve4v.csv icelandv.csv $scad

mv italymixedv.csv austriav.csv bri_frv.csv chv.csv cr_frv.csv D1_D2v.csv D3v.csv D4v.csv Dotherv.csv Fgal1v.csv fr1_5v.csv fr4v.csv fr6v.csv fr_Pv.csv galiber_frv.csv ganon_frv.csv gv_frv.csv polandv.csv vallon_frv.csv $alps

mv SCv.csv spainv.csv alpsv.csv norwayv.csv swedenv.csv scadv.csv $ind

echo "all jobs sended"

#xlsx2csv af_ref.xlsx af_ref.csv
#xlsx2csv af_alt.xlsx af_alt.csv
#xlsx2csv gt_het.xlsx gt_het.csv

#awk '{print $1 "_" $2}' chrom_pos.txt > temp.txt
#awk 'BEGIN {FS=","; OFS=","} NR==FNR{a[NR]=$0; next} {print a[FNR], $0}' temp.txt af_ref.csv > af_ref_temp.csv
#awk 'BEGIN{FS=OFS=","} {sum = ($2 + $3 + $4) / 3; print $0, sum}' af_ref_temp.csv > af_ref.csv

#awk 'BEGIN {FS=","; OFS=","} NR==FNR{a[NR]=$0; next} {print a[FNR], $0}' temp.txt af_alt.csv > af_alt_temp.csv
#awk 'BEGIN{FS=OFS=","} {sum = ($2 + $3 + $4) / 3; print $0, sum}' af_alt_temp.csv > af_alt.csv

#awk 'BEGIN {FS=","; OFS=","} NR==FNR{a[NR]=$0; next} {print a[FNR], $0}' temp.txt gt_het.csv > gt_het_temp.csv
#awk 'BEGIN{FS=OFS=","} {sum = ($2 + $3 + $4) / 3; print $0, sum}' gt_het_temp.csv > gt_het.csv

#python $af1

rm af_ref.xlsx af_alt.xlsx temp.txt af_ref_temp.csv af_alt_temp.csv gt_het.xlsx gt_het_temp.csv 
rm gt_het.csv af_alt.csv af_ref.csv
rm -r $spain
rm -r $ind
rm -r $alps
rm -r $scad

mkdir "$txt"

mv alpsv.txt spainv.txt scadv.txt SCv.txt norwayv.txt swedenv.txt austriav.txt bran_spv.txt bri_frv.txt chv.txt cr_frv.txt D1_D2v.txt D3v.txt D4v.txt Dotherv.txt E10_exp2v.txt E17_e15v.txt E1_E6_exp2v.txt E20v.txt E21v.txt E22v.txt E2v.txt E3_exp2v.txt E4_exp2_spv.txt E5v.txt est_E18v.txt Fgal1v.txt fr1_5v.txt fr4v.txt fr6v.txt fr_Pv.txt galiber_frv.txt ganon_frv.txt gv_frv.txt leit_spainv.txt ord_e19_spv.txt pajaresv.txt par_nor6v.txt pena_e9v.txt polandv.txt PORTv.txt s1v.txt s2v.txt S3v.txt S5v.txt spain_pyrenesv.txt vallon_frv.txt vega_spav.txt yeclav.txt nor1v.txt nor2v.txt nor4v.txt nor5v.txt nor6v.txt nor7v.txt sve2v.txt sve3v.txt sve4v.txt italymixedv.txt icelandv.txt $txt

rm -r $txt

rm chrom_pos_chr.txt chr_pthreshold.txt

echo "good job"

#python os.py
#this bash script calculates af, gt, maf, major, hf, hetaf using af.py, af1.py and then we get the csv file with all allele freq info
