#!/bin/bash
for i in {2..8};do
    cd chr$i
    mv _chr$i.vcf.gz chr$i.vcf.gz
    plink --vcf chr$i.vcf.gz --make-bed --out chr$i 
    plink --bfile chr$i --freq --out chr$i 
    plink --vcf chr$i.vcf.gz --recode --out chr$i
    cd ../
done
