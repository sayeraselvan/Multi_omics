#!/bin/bash

plink --vcf SC_933.vcf.gz --make-bed --out SC_933
plink --bfile SC_933 --freq --out SC_933
plink --vcf SC_933 --recode --out SC_933
