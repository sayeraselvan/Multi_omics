1. Connect to vetlinux05: ssh vetlinux05@pgnsrv043.vu-wien.ac.at
password: HA09B2205

check all available environments using 'conda info --envs'
then create your own environment: 'conda create --name Neda_slim'
to use or activate: 'conda activate Neda_slim'
to deactivate: 'conda deactivate'
#install a package in current environment: 'conda install python=3.8'
#install slim : 'conda install -c conda-forge slim'
#install numpy: 'conda install numpy'
check the list of stuff installed in your environment: 'conda list'

#2. install cutadapt: 'conda install cutadapt'
#check the version 'cutadapt --version'
#transfer files to vetgrid05:scp /Volumes/Temp2/SLIM/scripts/wrapper_burnin.py  vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/scripts

2. install trimmomatic 'conda install -c bioconda trimmomatic'
the adaptors used for DNA library preparations by Novogene are 'PCR_Primer1_rc' and 'PCR_Primer2_rc' in TruSeq2_PE.fa in Trimmomatic

3. transfer sequencing files to server
scp /Volumes/Temp2/bigcages/Pool-seq/raw_data/X204SC22122682-Z01-F001/01.RawData/DSIM_98/DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_1.fq.gz vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/DSIM_98

scp /Volumes/Temp2/bigcages/Pool-seq/raw_data/X204SC22122682-Z01-F001/01.RawData/DSIM_98/DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_2.fq.gz vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/DSIM_98

4. Trim reads using trimmomatic
#https://github.com/usadellab/Trimmomatic

4.1 TruSeq2-PE.fa was downloaded from https://github.com/usadellab/Trimmomatic and transferre to Vetlinux05: scp /Volumes/Temp2/bigcages/software/TruSeq2-PE.fa vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/

trimmomatic PE DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_1.fq.gz DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_2.fq.gz DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_1_paired.fq.gz DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_1_unpaired.fq.gz DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_2_paired.fq.gz DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_2_unpaired.fq.gz ILLUMINACLIP:/home/vetlinux05/Neda/DNA-seq/TruSeq2-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:15 MINLEN:35

5. reference genome was downloaded from here:
 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000754195.2/
 (https://www.ncbi.nlm.nih.gov/sra/SRX159097[accn])
We use this referene because Signor et al., 2018 mapped the inbred lines against this assembly

5.1 the names of the major chromosomes are changed
CM002914.1    20829647    18704147    0 #X
CM002910.1    23539531    24171762    0 #2L
CM002911.1    21544594    21566760    0 #2R
CM002912.1    24153973    25303481    0 #3L
CM002913.1    27160941    29984108    0 #3R
CM002916.1    1026345    906958    0 #4

#transfer the renamed reference genome to vetlinus05
5.2
scp /Volumes/Temp2/bigcages/ref_genome/ncbi_dataset/ncbi_dataset/data/GCA_000754195.3/GCA_000754195.3_ASM75419v3_genomic_renamed.fna vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/ref_genome

6. install bwa-mem2
conda install -c bioconda bwa-mem2

6.1 index genome
bwa-mem2 index -p dsim_Hu_v3_renamed GCA_000754195.3_ASM75419v3_genomic_renamed.fna

6.2 map
bwa-mem2 mem -t 20 /home/vetlinux05/Neda/DNA-seq/ref_genome/dsim_Hu_v3_renamed DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_1_paired.fq.gz DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_2_paired.fq.gz > DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2.sam

7. install picard
'conda install -c bioconda picard'

7.1 sort sam file and convert to bam
picard SortSam I=DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2.sam O=DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_sorted.bam SORT_ORDER=coordinate
      
# You can then remove duplicates using Picard MarkDuplicates
picard -Xmx20g MarkDuplicates REMOVE_DUPLICATES=true I=DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_sorted.bam O=DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_sorted_rmdup.bam M=rmdup_metrics.txt VALIDATION_STRINGENCY=SILENT

# then filter reads that mapped with low quality
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_sorted_rmdup.bam >> DSIM_98_EKDN220048412-1A_HLTJ3DSX5_L2_sorted_rmdup.bam.filter.bam

#-q 20: filters out reads having a low mapping quality (<20).
#-f 0x0002: filters reads that are mapped in a proper pair.
#-F 0x0004: removes unmapped reads.
#-F 0x0008: removes reads where the mate is unmapped.

8. download teh raw sequencing reads from Signor et al. 2018
 8.1 install sratoolkit
 
tar -vxzf /Volumes/Temp3/D_sim_Signor/tools/sratoolkit.3.0.10-mac-x86_64.tar.gz
export PATH=$PATH:/Volumes/Temp3/D_sim_Signor/tools/sratoolkit.3.0.10-mac-x86_64/bin

8.2 downloaded the info regarding the sample name and run name from NCBI here:/Volumes/Temp3/D_sim_Signor/SraRunInfo.csv

using this guideline: https://erilu.github.io/python-fastq-downloader/
8.3 downloaded the raw reads using /Volumes/Temp3/D_sim_Signor/scripts/download_sra.py

9. trim the downloaded files

###############################
# Pool-seq data
###############################
1. transfer files to vetlinux

scp /Volumes/Temp2/bigcages/Pool-seq/raw_data/X204SC22122682-Z01-F001/01.RawData/TS_R4_F40/TS_R4_F40_EKDN220048410-1A_HLTJ3DSX5_L2_1.fq.gz vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/TS_R4_F40

scp /Volumes/Temp2/bigcages/Pool-seq/raw_data/X204SC22122682-Z01-F001/01.RawData/TS_R4_F40/TS_R4_F40_EKDN220048410-1A_HLTJ3DSX5_L2_2.fq.gz vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/TS_R4_F40

 2. remove adapters and trim reads using cutadapt
 
 #sequencing adapters:
 #5' adapter
 #5'-AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT-3
 #3' adapter
 #5'-GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG-3'
 
 #-a is for 3' adapter and -g is for 5' adapter

 2.1 install cutadapt
 conda install -n Neda_slim cutadapt
 cutadapt --version
 
2.1 remove adaptors, remove low quality bases from the end of reads, and short reads

cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -A GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG -o TS_R4_F40_trimmed_1.fq.gz -p TS_R4_F40_trimmed_2.fq.gz TS_R4_F40_EKDN220048410-1A_HLTJ3DSX5_L2_1.fq.gz TS_R4_F40_EKDN220048410-1A_HLTJ3DSX5_L2_2.fq.gz --cores=12 --quality-cutoff 18 --minimum-length 35 --no-indels --pair-filter=any --action=trim

3. map the reads to the genome (genome has been indexed before, see above)

3.1
bwa-mem2 mem -t 20 /home/vetlinux05/Neda/DNA-seq/ref_genome/dsim_Hu_v3_renamed TS_R4_F40_trimmed_1.fq.gz TS_R4_F40_trimmed_2.fq.gz > TS_R4_F40.sam

4. sort sam file and convert to bam
picard SortSam I=TS_R4_F40.sam O=TS_R4_F40_sorted.bam SORT_ORDER=coordinate
      
5. Remove duplicates using Picard MarkDuplicates
picard -Xmx20g MarkDuplicates REMOVE_DUPLICATES=true I=TS_R4_F40_sorted.bam O=TS_R4_F40_sorted_rmdup.bam M=rmdup_metrics.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT

6. Filter reads that mapped with low quality

6.1 install samtools
conda install -n Neda_slim bioconda::samtools

6.2 filter the reads to keep only properly mapped paired end reads
samtools view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b TS_R4_F40_sorted_rmdup.bam >> TS_R4_F40_sorted_rmdup_filter.bam

7. compute coverage
7.1 install samtools

conda activate Neda_slim
conda config --add channels bioconda
conda config --add channels conda-forge
conda install -c bioconda samtools

7.2 compute coverage
# https://www.htslib.org/doc/samtools-coverage.html
#get a list of chromosome names

The tabulated form uses the following headings.

#rname    Reference name / chromosome
#startpos    Start position
#endpos    End position (or sequence length)
#numreads    Number reads aligned to the region (after filtering)
#covbases    Number of covered bases with depth >= 1
#coverage    Percentage of covered bases [0..100]
#meandepth    Mean depth of coverage
#meanbaseq    Mean baseQ in covered region
#meanmapq    Mean mapQ of selected reads

#first index the bam file
samtools index TS_R4_F40_sorted_rmdup_filter.bam
samtools index TS_R4_F40_sorted_rmdup.bam
samtools index TS_R4_F40_sorted.bam

#then compute coverage
samtools coverage -r X TS_R4_F40_sorted.bam
samtools coverage -r 2L TS_R4_F40_sorted.bam
samtools coverage -r 2R TS_R4_F40_sorted.bam
samtools coverage -r 3L TS_R4_F40_sorted.bam
samtools coverage -r 3R TS_R4_F40_sorted.bam
samtools coverage -r 4 TS_R4_F40_sorted.bam

samtools coverage -r X TS_R4_F40_sorted.bam #average coverage 177.088
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
X    1    20829647    28632514    19382221    93.0511    166.539    35.4    49.7
samtools coverage -r 2L TS_R4_F40_sorted.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
2L    1    23539531    33420584    21655582    91.9967    184.752    35.5    53.5
samtools coverage -r 2R TS_R4_F40_sorted.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
2R    1    21544594    29852738    19704439    91.4589    181.77    35.5    53.5
samtools coverage -r 3L TS_R4_F40_sorted.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
3L    1    24153973    33743512    22231071    92.039    184.217    35.5    55.4
samtools coverage -r 3R TS_R4_F40_sorted.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
3R    1    27160941    37816578    25926798    95.4562    187.308    35.5    58.5
samtools coverage -r 4 TS_R4_F40_sorted.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
4    1    1026345    1475889    814580    79.3671    157.942    35.5    48.9

#then compute coverage
samtools coverage -r X TS_R4_F40_sorted_rmdup.bam
samtools coverage -r 2L TS_R4_F40_sorted_rmdup.bam
samtools coverage -r 2R TS_R4_F40_sorted_rmdup.bam
samtools coverage -r 3L TS_R4_F40_sorted_rmdup.bam
samtools coverage -r 3R TS_R4_F40_sorted_rmdup.bam
samtools coverage -r 4 TS_R4_F40_sorted_rmdup.bam

samtools coverage -r X TS_R4_F40_sorted_rmdup.bam #average coverage 143.7
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
X    1    20829647    21992566    19375545    93.0191    129.282    35.5    51.7
samtools coverage -r 2L TS_R4_F40_sorted_rmdup.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
2L    1    23539531    27409530    21649629    91.9714    151.841    35.6    53.6
samtools coverage -r 2R TS_R4_F40_sorted_rmdup.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
2R    1    21544594    24565411    19699081    91.434    149.419    35.6    53.4
samtools coverage -r 3L TS_R4_F40_sorted_rmdup.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
3L    1    24153973    27794555    22225009    92.0139    151.46    35.6    55.3
samtools coverage -r 3R TS_R4_F40_sorted_rmdup.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
3R    1    27160941    31004533    25919996    95.4311    153.682    35.6    58.6
samtools coverage -r 4 TS_R4_F40_sorted_rmdup.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
4    1    1026345    1170870    814029    79.3134    126.518    35.6    48.6

samtools coverage -r X TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 2L TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 2R TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 3L TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 3R TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 4 TS_R4_F40_sorted_rmdup_filter.bam

samtools coverage -r X TS_R4_F40_sorted_rmdup_filter.bam #average coverage 130.75
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
X    1    20829647    18704147    19142855    91.902    114.817    35.5    58.6
samtools coverage -r 2L TS_R4_F40_sorted_rmdup_filter.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
2L    1    23539531    24171762    20555724    87.3243    138.092    35.6    59
samtools coverage -r 2R TS_R4_F40_sorted_rmdup_filter.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
2R    1    21544594    21566760    18658900    86.6059    135.478    35.6    59
samtools coverage -r 3L TS_R4_F40_sorted_rmdup_filter.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
3L    1    24153973    25303481    21418064    88.673    140.955    35.6    59.2
samtools coverage -r 3R TS_R4_F40_sorted_rmdup_filter.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
3R    1    27160941    29984108    25865954    95.2322    150.2    35.6    59.5
samtools coverage -r 4 TS_R4_F40_sorted_rmdup_filter.bam
#rname    startpos    endpos    numreads    covbases    coverage    meandepth    meanbaseq    meanmapq
4    1    1026345    906958    776060    75.614    104.985    35.7    57.7

8. filter VCF file

8.1 transfer CVF file to vetgrid05
scp /Volumes/Temp3/D_sim_Signor/simulans_multisamp_all_chr.vcf vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/vcf_haplo/

8.2 install vcftools
conda install -n Neda_slim vcftools

8.3 check heterozygosity
vcftools --vcf simulans_multisamp_all_chr.vcf --het --out simulans_het_per_in
#After filtering, kept 195 out of 195 Individuals
#Outputting Individual Heterozygosity
#    Individual Heterozygosity: Only using biallelic SNPs

8.4 amount of missing data per individual
vcftools --vcf simulans_multisamp_all_chr.vcf --missing-indv --out simulans_het_per_in
#After filtering, kept 195 out of 195 Individuals
#Outputting Individual Missingness
#After filtering, kept 7305030 out of a possible 7305030 Sites
#Run Time = 164.00 seconds

8.5 vcftools
follow the basic filtering from Signor paper: https://github.com/signor-molevol/signor_popgen_2016/blob/master/documentation/Documentation.txt

A. #this command removes removes SNPs that are not on major chromosomes, that are not bi allelic, removes indels, lines that have >25% missing data or heterozygosity, and loci with >10% missing data
vcftools --vcf simulans_multisamp_all_chr.vcf --chr X --chr 2L --chr 2R --chr 3L --chr 3R --chr 4 --remove-indv Sz116 --remove-indv Sz12 --remove-indv Sz158 --remove-indv Sz21 --remove-indv Sz237 --remove-indv Sz286 --remove-indv Sz76 --remove-indv Sz176 --remove-indv md106ts --remove-indv nc485 --remove-indv newc --remove-indv w501 --remove-indv 128 --remove-indv 132 --remove-indv 176 --remove-indv 74 --remove-indv 84 --remove-indv C167.4 --remove-indv MD106 --remove-indv MD199 --remove-indv STE  --remove-indv Sz57 --remove-indv Sz210 --max-missing-count 17 --remove-indels --min-alleles 2 --max-alleles 2 --recode --out simulans_multisamp_majorchr_filtered_10missingData

'''
Warning: Expected at least 2 parts in FORMAT entry: ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
Warning: Expected at least 2 parts in INFO entry: ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=MLEAC,Number=A,Type=Integer,Description="Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
Warning: Expected at least 2 parts in INFO entry: ID=MLEAF,Number=A,Type=Float,Description="Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed">
Excluding individuals in 'exclude' list
After filtering, kept 172 out of 195 Individuals
Outputting VCF file...
After filtering, kept 5166011 out of a possible 7305030 Sites
Run Time = 1867.00 seconds
'''

8.6 vcftools
#this command removes removes SNPs that are not on major chromosomes, that are not bi allelic, removes indels, lines that have >25% missing data or heterozygosity, and loci with ANY missing data

vcftools --vcf simulans_multisamp_all_chr.vcf --chr X --chr 2L --chr 2R --chr 3L --chr 3R --chr 4 --remove-indv Sz116 --remove-indv Sz12 --remove-indv Sz158 --remove-indv Sz21 --remove-indv Sz237 --remove-indv Sz286 --remove-indv Sz76 --remove-indv Sz176 --remove-indv md106ts --remove-indv nc485 --remove-indv newc --remove-indv w501 --remove-indv 128 --remove-indv 132 --remove-indv 176 --remove-indv 74 --remove-indv 84 --remove-indv C167.4 --remove-indv MD106 --remove-indv MD199 --remove-indv STE  --remove-indv Sz57 --remove-indv Sz210 --max-missing 1 --remove-indels --min-alleles 2 --max-alleles 2 --recode --out simulans_multisamp_majorchr_filtered_NOmissingData

'''Excluding individuals in 'exclude' list
After filtering, kept 172 out of 195 Individuals
Outputting VCF file...
After filtering, kept 29095 out of a possible 7305030 Sites
Run Time = 145.00 seconds'''

9. downsample the bam file
picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.73.bam PROBABILITY=0.73 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.59.bam PROBABILITY=0.59 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.44.bam PROBABILITY=0.44 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.29.bam PROBABILITY=0.29 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.15.bam PROBABILITY=0.15 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.07.bam PROBABILITY=0.07 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.04.bam PROBABILITY=0.04 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

picard -Xmx20g DownsampleSam I=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam O=/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/downsample/TS_R4_F40_sorted_rmdup_filter_cov0.02.bam PROBABILITY=0.02 RANDOM_SEED=1 VALIDATION_STRINGENCY=SILENT

#index the bam files
samtools index TS_R4_F40_sorted_rmdup_filter_cov0.73.bam

samtools coverage -r X TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 2L TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 2R TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 3L TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 3R TS_R4_F40_sorted_rmdup_filter.bam
samtools coverage -r 4 TS_R4_F40_sorted_rmdup_filter.bam

10. GATK
scp /Volumes/Temp2/bigcages/software/gatk-4.5.0.0.zip  vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/tools

conda env create -n gatk -f gatkcondaenv.yml
conda activate gatk
conda list #this should be among the installed packages: gatkpythonpackages


#remove SNPs that are not on chromosomes

11. Install harp
11.1 on vetlinux05
scp /Volumes/Temp2/bigcages/software/harp_linux_140925_103521.zip vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/

unzip harp_linux_140925_103521.zip

11.2 it seems it neds boost
#first try
conda config --add channels conda-forge
conda config --set channel_priority strict
conda install boost boost-cpp libboost libboost-devel libboost-headers libboost-python libboost-python-devel

The following NEW packages will be INSTALLED:

  boost              conda-forge/linux-64::boost-1.84.0-hb563948_2
  boost-cpp          conda-forge/linux-64::boost-cpp-1.84.0-h44aadfe_2
  libboost           conda-forge/linux-64::libboost-1.84.0-h8013b2b_2
  libboost-devel     conda-forge/linux-64::libboost-devel-1.84.0-h00ab1b0_2
  libboost-headers   conda-forge/linux-64::libboost-headers-1.84.0-ha770c72_2
  libboost-python    conda-forge/linux-64::libboost-python-1.84.0-py38hae673b5_2
  libboost-python-d~ conda-forge/linux-64::libboost-python-devel-1.84.0-py38hb563948_2

#second try
conda install conda-forge::boost

#third try
scp /Volumes/Temp2/bigcages/software/b2-5.1.0.tar.gz vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/tools

tar -vxzf b2-5.1.0.tar.gz

./bootstrap.sh
./b2 install --prefix=/home/vetlinux05/Neda/DNA-seq/tools/b2-5.1.0

export PATH=$PATH:/home/vetlinux05/Neda/DNA-seq/tools/b2-5.1.0/bin/b2

#run in /home/vetlinux05/Neda/DNA-seq/tools/dkessner-harp-a9af1a12b391/src
/home/vetlinux05/Neda/DNA-seq/tools/b2-5.1.0/bin/b2 Jamroot #this didn't work

export PATH=$PATH:/home/vetlinux05/Neda/DNA-seq/tools/harp_linux_140925_103521/bin

11.3 I will run harp on mac
#new trial, let's do it on mac
touch ~/.bash_profile; open ~/.bash_profile

export PATH=/Volumes/Temp2/bigcages/software/harp_osx_140925_100959/bin:$PATH
echo $PATH
#touch ~/.zshrc; open ~/.zshrc

12. troubleshooting HAF

#it just gives error after error so I made the following changes:
#error 1
#in HAFpipe-line-master/scripts/make_SNPtable_from_vcf.sh, I changed 'grep 1000' to 'grep 10000' in line 89 because the vcf file has many fragmented contigs, so the first line of header is not in the top 1000 lines.

#error 2
#I also removed '-P' in front of grep on line 93

12.1
#first transfer the filtered VCF file from Vetlinux05
scp vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/vcf_haplo/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf /Volumes/Temp3/D_sim_Signor

12.1.1 Step 1
#HAF: function 1, which is making SNP table from the VCF file
#convert VCF file to SNP table
 /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/2L -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 2L -s /Volumes/Temp3/D_sim_Signor/SNP_table/2L/snp_table_2L -k

 /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/2R -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 2R -s /Volumes/Temp3/D_sim_Signor/SNP_table/2R/snp_table_2R -k
 
  /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3L -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 3L -s /Volumes/Temp3/D_sim_Signor/SNP_table/3L/snp_table_3L -k

 /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3R -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 3R -s /Volumes/Temp3/D_sim_Signor/SNP_table/3R/snp_table_3R -k
 
  /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/4 -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 4 -s /Volumes/Temp3/D_sim_Signor/SNP_table/4/snp_table_4 -k
 
12.1.2 Step 2
#HAF: function 2, which is imputing missing SNP

/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/2L -s  /Volumes/Temp3/D_sim_Signor/SNP_table/2L/snp_table_2L -i npute

/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/2R -s  /Volumes/Temp3/D_sim_Signor/SNP_table/2R/snp_table_2R -i npute

/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3L -s  /Volumes/Temp3/D_sim_Signor/SNP_table/3L/snp_table_3L -i npute

/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3R -s  /Volumes/Temp3/D_sim_Signor/SNP_table/3R/snp_table_3R -i npute

/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/4 -s  /Volumes/Temp3/D_sim_Signor/SNP_table/4/snp_table_4 -i npute

12.1.3 Step 3 and 4

scp vetlinux05@pgnsrv043.vu-wien.ac.at:/home/vetlinux05/Neda/DNA-seq/TS_R4_F40/TS_R4_F40_sorted_rmdup_filter.bam /Volumes/Temp2/bigcages/Pool-seq/bam_files

#HAF: function 3, which is inferring haplotype frequency and function 4 is estimating the allele frequencies
# -a .0000000304 (which is 3.04 cM/Mb) is from https://www.biorxiv.org/content/10.1101/2022.09.12.507595v1.full

#first index the bam file
samtools index /Volumes/Temp2/bigcages/Pool-seq/bam_files/TS_R4_F40_sorted_rmdup_filter.bam

#error 3
#the index_snp_table in HAFpipe-line-master folder is a binary and cannot be run on mac. So I replaced it with the index_snp_table from harp_osx_140925_100959

#error 4
#Also in infer_haplotype_freqs.sh the start of window was in scientific notion so I added this to line 129-130 to convert it to integer.

    ### Convert start to integer using awk
    #start=$(awk -v start="${start}" 'BEGIN { printf "%.0f", start }')
    
#error 5
#it seems that harp will throw an error if the ref_genome is a mix of small and capital letters
#so, first convert the genome to all caps

convert lowercase to uppercase in reference genome using /Volumes/Temp3/D_sim_Signor/scripts/convert_fasta_uppercase.py

13. #finally run HAP-pipe
/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 3,4 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/2L -s  /Volumes/Temp3/D_sim_Signor/SNP_table/2L/snp_table_2L.npute -c 2L -r /Volumes/Temp2/bigcages/ref_genome/ncbi_dataset/ncbi_dataset/data/GCA_000754195.3/GCA_000754195.3_ASM75419v3_genomic_renamed_uppercase.fna -e sanger -b /Volumes/Temp2/bigcages/Pool-seq/bam_files/TS_R4_F40_sorted_rmdup_filter.bam -g 40 -a .0000000304

#error 6
#some genomic ranges were not valid (error by HARP).
#I checked the positions in SNP_table (estimated from VCF) and compared them with the genomic positions in the reference genome using /Volumes/Temp3/D_sim_Signor/scripts/compare_snptable_ref.py. There are some positions that are not present in the reference. It seems that the reference that Sara Signor mapped the data to (and hence made VCF) is different from the one in genbank. But chromosome 3R is the same. So I will only use this choromosome for the benchmarking.

/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 3,4 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3R -s  /Volumes/Temp3/D_sim_Signor/SNP_table/3R/snp_table_3R.npute -c 3R -r /Volumes/Temp2/bigcages/ref_genome/ncbi_dataset/ncbi_dataset/data/GCA_000754195.3/GCA_000754195.3_ASM75419v3_genomic_renamed_uppercase.fna -e sanger -b /Volumes/Temp2/bigcages/Pool-seq/bam_files/TS_R4_F40_sorted_rmdup_filter.bam -g 40 -a .0000000304

13.1 #run HAP-pipe for the full overage bam files with a SNP_table that heterozygoutes are converted to N
#convert VCF file to SNP table
 /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1,2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3R/3R_hetero_N -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 3R -s /Volumes/Temp3/D_sim_Signor/SNP_table/3R/3R_hetero_N/snp_table_3R_hetero_N -i npute
 
/Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 3,4 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3R/3R_hetero_N -s  /Volumes/Temp3/D_sim_Signor/SNP_table/3R/3R_hetero_N/snp_table_3R_hetero_N -c 3R -r /Volumes/Temp2/bigcages/ref_genome/ncbi_dataset/ncbi_dataset/data/GCA_000754195.3/GCA_000754195.3_ASM75419v3_genomic_renamed_uppercase.fna -e sanger -b /Volumes/Temp2/bigcages/Pool-seq/bam_files/TS_R4_F40_sorted_rmdup_filter.bam -g 40 -a .0000000304

13.2 #make the SNPTable using only lines that were used in the crosses. The list is in /Volumes/Temp2/bigcages/genotype/inbred_lines_RR.xlsx, sheet:  lines_VCF_filtered, column G

 /Volumes/Temp2/bigcages/software/HAFpipe-line-master/HAFpipe_wrapper.sh -t 1,2 -d /Volumes/Temp2/bigcages/software/HAFpipe-line-master -o /Volumes/Temp3/D_sim_Signor/SNP_table/3R/3R_hetero_N_cross_lines -v /Volumes/Temp3/D_sim_Signor/simulans_multisamp_majorchr_filtered_10missingData.recode.vcf -c 3R -s /Volumes/Temp3/D_sim_Signor/SNP_table/3R/3R_hetero_N_cross_lines/snp_table_3R_hetero_N_cross_lines -i npute --subsetlist /Volumes/Temp3/D_sim_Signor/SNP_table/keep_lines

#I got an error with the above command so I used -i simpute instead. But the output was empty. I ran it without '-i npute', still empty output. I think there is an error in using '--subsetlist' 

#

14. #compute allele frequency based on read count
samtools mpileup --max-depth 1000000 -B -Q 0 -f /Volumes/Temp2/bigcages/ref_genome/ncbi_dataset/ncbi_dataset/data/GCA_000754195.3/GCA_000754195.3_ASM75419v3_genomic_renamed.fna /Volumes/Temp2/bigcages/Pool-seq/bam_files/TS_R4_F40_sorted_rmdup_filter.bam | java -Xmx20g -jar /Volumes/Temp3/programs/popoolation2_1201/mpileup2sync.jar --input /dev/stdin --min-qual 20 --output TS_R4_F40_sorted_rmdup_filter.sync --fastq-type sanger --threads 10


15. HARP # I will run harp to see how the estimated haplotype freqeuncy changes

index_snp_table /Volumes/Temp3/D_sim_Signor/SNP_table/2L/snp_table_2L.npute 50000

harp like -b /Volumes/Temp2/bigcages/Pool-seq/bam_files/TS_R4_F40_sorted_rmdup_filter.bam --refseq /Volumes/Temp2/bigcages/ref_genome/ncbi_dataset/ncbi_dataset/data/GCA_000754195.3/chr2L_uppercase.fna.txt -r 2L:1-1580001 --stem /Volumes/Temp3/D_sim_Signor/SNP_table/2L/TS_R4_F40_sorted_rmdup_filter.bam.2L-1-1580001 --snps /Volumes/Temp3/D_sim_Signor/SNP_table/2L/snp_table_2L.npute

harp freq --hlk /Volumes/Temp3/D_sim_Signor/SNP_table/2L/TS_R4_F40_sorted_rmdup_filter.bam.2L-1-1580001.output/TS_R4_F40_sorted_rmdup_filter.bam.2L-1-1580001.hlk -r 2L:1-1580001


