#!/bin/bash
#$ -q broad
#$ -l h_vmem=64g
#$ -l h_rt=400:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002
  reuse Tabix
  reuse Bcftools
use UGER
gcta=/seq/vgb/software/gcta/current
plink=/seq/vgb/software/plink2/current/plink

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# PURPOSE:
# The purpose of this script is to list sample barcode ~ dog ID pairs that exist in the Darwin's Ark database;
# then, to construct the GWAS data set according to proper exclusions and genotype data filtering practices.

cd /seq/vgb/dd/gwas/geno/construction/

while read m
do
tabix -p vcf ${m}
done < ./merge/${1}

bcftools merge -l ./merge/${1} --threads 4 -o ./merge/${1}.vcf.gz -Oz
tabix -p vcf ./merge/${1}.vcf.gz
$plink --dog --double-id --vcf ./merge/${1}.vcf.gz --make-bed --out ./merge/${1}
