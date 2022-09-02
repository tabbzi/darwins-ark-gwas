#!/bin/bash
#$ -q broad
#$ -l h_vmem=10g
#$ -l h_rt=24:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002

gcta=/seq/vgb/software/gcta/current
cd /seq/vgb/dd/gwas/
# $gcta --autosome-num 38 --autosome --bfile /seq/vgb/dd/gwas/geno/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt --ld-score-region 200 --out /seq/vgb/dd/gwas/ld/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_ld_200kb

$gcta --autosome --autosome-num 38 --chr ${CHR} --bfile ${GENO} --ld-score-region ${KB} --out ${OUT}
