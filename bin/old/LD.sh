#!/bin/bash
#$ -q broad
#$ -l h_vmem=24g
#$ -l h_rt=48:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002
plink=/seq/vgb/software/plink2/current/plink
cd /seq/vgb/dd/gwas/

$plink --dog --bfile ./geno/${1} --r2 --out ./ld/byChrLong/${1}_chr_${2} --chr ${2}
