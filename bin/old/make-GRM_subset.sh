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

$gcta --make-grm --autosome-num 38 --autosome --bfile ${GENO} --out ${OUT} --extract ${SNPs}
