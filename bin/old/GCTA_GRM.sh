#!/bin/bash
#$ -q broad
#$ -l h_vmem=32g
#$ -l h_rt=24:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002

gcta=/seq/vgb/software/gcta/current
cd /seq/vgb/dd/gwas/

# submit jobs to generate LD scores
qsub -N anc_run_${DATE} -hold_jid anc_preproc_${DATE} 

$gcta --make-grm --autosome-num 38 --autosome --bfile ./geno/${1} --out ./grm/${1} 

$gcta --make-grm --autosome-num 38 --autosome --bfile ${GENO} --out ${OUT} --extract ${SNPs}
