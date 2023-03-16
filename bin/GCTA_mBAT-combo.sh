#!/bin/bash
#$ -q broad
#$ -l h_vmem=20g
#$ -l h_rt=12:00:00
#$ -o /seq/vgb/dd/gwas/logs/
#$ -e /seq/vgb/dd/gwas/logs/
#$ -M kmorrill@broadinstitute.org
#$ -m e
source /broad/software/scripts/useuse
umask 002

use GCC-5.2
use .htslib-1.8
use R-3.5
use BEDTools
use Parallel

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/env/

gcta=/seq/vgb/software/gcta/current

summ_stats=${summ_stats:-"/seq/vgb/dd/gwas/annotate/gwa/OCD/OCD-TS.NORDiC.ma"}
ld_ref=${ld_ref:-"/seq/vgb/dd/gwas/annotate/gwa/OCD/OCD-TS.NORDiC.1000G"}
region_bed=${region_bed:-"/seq/vgb/dd/gwas/annotate/regions/hg19/orthologous.map.RoCCs.hum.hg19.bed"}
output=${output:-"/seq/vgb/dd/gwas/mbat/hum/OCD-TS.NORDiC.RoCCs"}

${gcta} --bfile ${ld_ref} \
        --mBAT-combo ${summ_stats} \
        --mBAT-gene-list ${region_bed} \
        --mBAT-print-all \
        --out ${output}
