#!/bin/bash
#$ -q broad
#$ -l h_vmem=32g
#$ -l h_rt=8:00:00
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

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/env/

# software:
gemma=/seq/vgb/software/gemma/current
gcta=/seq/vgb/software/gcta/current
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

PATH=$PATH:/seq/vgb/software/anaconda3/current
source activate /seq/vgb/dd/gwas/env/Renv
Rscript /seq/vgb/dd/gwas/bin/manhattan-plot-levels.R ${MLMA} ${CLUMPS} <(echo ${TITLE}) 1e-6
