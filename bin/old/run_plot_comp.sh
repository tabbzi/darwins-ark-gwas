#!/bin/bash
#$ -q broad
#$ -l h_vmem=8g
#$ -l h_rt=2:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002

use GCC-5.2
use .htslib-1.8
use R-3.5

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/seq/vgb/dd/ancestry/env/
PATH=$PATH:/seq/vgb/software/anaconda3/current
source activate /seq/vgb/dd/gwas/env/Renv-2021

gwasA=${1}
gwasB=${2}

Rscript /seq/vgb/dd/gwas/bin/PC_qq.r ${gwasA} ${gwasB}
