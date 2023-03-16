#!/bin/bash
#$ -q broad
#$ -l h_vmem=32g
#$ -l h_rt=36:00:00
#$ -o /seq/vgb/dd/gwas/logs/
#$ -e /seq/vgb/dd/gwas/logs/
#$ -M kmorrill@broadinstitute.org
#$ -m e
#$ -pe smp 4
#$ -binding linear:4
#$ -R y
source /broad/software/scripts/useuse
umask 002

use GCC-5.2
use .htslib-1.8
use R-3.5
use BEDTools
use Parallel

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/env/

source activate ldsc
