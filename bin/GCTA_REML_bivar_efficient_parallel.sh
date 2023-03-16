#!/bin/bash
#$ -q broad
#$ -l h_vmem=20g
#$ -l h_rt=2:00:00
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

set -a

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/env/

# software:
gemma=/seq/vgb/software/gemma/current
gcta=/seq/vgb/software/gcta/current
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

# prepare:
DIR='/seq/vgb/dd/gwas'
DATE=${DATE:-'2022-11-20'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

SCRIPT=${DIR}'/bin/GCTA_REML_bivar_efficient_parallel_run.sh'

mkdir ${DIR}'/cor/'${DATE}

# genetic dataset:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465'}

# phenotypes:
PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
PHE=${PHE:-'all'}
PHEFILE=${PHEFILE:-${PHEDIR}'/'${PHE}'.tsv'}

# selections:
KEY=${KEY:-${PHEDIR}'/'${PHE}'.key.txt'}
PAIRS=${PAIRS:-${PHEDIR}'/'${PHE}'.pairs.tsv'}

# discrete covariates:
DCOVDIR=${DCOVDIR:-${DIR}'/covar/'${DATE}'/dcov'}
DCOV=${DCOV:-'datatype'}
DCOVFILE=${DCOVFILE:-${DCOVDIR}'/'${DCOV}'.tsv'}
if [ ${DCOV} == 'NA' ] ; then
  echo "Not using any discrete covariates"
  ARG_DCOVAR=''
else
  ARG_DCOVAR='--covar '${DCOVFILE}
fi

# quantitative covariates:
QCOVDIR=${QCOVDIR:-${DIR}'/covar/'${DATE}'/qcov'}
QCOV=${QCOV:-"NA"}
QCOVFILE=${QCOVFILE:-${QCOVDIR}'/'${QCOV}'.tsv'}
if [ ${QCOV} == 'NA' ] ; then
  echo "Not using any quantitative covariates"
  ARG_QCOVAR=''
else
  ARG_QCOVAR='--qcovar '${QCOVFILE}
fi

# relatedness cutoff:
RELCUT=${RELCUT:-"0.75"}
if [ ${RELCUT} == 'NA' ] ; then
  echo "Not using a relatedness cutoff"
  useRELCUT='F'
else
  echo "Using relatedness cutoff"
  useRELCUT='T'
fi

# grm:
if [ ${useRELCUT} == 'T' ] ; then
  INPUT_GRM_ALL=${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT}
  INPUT_MGRM_LIST=${DIR}'/grm/'${GENO}'.score.ld.rel-cutoff-'${RELCUT}'.list.txt'
else
  INPUT_GRM_ALL=${DIR}'/grm/'${GENO}
  INPUT_MGRM_LIST=${DIR}'/grm/'${GENO}'.score.ld.list.txt'
fi

parallel --colsep '\t' --jobs 4 --env DIR,DATE,INPUT_GRM_ALL,INPUT_MGRM_LIST,PHEFILE,ARG_DCOVAR,ARG_QCOVAR,RELCUT,useRELCUT,GENO,PHE,DCOV,QCOV bash ${SCRIPT} :::: ${PAIRS}
