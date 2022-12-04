#!/bin/bash
#$ -q broad
#$ -l h_vmem=16g
#$ -l h_rt=2:00:00
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
magma=/seq/vgb/software/magma/current

# prepare:
DIR='/seq/vgb/dd/gwas'
DATE=${DATE:-'2022-07-15'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

mkdir ${DIR}/magma/${DATE}

# genetic dataset:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.005_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3104'}

# tests:
ANNODIR=${ANNODIR:-${DIR}'/magma/anno/'${DATE}}
ANNO=${ANNO:-'magma.zoo-RoCCs'}
ANNOFILE=${ANNOFILE:-${ANNODIR}'/'${GENO}'_'${ANNO}'.genes.annot'}

# phenotypes:
PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
PHE=${PHE:-'fa.04'}
PHEFILE=${PHEFILE:-${PHEDIR}'/'${PHE}'.tsv'}
PN=${PN:-'1'}

# discrete covariates:
DCOVDIR=${DCOVDIR:-${DIR}'/covar/'${DATE}'/dcov'}
DCOV=${DCOV:-'datatype'}
DCOVFILE=${DCOVFILE:-${DCOVDIR}'/'${DCOV}'.tsv'}

# quantitative covariates:
QCOVDIR=${QCOVDIR:-${DIR}'/covar/'${DATE}'/qcov'}
QCOV=${QCOV:-'age.hgt.'${PHE}}
if [ ${QCOV} == 'NA' ] ; then
  echo "Not using a quantitative covariate"
  useQCOV=F
else
  useQCOV=T
fi
QCOVFILE=${QCOVFILE:-${QCOVDIR}'/'${QCOV}'.tsv'}

# output:
OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}}

# n:
if [ ${useQCOV} == T ] ; then
  N=`grep -wf <(grep -wf <(awk '$3!="NA" {print $1,$2}' ${DCOVFILE}) <(awk '$3!="NA" {print $1,$2}' ${QCOVFILE})) <(awk '$3!="NA" {print $1,$2}' ${PHEFILE}) | sed -n '$='`
else
  N=`grep -wf <(awk '$3!="NA" {print $1,$2}' ${DCOVFILE}) <(awk '$3!="NA" {print $1,$2}' ${PHEFILE}) | sed -n '$='`
fi

# GENE / REGION BASED
${magma} --bfile ${DIR}'/geno/'${GENO} \
         --pval ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' N=${N} \
         --gene-annot ${ANNOFILE} \
         --out ${DIR}'/magma/'${DATE}'/'${OUTPUT}'.'${ANNO} \
         --gene-settings fixed-permp=5000
