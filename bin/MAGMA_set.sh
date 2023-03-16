#!/bin/bash
#$ -q broad
#$ -l h_vmem=4g
#$ -l h_rt=1:00:00
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
DATE=${DATE:-'2022-11-20'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

mkdir ${DIR}'/magma/'${DATE}

# genetic dataset:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465'}

# tests:
ANNODIR=${ANNODIR:-${DIR}'/magma/anno/'${DATE}}
ANNO=${ANNO:-'magma.orthologous.map.genes'}
ANNOFILE=${ANNOFILE:-${ANNODIR}'/'${GENO}'_'${ANNO}'.genes.annot'}

# sets:
SETDIR=${SETDIR:-${DIR}'/magma/sets'}
SET=${SET:-'gtex.set'}
SETFILE=${SETFILE:-${SETDIR}'/'${SET}}

# phenotypes:
PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
PHE=${PHE:-'fa.01-filled-by-mean'}
PHEFILE=${PHEFILE:-${PHEDIR}'/'${PHE}'.tsv'}
PN=${PN:-'1'}

# environmental factors:
#GXEDIR=${GXEDIR:-${DIR}'/covar'/${DATE}}
#GXE=${GXE:-'sex.psychmed'}
#GXEFILE=${GXEFILE:-${GXEDIR}'/'${GXE}'.tsv'}

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
QCOV=${QCOV:-'age.hgt.fa.01'}
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

# output:
if [ ${useRELCUT} == 'T' ] ; then
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}'_rel-cutoff-'${RELCUT}}
else
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}}
fi

# n:
n=`grep ^"n" ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.'*'.hsq' | awk '{print $2}' | sort | uniq`

echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN" > ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'
tail -n+2 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' | awk -F "\t" '$8!="inf"{print $0}' | awk -v n=${n} -F "\t" 'OFS=FS {print $2,$4,$5,$6,$7,$8,$9,n}' >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'

# SET BASED
${magma} --gene-results ${DIR}'/magma/'${DATE}'/'${OUTPUT}'.'${ANNO}'.genes.raw' \
         --set-annot ${SETFILE} \
         --out ${DIR}'/magma/'${DATE}'/'${OUTPUT}'.'${ANNO}'.'${SET}
