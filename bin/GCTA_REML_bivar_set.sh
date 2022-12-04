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
gemma=/seq/vgb/software/gemma/current
gcta=/seq/vgb/software/gcta/current
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

# prepare:
DIR='/seq/vgb/dd/gwas'
DATE=${DATE:-'2022-07-15'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

mkdir ${DIR}/cor/${DATE}

# genetic dataset:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.02_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3104'}

while read -r PA PB A B
do
  echo $PA
  echo $PB
  echo $A
  echo $B
  export PNA=${PA}
  export PNB=${PB}
  export A=${A}
  export B=${B}

  # phenotypes:
  PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
  PHE=${PHE:-'bq.all'}
  PHEFILE=${PHEFILE:-${PHEDIR}'/'${PHE}'.tsv'}

  # selections:
  KEY=${KEY:-${PHEDIR}'/'${PHE}'.key.txt'}
  PNA=${PNA:-'1'}
  PNB=${PNB:-'2'}

  # discrete covariates:
  DCOVDIR=${DCOVDIR:-${DIR}'/covar/'${DATE}'/dcov'}
  DCOV=${DCOV:-'datatype'}
  DCOVFILE=${DCOVFILE:-${DCOVDIR}'/'${DCOV}'.tsv'}

  # quantitative covariates:
  QCOVDIR=${QCOVDIR:-${DIR}'/covar/'${DATE}'/qcov'}
  QCOV=${QCOV:-'age'}
  if [ ${QCOV} == 'NA' ] ; then
    echo "Not using a quantitative covariate"
    useQCOV=F
  else
    useQCOV=T
  fi
  QCOVFILE=${QCOVFILE:-${QCOVDIR}'/'${QCOV}'.tsv'}

  # output:
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}'_'${A}'-vs-'${B}}

  # GCTA-GREML: restricted maximum likelihood to estimate SNP-based heritability
  echo "Performing bivariate GCTA-GREML analysis..."
  if [ ${useQCOV} == T ] ; then
    ${gcta} --grm ${DIR}'/grm/'${GENO} \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --qcovar ${QCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.no-lds'
  else
    ${gcta} --grm ${DIR}'/grm/'${GENO} \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.no-lds'
  fi

  echo "Performing bivariate GCTA-GREML without constraint..."
  if [ ${useQCOV} == T ] ; then
    ${gcta} --grm ${DIR}'/grm/'${GENO} \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --qcovar ${QCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --reml-bivar-no-constrain \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.no-lds.no-constraint'
  else
    ${gcta} --grm ${DIR}'/grm/'${GENO} \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --reml-bivar-no-constrain \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.no-lds.no-constraint'
  fi

  echo "Performing bivariate GCTA-GREML-LDMS..."
  if [ ${useQCOV} == T ] ; then
    ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.list.txt' \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --qcovar ${QCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.lds'
  else
    ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.list.txt' \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.lds'
  fi

  echo "Performing bivariate GCTA-GREML-LDMS..."
  if [ ${useQCOV} == T ] ; then
    ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.list.txt' \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --qcovar ${QCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --reml-bivar-no-constrain \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.lds.no-constraint'
  else
    ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.list.txt' \
            --pheno ${PHEFILE} \
            --covar ${DCOVFILE} \
            --reml-bivar ${PNA} ${PNB} \
            --reml-bivar-no-constrain \
            --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.REML.lds.no-constraint'
  fi
done < ${input}
