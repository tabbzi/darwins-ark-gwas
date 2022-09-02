#!/bin/bash
#$ -q broad
#$ -l h_vmem=20g
#$ -l h_rt=24:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
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

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/seq/vgb/dd/ancestry/env/

gemma=/seq/vgb/software/gemma/current
gcta=/seq/vgb/software/gcta/current
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

# USAGE: qsub GCTA_GWAS.job geno-prefix grm-prefix pheno-prefix covar-prefix qcovar-prefix assoc? clump? plot? directory reml?

DIR='/seq/vgb/dd/gwas'
GENO=${GENO:-'DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt'}
GRM=${GRM:-'DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt'}
PHENO=${PHENO:-'Q121A'}
PHENODIR=${PHENODIR:-${DIR}'/pheno'}
COVAR=${COVAR:-'DarwinsArk'}

QCOVAR=${QCOVAR:-'NA'}
if [ ${QCOVAR} == "NA" ] ; then
  echo "Not using a quantitative covariate"
  useQCOVAR=F
else
  useQCOVAR=T
fi
QCOVARDIR=${QCOVARDIR:-${DIR}'/covar'}

ASSOC=${ASSOC:-T}
CLUMP=${CLUMP:-T}
PLOT=${PLOT:-F}

DATE=${DATE:-'current'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

# run REML SNP-based heritability?
REML=${REML:-F}

# run COJO?
COJO=${COJO:-F}

# run --mlma-loco as well? will save to .loco.mlma
LOCO=${LOCO:-F}

# run --mlma-no-preadj-covar as well? will save to _no-preadj.mlma
PREADJ=${PREADJ:-F}

mkdir ${DIR}/assoc/${DATE}

# Overview Statistics
# Which individuals going into this GWAS with no missing values?
# use PLINK to get sample counts per variant
if [ ${useQCOVAR} == T ] ; then
  grep -wf <(grep -wf <(awk '$3!="NA" {print $1,$2}' ${DIR}/covar/${COVAR}.covar) <(awk '$3!="NA" {print $1,$2}' ${QCOVARDIR}/${QCOVAR}.qcovar)) <(awk '$3!="NA" {print $1,$2}' ${PHENODIR}/${PHENO}.pheno) > ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_keep.txt
else
  grep -wf <(awk '$3!="NA" {print $1,$2}' ${DIR}/covar/${COVAR}.covar) <(awk '$3!="NA" {print $1,$2}' ${PHENODIR}/${PHENO}.pheno) > ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_keep.txt
fi

# ${plink2} --dog --bfile ${DIR}/geno/${GENO} --keep ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_keep.txt --missing variant-only --out ${DIR}/assoc/${DATE}/${GENO}_${PHENO}

# GCTA-GREML: restricted maximum likelihood to estimate SNP-based heritability
if [ ${REML} == T ]; then
  echo "Performing GCTA-GREML to estimate variance explained for phenotype ${PHENO} by SNPs using relationship matrix of ${GRM}"

  mkdir ${DIR}/reml/${DATE}

  if [ ${useQCOVAR} == T ] ; then
    ${gcta} --thread-num 10 --grm ${DIR}/grm/${GRM} --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR} --covar ${DIR}/covar/${COVAR}.covar --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar
  else
    ${gcta} --thread-num 10 --grm ${DIR}/grm/${GRM} --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR} --covar ${DIR}/covar/${COVAR}.covar
  fi
fi

# GCTA-MLMA: mixed linear model-based association
if [ ${ASSOC} == T ]; then
  echo "Performing GCTA-MLMA for phenotype ${PHENO} across genotypes of ${GENO}"

  if [ $useQCOVAR == T ] ; then
    ${gcta} --bfile ${DIR}/geno/${GENO} \
            --thread-num 10 \
            --mlma \
            --autosome-num 38 \
            --autosome \
            --grm ${DIR}/grm/${GRM} \
            --pheno ${PHENODIR}/${PHENO}.pheno \
            --covar ${DIR}/covar/${COVAR}.covar \
            --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar \
            --out ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}
  else
    ${gcta} --bfile ${DIR}/geno/${GENO} \
            --thread-num 10 \
            --mlma \
            --autosome-num 38 \
            --autosome \
            --grm ${DIR}/grm/${GRM} \
            --pheno ${PHENODIR}/${PHENO}.pheno \
            --covar ${DIR}/covar/${COVAR}.covar \
            --out ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}
  fi

fi

# GCTA-MLMA-LOCO: needs far more memory to run
if [ ${LOCO} == T ]; then
  echo "Performing GCTA-MLMA for phenotype ${PHENO} across genotypes of ${GENO}"

  if [ $useQCOVAR == T ] ; then
    ${gcta} --bfile ${DIR}/geno/${GENO} \
            --thread-num 10 \
            --mlma-loco \
            --autosome-num 38 \
            --autosome \
            --grm ${DIR}/grm/${GRM} \
            --pheno ${PHENODIR}/${PHENO}.pheno \
            --covar ${DIR}/covar/${COVAR}.covar \
            --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar \
            --out ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}
  else
    ${gcta} --bfile ${DIR}/geno/${GENO} \
            --thread-num 10 \
            --mlma-loco \
            --autosome-num 38 \
            --autosome \
            --grm ${DIR}/grm/${GRM} \
            --pheno ${PHENODIR}/${PHENO}.pheno \
            --covar ${DIR}/covar/${COVAR}.covar \
            --out ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}
  fi
fi

if [ ${CLUMP} == T ]; then

# CLUMP association results with stringent significance threshold (5e-8)

echo "Performing PLINK clumping at stringent threshold (p = 5e-8)"

mkdir ${DIR}/clump/${DATE}

export DIR=${DIR}
export GENO=${GENO}
export DATE=${DATE}
export PHENO=${PHENO}
export OUT=${GENO}_${PHENO}_qcov-${QCOVAR}

export KB=1000
export Rsq=0.2
export p=5e-08
${DIR}/bin/clump.sh

echo "Performing PLINK clumping at loose threshold (p = 5e-6)"

export KB=1000
export Rsq=0.2
export p=1e-06
${DIR}/bin/clump.sh

rm ${DIR}/clump/${DATE}/*.log
rm ${DIR}/clump/${DATE}/*.nosex

fi

if [ ${PLOT} == T ]; then
  # use Anaconda
  PATH=$PATH:/seq/vgb/software/anaconda3/current
  source activate /seq/vgb/dd/gwas/env/Renv-2021
  mkdir ${DIR}/plots/${DATE}
  cd ${DIR}/plots/${DATE}
  if [ -f /seq/vgb/dd/gwas/clump/${DATE}/${GENO}_${PHENO}_p-5e-08_kb-1000_r2-0.2.clumped.ranges ]; then
    Rscript ${DIR}/plots/GWAS_plots.R ${DIR} ${DATE} ${GENO} ${PHENO} ${COVAR} ${QCOVAR} 5e-08 1000 0.2
  else
    Rscript ${DIR}/plots/GWAS_plots.R ${DIR} ${DATE} ${GENO} ${PHENO} ${COVAR} ${QCOVAR} 1e-06 1000 0.2
  fi

  if [ -f /seq/vgb/dd/gwas/clump/${DATE}/${GENO}_${PHENO}_p-5e-08_kb-1000_r2-0.2.loco.clumped.ranges ]; then
    Rscript ${DIR}/plots/GWAS_plots_loco.R ${DIR} ${DATE} ${GENO} ${PHENO} ${COVAR} ${QCOVAR} 5e-08 1000 0.2
  else
    Rscript ${DIR}/plots/GWAS_plots_loco.R ${DIR} ${DATE} ${GENO} ${PHENO} ${COVAR} ${QCOVAR} 1e-06 1000 0.2
  fi

fi

# GCTA-COJO: multi-SNP-based conditional & joint association analysis using GWAS summary data
if [ ${COJO} == T ]; then
  echo "Performing GCTA-COJO on results: ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.mlma"

  echo "SNP A1 A2 freq b se p N" > ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.ma
  sampcount=$(< "${DIR}/assoc/${DATE}/${GENO}_${PHENO}_keep.txt" wc -l)
  tail -n+2 ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.mlma | awk -v s=${sampcount} '{print $2,$4,$5,$6,$7,$8,$9,s}' >> ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_${Q.ma

  ${gcta} --bfile ${DIR}/geno/${GENO} \
          --cojo-file ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.ma \
          --cojo-slct \
          --autosome-num 38 \
          --autosome \
          --out ${DIR}/assoc/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}
fi
