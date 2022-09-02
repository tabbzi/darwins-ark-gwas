#!/bin/bash
#$ -q broad
#$ -l h_vmem=5g
#$ -l h_rt=24:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
#$ -M kmorrill@broadinstitute.org
#$ -m e
#$ -pe smp 5
#$ -binding linear:5
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

# USAGE: qsub GCTA_GWAS.job geno-prefix grm-prefix pheno-prefix covar-prefix qcovar-prefix assoc? clump? plot? directory reml?

DIR='/seq/vgb/dd/gwas/'
GENO=${1}
GRM=${2}
PHENO=${3}
COVAR=${4}
if [ $5 == "NA" ] ; then
  echo "Not using a quantitative covariate"
  useQCOVAR=F
  QCOVAR="NA"
else
  QCOVAR=${5}
  useQCOVAR=T
fi
ASSOC=${6:=T}
CLUMP=${7:=T}
PLOT=${8:=F}
# TODAY=`date +%Y-%m-%d`
# DATE=${9:=${TODAY}}
if [ -z "$9" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
else
  DATE=${9}
fi
REML=${10:=F}
echo ${REML}

if [ ${REML} == T ]; then
if [ ${useQCOVAR} == T ] ; then
${gcta} --thread-num 10 --grm /seq/vgb/dd/gwas/grm/${GENO} --pheno /seq/vgb/dd/gwas/pheno/${PHENO}.pheno --mpheno 1 --reml --out /seq/vgb/dd/gwas/reml/${PHENO} --covar /seq/vgb/dd/gwas/covar/${COVAR}.covar --qcovar /seq/vgb/dd/gwas/covar/${QCOVAR}.qcovar
else
${gcta} --thread-num 10 --grm /seq/vgb/dd/gwas/grm/${GENO} --pheno /seq/vgb/dd/gwas/pheno/${PHENO}.pheno --mpheno 1 --reml --out /seq/vgb/dd/gwas/reml/${PHENO} --covar /seq/vgb/dd/gwas/covar/${COVAR}.covar
fi
fi

cd ${DIR}

if [ ${ASSOC} == T ]; then

# PERFORM Mixed Linear Model-based Association in GCTA

  echo "Performing GCTA-MLMA for phenotype ${PHENO}.pheno across genotypes ${GENO}"

  mkdir ./assoc/${DATE}

  if [ $useQCOVAR == T ] ; then
    ${gcta} --bfile ./geno/${GENO} \
        --thread-num 10 \
        --mlma \
        --autosome-num 38 \
        --autosome \
        --grm ./grm/${GRM} \
        --pheno ./pheno/${PHENO}.pheno \
        --covar ./covar/${COVAR}.covar \
        --qcovar ./covar/${QCOVAR}.qcovar \
        --out ./assoc/${DATE}/${GENO}_${PHENO}
  else
    ${gcta} --bfile ./geno/${GENO} \
        --thread-num 10 \
        --mlma \
        --autosome-num 38 \
        --autosome \
        --grm ./grm/${GRM} \
        --pheno ./pheno/${PHENO}.pheno \
        --covar ./covar/${COVAR}.covar \
        --out ./assoc/${DATE}/${GENO}_${PHENO}
  fi

fi

if [ ${CLUMP} == T ]; then

# CLUMP association results with stringent significance threshold (5e-8)

echo "Performing PLINK clumping at stringent threshold (p = 5e-8)"

mkdir ./clump/${DATE}

# $plink --bfile ./geno/${GENO} \
#      --clump ./assoc/${DATE}/${GENO}_${PHENO}.mlma \
#      --dog \
#      --allow-no-sex \
#      --clump-field p \
#      --clump-p1 5e-8 \
#      --clump-kb 500 \
#      --clump-range ./annotate/Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.bed \
#      --out ./clump/${DATE}/${GENO}_${PHENO}_strng

export DIR=${DIR}
export GENO=${GENO}
export DATE=${DATE}
export PHENO=${PHENO}

export KB=1000
export Rsq=0.2
export p=5e-08
./bin/clump.sh 

echo "Performing PLINK clumping at loose threshold (p = 5e-6)"

# CLUMP	association results with loose significance threshold (5e-6)
# $plink --bfile ./geno/${GENO} \
#      --clump ./assoc/${DATE}/${GENO}_${PHENO}.mlma \
#      --dog \
#      --allow-no-sex \
#      --pheno ./pheno/${PHENO}.pheno \
#      --clump-field p \
#      --clump-p1 1e-6 \
#      --clump-kb 500 \
#      --clump-range ./annotate/Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.bed \
#      --out ./clump/${DATE}/${GENO}_${PHENO}_loose

export KB=1000
export Rsq=0.2
export p=1e-06
./bin/clump.sh 

rm ./clump/${DATE}/*.log
rm ./clump/${DATE}/*.nosex

#_p-${p}_kb-${KB}_r2-${Rsq}

#tr -s ' ' '\t' < ./clump/${DATE}/${GENO}_${PHENO}_strng.clumped > ./clump/${DATE}/${GENO}_${PHENO}_p-5e-08_kb-1000_r2-0.5.clumped.tab
#tr -s ' ' '\t' < ./clump/${DATE}/${GENO}_${PHENO}_loose.clumped > ./clump/${DATE}/${GENO}_${PHENO}_p-1e-06_kb-1000_r2-0.5.clumped.tab
#sed 's/(1)//g' < ./clump/${DATE}/${GENO}_${PHENO}_strng.clumped.tab > ./clump/${DATE}/${GENO}_${PHENO}_p-5e-08_kb-1000_r2-0.5.clumped.tab.rm
#sed 's/(1)//g' < ./clump/${DATE}/${GENO}_${PHENO}_loose.clumped.tab > ./clump/${DATE}/${GENO}_${PHENO}_p-1e-06_kb-1000_r2-0.5.clumped.tab.rm

#rm ${GENO}_${PHENO}*.clumped
#rm ${GENO}_${PHENO}*.clumped.tab

fi

if [ ${PLOT} == T ]; then
  # use Anaconda
  PATH=$PATH:/seq/vgb/software/anaconda3/current
  source activate /seq/vgb/dd/gwas/env/Renv-2021
  mkdir ./plots/${DATE}
  cd ./plots/${DATE}
  if [ -f /seq/vgb/dd/gwas/clump/${DATE}/${GENO}_${PHENO}_p-5e-08_kb-1000_r2-0.2.clumped.ranges ]; then
    Rscript /seq/vgb/dd/gwas/plots/GWAS_plots.R ${DIR} ${DATE} ${GENO} ${PHENO} ${COVAR} ${QCOVAR} 5e-08 1000 0.2
  else
    Rscript /seq/vgb/dd/gwas/plots/GWAS_plots.R ${DIR} ${DATE} ${GENO} ${PHENO} ${COVAR} ${QCOVAR} 1e-06 1000 0.2
  fi
fi


