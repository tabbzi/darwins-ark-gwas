#!/bin/bash
#$ -q broad
#$ -l h_vmem=32g
#$ -l h_rt=8:00:00
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

mkdir ${DIR}/reml/${DATE}
mkdir ${DIR}/assoc/${DATE}
mkdir ${DIR}/clump/${DATE}
mkdir ${DIR}/plots/${DATE}

# genetic dataset:
GENO=${GENO:-'DarwinsArk_GeneticData_Aug2022'}

# options:
REML=${REML:-'T'}       # run GCTA-REML?
ASSOC=${ASSOC:-'F'}     # run GCTA-MLMA?
LOCO=${LOCO:-'F'}       # run GCTA-MLMA-LOCO? efficiency -, power +
LSPLIT=${LSPLIT:-'T'}   # run GCTA-MLMA-LOCO by splitting GRM? efficiency+
PREADJ=${PREADJ:-'F'}   # pre-adjust by covariates? efficiency +, power -
PLOT=${PLOT:-'T'}       # generate plots?

# phenotypes:
PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
PHE=${PHE:-'fa.no-na.1'}
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

# quantitative covariates:
QCOVDIR=${QCOVDIR:-${DIR}'/covar/'${DATE}'/qcov'}
QCOV=${QCOV:-'age.'${PHE}}
if [ ${QCOV} == 'NA' ] ; then
  echo "Not using a quantitative covariate"
  useQCOV=F
else
  useQCOV=T
fi
QCOVFILE=${QCOVFILE:-${QCOVDIR}'/'${QCOV}'.tsv'}

# relatedness cutoff:
RELCUT=${RELCUT:-"NA"}
if [ ${RELCUT} == 'NA' ] ; then
  echo "Not using a relatedness cutoff"
  useRELCUT=F
else
  echo "Using relatedness cutoff"
  useRELCUT=T
fi

# output:
OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}}

# clump settings:
CLUMP=${CLUMP:-T}
CLUMPkb=${CLUMPkb:-"1000"}
CLUMPr2=${CLUMPr2:-"0.2"}

# metadata:
echo "Genome-wide Association Study" >  ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
echo ${GENO} >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
echo ${DATE} >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
echo "PHE (phenotype): ${PHE}" >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
echo "DCOV (discrete covariates): ${DCOV}" >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
echo "QCOV (quantative covariates): ${QCOV}" >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
if [ ${PREADJ} == 'T' ] ; then
  echo "Phenotypes pre-adjusted by covariates for association tests." >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
else
  echo "Covariates fitted together with SNP in association tests." >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
fi
echo "# dogs included in analysis:" >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'

if [ ${useQCOV} == T ] ; then
  grep -wf <(grep -wf <(awk '$3!="NA" {print $1,$2}' ${DCOVFILE}) <(awk '$3!="NA" {print $1,$2}' ${QCOVFILE})) <(awk '$3!="NA" {print $1,$2}' ${PHEFILE}) | wc -l >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
else
  grep -wf <(awk '$3!="NA" {print $1,$2}' ${DCOVFILE}) <(awk '$3!="NA" {print $1,$2}' ${PHEFILE}) | wc -l >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
fi

if [ ${useRELCUT} == T ]; then
  # GCTA-GREML: restricted maximum likelihood to estimate SNP-based heritability
  if [ ${REML} == T ]; then
    echo "Performing GCTA-GREML..."
    if [ ${useQCOV} == T ] ; then
      ${gcta} --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.REML.no-lds'
              #--gxe ${GXEFILE} \
    else
      ${gcta} --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.REML.no-lds'
              #--gxe ${GXEFILE} \
    fi

    echo "Performing GCTA-GREML without constraint..."
    if [ ${useQCOV} == T ] ; then
      ${gcta} --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --reml \
              --reml-no-constrain \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.REML.no-lds.no-constraint'
              #--gxe ${GXEFILE} \
    else
      ${gcta} --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --reml \
              --reml-no-constrain \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.REML.no-lds.no-constraint'
              #--gxe ${GXEFILE} \
    fi

    echo "Performing GCTA-GREML-LDMS..."
    if [ ${useQCOV} == T ] ; then
      ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.rel-cutoff-'${RELCUT}'list.txt' \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.REML.lds'
              #--gxe ${GXEFILE} \
    else
      ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.rel-cutoff-'${RELCUT}'list.txt' \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.REML.lds'
              #--gxe ${GXEFILE} \
    fi
  fi

  # GCTA-MLMA: mixed linear model-based association
  if [ ${ASSOC} == T ]; then
    echo "Performing GCTA-MLMA..."

    if [ $useQCOV == T ] ; then
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --mlma \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}
    else
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --mlma \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}
    fi
  fi

  if [ ${LSPLIT} == T ]; then
    echo "Performing GCTA-MLMA-LOCO (chr splitting)..."

    if [ $useQCOV == T ] ; then
      for CHR in `seq 1 38`
      do
        ${gcta} --bfile ${DIR}'/geno/'${GENO} \
                --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
                --mlma-subtract-grm ${DIR}'/grm/'${GENO}'_chr-'${CHR}'_rel-cutoff-'${RELCUT} \
                --chr ${CHR} \
                --pheno ${PHEFILE} \
                --mpheno ${PN} \
                --covar ${DCOVFILE} \
                --qcovar ${QCOVFILE} \
                --mlma \
                --mlma-no-preadj-covar \
                --autosome-num 38 \
                --autosome \
                --thread-num 4 \
                --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}'_rel-cutoff-'${RELCUT}
      done
    else
      for CHR in `seq 1 38`
      do
        ${gcta} --bfile ${DIR}'/geno/'${GENO} \
                --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
                --mlma-subtract-grm ${DIR}'/grm/'${GENO}'_chr-'${CHR}'_rel-cutoff-'${RELCUT} \
                --chr ${CHR} \
                --pheno ${PHEFILE} \
                --mpheno ${PN} \
                --covar ${DCOVFILE} \
                --mlma \
                --mlma-no-preadj-covar \
                --autosome-num 38 \
                --autosome \
                --thread-num 4 \
                --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}'_rel-cutoff-'${RELCUT}
      done
    fi
    cat ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-1''_rel-cutoff-'${RELCUT}'.mlma' > ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.loco.mlma'
    for CHR in `seq 2 38`
    do
      tail -n+2 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}'_rel-cutoff-'${RELCUT}'.mlma' >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}'.loco.mlma'
    done
    rm {DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'*'_rel-cutoff-'${RELCUT}'.mlma'
  fi

  if [ ${LOCO} == T ]; then
    echo "Performing GCTA-MLMA-LOCO (no splitting)..."

    if [ $useQCOV == T ] ; then
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --mlma-loco \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}
    else
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --grm ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT} \
              --covar ${DCOVFILE} \
              --mlma-loco \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_rel-cutoff-'${RELCUT}
    fi
  fi
else
  # GCTA-GREML: restricted maximum likelihood to estimate SNP-based heritability
  if [ ${REML} == T ]; then
    echo "Performing GCTA-GREML..."
    if [ ${useQCOV} == T ] ; then
      ${gcta} --grm ${DIR}'/grm/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.no-lds'
              #--gxe ${GXEFILE} \
    else
      ${gcta} --grm ${DIR}'/grm/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.no-lds'
              #--gxe ${GXEFILE} \
    fi

    echo "Performing GCTA-GREML without constraint..."
    if [ ${useQCOV} == T ] ; then
      ${gcta} --grm ${DIR}'/grm/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --reml \
              --reml-no-constrain \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.no-lds.no-constraint'
              #--gxe ${GXEFILE} \
    else
      ${gcta} --grm ${DIR}'/grm/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --reml \
              --reml-no-constrain \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.no-lds.no-constraint'
              #--gxe ${GXEFILE} \
    fi

    echo "Performing GCTA-GREML-LDMS..."
    if [ ${useQCOV} == T ] ; then
      ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.list.txt' \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.lds'
              #--gxe ${GXEFILE} \
    else
      ${gcta} --mgrm ${DIR}'/grm/'${GENO}'.score.ld.list.txt' \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --reml \
              --thread-num 4 \
              --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.lds'
              #--gxe ${GXEFILE} \
    fi
  fi

  # GCTA-MLMA: mixed linear model-based association
  if [ ${ASSOC} == T ]; then
    echo "Performing GCTA-MLMA..."

    if [ $useQCOV == T ] ; then
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --grm ${DIR}'/grm/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --mlma \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}
    else
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --grm ${DIR}'/grm/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --covar ${DCOVFILE} \
              --mlma \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}
    fi
  fi

  if [ ${LSPLIT} == T ]; then
    echo "Performing GCTA-MLMA-LOCO (chr splitting)..."

    if [ $useQCOV == T ] ; then
      for CHR in `seq 1 38`
      do
        ${gcta} --bfile ${DIR}'/geno/'${GENO} \
                --grm ${DIR}'/grm/'${GENO} \
                --mlma-subtract-grm ${DIR}'/grm/'${GENO}'_chr-'${CHR} \
                --chr ${CHR} \
                --pheno ${PHEFILE} \
                --mpheno ${PN} \
                --covar ${DCOVFILE} \
                --qcovar ${QCOVFILE} \
                --mlma \
                --mlma-no-preadj-covar \
                --autosome-num 38 \
                --autosome \
                --thread-num 4 \
                --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}
      done
    else
      for CHR in `seq 1 38`
      do
        ${gcta} --bfile ${DIR}'/geno/'${GENO} \
                --grm ${DIR}'/grm/'${GENO} \
                --mlma-subtract-grm ${DIR}'/grm/'${GENO}'_chr-'${CHR} \
                --chr ${CHR} \
                --pheno ${PHEFILE} \
                --mpheno ${PN} \
                --covar ${DCOVFILE} \
                --mlma \
                --mlma-no-preadj-covar \
                --autosome-num 38 \
                --autosome \
                --thread-num 4 \
                --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}
      done
    fi
    cat ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-1.mlma' > ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma'
    for CHR in `seq 2 38`
    do
      tail -n+2 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}'.mlma' >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma'
    done
    rm {DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'*'.mlma'
  fi

  if [ ${LOCO} == T ]; then
    echo "Performing GCTA-MLMA-LOCO (no splitting)..."

    if [ $useQCOV == T ] ; then
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --grm ${DIR}'/grm/'${GENO} \
              --covar ${DCOVFILE} \
              --qcovar ${QCOVFILE} \
              --mlma-loco \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}
    else
      ${gcta} --bfile ${DIR}'/geno/'${GENO} \
              --pheno ${PHEFILE} \
              --mpheno ${PN} \
              --grm ${DIR}'/grm/'${GENO} \
              --covar ${DCOVFILE} \
              --mlma-loco \
              --mlma-no-preadj-covar \
              --autosome-num 38 \
              --autosome \
              --thread-num 4 \
              --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}
    fi
  fi
fi

if [ ${CLUMP} == T ]; then

# CLUMP association results with stringent significance threshold (5e-8)

echo "Performing PLINK clumping at stringent threshold (p = 5e-8)"

export DIR=${DIR}
export GENO=${GENO}
export DATE=${DATE}
export PHEFILE=${PHEFILE}
export OUTPUT=${OUTPUT}

export KB=${CLUMPkb}
export Rsq=${CLUMPr2}
export p=5e-08
${DIR}/bin/PLINK_clump.sh

echo "Performing PLINK clumping at loose threshold (p = 1e-6)"

export KB=${CLUMPkb}
export Rsq=${CLUMPr2}
export p=1e-6
${DIR}/bin/PLINK_clump.sh

rm ${DIR}/clump/${DATE}/*.log
rm ${DIR}/clump/${DATE}/*.nosex

fi

if [ ${PLOT} == T ]; then
  # use Anaconda
  PATH=$PATH:/seq/vgb/software/anaconda3/current
  source activate /seq/vgb/dd/gwas/env/Renv
  Rscript ${DIR}'/bin/manhattan-plot.R' ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' ${DIR}'/clump/'${DATE}'/'${OUTPUT}'_p-'${p}'_kb-'${KB}'_r2-'${Rsq}'.loco.clumped' <(echo "${PHE}")
fi
