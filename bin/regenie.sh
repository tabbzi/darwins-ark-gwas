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

# software:
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

# regenie:
PATH=$PATH:/seq/vgb/software/anaconda3/current
conda activate /seq/vgb/software/regenie/env

# prepare:
DIR='/seq/vgb/dd/gwas'
DATE=${DATE:-'2022-11-20'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

mkdir ${DIR}'/reml/'${DATE}
mkdir ${DIR}'/assoc/'${DATE}
mkdir ${DIR}'/clump/'${DATE}
mkdir ${DIR}'/plots/'${DATE}

# genetic dataset:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465'}

# options:
REML=${REML:-'T'}       # run GCTA-REML?
MODEL=${LOCO:-'T'}       # run GCTA-MLMA-LOCO? efficiency -, power +
PREADJ=${PREADJ:-'T'}   # pre-adjust by covariates? efficiency +, power -
PLOT=${PLOT:-'T'}       # generate plots?
MBAT=${MBAT:-'T'}       # run mBAT-combo?

# chromosome range:
START_CHR=${START_CHR:-'1'}
END_CHR=${END_CHR:-'38'}

# phenotypes:
PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
PHE=${PHE:-'mq.121'}
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
QCOV=${QCOV:-'age.'${PHE}}
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

# clump settings:
CLUMP=${CLUMP:-'T'}
CLUMPkb=${CLUMPkb:-"1000"}
CLUMPr2=${CLUMPr2:-"0.2"}

# preadj:
if [ ${PREADJ} == 'T' ] ; then
  ARG_PREADJ=''
else
  ARG_PREADJ='--mlma-no-preadj-covar'
fi

# output:
if [ ${useRELCUT} == 'T' ] ; then
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}'_rel-cutoff-'${RELCUT}}
else
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}}
fi

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

if [ ${QCOV} == 'NA' ] ; then
  grep -wf <(awk '$3!="NA" {print $1,$2}' ${DCOVFILE}) <(awk '$3!="NA" {print $1,$2}' ${PHEFILE}) | wc -l >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
else
  grep -wf <(grep -wf <(awk '$3!="NA" {print $1,$2}' ${DCOVFILE}) <(awk '$3!="NA" {print $1,$2}' ${QCOVFILE})) <(awk '$3!="NA" {print $1,$2}' ${PHEFILE}) | wc -l >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.info.txt'
fi

# GCTA-GREML: restricted maximum likelihood to estimate SNP-based heritability
if [ ${REML} == 'T' ] ; then
  if [ ${useRELCUT} == 'T' ] ; then
    INPUT_GRM_ALL=${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT}
    INPUT_MGRM_LIST=${DIR}'/grm/'${GENO}'.score.ld.rel-cutoff-'${RELCUT}'.list.txt'
  else
    INPUT_GRM_ALL=${DIR}'/grm/'${GENO}
    INPUT_MGRM_LIST=${DIR}'/grm/'${GENO}'.score.ld.list.txt'
  fi

  echo "Performing GCTA-GREML with constraint..."
  ${gcta} --grm ${INPUT_GRM_ALL} \
          --reml \
          --pheno ${PHEFILE} \
          --mpheno ${PN} \
          ${ARG_DCOVAR} \
          ${ARG_QCOVAR} \
          --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.no-lds'

  echo "Performing GCTA-GREML without constraint..."
  ${gcta} --grm ${INPUT_GRM_ALL} \
          --reml \
          --reml-no-constrain \
          --pheno ${PHEFILE} \
          --mpheno ${PN} \
          ${ARG_DCOVAR} \
          ${ARG_QCOVAR} \
          --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.no-lds.no-constraint'

  echo "Performing GCTA-GREML-LDMS with constraint..."
  ${gcta} --mgrm ${INPUT_MGRM_LIST} \
          --reml \
          --pheno ${PHEFILE} \
          --mpheno ${PN} \
          ${ARG_DCOVAR} \
          ${ARG_QCOVAR} \
          --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.lds'

  echo "Performing GCTA-GREML-LDMS without constraint..."
  ${gcta} --mgrm ${INPUT_MGRM_LIST} \
          --reml \
          --reml-no-constrain \
          --pheno ${PHEFILE} \
          --mpheno ${PN} \
          ${ARG_DCOVAR} \
          ${ARG_QCOVAR} \
          --out ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.lds.no-constraint'
fi

#
if [ ${LOCO} == 'T' ] ; then
  echo "Performing GCTA-MLMA-LOCO..."

  # Generate range of chromosomes:
  CHR_RANGE=`seq ${START_CHR} ${END_CHR}`

  # Set the input GRM files based on whether the useRELCUT variable is set:
  if [ ${useRELCUT} == 'T' ] ; then
    INPUT_GRM_ALL=${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT}
    INPUT_GRM_SUB=${DIR}'/grm/'${GENO}'_chr-{1}_rel-cutoff-'${RELCUT}
  else
    INPUT_GRM_ALL=${DIR}'/grm/'${GENO}
    INPUT_GRM_SUB=${DIR}'/grm/'${GENO}'_chr-{1}_rel-cutoff-'${RELCUT}
  fi

  # Run the GCTA command in parallel, using all available cores
  echo "Running chromosomes ${CHR_RANGE}"
  parallel --jobs 4 ${gcta} \
          --bfile ${DIR}'/geno/'${GENO} \
          --grm ${INPUT_GRM_ALL} \
          --mlma-subtract-grm ${INPUT_GRM_SUB} \
          --mlma \
          --autosome \
          --autosome-num 38 \
          --chr {1} \
          --thread-num 4 \
          --pheno ${PHEFILE} \
          --mpheno ${PN} \
          ${ARG_DCOVAR} \
          ${ARG_QCOVAR} \
          ${ARG_PREADJ} \
          --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-{1}' ::: ${CHR_RANGE}

  # Wait for all the parallel jobs to finish:
  wait

  # Extract the first line (column headers) from the first file:
  RANGE=(${CHR_RANGE})
  head -n1 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${RANGE[0]}'.mlma' > ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma'

  # Concatenate the first line with the rest of the files
  for CHR in ${CHR_RANGE}
  do
    tail -n+2 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'${CHR}'.mlma' >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma'
  done

  # Remove the input files:
  rm ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'_chr-'*'.mlma'
fi

if [ ${CLUMP} == 'T' ] ; then

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
export p=1e-06
${DIR}/bin/PLINK_clump.sh

rm ${DIR}/clump/${DATE}/*.log
rm ${DIR}/clump/${DATE}/*.nosex

fi

if [ ${PLOT} == 'T' ] ; then
  # use Anaconda
  PATH=$PATH:/seq/vgb/software/anaconda3/current
  # source activate /seq/vgb/dd/gwas/env/Renv
  source activate /seq/vgb/dd/gwas/env/Renv-2023
  Rscript ${DIR}'/bin/manhattan-plot.R' --mlma ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' --clump ${DIR}'/clump/'${DATE}'/'${OUTPUT}'_p-'${p}'_kb-'${KB}'_r2-'${Rsq}'.loco.clumped' --title ${PHE}
fi

if [ ${MBAT} == 'T' ]; then

  # get sample size from REML
  n=`grep ^"n" ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.'*'.hsq' | awk '{print $2}' | sort | uniq`

  # make COJO summary statistics file
  # SNP A1 A2      freq     BETA       SE         P     N
  echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN" > ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'
  tail -n+2 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' | awk -F "\t" '$8!="inf"{print $0}' | awk -v n=${n} -F "\t" 'OFS=FS {print $2,$4,$5,$6,$7,$8,$9,n}' >> ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'

  # perform mBAT-combo on gene set
  ${gcta} --bfile ${DIR}'/geno/'${GENO} \
          --mBAT-combo ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
          --mBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.genes.dog.bed' \
          --mBAT-print-all-p \
          --thread-num 4 \
          --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.region-test.genes'

  # perform mBAT-combo on UNICORN set
  ${gcta} --bfile ${DIR}'/geno/'${GENO} \
          --mBAT-combo ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
          --mBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.UNICORNs.dog.bed' \
          --mBAT-print-all-p \
          --thread-num 4 \
          --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.region-test.UNICORNs'

  # perform mBAT-combo on RoCC set
  ${gcta} --bfile ${DIR}'/geno/'${GENO} \
          --mBAT-combo ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
          --mBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.RoCCs.dog.bed' \
          --mBAT-print-all-p \
          --thread-num 4 \
          --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.region-test.RoCCs'

  # perform mBAT-combo on cCRE set
  ${gcta} --bfile ${DIR}'/geno/'${GENO} \
          --mBAT-combo ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
          --mBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.cCREs.dog.bed' \
          --mBAT-print-all-p \
          --thread-num 4 \
          --out ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.region-test.cCREs'

fi
