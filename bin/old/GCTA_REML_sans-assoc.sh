#!/bin/bash
#$ -q broad
#$ -l h_vmem=15g
#$ -l h_rt=12:00:00
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

DIR='/seq/vgb/dd/gwas'
GENO=${GENO:-'DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt'}
PHENO=${PHENO:-'Q121A'}
PHENODIR=${PHENODIR:-${DIR}'/pheno'}
COVAR=${COVAR:-'DarwinsArk'}

p=${p:-"1e-6"}
KB=${KB:-"250"}
Rsq=${Rsq:-"0.2"}

QCOVAR=${QCOVAR:-'NA'}
if [ ${QCOVAR} == "NA" ] ; then
  echo "Not using a quantitative covariate"
  useQCOVAR=F
else
  useQCOVAR=T
fi
QCOVARDIR=${QCOVARDIR:-${DIR}'/covar'}

DATE=${DATE:-'mutt_paper'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

if [ ! -f ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXCLUDED.grm.N.bin ]; then
/seq/vgb/dd/gwas/bin/clumps2bed_noAnnot.sh ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges

${plink2} --dog --bfile ${DIR}/geno/${GENO} --extract bed1 ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed --write-snplist --out ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed

${gcta} --autosome-num 38 --autosome --bfile ${DIR}/geno/${GENO} --exclude ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.snplist --make-grm --out ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXCLUDED
${gcta} --autosome-num 38 --autosome --bfile ${DIR}/geno/${GENO} --extract ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.snplist --make-grm --out ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXTRACTED
fi

# ${gcta} --thread-num 10 --grm ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXCLUDED --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_REML_exclude-clumps --covar ${DIR}/covar/${COVAR}.covar --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar
# ${gcta} --thread-num 10 --grm ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXTRACTED --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_REML_extract-clumps --covar ${DIR}/covar/${COVAR}.covar --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar

# ${gcta} --thread-num 10 --grm ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXCLUDED --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_REML_exclude-clumps_highly-admixed --covar ${DIR}/covar/${COVAR}.covar --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar --keep /seq/vgb/dd/gwas/samples/highly-admixed.txt
# ${gcta} --thread-num 10 --grm ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXTRACTED --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_REML_extract-clumps_highly-admixed --covar ${DIR}/covar/${COVAR}.covar --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar --keep /seq/vgb/dd/gwas/samples/highly-admixed.txt

echo "${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXTRACTED" > ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.txt
echo "${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.EXCLUDED" >> ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.txt
if [ ${QCOVAR} == "NA" ] ; then
${gcta} --thread-num 10 --mgrm ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.txt --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_REML_by-assoc --covar ${DIR}/covar/${COVAR}.covar
else
${gcta} --thread-num 10 --mgrm ${DIR}/clump/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_p-${p}_kb-${KB}_r2-${Rsq}.loco.clumped.ranges.comp.noAnnot.bed.txt --pheno ${PHENODIR}/${PHENO}.pheno --mpheno 1 --reml --out ${DIR}/reml/${DATE}/${GENO}_${PHENO}_qcov-${QCOVAR}_REML_by-assoc --covar ${DIR}/covar/${COVAR}.covar --qcovar ${QCOVARDIR}/${QCOVAR}.qcovar
fi
