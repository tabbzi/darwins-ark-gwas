#!/bin/bash
#$ -q broad
#$ -l h_vmem=16g
#$ -l h_rt=2:00:00
#$ -o /seq/vgb/dd/gwas/logs/
#$ -e /seq/vgb/dd/gwas/logs/
#$ -M kmorrill@broadinstitute.org
#$ -m e

gcta=/seq/vgb/software/gcta/current

DIR='/seq/vgb/dd/gwas'
DATE=${DATE:-'2022-11-20'}
if [ ${DATE} == "NA" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
fi

mkdir ${DIR}'/fastBAT/'${DATE}

# genotypes:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465'}

# phenotypes:
PHEDIR=${PHEDIR:-${DIR}'/pheno/'${DATE}}
PHE=${PHE:-'mq.121'}
PHEFILE=${PHEFILE:-${PHEDIR}'/'${PHE}'.tsv'}

# covariates (discrete):
DCOVDIR=${DCOVDIR:-${DIR}'/covar/'${DATE}'/dcov'}
DCOV=${DCOV:-'datatype'}
DCOVFILE=${DCOVFILE:-${DCOVDIR}'/'${DCOV}'.tsv'}
if [ ${DCOV} == 'NA' ] ; then
  echo "Not using any discrete covariates"
  ARG_DCOVAR=''
else
  ARG_DCOVAR='--covar '${DCOVFILE}
fi

# covariates (quantitative):
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

# output:
if [ ${useRELCUT} == 'T' ] ; then
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}'_rel-cutoff-'${RELCUT}}
else
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}}
fi

# get sample size from REML
n=`grep ^"n" ${DIR}'/reml/'${DATE}'/'${OUTPUT}'.REML.'*'.hsq' | awk '{print $2}' | sort | uniq`

# make COJO summary statistics file
echo -e "SNP\tA1\tA2\tfreq\tBETA\tSE\tP\tN" > ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'
tail -n+2 ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' | awk -F "\t" '$8!="inf"{print $0}' | awk -v n=${n} -F "\t" 'OFS=FS {print $2,$4,$5,$6,$7,$8,$9,n}' >> ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'

# perform fastBAT on gene set
${gcta} --bfile ${DIR}'/geno/'${GENO} \
        --fastBAT ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
        --fastBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.genes.dog.bed' \
        --out ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.genes'
paste <(echo -e "PHE\tQCOV\tDCOV\tSET") <(head -n1 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.genes.gene.fastbat') > ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.genes.gene.fastbat.tsv'
tail -n+2 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.genes.gene.fastbat' | awk -F "\t" -v P=${PHE} -v Q=${QCOV} -v D=${DCOV} 'OFS=FS {print P,Q,D,"gene",$0}' >> ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.genes.gene.fastbat.tsv'

# perform fastBAT on UNICORN set
${gcta} --bfile ${DIR}'/geno/'${GENO} \
        --fastBAT ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
        --fastBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.UNICORNs.dog.bed' \
        --out ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.UNICORNs'
paste <(echo -e "PHE\tQCOV\tDCOV\tSET") <(head -n1 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.UNICORNs.gene.fastbat') > ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.UNICORNs.gene.fastbat.tsv'
tail -n+2 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.UNICORNs.gene.fastbat' | awk -F "\t" -v P=${PHE} -v Q=${QCOV} -v D=${DCOV} 'OFS=FS {print P,Q,D,"UNICORN",$0}' >> ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.UNICORNs.gene.fastbat.tsv'

# perform fastBAT on RoCC set
${gcta} --bfile ${DIR}'/geno/'${GENO} \
        --fastBAT ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
        --fastBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.RoCCs.dog.bed' \
        --out ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.RoCCs'
paste <(echo -e "PHE\tQCOV\tDCOV\tSET") <(head -n1 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.RoCCs.gene.fastbat') > ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.RoCCs.gene.fastbat.tsv'
tail -n+2 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.RoCCs.gene.fastbat' | awk -F "\t" -v P=${PHE} -v Q=${QCOV} -v D=${DCOV} 'OFS=FS {print P,Q,D,"RoCC",$0}' >> ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.RoCCs.gene.fastbat.tsv'

# perform fastBAT on cCRE set
${gcta} --bfile ${DIR}'/geno/'${GENO} \
        --fastBAT ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma' \
        --fastBAT-gene-list '/seq/vgb/dd/gwas/annotate/regions/cf3/orthologous.map.cCREs.dog.bed' \
        --out ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.cCREs'
paste <(echo -e "PHE\tQCOV\tDCOV\tSET") <(head -n1 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.cCREs.gene.fastbat') > ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.cCREs.gene.fastbat.tsv'
tail -n+2 ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.cCREs.gene.fastbat' | awk -F "\t" -v P=${PHE} -v Q=${QCOV} -v D=${DCOV} 'OFS=FS {print P,Q,D,"cCRE",$0}' >> ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.region-test.cCREs.gene.fastbat.tsv'


rm ${DIR}'/fastBAT/'${DATE}'/'${OUTPUT}'.loco.mlma.ma'
