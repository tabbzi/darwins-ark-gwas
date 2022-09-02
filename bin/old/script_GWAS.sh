# USAGE: ./bin/script_GWAS.sh geno-prefix grm-prefix pheno-prefix covar-prefix qcovar-prefix assoc? clump? plot?

gcta=/seq/vgb/kmorrill/bin/gcta_1.91.7beta/gcta64

GENO=${1}
GRM=${2}
PHENO=${3}
COVAR=${4}
if [ $5 == "NA" ] ; then
  echo "Not using a quantitative covariate"
  useQCOVAR=F
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

if [ ${ASSOC} == T ]; then

# PERFORM Mixed Linear Model-based Association in GCTA

  echo "Performing GCTA-MLMA for phenotype ${PHENO}.pheno across genotypes ${GENO}"

  mkdir ./assoc/${DATE}

  if [ $useQCOVAR == T ] ; then
    ${gcta} --bfile ./geno/${GENO} \
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

# CLUMP association results with stringent significance threshold (6e-8)

echo "Performing PLINK clumping at stringent threshold (p = 6e-8)"

mkdir ./clump/${DATE}

plink --bfile ./geno/${GENO} \
      --clump ./assoc/${DATE}/${GENO}_${PHENO}.mlma \
      --dog \
      --allow-no-sex \
      --clump-field p \
      --clump-p1 6e-8 \
      --clump-range ./annotate/dog_genes_from_ensembl.bed4 \
      --out ./clump/${DATE}/${GENO}_${PHENO}_strng

echo "Performing PLINK clumping at loose threshold (p = 5e-6)"

# CLUMP	association results with loose significance threshold (5e-6)
plink --bfile ./geno/${GENO} \
      --clump ./assoc/${DATE}/${GENO}_${PHENO}.mlma \
      --dog \
      --allow-no-sex \
      --pheno ./pheno/${PHENO}.pheno \
      --clump-field p \
      --clump-p1 5e-6 \
      --clump-range ./annotate/dog_genes_from_ensembl.bed4 \
      --out ./clump/${DATE}/${GENO}_${PHENO}_loose

cd ./clump/${DATE}/

rm *.log
rm *.nosex

tr -s ' ' '\t' < ${GENO}_${PHENO}_strng.clumped > ${GENO}_${PHENO}_strng.clumped.tab
tr -s ' ' '\t' < ${GENO}_${PHENO}_loose.clumped > ${GENO}_${PHENO}_loose.clumped.tab
sed 's/(1)//g' < ${GENO}_${PHENO}_strng.clumped.tab > ${GENO}_${PHENO}_strng.clumped.tab.rm
sed 's/(1)//g' < ${GENO}_${PHENO}_loose.clumped.tab > ${GENO}_${PHENO}_loose.clumped.tab.rm

rm ${GENO}_${PHENO}*.clumped
rm ${GENO}_${PHENO}*.clumped.tab
fi

if [ ${PLOT} == T ]; then

  cd /seq/vgb/dd/gwas/plots/

  mkdir ./${DATE}

  i=1

  until [ -e ${GENO}_${PHENO}_noclump.qq.png ] ; do

    if [[ -f ../clump/${DATE}/${GENO}_${PHENO}_strng.clumped.ranges ]]; then
      echo "Generating R plots with stringent clumps colored and annotated"
      Rscript ./GWAS-plots_standard.r ../assoc/${DATE}/${GENO}_${PHENO}.mlma ../clump/${DATE}/${GENO}_${PHENO}_strng.clumped.tab.rm ../clump/${DATE}/${GENO}_${PHENO}_strng.clumped.ranges ${GENO} ${PHENO}_strng
    fi

    echo "Generating R plots without clumps"
    Rscript ./GWAS-plots_standard_noclump.r ../assoc/${DATE}/${GENO}_${PHENO}.mlma ${GENO} ${PHENO}_noclump

    i=`expr $i + 1`

    if [ ${i} == 3 ]; then
      echo "Plots failed to generate"
      break 1
    else
      echo "Tried generating plots ${i} times so far"
    fi
  done

  mv ${GENO}*${PHENO}*.png ./${DATE}/.

fi
