#!/bin/bash
DIR=${DIR:-"/seq/vgb/dd/gwas"}
GENO=${GENO:-'DarwinsArk_GeneticData_Aug2022'}
DATE=${DATE:-'2022-07-15'}
KB=${KB:-250}
Rsq=${Rsq:-0.2}
p=${p:-"1e-6"}
ANNOT=${ANNOT:-"Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.withHuman.bed"}

plink=/seq/vgb/software/plink2/current/plink

echo "Starting clumping..."

if [ -f ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.mlma' ]; then
$plink --bfile ${DIR}'/geno/'${GENO} \
      --clump ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.mlma' \
      --dog \
      --allow-no-sex \
      --pheno ${PHEFILE} \
      --clump-field p \
      --clump-p1 ${p} \
      --clump-kb ${KB} \
      --clump-r2 ${Rsq} \
      --clump-range ${DIR}'/annotate/'${ANNOT} \
      --out ${DIR}'/clump/'${DATE}'/'${OUTPUT}'_p-'${p}'_kb-'${KB}'_r2-'${Rsq}
else
echo "No results file (non-loco)?"
fi

if [ -f ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' ]; then
  $plink --bfile ${DIR}'/geno/'${GENO} \
        --clump ${DIR}'/assoc/'${DATE}'/'${OUTPUT}'.loco.mlma' \
        --dog \
        --allow-no-sex \
        --pheno ${PHEFILE} \
        --clump-field p \
        --clump-p1 ${p} \
        --clump-kb ${KB} \
        --clump-r2 ${Rsq} \
        --clump-range ${DIR}'/annotate/'${ANNOT} \
        --out ${DIR}'/clump/'${DATE}'/'${OUTPUT}'_p-'${p}'_kb-'${KB}'_r2-'${Rsq}'.loco'
else
echo "No results file (loco)?"
fi
