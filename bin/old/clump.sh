DIR=${DIR:-"/seq/vgb/dd/gwas"}
GENO=${GENO:-"DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt"}
DATE=${DATE:-"mutt_paper"}
KB=${KB:-250}
Rsq=${Rsq:-0.2}
p=${p:-"1e-6"}
PHENO=${PHENO:-"Q121A"}
OUT=${OUT:-${GENO}_${PHENO}}
# ANNOT=${ANNOT:-"Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.expanded.bed"}
ANNOT=${ANNOT:-"Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.withHuman.bed"}

plink=/seq/vgb/software/plink2/current/plink

echo "Starting clumping..."

if [ -f ${DIR}/assoc/${DATE}/${OUT}.mlma ]; then
$plink --bfile ${DIR}/geno/${GENO} \
      --clump ${DIR}/assoc/${DATE}/${OUT}.mlma \
      --dog \
      --allow-no-sex \
      --pheno ${DIR}/pheno/${PHENO}.pheno \
      --clump-field p \
      --clump-p1 ${p} \
      --clump-kb ${KB} \
      --clump-r2 ${Rsq} \
      --clump-range ${DIR}/annotate/${ANNOT} \
      --out ${DIR}/clump/${DATE}/${OUT}_p-${p}_kb-${KB}_r2-${Rsq}
else
echo "No results file (non-loco)?"
fi

if [ -f ${DIR}/assoc/${DATE}/${OUT}.loco.mlma ]; then
$plink --bfile ${DIR}/geno/${GENO} \
      --clump ${DIR}/assoc/${DATE}/${OUT}.loco.mlma \
      --dog \
      --allow-no-sex \
      --pheno ${DIR}/pheno/${PHENO}.pheno \
      --clump-field p \
      --clump-p1 ${p} \
      --clump-kb ${KB} \
      --clump-r2 ${Rsq} \
      --clump-range ${DIR}/annotate/${ANNOT} \
      --out ${DIR}/clump/${DATE}/${OUT}_p-${p}_kb-${KB}_r2-${Rsq}.loco
else
echo "No results file (loco)?"
fi
