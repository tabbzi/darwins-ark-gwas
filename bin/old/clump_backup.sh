DIR=${DIR:-"/seq/vgb/dd/gwas"}
GENO=${GENO:-"DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155"}
DATE=${DATE:-"MuttPaper"}
KB=${KB:-1000}
Rsq=${Rsq:-0.2}
p=${p:-"1e-04"}
PHENO=${PHENO:-"Q121A"}

plink=/seq/vgb/software/plink2/current/plink

$plink --bfile ${DIR}/geno/${GENO} \
      --clump ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.mlma \
      --dog \
      --allow-no-sex \
      --pheno ${DIR}/pheno/${PHENO}.pheno \
      --clump-field p \
      --clump-p1 ${p} \
      --clump-kb ${KB} \
      --clump-r2 ${Rsq} \
      --clump-range ${DIR}/annotate/Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.bed \
      --out ${DIR}/clump/${DATE}/${GENO}_${PHENO}_p-${p}_kb-${KB}_r2-${Rsq}

if [ -f ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.loco.mlma ]; then
$plink --bfile ${DIR}/geno/${GENO} \
      --clump ${DIR}/assoc/${DATE}/${GENO}_${PHENO}.loco.mlma \
      --dog \
      --allow-no-sex \
      --pheno ${DIR}/pheno/${PHENO}.pheno \
      --clump-field p \
      --clump-p1 ${p} \
      --clump-kb ${KB} \
      --clump-r2 ${Rsq} \
      --clump-range ${DIR}/annotate/Canis_lupus_familiaris.CanFam3.1.ensembl.gene_annotations.bed \
      --out ${DIR}/clump/${DATE}/${GENO}_${PHENO}_p-${p}_kb-${KB}_r2-${Rsq}.loco
fi
