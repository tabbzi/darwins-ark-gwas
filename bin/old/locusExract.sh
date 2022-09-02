#!/bin/bash
#$ -q broad
#$ -l h_vmem=8g
#$ -l h_rt=1:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

DIR="/seq/vgb/dd/gwas/locus"
GENO="/seq/vgb/dd/gwas/geno/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt"
MLMA="/seq/vgb/dd/gwas/assoc/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_qcov-${Q}.loco.mlma"
OUT="/seq/vgb/dd/gwas/locus/${P}_qcov-${Q}_chr${chr}_start-${start}_end-${end}_snp-${snp}"
printf "${chr}\t${start}\t${end}" > ${OUT}_locus.bed
/seq/vgb/software/plink2/dev --dog --bfile ${GENO} --extract bed1 ${OUT}_locus.bed --make-bed --out ${OUT}

# linkage:
/seq/vgb/software/plink2/current/plink --dog --bfile ${OUT} --r2 --ld-snp ${snp} --ld-window-r2 0 --ld-window 10000 --ld-window-kb 10000 --out ${OUT}
# association:
awk -v c=${chr} -v a=${start} -v b=${end} -F "\t" '$1==c&&$3>=a&&$3<=b {print $2"\t"$3"\t"$7"\t"$9}' ${MLMA} > ${OUT}.ass.tsv
# effects:
awk -v c=${chr} -v a=${start} -v b=${end} -F "\t" '$1==c&&$2>=a&&$2<=b {print $0}' ${GENO}.annotation.tsv > ${OUT}.eff.tsv
# phyloP:
/seq/vgb/software/bedtools/bedtools intersect -a <(awk '{print "chr"$0}' ${OUT}_locus.bed) -b /seq/vgb/dd/gwas/zoo/dog-centered-Feb2021/chr${chr}.bed.gz -wb > ${OUT}.zoo.tsv


# clean up
rm ${OUT}.bim
rm ${OUT}.bed
rm ${OUT}.fam
rm ${OUT}.log
rm ${OUT}_locus.bed
