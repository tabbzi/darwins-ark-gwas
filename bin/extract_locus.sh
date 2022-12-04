#!/bin/bash
#$ -q broad
#$ -l h_vmem=32g
#$ -l h_rt=1:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002
plink=/seq/vgb/software/plink2/current/plink
plink2=/seq/vgb/software/plink2/dev

# snp
# Q
# D
# P
# chr
# start
# end

P=${P:-'bq.017'}
D=${D:-'datatype'}
Q=${Q:-'age.hgt'}
chr=${chr:-'30'}
start=${start:-'33522500'}
end=${end:-'33529500'}
snp=${snp:-'30:33526686:G:A'}

DIR=${DIR:-"/seq/vgb/dd/gwas"}
DATE=${DATE:-"2022-07-15"}

mkdir ${DIR}'/locus/'${DATE}

GENO=${GENO:-"DarwinsArk_gp-0.70_biallelic-snps_maf-0.02_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3104"}
MLMA=${MLMA:-${DIR}'/assoc/'${DATE}'/'${GENO}'_phe-'${P}'_dcov-'${D}'_qcov-'${Q}'.loco.mlma'}
OUT=${OUT:-${DIR}'/locus/'${DATE}'/'${GENO}'_phe-'${P}'_dcov-'${D}'_qcov-'${Q}'_chr-'${chr}'_start-'${start}'_end-'${end}'_snp-'${snp}}

echo ${P}
echo ${Q}
echo ${D}
echo ${chr}
echo ${start}
echo ${end}
echo ${snp}
echo ${MLMA}
echo ${OUT}

printf "${chr}\t${start}\t${end}" > ${OUT}_locus.bed
/seq/vgb/software/plink2/dev --dog --bfile ${DIR}'/geno/'${GENO} --extract bed1 ${OUT}'_locus.bed' --make-bed --out ${OUT}

# linkage:
/seq/vgb/software/plink2/current/plink --dog --bfile ${OUT} --r2 --ld-snp ${snp} --ld-window-r2 0 --ld-window 10000 --ld-window-kb 10000 --out ${OUT}
# association:
awk -v c=${chr} -v a=${start} -v b=${end} -F "\t" '$1==c&&$3>=a&&$3<=b {print $2"\t"$3"\t"$7"\t"$9}' ${MLMA} > ${OUT}.ass.tsv
# effects:
awk -v c=${chr} -v a=${start} -v b=${end} -F "\t" '$1==c&&$2>=a&&$2<=b {print $0}' ${DIR}'/geno/'${GENO}.annotated.tsv > ${OUT}.eff.tsv
# phyloP:
/seq/vgb/software/bedtools/bedtools intersect -a <(awk '{print "chr"$0}' ${OUT}_locus.bed) -b /seq/vgb/dd/gwas/zoo/dog-centered-Feb2021/chr${chr}.bed.gz -wb > ${OUT}.zoo.tsv


# clean up
rm ${OUT}.bim
rm ${OUT}.bed
rm ${OUT}.fam
rm ${OUT}.log
rm ${OUT}_locus.bed
