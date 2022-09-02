# USAGE: ./mlma-to-wig.sh chr2 /seq/vgb/dd/gwas/assoc/Q125A_nocrop.mlma output-prefix

# CHR=${1}
# TARGET=${2}
# OUT=${3}
# TITLE=${4}

# echo "track type=wiggle_0 name=\"${TITLE}\" visibility=full autoScale=off viewLimits=0.0:5.2 color=51,160,44" > ${OUT}.wig
# echo "variableStep chrom=${CHR}" >> ${OUT}.wig
# awk -v CHR="${CHR}" '$2 ~ CHR":" { print $3,(-log($9)/log(10)); }' ${TARGET} >> ${OUT}.wig

TARGET=${TARGET:-${1}}
TITLE=${TITLE:-${2}}

echo "track type=wiggle_0 name=\"${TITLE}\" visibility=full autoScale=off viewLimits=0.0:5.2 color=51,160,44" > ${TARGET}.wig

