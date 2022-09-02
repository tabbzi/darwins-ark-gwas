# USAGE: ./liftOver_clumps.sh FILE.clumped.ranges

paste <(grep -vf <(head -n 1 ${1}.clumped.ranges) ${1}.clumped.ranges | awk '{print $5}' | sed s/"\.\."/"\t"/g | sed s/":"/"\t"/g) <(grep -vf <(head -n 1 ${1}.clumped.ranges) ${1}.clumped.ranges | awk '{print $2,$3}') | sed s/" "/":pval:"/g > ${1}_temp.bed

./annotate/liftOver/liftOver -minMatch=0.1 ${1}_temp.bed ./annotate/liftOver/canFam3ToHg38.over.chain ${1}_lifted_hg38.bed ${1}_unlifted.bed

printf "index\tdogP\tlifted_chr\tlifted_pos\ttrait\thumanP\teffect\tgenes\tmapped_gene\n" > ${1}_annotations.txt
bedtools intersect -a ${1}_lifted_hg38.bed -b ./annotate/GWAS_Catalog/GWAS_Catalog_fixed.bed -wb | awk -F "\t" 'BEGIN {OFS=FS} {print $8,$4}' | sort | join -t $'\t' - <(sort ./annotate/GWAS_Catalog/GWAS_Catalog_fixed.tsv) | awk -F "\t" 'BEGIN {OFS=FS} {print $2,$3,$6,$7,$8,$9,$10,$11}' | sed s/":pval:"/"\t"/g >> ${1}_annotations.txt

rm ${1}_temp.bed
rm ${1}_unlifted.bed
