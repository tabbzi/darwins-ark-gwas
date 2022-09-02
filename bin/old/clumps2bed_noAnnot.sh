tail -n+2 ${1} | awk '{print $1,$5}' | sed s/"\.\."/" "/g | tr ":" " " | awk '{print $1,$3,$4}' | tr " " "\t" > ${1}.bed.tosort

export in=${1}.bed.tosort
export out=${1}.bed
/seq/vgb/dd/gwas/bin/sort_bed.sh

/seq/vgb/software/bedtools/bedtools merge -i ${1}.bed > ${1}.comp.noAnnot.bed
