tail -n+2 ${1} | awk '{print $1,$3,$5}' | sed s/"\.\."/" "/g | tr ":" " " | awk '{print $1,$4,$5,$2}' | tr " " "\t" > ${1}.bed.tosort

export in=${1}.bed.tosort
export out=${1}.bed
/seq/vgb/dd/gwas/bin/sort_bed.sh

/seq/vgb/software/bedtools/bedtools merge -i ${1}.bed | /seq/vgb/software/bedtools/bedtools map -a - -b ${1}.bed -c 4 -o min | awk -F "\t" 'OFS=FS {print $0,NR}' > ${1}.comp.bed
