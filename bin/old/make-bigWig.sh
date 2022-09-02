# USAGE: ./make_bedGraph.sh /seq/vgb/dd/gwas/assoc/current/... "my_title"

TARGET=${TARGET:-${1}}
TITLE=${TITLE:-${2}}

echo "track type=bedGraph name=\"${TITLE}\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 autoScale=off" > ${TARGET}.bedGraph 
awk -F "\t" 'OFS=FS { print "chr"$1,$3-1,$3,(-log($9)/log(10)); }' ${TARGET} | tail -n+2 | sort -k1,1 -k2,2n >> ${TARGET}.bedGraph
/seq/vgb/dd/gwas/bin/bedGraphToBigWig ${TARGET}.bedGraph /seq/vgb/dd/gwas/tracks/canFam3.chrom.sizes ${TARGET}.bw
