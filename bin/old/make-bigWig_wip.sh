#!/bin/bash
#$ -q broad
#$ -l h_vmem=5g
#$ -l h_rt=24:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
#$ -M kmorrill@broadinstitute.org
#$ -m e
#$ -pe smp 2
#$ -binding linear:2
#$ -R y
source /broad/software/scripts/useuse
umask 002

# USAGE: qsub -v TITLE="Size",OUT="/seq/vgb/dd/gwas/assoc/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_Q121A_qcov-Q121A" /seq/vgb/dd/gwas/bin/make_bedGraph.sh /seq/vgb/dd/gwas/assoc/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_Q121A_qcov-Q121A.loco.mlma

# add track settings for log10p and effect
echo "track type=bedGraph name=\"${TITLE}\" description=\"Association p-values (-log10) from GWAS for ${TITLE} in diverse cohort of pet dogs from Darwin's Ark\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 autoScale=off" > ${OUTPUT}_log10p.bedGraph 
echo "track type=bedGraph name=\"${TITLE}\" description=\"Effect sizes (b) from GWAS for ${TITLE} in diverse cohort of pet dogs from Darwin's Ark\" visibility=full color=200,100,0 altColor=0,100,200 priority=20 autoScale=off" > ${OUTPUT}_log10p.bedGraph 

# get statistics of interest and convert from 1-base to 0-base positions, make bedGraph files
awk -F "\t" 'OFS=FS { print "chr"$1,$3-1,$3,(-log($9)/log(10)); }' ${TARGET} | tail -n+2 | sort -k1,1 -k2,2n >> ${TARGET}.bedGraph
/seq/vgb/dd/gwas/bin/bedGraphToBigWig ${TARGET}.bedGraph /seq/vgb/dd/gwas/tracks/canFam3.chrom.sizes ${TARGET}.bigWig
