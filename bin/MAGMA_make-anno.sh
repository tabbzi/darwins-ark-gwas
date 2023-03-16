#!/bin/bash
#$ -q broad
#$ -l h_vmem=16g
#$ -l h_rt=2:00:00
#$ -o /seq/vgb/dd/gwas/logs/
#$ -e /seq/vgb/dd/gwas/logs/
#$ -M kmorrill@broadinstitute.org
#$ -m e
source /broad/software/scripts/useuse
umask 002

use GCC-5.2
use .htslib-1.8
use R-3.5
use BEDTools

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/env/

# software:
magma=/seq/vgb/software/magma/current

# prepare:
DIR='/seq/vgb/dd/gwas'

mkdir ${DIR}/magma/anno/${DATE}

# snp dataset:
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.005_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3104'}

# gene dataset:
GENE=${GENE:-'/seq/vgb/dd/gwas/zoo/zooUNICORNs_refined_dog_genome.ann.bed'}

# outputs:
ANNODIR=${ANNODIR:-${DIR}'/magma/anno/'${DATE}}
ANNO=${ANNO:-'magma.zoo-UNICORNs'}
ANNOFILE=${ANNOFILE:-${ANNODIR}'/'${ANNO}'.genes.annot'}

${magma} --annotate nonhuman \
         --snp-loc ${DIR}'/geno/'${GENO}'.bim' \
         --gene-loc ${GENE} \
         --out ${ANNODIR}'/'${GENO}'_'${ANNO}
