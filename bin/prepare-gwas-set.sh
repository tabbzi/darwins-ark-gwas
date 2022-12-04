#!/bin/bash
DIR='/seq/vgb/dd/gwas'
GENO=${GENO:-'DarwinsArk_gp-0.70_biallelic-snps_maf-0.02_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3104'}

mkdir ${DIR}/grm/
mkdir ${DIR}/lds/

# generate GRM across all variants
qsub -q broad \
     -l h_vmem=12g \
     -l h_rt=4:00:00 \
     -N ${GENO}'_GRM' \
     -v DIR=${DIR},GENO=${GENO},RELCUT=${RELCUT} \
     ${DIR}/bin/GCTA_GRM.sh

# submit jobs for GRMs per chromosome and jobs for LD scoring
for CHR in `seq 1 38`
do
  KB=${KB:-'250'}
  qsub -q broad \
       -hold_jid ${GENO}'_GRM' \
       -l h_vmem=12g \
       -l h_rt=4:00:00 \
       -N ${GENO}'_LD-score_chr-'${CHR} \
       -v DIR=${DIR},CHR=${CHR},GENO=${GENO},KB=${KB} \
       ${DIR}/bin/GCTA_LD-score.sh

  qsub -q broad \
      -hold_jid ${GENO}'_GRM' \
      -l h_vmem=12g \
      -l h_rt=4:00:00 \
      -N ${GENO}'_GRM_chr-'${CHR} \
      -v DIR=${DIR},CHR=${CHR},GENO=${GENO},RELCUT=${RELCUT} \
      ${DIR}/bin/GCTA_GRM_by-chr.sh
done

# submit job to split variants by LD scores and generate GRMs
qsub -q broad \
     -hold_jid ${GENO}_LD-score_chr-* \
     -l h_vmem=24g \
     -l h_rt=4:00:00 \
     -N ${GENO}'_GRM_LD-strat' \
     -v DIR=${DIR},GENO=${GENO},RELCUT=${RELCUT} \
     ${DIR}/bin/GCTA_GRM_LD-strat.sh

# submit job to annotate variant effects using SnpEff
qsub -v GENO=${GENO} ${DIR}/bin/snpEff.sh

# submit job to estimate taggings for genetic correlations using LDAK

# make MAGMA annotation files
qsub -v GENO=${GENO},GENE='/seq/vgb/dd/gwas/zoo/zooUNICORNs_refined_dog_genome.ann.bed',ANNO='magma.zoo-UNICORNs' ${DIR}/bin/MAGMA_make-anno.sh
qsub -v GENO=${GENO},GENE='/seq/vgb/dd/gwas/zoo/zooRoCCs_GrEq20_dog_genome.ann.bed',ANNO='magma.zoo-RoCCs' ${DIR}/bin/MAGMA_make-anno.sh
qsub -v GENO=${GENO},GENE='/seq/vgb/dd/gwas/zoo/zooUCEs_gt19bp_dog_genome.ann.bed',ANNO='magma.zoo-UCEs' ${DIR}/bin/MAGMA_make-anno.sh
