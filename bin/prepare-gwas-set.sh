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
     -v DIR=${DIR},GENO=${GENO} \
     ${DIR}/bin/GCTA_GRM.sh

# submit jobs for GRMs per chromosome and jobs for LD scoring
for CHR in `seq 1 38`
do
  KB=${KB:-'250'}
  qsub -q broad \
       -l h_vmem=12g \
       -l h_rt=4:00:00 \
       -N ${GENO}'_LD-score_chr-'${CHR} \
       -v DIR=${DIR},CHR=${CHR},GENO=${GENO},KB=${KB} \
       ${DIR}/bin/GCTA_LD-score.sh

  qsub -q broad \
      -l h_vmem=12g \
      -l h_rt=4:00:00 \
      -N ${GENO}'_GRM_chr-'${CHR} \
      -v DIR=${DIR},CHR=${CHR},GENO=${GENO} \
      ${DIR}/bin/GCTA_GRM_by-chr.sh
done

# submit job to split variants by LD scores and generate GRMs
qsub -q broad \
     -hold_jid ${GENO}_LD-score_chr-* \
     -l h_vmem=16g \
     -l h_rt=4:00:00 \
     -N ${GENO}'_GRM_LD-strat' \
     -v DIR=${DIR},GENO=${GENO} \
     ${DIR}/bin/GCTA_GRM_LD-strat.sh
