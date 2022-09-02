#!/bin/bash
source /broad/software/scripts/useuse
umask 002
use GCC-5.2
use .htslib-1.8
use R-3.5
PATH=$PATH:/seq/vgb/software/anaconda3/current
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/env/
source activate ${DIR}/env/Renv
gcta='/seq/vgb/software/gcta/current'

Rscript ${DIR}'/bin/ld_score_quartiles.R' ${DIR}'/lds' ${GENO}

$gcta --bfile ${DIR}'/geno/'${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}'/lds/'${GENO}'.score.ld.Q000-Q025.txt' \
      --out ${DIR}'/grm/'${GENO}'.score.ld.Q000-Q025'
echo ${DIR}'/grm/'${GENO}'.score.ld.Q000-Q025' > ${DIR}'/grm/'${GENO}'.score.ld.list.txt'

$gcta --bfile ${DIR}'/geno/'${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}'/lds/'${GENO}'.score.ld.Q025-Q050.txt' \
      --out ${DIR}'/grm/'${GENO}'.score.ld.Q025-Q050'
echo ${DIR}'/grm/'${GENO}'.score.ld.Q025-Q050' >> ${DIR}'/grm/'${GENO}'.score.ld.list.txt'

$gcta --bfile ${DIR}'/geno/'${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}'/lds/'${GENO}'.score.ld.Q050-Q075.txt' \
      --out ${DIR}'/grm/'${GENO}'.score.ld.Q050-Q075'
echo ${DIR}'/grm/'${GENO}'.score.ld.Q050-Q075' >> ${DIR}'/grm/'${GENO}'.score.ld.list.txt'

$gcta --bfile ${DIR}'/geno/'${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}'/lds/'${GENO}'.score.ld.Q075-Q100.txt' \
      --out ${DIR}'/grm/'${GENO}'.score.ld.Q075-Q100'
echo ${DIR}'/grm/'${GENO}'.score.ld.Q075-Q100' >> ${DIR}'/grm/'${GENO}'.score.ld.list.txt'
