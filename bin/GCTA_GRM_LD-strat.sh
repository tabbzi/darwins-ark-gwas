#!/bin/bash
source activate ${DIR}/env/Renv
gcta=/seq/vgb/software/gcta/current

Rscript ${DIR}/bin/ld_score_quartiles.R ${DIR}/lds/${GENO}

$gcta --bfile ${DIR}/geno/${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}/lds/${GENO}'.score.ld.Q000-Q025.txt' \
      --out ${DIR}/grm/${GENO}'.score.ld.Q000-Q025'

$gcta --bfile ${DIR}/geno/${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}/lds/${GENO}'.score.ld.Q025-Q050.txt' \
      --out ${DIR}/grm/${GENO}'.score.ld.Q025-Q050'

$gcta --bfile ${DIR}/geno/${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}/lds/${GENO}'.score.ld.Q050-Q075.txt' \
      --out ${DIR}/grm/${GENO}'.score.ld.Q050-Q075'

$gcta --bfile ${DIR}/geno/${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --extract ${DIR}/lds/${GENO}'.score.ld.Q075-Q100.txt' \
      --out ${DIR}/grm/${GENO}'.score.ld.Q075-Q100'
