#!/bin/bash
gcta=/seq/vgb/software/gcta/current
$gcta --bfile ${DIR}/geno/${GENO} \
      --ld-score-region ${KB} \
      --chr ${CHR}  \
      --autosome \
      --autosome-num 38 \
      --out ${DIR}/lds/${GENO}'_kb-'${KB}'_chr-'${CHR}
