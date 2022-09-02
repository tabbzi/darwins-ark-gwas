#!/bin/bash
gcta=/seq/vgb/software/gcta/current
$gcta --bfile ${DIR}/geno/${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --chr ${CHR} \
      --out ${DIR}/grm/${GENO}_chr-${CHR}
