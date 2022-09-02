#!/bin/bash
gcta=/seq/vgb/software/gcta/current
$gcta --bfile ${DIR}/geno/${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --out ${DIR}/grm/${GENO}
