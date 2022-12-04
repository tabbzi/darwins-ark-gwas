#!/bin/bash
gcta=/seq/vgb/software/gcta/current
$gcta --bfile ${DIR}'/geno/'${GENO} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --chr ${CHR} \
      --out ${DIR}'/grm/'${GENO}'_chr-'${CHR}

$gcta --grm ${DIR}'/grm/'${GENO}'_chr-'${CHR} \
      --make-grm \
      --autosome \
      --autosome-num 38 \
      --keep ${DIR}'/grm/'${GENO}'_rel-cutoff-'${RELCUT}'.grm.id' \
      --out ${DIR}'/grm/'${GENO}'_chr-'${CHR}'_rel-cutoff-'${RELCUT}
