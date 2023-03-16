#!/bin/bash
A=${1}
B=${2}
PNA=${3}
PNB=${4}

printenv

echo ${DIR}
echo ${A}
echo ${B}

if [ ${useRELCUT} == 'T' ] ; then
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}'_rel-cutoff-'${RELCUT}'_'${A}'-vs-'${B}}
else
  OUTPUT=${OUTPUT:-${GENO}'_phe-'${PHE}'_dcov-'${DCOV}'_qcov-'${QCOV}'_'${A}'-vs-'${B}}
fi

echo ${OUTPUT}

echo "Performing bivariate GCTA-GREML with constraint between ${A} and ${B}..."
${gcta} --grm ${INPUT_GRM_ALL} \
        --reml-bivar ${PNA} ${PNB} \
        --pheno ${PHEFILE} \
        ${ARG_DCOVAR} \
        ${ARG_QCOVAR} \
        --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.bivar-REML.no-lds'

echo "Performing bivariate GCTA-GREML without constraint between ${A} and ${B}..."
${gcta} --grm ${INPUT_GRM_ALL} \
        --reml-bivar ${PNA} ${PNB} \
        --reml-bivar-no-constrain \
        --pheno ${PHEFILE} \
        ${ARG_DCOVAR} \
        ${ARG_QCOVAR} \
        --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.bivar-REML.no-lds.no-constraint'

echo "Performing bivariate GCTA-GREML-LDMS with constraint between ${A} and ${B}..."
${gcta} --mgrm ${INPUT_MGRM_LIST} \
        --reml-bivar ${PNA} ${PNB} \
        --pheno ${PHEFILE} \
        ${ARG_DCOVAR} \
        ${ARG_QCOVAR} \
        --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.bivar-REML.lds'

echo "Performing bivariate GCTA-GREML-LDMS without constraint between ${A} and ${B}..."
${gcta} --mgrm ${INPUT_MGRM_LIST} \
        --reml-bivar ${PNA} ${PNB} \
        --reml-bivar-no-constrain \
        --pheno ${PHEFILE} \
        ${ARG_DCOVAR} \
        ${ARG_QCOVAR} \
        --out ${DIR}'/cor/'${DATE}'/'${OUTPUT}'.bivar-REML.lds.no-constraint'
