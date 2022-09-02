#!/bin/bash
#$ -q broad
#$ -l h_vmem=32g
#$ -l h_rt=24:00:00
#$ -o /seq/vgb/dd/gwas/logs/out/
#$ -e /seq/vgb/dd/gwas/logs/err/
source /broad/software/scripts/useuse
umask 002

use GCC-5.2
use .htslib-1.8
use R-3.5
use BEDTools

conda activate RSAIGE
FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"

gemma=/seq/vgb/software/gemma/current
gcta=/seq/vgb/software/gcta/current
plink=/seq/vgb/software/plink2/current/plink

# USAGE: qsub GCTA_GWAS.job geno-prefix grm-prefix pheno-prefix covar-prefix qcovar-prefix assoc? clump? plot?

DIR='/seq/vgb/dd/gwas/'
GENO=${1}
GRM=${2}
PHENO=${3}
COVAR=${4}
if [ $5 == "NA" ] ; then
  echo "Not using a quantitative covariate"
  useQCOVAR=F
  QCOVAR="NA"
else
  QCOVAR=${5}
  useQCOVAR=T
fi
ASSOC=${6:=T}
CLUMP=${7:=T}
PLOT=${8:=F}
# TODAY=`date +%Y-%m-%d`
# DATE=${9:=${TODAY}}
if [ -z "$9" ] ; then
  DATE=`date +%Y-%m-%d`
  echo "Setting date to current date ${DATE}"
else
  DATE=${9}
fi

cd ${DIR}

if [ ${ASSOC} == T ]; then

#step 1: fit the NULL glmm. For quantitative traits, if not normally distributed, inverse normalization needs to be specified to be TRUE --invNormalize=TRUE

Rscript step1_fitNULLGLMM.R     \
        --plinkFile=${GENO} \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
	--invNormalize=TRUE	\
        --outputPrefix=./output/example_quantitative \
        --nThreads=4 \
        --LOCO=FALSE	\
	--tauInit=1,0
	


Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=4 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.vcf.genotype.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE    


########sav file
Rscript step2_SPAtests.R	\
	--savFile=./input/dosage_10markers.sav	\
	--savFileIndex=./input/dosage_10markers.sav.s1r	\
	--minMAF=0.0001 \
        --minMAC=4 \
	--vcfField=DS \
	--chrom=1 \
        --sampleFile=./input/sampleIDindosage.txt \
        --GMMATmodelFile=./output/example.rda \
        --varianceRatioFile=./output/example.varianceRatio.txt \
        --SAIGEOutputFile=./output/example.SAIGE.sav.txt \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE


Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
        --minMAF=0.0001 \
        --minMAC=1 \
        --sampleFile=./input/sampleIDindosage.txt \
	--GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.vcf.genotype.txt_cond \
        --numLinesOutput=2 \
        --IsOutputAFinCaseCtrl=TRUE     \
        --condition=1:4_1/2  #conditional analysis can be performed if a conditioning genetic marker is specified (chr:pos_ref/alt)


#For gene-based test
#step 1: fit the NULL glmm
#step 1 model result from the single-variant assoc test can be re-used, except that for gene-based tests, variance ratios for multiple MAC categories and a sparse GRM need to be used. If IsSparseKin=TRUE and no sparseSigmaFile and sparseSigmaSampleIDFile are specified, a sparse GRM will be created based on the relatednessCutoff. sparseSigmaFile and sparseSigmaSampleIDFile can be used to specify a pre-calcuated sparse GRM and the sample ids for the sparse GRM. Tested samples would be a subset of samples in the pre-calcuated GRM.

#To activate the variance ratio estimation based multiple MAC categories, --isCateVarianceRatio=TRUE
#cateVarRatioMinMACVecExclude and cateVarRatioMaxMACVecInclude are used to specify the MAC categories
#by default --cateVarRatioMinMACVecExclude=0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5
#--cateVarRatioMaxMACVecInclude=1.5,2.5,3.5,4.5,5.5,10.5,20.5
#corresponding to
#0.5 < MAC <=  1.5
#1.5 < MAC <=  2.5
#2.5 < MAC <=  3.5
#3.5 < MAC <=  4.5
#4.5 < MAC <=  5.5
#5.5 < MAC <=  10.5
#10.5 < MAC <=  20.5
#20.5 < MAC

###The sparse GRM that has been previously calculated in the job for quantitative traits can be used. 
###Note: a sparse GRM can be calculated once for each data set and re-used for all phenotypes for that data set as long as all samples tested are included in the sparse GRM.   
###--skipModelFitting=TRUE becuase the null GLMM has already been fit when performing single-variant assoc test above
###--isCateVarianceRatio=TRUE categorical variance ratios for different MAC categories need to be calculated 
###NOTE:Please store the single variance ratio for single-variant assoc test before this step. e.g rename the file, since the variance ratio file will contain categorical variance ratios if --isCateVarianceRatio=TRUE

##quantitative traits
Rscript step1_fitNULLGLMM.R     \
	--plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_quantitative \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=quantitative       \
        --invNormalize=TRUE     \
        --outputPrefix=./output/example_quantitative \
	--outputPrefix_varRatio=./output/example_quantitative_cate	\
	--sparseGRMFile=./output/example_binary_cate.varianceRatio.txt.sparseGRM.mtx    \
        --sparseGRMSampleIDFile=./output/example_binary.varianceRatio.txt.sparseGRM.mtx.sample  \
        --nThreads=4 \
        --LOCO=FALSE	\
	--skipModelFitting=FALSE \
        --IsSparseKin=TRUE      \
        --isCateVarianceRatio=TRUE	


##binay traits
Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary       \
        --invNormalize=TRUE     \
        --outputPrefix=./output/example_binary \
        --outputPrefix_varRatio=./output/example_binary_cate      \
        --sparseGRMFile=./output/example_binary_cate.varianceRatio.txt.sparseGRM.mtx    \
        --sparseGRMSampleIDFile=./output/example_binary.varianceRatio.txt.sparseGRM.mtx.sample  \
        --nThreads=4 \
        --LOCO=FALSE    \
        --skipModelFitting=FALSE \
        --IsSparseKin=TRUE      \
        --isCateVarianceRatio=TRUE



#Perform gene-based/region-based tests according to the group file specified in groupFile
#IsSingleVarinGroupTest=TRUE is to perform single-variant assoc tests as well for markers included in the gene-based tests
#only vcf, sav, and bgen dosage file formats can be used for gene-based tests
#to perform gene-based tests, --groupFile is used to specify a group file
#Each line is for one gene/set of
#          variants. The first element is for gene/set name. The rest of
#          the line is for variant ids included in this gene/set. For
#          vcf/sav, the genetic marker ids are in the format
#          chr:pos_ref/alt. For begen, the genetic marker ids should
#          match the ids in the bgen file. Each element in the line is
#          seperated by tab.

Rscript step2_SPAtests.R \
        --vcfFile=./input/genotype_10markers.vcf.gz \
        --vcfFileIndex=./input/genotype_10markers.vcf.gz.tbi \
        --vcfField=GT \
        --chrom=1 \
	--minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_quantitative.SAIGE.gene.txt \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest_simulation.txt    \
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsSingleVarinGroupTest=TRUE
        --condition=1:4_1/2


##another example, conditional analysis for gene-based tests
Rscript step2_SPAtests.R \
        --vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
        --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
        --vcfField=DS \
        --chrom=chr1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_quantitative.rda \
        --varianceRatioFile=./output/example_quantitative_cate.varianceRatio.txt \
	--SAIGEOutputFile=./output/example_quantitative.SAIGE.gene_conditional.txt \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest.txt    \
        --sparseSigmaFile=./output/example_quantitative_cate.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsSingleVarinGroupTest=TRUE   \
	--condition=chr1:32302_A/C 

##binary traits
Rscript step2_SPAtests.R \
        --vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
        --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
        --vcfField=DS \
        --chrom=chr1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary_cate_v2.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.gene.txt \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest.txt    \
        --sparseSigmaFile=./output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsSingleVarinGroupTest=TRUE   \
	--IsOutputPvalueNAinGroupTestforBinary=TRUE	\
	--IsAccountforCasecontrolImbalanceinGroupTest=TRUE

Rscript step2_SPAtests.R \
        --vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
        --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
        --vcfField=DS \
        --chrom=chr1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary_cate_v2.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.gene_conditional.txt_0.37 \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest.txt    \
        --sparseSigmaFile=./output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsSingleVarinGroupTest=TRUE   \
        --IsOutputPvalueNAinGroupTestforBinary=TRUE     \
        --IsAccountforCasecontrolImbalanceinGroupTest=TRUE	\
	--condition=chr1:32302_A/C



####Specify customized weights for markers in the gene- or region-based tests
#weightsIncludeinGroupFile logical. Whether to specify customized weight for makers in gene- or region-based tests. If TRUE, weights are included in the group file. For vcf/sav, the genetic marker ids and weights are in the format chr:pos_ref/alt;weight. For bgen, the genetic marker ids should match the ids in the bgen filE, e.g. SNPID;weight. Each element in the line is seperated by tab. By default, FALSE
#weights_for_G2_cond vector of float. weights for conditioning markers for gene- or region-based tests. The length equals to the number of conditioning markers, delimited by comma.
Rscript step2_SPAtests.R \
        --vcfFile=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav \
        --vcfFileIndex=./input/seedNumLow_126001_seedNumHigh_127000_nfam_1000_nindep_0.sav.s1r \
        --vcfField=DS \
        --chrom=chr1 \
        --minMAF=0 \
        --minMAC=0.5 \
        --maxMAFforGroupTest=0.01       \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary_cate_v2.varianceRatio.txt \
        --SAIGEOutputFile=./output/example_binary.SAIGE.gene_conditional_withspecifiedWeights.txt.0.37 \
        --numLinesOutput=1 \
        --groupFile=./input/groupFile_geneBasedtest_withWeights.txt    \
        --sparseSigmaFile=./output/example_binary_cate_v2.varianceRatio.txt_relatednessCutoff_0.125.sparseSigma.mtx       \
        --IsOutputAFinCaseCtrl=TRUE     \
        --IsSingleVarinGroupTest=TRUE   \
        --IsOutputPvalueNAinGroupTestforBinary=TRUE     \
        --IsAccountforCasecontrolImbalanceinGroupTest=TRUE      \
	--weightsIncludeinGroupFile=TRUE	\
	--weights_for_G2_cond=3,1	\
	--condition=chr1:32302_A/C,chr1:32304_A/C	

#--sampleFile=./input/samplelist.txt \
