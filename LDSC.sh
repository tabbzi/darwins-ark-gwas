./ldsc.py \
--out test \
--rg phe-bq.017_cov-age.hgt.sumstats.gz,phe-bq.018_cov-age.hgt.sumstats.gz,phe-lineage.int.01_cov-NA.sumstats.gz \
--w-ld ../lds/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465 \
--ref-ld ../lds/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465 \
--return-silly-things \
--intercept-gencov 1,1,1
