in=${in:-"/seq/vgb/dd/gwas/bed/permBeds/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_coat-pattern_ticking_piebald-only_qcov-NA_p-1e-6_kb-250_r2-0.2.loco.clumped.ranges.bed.tosort"}
out=${out:-"/seq/vgb/dd/gwas/bed/permBeds/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_coat-pattern_ticking_piebald-only_qcov-NA_p-1e-6_kb-250_r2-0.2.loco.clumped.ranges.bed"}

/seq/vgb/software/bedtools/bedtools sort -g /seq/vgb/resources/dog/bed/canFam3.autosome.genome -i ${in} > ${out}
