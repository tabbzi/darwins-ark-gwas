#!/bin/bash
#$ -q broad
#$ -l h_vmem=6g
#$ -l h_rt=1:00:00
#$ -e /seq/vgb/dd/gwas/logs/
#$ -o /seq/vgb/dd/gwas/logs/
umask 002
source /broad/software/scripts/useuse
PATH=$PATH:/seq/vgb/software/anaconda3/current
source activate /seq/vgb/dd/gwas/pbs/env

P=${P:-"Q121A"}
Q=${Q:-"Q121A"}



python3 /seq/vgb/dd/gwas/bin/join_clmp_to_mlma.py -Q ${Q} -D "all dogs" -C /seq/vgb/dd/gwas/clump/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_qcov-${Q}_p-1e-6_kb-250_r2-0.2.loco.clumped.ranges -P ${P} -S /seq/vgb/dd/gwas/assoc/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_qcov-${Q}.loco.mlma -O /seq/vgb/dd/gwas/summ/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_qcov-${Q}_p-1e-6_kb-250_r2-0.2.loco.tsv

# python3 /seq/vgb/dd/gwas/bin/join_clmp_to_mlma.py -Q ${Q} -D "highly admixed" -C /seq/vgb/dd/gwas/clump/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_non-purebred_highly-admixed_qcov-${Q}_p-1e-6_kb-250_r2-0.2.loco.clumped.ranges -P ${P} -S /seq/vgb/dd/gwas/assoc/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_non-purebred_highly-admixed_qcov-${Q}.loco.mlma -O /seq/vgb/dd/gwas/summ/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_non-purebred_highly-admixed_qcov-${Q}_p-1e-6_kb-250_r2-0.2.loco.tsv

# python3 /seq/vgb/dd/gwas/bin/join_clmp_to_mlma.py -Q ${Q} -D "random sample" -C /seq/vgb/dd/gwas/clump/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_matched_N_qcov-${Q}_p-1e-6_kb-250_r2-0.2.loco.clumped.ranges -P ${P} -S /seq/vgb/dd/gwas/assoc/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_matched_N_qcov-${Q}.loco.mlma -O /seq/vgb/dd/gwas/summ/mutt_paper/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_${P}_matched_N_qcov-${Q}_p-1e-6_kb-250_r2-0.2.loco.tsv
