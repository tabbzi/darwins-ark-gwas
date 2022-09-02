join <(sort pheno/${1}.pheno ) <(sort covar/sex_status_type.covar) |\
join - <(sort covar/height.qcovar ) |\
join - <(sort geno/darwins_dogs_merged_maf-0.01_gen-0.10.fam ) |\
join - <(sort grm/darwins_dogs_merged_maf-0.01_gen-0.10_grm-cutoff_0.5.grm.id ) |\
grep -vf geno/remove_younger_10mo.txt | awk '{print $1}' | sort | join - <(sort pheno/${1}.pheno) | awk '{print $3}' | sort | uniq -c

join <(sort pheno/${1}.pheno ) <(sort covar/sex_status_type.covar) |\
join - <(sort covar/height.qcovar ) |\
join - <(sort geno/darwins_dogs_merged_maf-0.01_gen-0.10.fam ) |\
join - <(sort grm/darwins_dogs_merged_maf-0.01_gen-0.10_grm-cutoff_0.5.grm.id ) |\
grep -vf geno/remove_younger_10mo.txt | awk '{print $1}' | sort | join - <(sort pheno/${1}.pheno) | wc -l
