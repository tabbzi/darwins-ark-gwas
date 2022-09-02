library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly=T)
file_dir <- args[1]
file_prefix <- args[2]

file_names <- list.files(path = file_dir, pattern = paste("(",file_prefix, ")*.score.ld$", sep = ""), full.names=TRUE)
ld_scores <- lapply(file_names, read_table2)
ld_scores = rbindlist(ld_scores)

Q25 = quantile(ld_scores$ldscore_SNP, 0.25)
print(Q25)
Q50 = quantile(ld_scores$ldscore_SNP, 0.50)
print(Q50)
Q75 = quantile(ld_scores$ldscore_SNP, 0.75)
print(Q75)

ld1_snp = ld_scores %>% filter(ldscore_SNP <= Q25) %>% pull(SNP)
ld2_snp = ld_scores %>% filter(ldscore_SNP > Q25 & ldscore_SNP <= Q50) %>% pull(SNP)
ld3_snp = ld_scores %>% filter(ldscore_SNP > Q50 & ldscore_SNP <= Q75) %>% pull(SNP)
ld4_snp = ld_scores %>% filter(ldscore_SNP > Q75) %>% pull(SNP)

write_lines(ld1_snp, paste(file_dir,"/",file_prefix, ".score.ld.Q000-Q025.txt", sep = ""))
write_lines(ld2_snp, paste(file_dir,"/",file_prefix, ".score.ld.Q025-Q050.txt", sep = ""))
write_lines(ld3_snp, paste(file_dir,"/",file_prefix, ".score.ld.Q050-Q075.txt", sep = ""))
write_lines(ld4_snp, paste(file_dir,"/",file_prefix, ".score.ld.Q075-Q100.txt", sep = ""))
