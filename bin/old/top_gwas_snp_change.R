# This script checks the changes to the top SNPs after inclusion of top 10 PCs

# Libraries
library(tidyverse)
library(data.table)
library(ggpubr)

# Arguments
args <- commandArgs(trailingOnly=T)
mlmaFileA <- args[1]
mlmaFileB <- args[2]

# Defaults

summFile = "/seq/vgb/dd/gwas/summ/mutt_paper/DarwinsArk_20191115_ResultsGWAS_Updated_r2-0.8.tsv"

# read MLMA files containing GWAS summary statistics
gwasA = read_tsv(mlmaFileA) # original GWAS
gwasB = read_tsv(mlmaFileB) # new GWAS

# read summary file containing associated loci
summ = read_tsv(summFile)

# identify which GWAS by matching of top SNP + association p-value from summary file to MLMA summary statistics
loci = summ %>%
  merge(gwasA,
        by = c("SNP","p","b")) %>%
  select(Index, Class, Trait, Set, N, Chromosome, Start, End, POS, SNP, p, b)

# delete MLMA from memory
gwasA = NULL

delta = loci %>%
  merge((gwasB %>% select(SNP,p,b)), by = "SNP") %>%
  mutate(delta_p = p.y - p.x,
         fold_p = (p.y - p.x)/p.x,
         delta_logp = -log10(p.y) + log10(p.x),
         fold_logp = (-log10(p.y) + log10(p.x))/(-log10(p.x)),
         delta_b = b.y - b.x,
         fold_b = (b.y - b.x)/b.x)

write.table(x = delta, file = paste(mlmaFileB,".delta.csv", sep = ""), sep = ",", row.names = F)