# This script checks the changes to the top SNPs after inclusion of top 10 PCs

# Libraries
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggrepel)

# Arguments
args <- commandArgs(trailingOnly=T)
mlmaFileA <- args[1]
mlmaFileB <- args[2]
traitString = paste(readLines(args[3]), collapse=" ")

print("Comparing results from the following GWAS...")
print(traitString)
print("original GWAS A:")
print(mlmaFileA)
print("new GWAS B:")
print(mlmaFileB)

# Read MLMA files containing GWAS summary statistics
gwasA = read_tsv(mlmaFileA, na = c("nan","-nan","NA","inf","-inf","NULL")) # original GWAS
gwasB = read_tsv(mlmaFileB, na = c("nan","-nan","NA","inf","-inf","NULL")) # new GWAS

# Perform SNP-by-SNP comparison
comp = gwasA %>%
  select(SNP,p,b) %>%
  filter(!is.na(p)) %>%
  merge((gwasB %>% select(SNP,p,b) %>% filter(!is.na(p))), by = "SNP")

# c(0, 10^(-20:-3), 5e-3,0.01,0.05,1)
# c(0, 10^(-20:-10), sort(c(5*10^(-9:-4),10^(-9:-4))), 0.001, 0.01, 1)

gwas_survival = comp %>%
  select(SNP,p.x,b.x,p.y,b.y) %>%
  mutate(loose = if_else(p.x < 1e-5, TRUE, FALSE),
         suggestive = if_else(p.x < 1e-6, TRUE, FALSE),
         significant = if_else(p.x < 5e-8, TRUE, FALSE),
         bin = cut(p.y, ordered_result = T, right = T,
                   breaks = c(0, 10^(-20:-10), sort(c(5*10^(-9:-4),10^(-9:-4))), 0.001, 0.01, 1))) %>%
  pivot_longer(cols = c("loose","suggestive","significant"),
               names_to = "threshold",
               values_to = "member") %>%
  group_by(threshold) %>%
  mutate(nTotalTop = sum(member == TRUE)) %>%
  mutate(threshold = factor(x = threshold,
                            levels = c("loose","suggestive","significant"))) %>%
  group_by(threshold,bin) %>%
  summarize(nBinTop = sum(member == TRUE),
            nTotalTop = unique(nTotalTop)) %>%
  group_by(threshold) %>%
  arrange(bin) %>%
  mutate(nBinTopCum = cumsum(nBinTop)) %>%
  mutate(propTopSurvive = (nTotalTop-nBinTopCum)/nTotalTop,
         propTopDead = (nBinTopCum)/nTotalTop) %>%
  arrange(threshold,bin) %>%
  mutate(trait = traitString)

write.table(x = gwas_survival,
            file = paste(mlmaFileB,".survival.csv", sep = ""),
            sep = ",",
            row.names = F)

p = gwas_survival %>%
  ggplot(aes(x = bin,
             y = propTopDead,
             group = threshold,
             color = threshold)) +
  coord_flip() +
  scale_y_reverse(limits = c(1,0)) +
  geom_vline(xintercept = "5e-8", alpha = 0.25, linetype = "dashed") +
  geom_vline(xintercept = "1e-6", alpha = 0.25, linetype = "dotted") +
  geom_line() +
  labs(x = "p-value bin\nin compared GWAS",
       y = "proportion of SNPs\nfrom original GWAS",
       color = "p-value of SNPs\nin original GWAS",
       subtitle = traitString) +
  theme_classic()

ggsave(filename = paste(mlmaFileB,".survival-curve.png", sep = ""),
       plot = p,
       device = "png",
       dpi = 300,
       unit = "in",
       width = 6,
       height = 6,
       bg = "white"
)