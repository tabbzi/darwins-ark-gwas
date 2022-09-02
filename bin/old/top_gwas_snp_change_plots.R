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
print("Comparing results from the following GWAS...")
print("original GWAS A:")
print(mlmaFileA)
print("new GWAS B:")
print(mlmaFileB)

# Defaults
summFile = "/seq/vgb/dd/gwas/summ/mutt_paper/DarwinsArk_20191115_ResultsGWAS_Updated_r2-0.8.tsv"

# Read MLMA files containing GWAS summary statistics
gwasA = read_tsv(mlmaFileA) # original GWAS
gwasB = read_tsv(mlmaFileB) # new GWAS

# Perform SNP-by-SNP comparison
comp = gwasA %>%
  select(SNP,p,b) %>%
  merge((gwasB %>% select(SNP,p,b)), by = "SNP") %>%
  mutate(delta_p = p.y - p.x,
         fold_p = (p.y - p.x)/p.x,
         delta_b = b.y - b.x,
         fold_b = (b.y - b.x)/b.x)

comp_summ = comp %>%
  summarize(mean_delta_p = mean(delta_p, na.rm = T),
            sd_delta_p = sd(delta_p, na.rm = T),
            mean_fold_p = mean(fold_p, na.rm = T),
            sd_fold_p = sd(fold_p, na.rm = T),
            mean_delta_b = mean(delta_b, na.rm = T),
            sd_delta_b = sd(delta_b, na.rm = T),
            mean_fold_b = mean(fold_b, na.rm = T),
            sd_fold_b = sd(fold_b, na.rm = T))
write.table(x = comp_summ,
            file = paste(mlmaFileB,".overall-delta.csv", sep = ""),
            sep = ",",
            row.names = F)


# Read summary file containing associated loci
summ = read_tsv(summFile)

# Identify which GWAS by matching of top SNP + association p-value from summary file to MLMA summary statistics
loci = summ %>%
  merge(gwasA,
        by = c("SNP","p","b")) %>%
  select(Index,
         Class,
         Trait,
         Set,
         N,
         Chromosome,
         Start,
         End,
         POS,
         SNP)

# Report statistics for top SNPs from original GWAS
delta = comp %>%
  filter(SNP %in% loci$SNP)

write.table(x = delta,
            file = paste(mlmaFileB,".top-delta.csv", sep = ""),
            sep = ",",
            row.names = F)

# Plots
# b effect size vs b effect size
p = comp %>%
  mutate(color = if_else(SNP %in% loci$SNP, "#972B28", "#000000"),
         alpha = if_else(SNP %in% loci$SNP, 1, 0.5),
         label = if_else(SNP %in% loci$SNP, SNP, NA_character_)) %>%
  ggplot(aes(x = b.x,
             y = b.y,
             color = color, 
             alpha = alpha,
             label = label)) +
  geom_abline(slope = 1, alpha = 0.5, linetype = "dashed") +
  geom_point(shape = 1) +
  geom_text_repel() +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(-1,1)) +
  scale_color_identity() +
  scale_alpha_identity() +
  labs(x = "effect size (b) in original GWAS",
       y = "effect size (b) in compared GWAS",
       title = unique(loci$Trait)) +
  theme_pubr()

ggsave(filename = paste(mlmaFileB,".overall-delta.b.png", sep = ""),
       plot = p,
       device = "png",
       dpi = 300,
       unit = "in",
       width = 4,
       height = 4,
       bg = "white"
)

# -log10(p) vs -log10(p)
maxlogp = max(max(-log10(comp$p.x)),
              max(-log10(comp$p.y)))
p = comp %>%
  mutate(color = if_else(SNP %in% loci$SNP, "#972B28", "#000000"),
         alpha = if_else(SNP %in% loci$SNP, 1, 0.5),
         label = if_else(SNP %in% loci$SNP, SNP, NA_character_)) %>%
  ggplot(aes(x = -log10(p.x),
             y = -log10(p.y),
             color = color,
             alpha = alpha,
             label = label)) +
  geom_abline(slope = 1, alpha = 0.5, linetype = "dashed") +
  geom_point(shape = 1) +
  geom_text_repel() +
  scale_x_continuous(limits = c(0,maxlogp+1)) +
  scale_y_continuous(limits = c(0,maxlogp+1)) +
  scale_color_identity() +
  scale_alpha_identity() +
  labs(x = "-log10(p) in original GWAS",
       y = "-log10(p) in compared GWAS",
       title = unique(loci$Trait)) +
  theme_pubr()

ggsave(filename = paste(mlmaFileB,".overall-delta.p.png", sep = ""),
       plot = p,
       device = "png",
       dpi = 300,
       unit = "in",
       width = 4,
       height = 4,
       bg = "white"
)