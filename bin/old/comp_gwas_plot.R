library(tidyverse)
library(ggplot2)
library(data.table)
library(ggpubr)

comp = read.delim("~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/comp/Q21A_uncond_vs_cond-Q126A.csv",
                  sep = ",",
                  header = T) %>% as.data.table()
pheno = "Q21_cond-Q127A"
titleGWAS = "Q21 Focus in distracting situation"
labX = "-log10(p) unconditioned"
labY = "-log10(p) conditioned on\ncoat length phenotype"

comp = read.delim("~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/comp/Q21A_uncond_vs_cond-SNP.csv",
                  sep = ",",
                  header = T) %>% as.data.table()

pheno = "Q21_cond-SNP"
titleGWAS = "Q21 Focus in distracting situation"
labX = "-log10(p) unconditioned"
labY = "-log10(p) conditioned on\ncoat length SNP"

library(RColorBrewer)
n = 38
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

thisPalette = sample(col_vector, n)
pie(rep(1,n), col=thisPalette)

comp[, SNPlabel := NA]
comp[Chr == 32 & log10_p_x > 6, SNPlabel := SNP]

#prob = 0.001
#comp[(p_x <= quantile(comp$p_x, probs =  prob) & p_y > quantile(comp$p_x, probs =  prob)) | (p_x > quantile(comp$p_y, probs =  prob) & p_y <= quantile(comp$p_y, probs =  prob)), SNPlabel := SNP]

p = ggplot(comp[Chr == 32], aes(x = log10_p_x,
                 y = log10_p_y,
                 fill = as.factor(Chr),
                 label = SNPlabel,
                 shape = as.factor(Chr))) +
  geom_point() +
  #stat_cor(data = comp, aes(x = log10_p_x, y = log10_p_y), method="pearson", inherit.aes = F, p.accuracy = 1e-8) +
  stat_cor(method="pearson", p.accuracy = 1e-8) +
  geom_text_repel() +
  coord_equal() +
  guides(fill = F, shape = F) +
  scale_shape_manual(values = rep(c(21,22,23),n)) +
  scale_fill_manual(values = thisPalette) +
  scale_x_continuous(limits = c(3,max(comp$log10_p_x,comp$log10_p_y))) +
  scale_y_continuous(limits = c(3,max(comp$log10_p_x,comp$log10_p_y))) +
  labs(x = labX,
       y = labY,
       title = titleGWAS) +
  theme_pubclean()

ggsave(filename = paste("~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/comp/",pheno,"_chr32",".pdf",sep=""), plot = p, device = "pdf", width = 6, height = 6, units = "in")
