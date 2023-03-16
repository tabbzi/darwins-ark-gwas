# GWAS Plots for Publication

# load libraries
#install.packages("svglite")
library(svglite)
library(readr)
library(tidyverse)
library(data.table)
library(ggtext)
library(stringr)
library(dplyr)

# set directory
setwd("~/Dropbox (UMass Medical School)/Papers/cross-species/dat/gwa/")

# set input
mlma="mlm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-bq.017_dcov-datatype_qcov-age.hgt.bq.017_rel-cutoff-0.75.loco.mlma"
mlma="mlm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-bq.017_dcov-datatype_qcov-age.hgt.bq.017_rel-cutoff-0.75"

# load data
data = read_tsv(mlma) %>% as.data.table()
colnames(data) = c("chr","snp","pos","A1","A2","frq","b","se","p")

#data = data[order(chr, pos)]

clump = read_table("mlm/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-bq.017_dcov-datatype_qcov-age.hgt.bq.017_rel-cutoff-0.75_p-1e-06_kb-1000_r2-0.2.loco.clumped") %>% as.data.table()
clump$clump = lapply(gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", clump$SP2, perl=TRUE)
                     , function(x) strsplit(x, ",")[[1]])
data.clumps = data.frame(index = with(clump, rep(SNP, lengths(clump))),
snps = unlist(clump$clump))

# bins
sig.high <- 5e-8
sig.low <- 1e-6
#data$bin = cut(x = data$p, breaks = c(0,5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.1,1))

data = merge(data, data.clumps, by.x = "snp", by.y = "snps", all.x = T) %>% as.data.table()
data[snp %in% clump$SNP, top := snp]

# prune insignificant SNPs at 10% from each CHR for plotting
data.sig = data %>% subset(!is.na(index) | !is.na(top))
data.not = data %>%
  subset(p >= 1e-6 & is.na(index)) %>%
  group_by(chr) %>%
  sample_frac(0.10)
data = bind_rows(data.sig,data.not) %>% as.data.table()

# add cumulative pos
data$pos = as.numeric(data$pos)
data.cum <- data %>% 
  group_by(chr) %>% 
  summarise(max.pos = max(pos)) %>% 
  mutate(pos.add = lag(cumsum(max.pos), default = 0)) %>% 
  select(chr, pos.add)

data.plot <- data %>% 
  inner_join(data.cum, by = "chr") %>% 
  mutate(pos.cum = pos + pos.add)

# add axis info
axis.set <- data.plot %>% 
  group_by(chr) %>% 
  summarize(center = mean(pos.cum))

ylim <- data.plot %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

# set color, fill, shape, alpha
data.plot[p >= sig.low & chr %% 2 == 0, color := "#707070"]
data.plot[p >= sig.low & chr %% 2 != 0, color := "#8f8f8f"]
data.plot[p >= sig.low, alpha := 0.75]
data.plot[p >= sig.low, shape := 16]
data.plot[p >= sig.low, size := 0.5]
data.plot[p >= sig.low, plot.group := "A. insignif"]

# index SNPs
data.plot[top == snp, alpha := 1.00]
data.plot[top == snp, shape := 23]
data.plot[top == snp, size := 2]
data.plot[top == snp, color := "#000000"]
data.plot[top == snp, fill := "#972B28"]
data.plot[top == snp, label := snp]
data.plot[top == snp, plot.group := "C. index SNP"]

# clumped SNPs
data.plot[!is.na(index), alpha := 0.75]
data.plot[!is.na(index), shape := 18]
data.plot[!is.na(index), size := 1]
data.plot[!is.na(index), color := "#a12d2b"]
data.plot[!is.na(index), plot.group := "B. clumped SNP"]

# plot!
# full 
p = ggplot(data.plot %>% arrange(plot.group), aes(x = pos.cum,
                      y = -log10(p),
                      size = size,
                      alpha = alpha,
                      shape = shape,
                      color = color,
                      fill = fill,
                      label = label)) +
  geom_hline(yintercept = -log10(sig.low), color = "grey40", linetype = "dotted") + 
  geom_hline(yintercept = -log10(sig.high), color = "grey40", linetype = "dashed") + 
  geom_point() +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_size_identity() +
  #scale_size_continuous(range = c(0.5,2)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

# save full as tiff
ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/cross-species/fig/Q242A_manhattan_full_p-5e-08_under-pruned-0.10.tiff",
       plot = p,
       device = "tiff",
       units = "in",
       dpi = 300,
       width = 9,
       height = 5)

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/cross-species/fig/Q242A_manhattan_full_p-5e-08_figure.tiff",
       plot = p,
       device = "tiff",
       units = "in",
       dpi = 300,
       width = 4.75*2,
       height = 1.5*2)

# save unclumped SNPs as tiff
p = ggplot(data.plot %>% arrange(plot.group) %>% subset(is.na(index) & is.na(top)), aes(x = pos.cum,
                                                  y = -log10(p),
                                                  alpha = alpha,
                                                  size = size,
                                                  shape = shape,
                                                  color = color,
                                                  fill = fill,
                                                  label = label)) +
  geom_point(size = 0.5) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/cross-species/fig/Q242A_manhattan_low_p-5e-08_under-pruned-0.10.tiff",
       plot = p,
       device = "tiff",
       units = "in",
       dpi = 300,
       width = 9,
       height = 5)

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/cross-species/fig/Q242A_manhattan_low_p-5e-08_under-pruned-0.10_figure.tiff",
       plot = p,
       device = "tiff",
       units = "in",
       dpi = 300,
       width = 4.75*2,
       height = 1.5*2)

# save clumped SNPs as vector
data.plot.subset = copy(data.plot)
data.plot.subset[is.na(index) & is.na(top), p := NA]

p = ggplot(data.plot.subset %>% arrange(plot.group), aes(x = pos.cum,
                                                  y = -log10(p),
                                                  size = size,
                                                  alpha = alpha,
                                                  shape = shape,
                                                  color = color,
                                                  fill = fill,
                                                  label = label)) +
  geom_hline(yintercept = -log10(sig.low),
             color = "grey40",
             linetype = "dotted") + 
  geom_hline(yintercept = -log10(sig.high),
             color = "grey40",
             linetype = "dashed") + 
  geom_point() +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center, expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  #scale_size_continuous(range = c(0.5,2)) +
  scale_size_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/cross-species/fig/Q242A_manhattan_top_p-5e-08.svg",
       plot = p,
       device = "svg",
       units = "in",
       width = 9,
       height = 5)

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/cross-species/fig/Q242A_manhattan_top_p-5e-08_figure.svg",
       plot = p,
       device = "svg",
       units = "in",
       width = 4.75*2,
       height = 1.5*2)

#### other GWAS ####
filepath = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/plots/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_CCD_RRBs_summedFreq_norm_qcov-CCD_RRBs_summedFreq_norm"

filepath = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/plots/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_CCD_RRBs_sums_qcov-CCD_RRBs_sums"

filepath = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/plots/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_mixed_breed-labrador_retriever_Q54A_qcov-Q54A"

filepath = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/plots/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_Q52A_qcov-Q52A"

filepath = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/plots/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_Q144A_freq_qcov-Q144A"

filepath = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/data/gwas/plots/DarwinsArk_gp-0.70_snps-only_maf-0.02_geno-0.20_hwe-midp-1e-20_het-0.25-1.00_N-2155_chr-pos-ref-alt_CCD_RRBs_CaseCont_144_qcov-Q144A"

# load data
data = read.table(paste(filepath,".loco.mlma", sep = ""),
                  header = T,
                  col.names = c("chr","snp","pos","A1","A2","frq","b","se","p")) %>% as.data.table()

clump = read.table(paste(filepath,"_p-1e-06_kb-250_r2-0.5.loco.clumped", sep = ""),
                   header = T) %>% as.data.table()
clump$clump = lapply(gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))", "", clump$SP2, perl=TRUE)
                     , function(x) strsplit(x, ",")[[1]])
data.clumps = data.frame(index = with(clump, rep(SNP, lengths(clump))),
                         snps = unlist(clump$clump))

# bins
sig.high <- 5e-8
sig.low <- 1e-6
#data$bin = cut(x = data$p, breaks = c(0,5e-8,1e-6,1e-5,1e-4,1e-3,0.01,0.1,1))

data = merge(data, data.clumps, by.x = "snp", by.y = "snps", all.x = T) %>% as.data.table()
data[snp %in% clump$SNP, top := snp]

# or, if no sig clumps:
#data$index = NA
#data$top = NA

# prune insignificant SNPs at 10% from each CHR for plotting
data.sig = data %>% subset(!is.na(index) | !is.na(top))
data.not = data %>%
  subset(p >= 1e-6 & is.na(index)) %>%
  group_by(chr) %>%
  sample_frac(0.10)
data.plot = bind_rows(data.sig,data.not) %>% as.data.table()

# add cumulative pos
data.plot$pos = as.numeric(data.plot$pos)
data.cum <- data.plot %>% 
  group_by(chr) %>% 
  summarise(max.pos = max(pos)) %>% 
  mutate(pos.add = lag(cumsum(max.pos), default = 0)) %>% 
  select(chr, pos.add)

data.plot <- data.plot %>% 
  inner_join(data.cum, by = "chr") %>% 
  mutate(pos.cum = pos + pos.add)

# add axis info
axis.set <- data.plot %>% 
  group_by(chr) %>% 
  summarize(center = mean(pos.cum))

ylim <- data.plot %>% 
  filter(p == min(p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

# set color, fill, shape, alpha
data.plot[p >= sig.low & chr %% 2 == 0, color := "#707070"]
data.plot[p >= sig.low & chr %% 2 != 0, color := "#8f8f8f"]
data.plot[p >= sig.low, alpha := 0.75]
data.plot[p >= sig.low, shape := 16]
data.plot[p >= sig.low, size := 0.5]
data.plot[p >= sig.low, plot.group := "A. insignif"]

# index SNPs
data.plot[top == snp, alpha := 1.00]
data.plot[top == snp, shape := 23]
data.plot[top == snp, size := 2]
data.plot[top == snp, color := "#000000"]
data.plot[top == snp, fill := "#972B28"]
data.plot[top == snp, label := snp]
data.plot[top == snp, plot.group := "C. index SNP"]

# clumped SNPs
data.plot[!is.na(index), alpha := 0.75]
data.plot[!is.na(index), shape := 18]
data.plot[!is.na(index), size := 1]
data.plot[!is.na(index), color := "#a12d2b"]
data.plot[!is.na(index), plot.group := "B. clumped SNP"]

# plot!
# full 
p = ggplot(data.plot %>% arrange(plot.group), aes(x = pos.cum,
                                                  y = -log10(p),
                                                  size = size,
                                                  alpha = alpha,
                                                  shape = shape,
                                                  color = color,
                                                  fill = fill,
                                                  label = label)) +
  geom_hline(yintercept = -log10(sig.low), color = "grey40", linetype = "dotted") + 
  geom_hline(yintercept = -log10(sig.high), color = "grey40", linetype = "dashed") + 
  geom_point() +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_size_identity() +
  #scale_size_continuous(range = c(0.5,2)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

# save full as tiff
ggsave(filename = paste(filepath,".tiff",sep=""),
       plot = p,
       device = "tiff",
       units = "in",
       dpi = 300,
       width = 10,
       height = 5)

library(ggpubr)
library(devtools)
devtools::install_github("aloy/qqplotr")
library(qqplotr)

gg <- ggplot(data = data, mapping = aes(sample = p)) +
  stat_qq_band(detrend = T) +
  stat_qq_line(detrend = T) +
  stat_qq_point(detrend = T) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")

ggplot(data = data.plot, aes(sample = p)) +
  coord_flip(xlim = c()) +
  stat_qq() +
  stat_qq_line() +
  theme_pubclean()

# save unclumped SNPs as tiff
p = ggplot(data.plot %>% arrange(plot.group) %>% subset(is.na(index) & is.na(top)), aes(x = pos.cum,
                                                                                        y = -log10(p),
                                                                                        alpha = alpha,
                                                                                        size = size,
                                                                                        shape = shape,
                                                                                        color = color,
                                                                                        fill = fill,
                                                                                        label = label)) +
  geom_point(size = 0.5) +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  scale_size_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/figures/Q121A_all-dogs_manhattan_low_p-1e-6_under-pruned-0.10.tiff",
       plot = p,
       device = "tiff",
       units = "in",
       dpi = 300,
       width = 9,
       height = 5)

# save clumped SNPs as vector
data.plot.subset = copy(data.plot)
data.plot.subset[is.na(index) & is.na(top), p := NA]

p = ggplot(data.plot.subset %>% arrange(plot.group), aes(x = pos.cum,
                                                         y = -log10(p),
                                                         size = size,
                                                         alpha = alpha,
                                                         shape = shape,
                                                         color = color,
                                                         fill = fill,
                                                         label = label)) +
  geom_hline(yintercept = -log10(sig.low),
             color = "grey40",
             linetype = "dotted") + 
  geom_hline(yintercept = -log10(sig.high),
             color = "grey40",
             linetype = "dashed") + 
  geom_point() +
  scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  #scale_size_continuous(range = c(0.5,2)) +
  scale_size_identity() +
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  labs(x = NULL, 
       y = "-log<sub>10</sub>(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

ggsave(filename = "~/Dropbox (UMass Medical School)/Papers/mutt-paper/figures/Q121A_all-dogs_manhattan.svg",
       plot = p,
       device = "svg",
       units = "in",
       width = 9,
       height = 5)