library(svglite)
library(tidyverse)
library(data.table)
library(ggtext)
library(stringr)
library(dplyr)

# arguments
args <- commandArgs(trailingOnly=T)
mlmaFile <- args[1]
clumpFile <- args[2]
title <- paste(readLines(args[3]), collapse=" ")

# load data
data = read.table(mlmaFile,
                  header = T,
                  col.names = c("chr","snp","pos","A1","A2","frq","b","se","p")) %>% as.data.table()

clump = read.table(clumpFile,
                  header = T) %>% as.data.table()
clumpSNPs = lapply(gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))",
                          "",
                          clump$SP2,
                          perl=TRUE), function(x) strsplit(x, ",")[[1]])
data.clumps = data.frame(index = with(clump, rep(SNP, lengths(clumpSNPs))),
snps = unlist(clumpSNPs))

# bins
sig.high <- 5e-8
sig.low <- 1e-6

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

# plot
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
  scale_color_identity() +
  scale_fill_identity() +
  scale_shape_identity() +
  scale_alpha_identity() +
  labs(x = NULL,
       y = "-log<sub>10</sub>(p)",
       title = as.character(title)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.y = element_markdown(),
    axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
  )

# save full as png
ggsave(filename = paste(mlmaFile, ".manhattan-plot_full_under-pruned-0.10.png", sep = ""),
       plot = p,
       device = "png",
       units = "in",
       dpi = 300,
       width = 4.75*2,
       height = 1.5*2)
