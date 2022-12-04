require(tidyverse)
require(ggtext)
require(data.table)
require(stringr)
require(dplyr)
require(svglite)
require(colorspace)

# arguments
args <- commandArgs(trailingOnly=T)
mlmaFile <- args[1]
clumpFile <- args[2]
title <- paste(readLines(args[3]), collapse=" ")

if (file.exists(clumpFile)){
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
  
  # snp bins
  sig.high <- 5e-8
  sig.low <- 1e-6
  
  # clump bins
  clump.high = clump %>% filter(P < sig.high) %>% pull(SNP)
  clump.low = clump %>% filter(P < sig.low & P > sig.high) %>% pull(SNP)
  
  data = merge(data, data.clumps, by.x = "snp", by.y = "snps", all.x = T) %>% as.data.table()
  data[snp %in% clump$SNP, top := snp]
  
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
  
  if (ylim < 10){
    ylim = 10
  }
  
  # set color, fill, shape, alpha
  data.plot = as.data.table(data.plot)
  data.plot[chr %% 2 == 0, color := "#707070"]
  data.plot[chr %% 2 != 0, color := "#8f8f8f"]
  data.plot[, alpha := 0.75]
  data.plot[, shape := 16]
  data.plot[, size := 0.5]
  data.plot[, plot.group := "A. insignif or unclumped"]
  
  
  # index SNPs
  data.plot[top == snp & p < sig.low & p > sig.high, alpha := 1.00]
  data.plot[top == snp & p < sig.low & p > sig.high, shape := 23]
  data.plot[top == snp & p < sig.low & p > sig.high, size := 2]
  data.plot[top == snp & p < sig.low & p > sig.high, color := "#000000"]
  data.plot[top == snp & p < sig.low & p > sig.high, fill := "#284697"]
  data.plot[top == snp & p < sig.low & p > sig.high, label := snp]
  data.plot[top == snp & p < sig.low & p > sig.high, plot.group := "D. index SNP low"]
  
  # index SNPs
  data.plot[top == snp & p < sig.high, alpha := 1.00]
  data.plot[top == snp & p < sig.high, shape := 23]
  data.plot[top == snp & p < sig.high, size := 2]
  data.plot[top == snp & p < sig.high, color := "#000000"]
  data.plot[top == snp & p < sig.high, fill := "#972B28"]
  data.plot[top == snp & p < sig.high, label := snp]
  data.plot[top == snp & p < sig.high, plot.group := "D. index SNP high"]
  
  # clumped SNPs level 1
  data.plot[!is.na(index) & index %in% clump.low, alpha := 0.75]
  data.plot[!is.na(index) & index %in% clump.low, shape := 18]
  data.plot[!is.na(index) & index %in% clump.low, size := 1]
  data.plot[!is.na(index) & index %in% clump.low, color := "#2b48a1"]
  data.plot[!is.na(index) & index %in% clump.low, plot.group := "B. clumped SNP low"]
  
  # clumped SNPs level 2
  data.plot[!is.na(index) & index %in% clump.high, alpha := 0.75]
  data.plot[!is.na(index) & index %in% clump.high, shape := 18]
  data.plot[!is.na(index) & index %in% clump.high, size := 1]
  data.plot[!is.na(index) & index %in% clump.high, color := "#a12d2b"]
  data.plot[!is.na(index) & index %in% clump.high, plot.group := "C. clumped SNP high"]
  
  
  print(data.plot[is.na(p)])
  print(data.plot[is.na(-log10(p))])
  print(data.plot[is.na(size)])
  print(data.plot[is.na(alpha)])
  print(data.plot[is.na(shape)])
  print(data.plot[is.na(color)])
  print(data.plot[is.na(fill)])
  
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
    scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
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
  ggsave(filename = paste(mlmaFile, ".manhattan-plot.png", sep = ""),
         plot = p,
         device = "png",
         units = "in",
         dpi = 300,
         width = 4.75*2,
         height = 3*2)
  
  # qq plot
  gg_qqplot <- function(ps, ci = 0.95) {
    N  <- length(ps)
    df <- data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(ppoints(N)),
      clower   = -log10(qbeta((1 - ci) / 2, 1:N, N - 1:N + 1)),
      cupper   = -log10(qbeta((1 + ci) / 2, 1:N, N - 1:N + 1))
    )
    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
      geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
      geom_line(aes(expected, cupper), linetype = 2) +
      geom_line(aes(expected, clower), linetype = 2) +
      xlab(log10Pe) +
      ylab(log10Po)
  }
  
  gg.qq = gg_qqplot(data$p) +
    coord_fixed(ratio = 1) +
    theme(axis.text=element_text(size=10,face="bold"), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(1,2,1,1), "lines"))
  
  ggsave(filename = paste(mlmaFile, ".qq-plot.png", sep = ""),
         plot = gg.qq,
         device = "png",
         units = "in",
         dpi = 300,
         width = 4.75*2,
         height = 4.75*2)
  
} else {
  # load data
  data = read.table(mlmaFile,
                    header = T,
                    col.names = c("chr","snp","pos","A1","A2","frq","b","se","p")) %>% as.data.table()
  
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
  
  if (ylim < 10){
    ylim = 10
  }
  
  # set color, fill, shape, alpha
  data.plot = as.data.table(data.plot)
  data.plot[chr %% 2 == 0, color := "#707070"]
  data.plot[chr %% 2 != 0, color := "#8f8f8f"]
  data.plot[, alpha := 0.75]
  data.plot[, shape := 16]
  data.plot[, size := 0.5]
  data.plot[, plot.group := "A. insignif or unclumped"]
  
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
    scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
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
  ggsave(filename = paste(mlmaFile, ".manhattan-plot.png", sep = ""),
         plot = p,
         device = "png",
         units = "in",
         dpi = 300,
         width = 4.75*2,
         height = 3*2)
  
  # qq plot
  gg_qqplot <- function(ps, ci = 0.95) {
    N  <- length(ps)
    df <- data.frame(
      observed = -log10(sort(ps)),
      expected = -log10(ppoints(N)),
      clower   = -log10(qbeta((1 - ci) / 2, 1:N, N - 1:N + 1)),
      cupper   = -log10(qbeta((1 + ci) / 2, 1:N, N - 1:N + 1))
    )
    log10Pe <- expression(paste("Expected -log"[10], plain(P)))
    log10Po <- expression(paste("Observed -log"[10], plain(P)))
    ggplot(df) +
      geom_point(aes(expected, observed), shape = 1, size = 3) +
      geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
      geom_line(aes(expected, cupper), linetype = 2) +
      geom_line(aes(expected, clower), linetype = 2) +
      xlab(log10Pe) +
      ylab(log10Po)
  }
  
  gg.qq = gg_qqplot(data$p) +
    coord_fixed(ratio = 1) +
    theme(axis.text=element_text(size=10,face="bold"), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(1,2,1,1), "lines"))
  
  ggsave(filename = paste(mlmaFile, ".qq-plot.png", sep = ""),
         plot = gg.qq,
         device = "png",
         units = "in",
         dpi = 300,
         width = 4.75*2,
         height = 4.75*2)
  
}
