# Required libraries
require(tidyverse)
require(ggtext)
require(ggpubr)
require(data.table)
require(stringr)
require(dplyr)
require(svglite)
require(colorspace)
require(argparse)

# Define the arguments that the script accepts
parser <- ArgumentParser()
parser$add_argument("--mlma", type = str, help = "Path to the mlma file", default = "/seq/vgb/dd/gwas/assoc/2022-11-20/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-mq.121_dcov-datatype_qcov-age_rel-cutoff-0.75.mlma")
parser$add_argument("--clump", type = str, help = "Path to the clump file", default = "/seq/vgb/dd/gwas/clump/2022-11-20/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-mq.121_dcov-datatype_qcov-age_rel-cutoff-0.75_p-1e-6_kb-1000_r2-0.2.loco.clumped")
parser$add_argument("--title", type = str, help = "Title of the plot", default = "Q121")
parser$add_argument("--width", type = float, help = "Width of the output plot", default = 4.75 * 2)
parser$add_argument("--height", type = float, help = "Height of the output plot", default = 3 * 2)

# Parse the arguments passed to the script
args <- parser$parse_args()

# Use the arguments in your script
mlmaFile <- args$mlma
clumpFile <- args$clump
title <- args$title
width <- args$width
height <- args$height

if (file.exists(clumpFile)){
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
} else
{
  # Define the column names for the empty data frame
  column_names <- c("CHR", "F", "SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2")
  
  # Create the empty data frame
  empty_data_frame <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(empty_data_frame) <- column_names
  
  # Use the empty data frame in place of clumpFile
  clump = empty_data_frame
  clumpSNPs = lapply(gsub("\\s*(\\([^()]*(?:(?1)[^()]*)*\\))",
                          "",
                          clump$SP2,
                          perl=TRUE), function(x) strsplit(x, ",")[[1]])
  data.clumps = data.frame(index = with(clump, rep(SNP, lengths(clumpSNPs))),
                           snps = unlist(clumpSNPs))
}

data = read.table(mlmaFile,
                  header = T,
                  col.names = c("chr","snp","pos","A1","A2","frq","b","se","p")) %>% as.data.table()

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

ylim <- abs(floor(log10(min(data.plot$p, na.rm = T)))) + 2

if (ylim < 10){
  ylim = 10
}

# set color, fill, shape, alpha
data.plot = as.data.table(data.plot)

# Set color for even chromosomes
data.plot[chr %% 2 == 0,
          color := "#707070"
]

# Set color for odd chromosomes
data.plot[chr %% 2 != 0,
          color := "#8f8f8f"
]

# Set default alpha, shape, size, and plot group
data.plot[,
          alpha := 0.75,
          shape := 16,
          size := 0.5,
          plot.group := "A. insignif or unclumped"
]

# index SNPs with p-values between sig.low and sig.high
data.plot[top == snp & p < sig.low & p > sig.high,
          alpha := 1.00,
          shape := 23,
          size := 2,
          color := "#000000",
          fill := "#284697",
          label := snp,
          plot.group := "D. index SNP low"
]

# index SNPs with p-values less than sig.high
data.plot[top == snp & p < sig.high,
          alpha := 1.00,
          shape := 23,
          size := 2,
          color := "#000000",
          fill := "#972B28",
          label := snp,
          plot.group := "D. index SNP high"
]

# clumped SNPs with p-values less than sig.low
data.plot[index %in% clump.low,
          alpha := 0.75,
          shape := 18,
          size := 1,
          color := "#2b48a1",
          plot.group := "B. clumped SNP low"
]

# clumped SNPs with p-values less than sig.high
data.plot[index %in% clump.high,
          alpha := 0.75,
          shape := 18,
          size := 1,
          color := "#a12d2b",
          plot.group := "C. clumped SNP high"
]

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
  # Calculate the number of non-NA elements in ps
  N <- length(ps[!is.na(ps)])
  
  # Create a data frame with the observed, expected, and confidence interval values
  df <- data.frame(
    observed = -log10(sort(ps[!is.na(ps)])),
    expected = -log10(ppoints(N)),
    clower   = -log10(qbeta((1 - ci) / 2, 1:N, N - 1:N + 1)),
    cupper   = -log10(qbeta((1 + ci) / 2, 1:N, N - 1:N + 1))
  )
  
  # Calculate the genomic inflation factor
  lambda <- median(df$observed/df$expected)
  
  # Create labels for the x and y axes
  log10Pe <- expression(paste("Expected -log"[10], plain(P)))
  log10Po <- expression(paste("Observed -log"[10], plain(P)))
  
  # Create the ggplot object
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 3) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    xlab(log10Pe) +
    ylab(log10Po) +
    labs(caption = paste("Genomic inflation factor =", lambda))
}

gg.qq = gg_qqplot(data$p) +
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=10,face="bold"), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,2,1,1), "lines")) +
  labs(title = as.character(title)) +

ggsave(filename = paste(mlmaFile, ".qq-plot.png", sep = ""),
       plot = gg.qq,
       device = "png",
       units = "in",
       dpi = 300,
       width = 4.75*2,
       height = 4.75*2)
