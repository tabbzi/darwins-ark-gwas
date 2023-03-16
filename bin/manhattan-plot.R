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
parser$add_argument("--mlma", 
                    help = "Path to the mlma file", 
                    default = "/seq/vgb/dd/gwas/assoc/2022-11-20/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-mq.121_dcov-datatype_qcov-age_rel-cutoff-0.75.loco.mlma")
parser$add_argument("--clump", 
                    help = "Path to the clump file", 
                    default = "/seq/vgb/dd/gwas/clump/2022-11-20/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-mq.121_dcov-datatype_qcov-age_rel-cutoff-0.75_p-1e-6_kb-1000_r2-0.2.loco.clumped")
parser$add_argument("--title", 
                    help = "Title of the plot", 
                    default = "Q121")
parser$add_argument("--width", 
                    help = "Width of the output plot", 
                    default = 4.75 * 2)
parser$add_argument("--height", 
                    help = "Height of the output plot", 
                    default = 3 * 2)
parser$add_argument("--manhattan", 
                    help = "Generate Manhattan plot?", 
                    default = TRUE)
parser$add_argument("--qq", 
                    help = "Generate Q-Q plot?", 
                    default = TRUE)
parser$add_argument("--inflation", 
                    help = "Generate genomic inflation and variant frequency plot?", 
                    default = TRUE)

# Parse the arguments passed to the script
args <- parser$parse_args()
print(args)

# Use the arguments in your script
mlmaFile <- args$mlma
clumpFile <- args$clump
title <- args$title
width <- args$width
height <- args$height

if (file.exists(clumpFile)){
  #clump = read.table(clumpFile,
  #                   header = T) %>% as.data.table()
  clump = read_table(clumpFile) %>% as.data.table()
  
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

#data = read.table(mlmaFile,
#                  header = T,
#                  col.names = c("chr","snp","pos","A1","A2","frq","b","se","p")) %>% as.data.table()

data = read_tsv(mlmaFile) %>% as.data.table()
colnames(data) = c("chr","snp","pos","A1","A2","frq","b","se","p")

# Manhattan Plot ----
if (args$manhattan == TRUE) {
  print("Running Manhattan plot...")
  data = merge(data, 
               data.clumps, 
               by.x = "snp", 
               by.y = "snps", 
               all.x = T) %>% as.data.table()
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
  # Set color for even chromosomes
  data.plot[chr %% 2 == 0,
            c("color") := list("#707070")
  ]
  
  # Set color for odd chromosomes
  data.plot[chr %% 2 != 0,
            c("color") := list("#8f8f8f")
  ]
  
  # Set default alpha, shape, size, and plot group
  data.plot[,
            c("alpha", "shape", "size", "plot.group") := list(0.75, 16, 0.5, "A. insignif or unclumped")
  ]
  
  # index SNPs with p-values between sig.low and sig.high
  data.plot[top == snp & p < sig.low & p > sig.high,
            c("alpha", "shape", "size", "color", "fill", "label", "plot.group") := list(1.00, 23, 2, "#000000", "#284697", snp, "D. index SNP low")
  ]
  
  # index SNPs with p-values less than sig.high
  data.plot[top == snp & p < sig.high,
            c("alpha", "shape", "size", "color", "fill", "label", "plot.group") := list(1.00, 23, 2, "#000000", "#972B28", snp, "D. index SNP high")
  ]
  
  # clumped SNPs with p-values less than sig.low
  data.plot[index %in% clump.low,
            c("alpha", "shape", "size", "color", "plot.group") := list(0.75, 18, 1, "#2b48a1", "B. clumped SNP low")
  ]
  
  # clumped SNPs with p-values less than sig.high
  data.plot[index %in% clump.high,
            c("alpha", "shape", "size", "color", "plot.group") := list(0.75, 18, 1, "#a12d2b", "C. clumped SNP high")
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
         width = 6*2,
         height = 3*2)
}

# Q-Q Plot ----

if (args$qq == TRUE){
  print("Running Q-Q plot...")
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
      labs(caption = paste("Genomic inflation factor =", round(lambda, digits = 2)))
  }
  
  gg.qq = gg_qqplot(data$p) +
    coord_fixed(ratio = 1) +
    theme(axis.text=element_text(size=10,face="bold"), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          plot.margin = unit(c(1,2,1,1), "lines")) +
    labs(title = as.character(title))
  
  ggsave(filename = paste(mlmaFile, ".qq-plot.png", sep = ""),
         plot = gg.qq,
         device = "png",
         units = "in",
         dpi = 300,
         width = 4.75*2,
         height = 4.75*2)
}

# Genomic Inflation Plot ----
if (args$inflation == TRUE) {
  print("Running genomic inflation plot...")
  gwa.plot = data %>% 
    select(frq,p) %>%
    mutate(bin = cut(frq, 
                     breaks = 10^(seq(from = -3.0, to = -0.3, by = 0.1)))) %>% 
    group_by(bin) %>% 
    mutate(lambda = qchisq(mean(p), df = 1) / qchisq(0.5, df = 1),
           maxfreq = max(frq)) %>%
    ungroup() %>%
    group_by(bin,lambda,maxfreq) %>%
    summarize(n = n())
  
  min_lambda = gwa.plot %>%
    ungroup() %>%
    pull(lambda) %>%
    min(na.rm = T)
  
  min_lambda = if (min_lambda < 0.80) {
    round(min_lambda-0.05, digits = 2)
  } else {
    0.80
  }
  
  max_lambda = gwa.plot %>%
    ungroup() %>%
    pull(lambda) %>%
    max(na.rm = T)
  
  max_lambda = if (max_lambda < 1.20) {
    round(max_lambda+0.05, digits = 2)
  } else {
    1.20
  }

  
  p1 = ggplot(gwa.plot,
              aes(x = -log10(maxfreq),
                  y = lambda)) +
    geom_hline(yintercept = 1, color = "red", linetype = "twodash") +
    geom_line() +
    labs(y = "genomic inflation") +
    scale_y_continuous(limits = c(min_lambda,max_lambda)) +
    theme_pubr() +
    theme(plot.margin = margin(0,0,0,0, "cm"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          axis.line.x=element_blank())
  
  p2 = ggplot(gwa.plot,
              aes(x = bin,
                  y = n)) +
    geom_col(fill = "white", color = "black") +
    geom_hline(yintercept = (gwa.plot %>% ungroup() %>% mutate(dist=1-lambda) %>% slice_min(order_by = dist, n= 1) %>% pull(n)), color = "red", linetype = "twodash") +
    labs(y = "# variants") +
    theme_pubr() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x=element_blank(),
          plot.margin = margin(0,0,0,0, "cm"))
  
  
  p = ggarrange(p1,p2,
                ncol=1,
                nrow=2,
                heights = c(1,2),
                align = "v") +
    theme(plot.margin = margin(0,0,0,0, "cm")) 
  
  ggsave(filename=paste(mlmaFile, ".inflation.png", sep = ""),
         plot=p,
         device="png",
         width=3*2,
         height=3*2,
         units = 'in',
         dpi = 300)
  
  ggsave(filename=paste(mlmaFile, ".inflation.pdf", sep = ""),
         plot=p,
         device="pdf",
         width=3*2,
         height=3*2,
         units = 'in',
         dpi = 300)
}