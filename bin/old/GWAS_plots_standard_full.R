# Libraries
require(ggplot2)
require(ggman)
require(ggplot2)
require(ggpubr)
require(colorspace)
require(data.table)

# Arguments
args = commandArgs(trailingOnly = T)
# args[1] = working directory (/seq/vgb/dd/gwas/)
# args[2] = subdirectory (2019-08-12)
# args[3] = geno prefix
# args[4] = pheno prefix
# args[5] = covar prefix
# args[6] = qcovar prefix
assoc.file = paste(args[1], "assoc/", args[2],"/", args[3], "_", args[4], ".mlma", sep = "")
clump.file = paste(args[1], "clump/", args[2],"/", args[3], "_", args[4], "_strng.clumped.tab.rm", sep = "")
genes.file = paste(args[1], "clump/", args[2],"/", args[3], "_", args[4], "_strng.clumped.ranges", sep = "")
covar.file = paste(args[1], "covar/", args[5], ".covar", sep = "")
if (args[6] != "NA"){
qcovar.file = paste(args[1], "covar/", args[6], ".qcovar", sep = "")
}
pheno.file = paste(args[1], "pheno/", args[4], ".pheno", sep = "")
string.file = paste(args[1], "pheno/", args[4], "_string.txt", sep = "")

# Load Data
# Association (GCTA MLMA)
assoc = fread(paste(assoc.file))
assoc = na.omit(assoc)

# Clump (PLINK --clump)
head = c("SNP","SP2","RANGES","GROUPS")
clumped.ranges = fread(genes.file)
clumped = fread(clump.file)
clump = data.table(clumped$SNP, clumped$SP2, clumped.ranges$RANGES, rownames(clumped))
colnames(clump) = head
clump$RANGES = lapply(clump$RANGES, function(x) gsub("\\[|\\]", "", x))
for (i in 1:nrow(clump)){
  if (nchar(clump[i,3]) >= 14){
    clump[i,3] = "multigenic"
  }
}


# Phenotypes
covar = na.omit(fread(covar.file))
pheno = na.omit(fread(pheno.file))
common = intersect(as.character(covar$V1), as.character(pheno$V1))

if (args[6] != "NA"){
  qcovar = na.omit(fread(qcovar.file))
  common = intersect(intersect(as.character(qcovar$V1),as.character(covar$V1)), as.character(pheno$V1))
}


pheno = pheno[as.character(pheno$V1) %in% common]

colnames(pheno) = c("FID","IID","phenotype")

pheno.string = readLines(string.file)

# Functions
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

# Q-Q PLOT
gg.qq = gg_qqplot(assoc$p) + 
  coord_fixed(ratio = 1) +
  theme(axis.text=element_text(size=10,face="bold"), aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,2,1,1), "lines"))

# DISTRIBUTION PLOT
gg.dist.invert = ggplot(pheno, aes(x = phenotype, fill = as.character(phenotype))) +
  geom_bar(color = "black", size = 0.3) +
  geom_text(stat='count', aes(label=..count..), position=position_dodge(width=0), hjust=-0.5, vjust=0) +
  labs(x = "Phenotype", y = "Number of Dogs") +
  guides(fill = guide_legend(title = "")) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10,face="bold"),
        plot.margin = unit(c(1,2,1,1), "lines")) +
  scale_fill_brewer(palette = "GnBu") +
  coord_flip(clip = "off") +
  guides(fill = F)

if (grepl(x = args[4], pattern = "fac") == T){
gg.dist.invert = ggplot(pheno, aes(x = phenotype)) +
  geom_density(color = "black", size = 0.3) +
  labs(x = "Phenotype", y = "Number of Dogs") +
  guides(fill = guide_legend(title = "")) +
  theme(rect = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y = element_blank(),
        axis.title.x=element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_text(size=10,face="bold"),
        plot.margin = unit(c(1,2,1,1), "lines")) +
  scale_fill_brewer(palette = "GnBu") +
  coord_flip(clip = "off") +
  guides(fill = F)
}

# Now without labels!
gg.clumps = ggmanClumps(clump, index.snp.column = "SNP", clumps.column = "SP2", group.column = "GROUPS") 
gg.manhattan = ggman(assoc, clumps = gg.clumps, chrom = "Chr", snp = "SNP", bp = "bp", pvalue = "p", relative.positions = T, pointSize = 0.5, sigLine = NA)

maxlim = 8

if (-log10(min(assoc$p)) >= 8){
  maxlim = -log10(min(assoc$p))
}

gg.manhattan = gg.manhattan[[1]] +
  geom_hline(yintercept = 7.22, colour="#990000", linetype="dashed") +
  labs(x = "Chromosome: Position", y = "-log10(P)") +
  guides(fill = F) +
  ggtitle("") +
  scale_y_continuous(limits = c(2,maxlim), breaks = seq(2,maxlim,1)) +
  theme(axis.text.y=element_text(size=10,face="bold"),
        axis.text.x=element_text(size=8,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(1,1,1,1), "lines"))

# Compiled Figure
gg.GWAS.arrange = annotate_figure(ggarrange(ncol = 2, widths = c(4,1),
                                            gg.manhattan,
                                            ggarrange(nrow = 2, heights = c(2,1),
                                                      gg.dist.invert, gg.qq)),
                                  fig.lab = pheno.string,
                                  fig.lab.pos	= "top.left",
                            fig.lab.face = "bold")

ggsave(
  plot = gg.GWAS.arrange,
  filename = paste(args[1], "plots/", args[2], "/", args[3], "_", args[4], "_unlabeled.png", sep=""),
  device = "png",
  dpi = 300,
  unit = "in",
  width = 12,
  height = 6,
  bg = "white"
 )
