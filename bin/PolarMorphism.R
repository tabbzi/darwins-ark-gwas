# Required libraries:
require(tidyverse)
require(ggtext)
require(ggpubr)
require(data.table)
require(stringr)
require(dplyr)
require(svglite)
require(colorspace)
require(argparse)
require(PolarMorphism)
require(whitening)
require(qvalue)

# Define the arguments that the script accepts
parser <- ArgumentParser()
parser$add_argument("--in",
                    help = "Input .tsv containing phenotype codes and paths to summary statistics",
                    default = "/seq/vgb/dd/gwas/polar/2022-11-20/fa.all.list.tsv")
parser$add_argument("--out", 
                    help = "Output path and prefix",
                    default = "/seq/vgb/dd/gwas/polar/2022-11-20/fa.all")

# Parse the arguments passed to the script:
args <- parser$parse_args()
print(args)

# Read inputs table:
inputs = read_tsv(args$in)

# Load all summary statistics:
for (id in inputs$index) {
  tbl <- read_tsv((inputs %>% 
                     filter(index==id) %>% 
                     pull(path)), 
                  col_names = T) %>%
    rename(snpid = SNP,
           a1 = A1,
           a2 = A2,
           freq = Freq,
           beta = b,
           pval = p) %>%
    mutate(z = beta/se)
  assign((inputs %>%
            filter(index==id) %>%
            pull(pheno)),
         tbl)
}
rm(tbl)

gwas_polar <- ConvertToPolar(dfnames = inputs$pheno,
                             snpid = "snpid",
                             whiten = T,
                             LDcorrect = F)

p = gwas_polar %>%
  filter(r > 4) %>%
  ggplot(aes(x = z.whitened.1,
             y = z.whitened.2)) +
  theme(aspect.ratio = 1) +
  geom_point() +
  theme_pubclean()
ggsave(plot = p,
       filename = "fa.test.whiten.png")

gwas_polar$r.pval <- PvalueForR(r = gwas_polar$r, p = 2)
gwas_polar$r.qval <- qvalue(p = gwas_polar$r.pval)$qvalues

gwas_polar_filt <- gwas_polar[gwas_polar$r.qval < 0.05,]
gwas_polar_filt$theta.pval <- PvalueForAngle(angle.trans = gwas_polar_filt$angle,
                                             r = gwas_polar_filt$r,
                                             kappa.file = "/seq/vgb/dd/gwas/polar/PolarMorphism/R/kappas.4foldtransform.Rda")
gwas_polar_filt$theta.qval <- qvalue(p = gwas_polar_filt$theta.pval)$qvalues

p = gwas_polar_filt %>%
  ggplot(aes(x = abs(angle), y = r, color = theta.qval < 0.05)) +
  geom_point()

ggsave(plot = p,
       filename = "fa.test.angle.png")

write_csv(x = gwas_polar_filt, file = "fa.test.csv")