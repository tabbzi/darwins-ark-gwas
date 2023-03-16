# PREAMBLE ----

## LIBRARIES ----
require(tidyverse)
require(ggtext)
require(ggpubr)
require(data.table)
require(stringr)
require(dplyr)
require(svglite)
require(colorspace)
require(argparse)
require(corrplot)

## ARGUMENTS ----
parser <- ArgumentParser()
parser$add_argument("--gen",
                    help = "Input .tsv with genetic correlations",
                    default = "~/Dropbox (UMass Medical School)/Papers/cross-species/dat/gen/all.fa-vs-fa.no-lds.hsq.tsv")
parser$add_argument("--cor",
                    help = "Input .tsv with survey correlations",
                    default = "~/Dropbox (UMass Medical School)/Projects/DarwinsArk/git/darwins-ark-data/dat/DarwinsArk_20221120_factor-analysis_discovery-no-na_qn-182_dn-7883_fn-25.factor-corr.tsv")
parser$add_argument("--out", 
                    help = "Output for plot",
                    default = "")

args <- parser$parse_args()
print(args)

# LOAD DATA ----
cx = read_tsv(args$cor)
gx = read_tsv(args$gen)

## Set phenotypes ----
phe = sort(unique(c(cx$phe1,cx$phe2)))

## Complete pairs of phenotypes ----
gx = gx %>%
  mutate(phe1 = factor(x = phe1,
                          levels = phe),
            phe2 = factor(x = phe2,
                          levels = phe)) %>%
  complete(phe1,phe2) %>%
  arrange(phe1,phe2) %>%
  pivot_wider(id_cols = "phe1",
              names_from = "phe2",
              values_from = "corr") %>%
  select(-phe1) %>%
  as.matrix()

colnames(gx) = phe
rownames(gx) = phe

cx = cx %>%
  mutate(phe1 = factor(x = phe1,
                       levels = phe),
         phe2 = factor(x = phe2,
                       levels = phe)) %>%
  complete(phe1,phe2) %>%
  arrange(phe1,phe2) %>%
  pivot_wider(id_cols = "phe1",
              names_from = "phe2",
              values_from = "corr") %>%
  select(-phe1) %>%
  as.matrix()

colnames(cx) = phe
rownames(cx) = phe

## Combine matrices ----
dx = gx
dx[lower.tri(dx, diag = FALSE)] = cx[lower.tri(cx, diag = FALSE)]
diag(dx) = 1

## Plot ----
pdf(out)
corrplot.mixed(dx,
               upper = "circle",
               lower = "circle",
               is.corr = F)
dev.off()
