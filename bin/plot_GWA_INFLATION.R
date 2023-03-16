gwa = read_tsv("/seq/vgb/dd/gwas/assoc/2022-11-20/DarwinsArk_gp-0.70_biallelic-snps_maf-0.001_geno-0.05_hwe-1e-20-midp-keep-fewhet_N-3465_phe-bq.146_dcov-datatype_qcov-age.hgt.bq.146_rel-cutoff-0.75.loco.mlma")
gwa.plot = gwa %>% 
  select(Freq,p) %>%
  mutate(bin = cut(Freq, 
                   breaks = 10^(seq(from = -3.0, to = -0.3, by = 0.1)))) %>% 
  group_by(bin) %>% 
  mutate(lambda = qchisq(mean(p), df = 1) / qchisq(0.5, df = 1),
         maxfreq = max(Freq)) %>%
  ungroup() %>%
  group_by(bin,lambda,maxfreq) %>%
  summarize(n = n())


# Load the tidyverse and ggplot2 packages
library(tidyverse)
library(ggplot2)

# Read in summary statistics file as a data frame
df <- read_tsv("summary_stats.tsv")

# Bin by -log10(maf) using the cut() function
df <- df %>% mutate(bin = cut(-log10(maf), breaks = c(0, 0.5, 1, 1.5, 2)))

# Calculate genomic inflation factor using the qchisq() function
df <- df %>% 
  group_by(bin) %>% 
  mutate(lambda = qchisq(mean(p), df = 1) / qchisq(0.5, df = 1))

# Calculate number of variants by bin using the n() function
df <- df %>% group_by(bin) %>% summarize(n = n())

# Create a ggplot with separate panels for variant count (n) and genomic inflation (lambda)
p = gwa.plot %>%
  mutate(bin = as.factor(bin)) %>%
  ggplot(aes(x = bin, group = 1)) +
  geom_col(aes(y = n), 
           color = "blue") +
  geom_line(aes(y = lambda), 
            color = "red") +
  scale_y_continuous(name = "Variant count (n)",
                     sec.axis = sec_axis(name = "Genomic inflation (lambda)")) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ylim.prim <- c(0, max(gwa.plot$n))   # in this example, precipitation
ylim.sec <- c(min(gwa.plot$lambda, na.rm = T), max(gwa.plot$lambda, na.rm = T))    # in this example, temperature

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

p = ggplot(gwa.plot, aes(x = -log10(maxfreq))) +
  geom_col(aes(y = n), position = "dodge") +
  geom_line(aes(y = lambda, group = 1), color = "red", size = 1) +
  scale_y_continuous(name = "n",
                sec.axis = sec_axis(~ (. - a)/b, name = "lambda"))

library(cowplot)



p1 = ggplot(gwa.plot,
            aes(x = -log10(maxfreq),
                y = lambda)) +
  geom_hline(yintercept = 1, color = "red", linetype = "twodash") +
  geom_line() +
  labs(y = "genomic inflation") +
  theme_pubr()

p2 = ggplot(gwa.plot,
            aes(x = bin,
                y = n)) +
  geom_hline(yintercept = (gwa.plot %>% ungroup() %>% mutate(dist=1-lambda) %>% slice_min(order_by = dist, n= 1) %>% pull(n)), color = "red", linetype = "twodash") +
  geom_col() +
  labs(y = "# variants") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p = ggarrange(p1,p2,
              ncol=1,
              nrow=2,
              heights = c(1,2),
              align = "h")

ggsave(filename="/seq/vgb/dd/gwas/inflat.plot.png",plot=p,device="png",width=3*2,height=5*2,units = 'in',dpi = 300)