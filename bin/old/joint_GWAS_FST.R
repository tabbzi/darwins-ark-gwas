# Libraries
require(ggplot2)
require(ggman)
require(ggplot2)
require(ggpubr)
require(colorspace)
require(data.table)
require(tidyr)
require(dplyr)

# Arguments
args = commandArgs(trailingOnly = T)
# args[1] = working directory
# args[2] = .mlma
# args[3] = .fst

setwd(args[1])

# Load Data
gwas = as.data.table(read.delim(args[2], header = T, sep = "\t")) %>%
       setnames("Chr","CHR") %>%
       setnames("bp","POS")
gwas$log.p = -log10(gwas$p)
#print(head(gwas))
fst = as.data.table(read.delim(args[3], header = T, sep = "\t"))
#print(head(fst))
joint = merge(gwas, fst, by = c("CHR","SNP","POS"))
print(head(joint))
setkey(joint, CHR, POS)

joint$POS = as.numeric(joint$POS)

joint$map.x = cumsum(joint$POS)

axis.set = joint %>% 
  group_by(CHR) %>% 
  summarize(center = (max(map.x) + min(map.x)) / 2)
ylim <- abs(floor(log10(min(joint$p)))) + 2 
sig <- 5e-8

nCHR <- length(unique(joint$CHR))

joint.plot.gwas <- ggplot(joint[log.p > 2], aes(x = map.x, 
                                 color = as.factor(CHR))) +
  geom_point(aes(y = -log10(p)), alpha = 0.75) +
  geom_hline(yintercept = -log10(sig), color = "grey40", linetype = "dashed") + 
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(2, ylim)) +
  scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "-log10(p)") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )

joint.plot.fst <- ggplot(joint[log.p > 2], aes(x = map.x, 
                                 color = as.factor(CHR))) +
  geom_point(aes(y = FST), alpha = 0.75) +
  scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_color_manual(values = rep(c("#f1935c", "#ba6b57"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  labs(x = NULL, 
       y = "Fst") + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank()
  )

joint.plot = ggarrange(joint.plot.gwas,joint.plot.fst,nrow=2,heights=c(2,1), align = "v")

ggsave(
  plot = joint.plot,
  filename = paste(args[3],".png",sep=""),
  device = "png",
  dpi = 300,
  unit = "in",
  width = 15,
  height = 6,
  bg = "white"
 )

