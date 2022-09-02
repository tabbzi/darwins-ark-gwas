library(tidyr)
library(readr)
args = commandArgs(trailingOnly = T)

data = read.table(args[1], header = T)
chisq = qchisq(1-data$p,1)
qclamb50 = median(chisq)/qchisq(0.5,1)
qclamb10 = median(chisq)/qchisq(0.1,1)
qclamb1 = median(chisq)/qchisq(0.01,1)
qclamb01 = median(chisq)/qchisq(0.001,1)

library(data.table)
results = data.table(percentile = c("50th","10th","1st","1/10th"), qclamb = c(qclamb50,qclamb10,qclamb1,qclamb01))
results[, gwas := args[1]]

write.table(x = results, file = args[2], quote = F, sep = ",", row.names=F)
