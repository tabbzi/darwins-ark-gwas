library(data.table)
library(tidyr)
library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly = T)
# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

grm = ReadGRMBin(args[1])

# Functions used in the functions above
var_vg_func <- function(N, var_pi=2e-5){
    return(2/(N^2*var_pi))
}

var_rg_func <- function(N1, N2, hsq1, hsq2, rg, rp, overlap=TRUE, var_pi=2e-5){
    if(overlap==T) var_rg=((1-rg*rp)^2+(rg-rp)^2)/(hsq1*hsq2*N1^2*var_pi)
    if(overlap==F) var_rg=(rg^2*(N1^2*hsq1^2+N2^2*hsq2^2)+2*hsq1*hsq2*N1*N2)/(2*hsq1^2*hsq2^2*N1^2*N2^2*var_pi)
    return(var_rg)
}

power_func <- function(ncp, alpha){
    pchisq(qchisq(alpha, df=1,lower.tail=F), ncp=ncp, df=1, lower.tail=F)
}

h2O_func <- function(ncase, ncontrol, K, h2L, var_pi=2e-5){
    n=ncase+ncontrol
    v=ncase/(ncase+ncontrol)
    z=dnorm(qnorm(K))
    c=(K*(1-K))^2/(v*(1-v)*z^2)
    h2O=h2L/c
    var_h2O=var_vg_func(n, var_pi)
    var_h2L=c^2*var_h2O
    return(list(h2L=h2L, var_h2L=var_h2L, h2O=h2O, var_h2O=var_h2O))
}

# Function for a quantitative trait
# n = sample size 
# hsq = variance explained by all SNPs
# alpha = significance level
# var_pi = variance of the off-diagonal elements of the GRM
# The output are: se (standard error), ncp (non-centrality parameter) and power
calcUniQt <- function(
    n     =1000, 
    hsq   =0.5, 
    alpha =0.05,
  var_pi=2e-5
){
    l <- list()
    var_vg <- var_vg_func(n, var_pi)
    l$se <- sqrt(var_vg)
    l$ncp <- hsq^2/var_vg;
    l$power <- power_func(l$ncp, alpha)
    return(l)
}

args[2]
if(is.na(args[2])){
  varpi = as.numeric(var(grm$off))
  print(varpi)
} else {
varpi = as.numeric(args[2])
}

allresults = list()
i = 1
for (npop in c(1000,1500,1750,2000,2500,3000,4000)){
for (h in as.numeric(c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))){

results = calcUniQt(n=npop, hsq = h, alpha = 0.05, var_pi = varpi) %>% as.data.table()
results$N = npop
results$h2 = h
allresults[[i]] = results
i = i+1
}

}

allresults = rbindlist(allresults)
print(allresults)

p = ggplot(allresults, aes(color = as.factor(N), y = power, x = h2)) +
      geom_line() +
      geom_pointrange(aes(ymin=power-se,ymax=power+se)) +
      labs(x = "True Heritability", y = "Power", color = "Sample Size") +
      scale_x_continuous(labels=scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels=scales::percent_format(accuracy = 1), limits = c(0,1)) +
      theme_pubr()

#ggsave(plot = p, paste(args[1],"_power.pdf",sep=""), device = "pdf", units = "in", height = 5, width = 8)
ggsave(plot = p, paste(args[1],"_power_varpi-",varpi,".pdf",sep=""), device = "pdf", units = "in", height = 4, width = 12)
