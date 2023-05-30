args <- commandArgs(trailingOnly = T)
file <- args[1]
library(tidyverse)
Data <- read_delim(file, delim="\t",col_names=c("Measurement","Levene","BF","Bartlett","FK","DGLM", "DRM", "CLS","SVLM","Zscore","QUAIL",
                                                   "TotalSamples","MAF","mean_effect","var_effect","env_effect","error_effect","ErrorDistribution",
                                                   "freq1","env1","freq2","env2","Gadd_mean","Gadd_sd","Gvar_mean","Gvar_sd","Env_mean","Env_sd","Error_mean","Error_sd")) %>% filter(Measurement=="unadjusted_p") %>% 
  select(BF, DRM, SVLM, TotalSamples, MAF, mean_effect, var_effect, ErrorDistribution)

CalculatePower <- function(samples, maf, g_mean, g_var, error, Data){
  Data1 <- Data %>% filter(MAF==maf & TotalSamples==samples & mean_effect==g_mean & var_effect==g_var & ErrorDistribution==error) %>%
    select(-TotalSamples, -MAF, -mean_effect, -var_effect, -ErrorDistribution)
  power <- apply(Data1, 2, function(x) (1-(length(which(x>0.05))/nrow(Data1)))) %>% t() %>% as.data.frame() %>%
    mutate(Sample=samples, MAF=maf,errorDistribution=error,Gmean=g_mean, Gvar=g_var)
  return(power)
}

parameters <- Data %>% group_by(TotalSamples, MAF, mean_effect, var_effect, ErrorDistribution) %>% 
  tally() %>% filter(var_effect>0, mean_effect==0)

powerRes <- list()
for (i in c(1:nrow(parameters))){
  powerRes[[i]] <- CalculatePower(parameters$TotalSamples[i], parameters$MAF[i], parameters$mean_effect[i], parameters$var_effect[i], parameters$ErrorDistribution[i], Data)
}

power <- bind_rows(powerRes)
write.table(power, gsub("SimResult","PowerCal",file), col.names = T, row.names=F, sep="\t", quote=F)

