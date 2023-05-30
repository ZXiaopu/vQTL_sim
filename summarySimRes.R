library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
file <- as.character(args[1])

SummaryP <- function(g_mean, g_var, error, Data){
  Res <- Data %>% filter(mean_effect==g_mean & var_effect==g_var & ErrorDistribution==error) %>%
    select(Levene, BF, Bartlett, FK, DGLM, DRM, CLS, SVLM, Zscore, QUAIL, Measurement)
  Res_unadjust <- Res %>% filter(Measurement=="unadjusted_p") %>% select(-Measurement)
  Res_RINT <- Res %>% filter(Measurement=="RINT_p") %>% select(-Measurement)
  
  OutP <- rbind(apply(Res_unadjust, 2, function(x) length(which(x<=0.05))/nrow(Res_unadjust)) %>% t() %>% as.data.frame() %>%
                  mutate(Measurement="unadjusted_p", g_mean=g_mean, g_var=g_var, errorDistribution=error),
                apply(Res_RINT, 2, function(x) length(which(x<=0.05))/nrow(Res_RINT)) %>% t() %>% as.data.frame() %>%
                  mutate(Measurement="RINT_p", g_mean=g_mean, g_var=g_var, errorDistribution=error)) %>%
    mutate(Gadd=mean(Data$Gadd_mean), Gvar=mean(Data$Gvar_mean), Env=mean(Data$Env_mean), Error=mean(Data$Error_mean))
  return(OutP)
}

SummaryEffect <- function(g_mean, g_var, error, Data){
  Res <- Data %>% filter(mean_effect==g_mean & var_effect==g_var & ErrorDistribution==error) %>%
    select(DGLM, DRM, SVLM, Zscore, QUAIL, Measurement, var_effect)

#  ResEffect_unadjust <- Res %>% filter(Measurement == "unadjusted_effect") %>% select(-Measurement)
#  ResEffect_unadjust1 <- apply(ResEffect_unadjust, 2, function(x) abs(x - (ResEffect_unadjust$var_effect)))
#  ResEffect_unadjust2 <- apply(ResEffect_unadjust1, 2, function(x) mean(x)) %>% t() %>% as.data.frame() %>%
#    mutate(Measurement="unadjust_effect_mean", g_mean=g_mean, g_var=g_var, errorDistribution=error) %>%
#    mutate(Gadd=mean(Data$Gadd_mean), Gvar=mean(Data$Gvar_mean), Env=mean(Data$Env_mean), Error=mean(Data$Error_mean))
  ResEffect_unadjust <- Res %>% filter(Measurement == "unadjusted_effect") %>% select(-Measurement)
  ResEffect_unadjust1 <- data.frame(t(apply(ResEffect_unadjust, 2, mean))) %>% mutate(type="mean")
  ResEffect_unadjust2 <- data.frame(t(apply(ResEffect_unadjust, 2, min))) %>% mutate(type="min")
  ResEffect_unadjust3 <- data.frame(t(apply(ResEffect_unadjust, 2, max))) %>% mutate(type="max")
  OutEffect <- rbind(ResEffect_unadjust1,ResEffect_unadjust2,ResEffect_unadjust3) %>%
    mutate(Measurement="unadjust", g_mean=g_mean, g_var=g_var, errorDistribution=error)
  
#  ResEffect_RINT <- Res %>% filter(Measurement=="RINT_effect") %>% select(-Measurement)  
#  ResEffect_RINT1 <- apply(ResEffect_RINT, 2, function(x) abs(x- (ResEffect_RINT$var_effect)))
#  ResEffect_RINT2 <- apply(ResEffect_RINT1, 2, function(x) mean(x)) %>% t() %>% as.data.frame() %>%
#    mutate(Measurement="RINT_effect_mean", g_mean=g_mean, g_var=g_var, errorDistribution=error) %>%
#    mutate(Gadd=mean(Data$Gadd_mean), Gvar=mean(Data$Gvar_mean), Env=mean(Data$Env_mean), Error=mean(Data$Error_mean))
#  OutEffect <- rbind(ResEffect_unadjust2, ResEffect_RINT2)
  return(OutEffect)
}

runSummary <- function(file){
  Data <- read_delim(file, delim="\t", col_names=c("Measurement","Levene","BF","Bartlett","FK","DGLM", "DRM", "CLS","SVLM","Zscore","QUAIL",
                                                   "TotalSamples","MAF","mean_effect","var_effect","env_effect","error_effect","ErrorDistribution",
                                                   "freq1","env1","freq2","env2","Gadd_mean","Gadd_sd","Gvar_mean","Gvar_sd","Env_mean","Env_sd","Error_mean","Error_sd"))
  parameters <- Data %>% group_by(mean_effect, var_effect, ErrorDistribution) %>% tally()
  
  summary.P <- list()
#  summary.Effect <- list()
  for (i in c(1:nrow(parameters))){
    summary.P[[i]] <- SummaryP(parameters$mean_effect[i], parameters$var_effect[i], parameters$ErrorDistribution[i], Data)
#    summary.Effect[[i]] <- SummaryEffect(parameters$mean_effect[i], parameters$var_effect[i], parameters$ErrorDistribution[i], Data)
  }
  
  summaryres.P <- bind_rows(summary.P)
#  summaryres.Effect <- bind_rows(summary.Effect)
  
  write.table(summaryres.P, gsub("SimResult","SummaryP",file), col.names = T, row.names=F, sep="\t", quote=F)
#  write.table(summaryres.Effect, gsub("SimResult","SummaryEffect.range",file), col.names = T, row.names=F, sep="\t", quote=F)
}

runSummary(file)
