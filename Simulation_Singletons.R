library(dplyr)
library(tidyr)
library(janitor)
library(car) # Levene's test, Brown-Forysthe test
library(stats) # Bartlett's test, FK test
library(dglm) # dglm
library(quantreg)

args <- commandArgs(trailingOnly = TRUE)
N <- as.numeric(args[1])
MAF <- as.numeric(args[2]) 
freq_env <- as.numeric(args[3])
env1 <- as.numeric(args[4])
g_mean <- as.numeric(args[5])
g_var <- as.numeric(args[6])
a_env <- as.numeric(args[7])
a_error <- as.numeric(args[8])
Dist <- args[9]
Exposure <- args[10] ## binary/uniform/normal
count <- args[11]

GenerateVariable <- function(N, MAF, freq_env, env1){
  ## minor homogeneity group contains 5 samples
  MHF <- 10
  Heter <- round(2*N*MAF - 20)
  MajorHomo <- N - Heter - MHF
  df <- data.frame(genotype=c(rep(0,MHF), rep(1,Heter), rep(2,MajorHomo)))
  
  if(Exposure=="binary"){
    freq_env1 <- freq_env
    freq_env2 <- 1-freq_env1
    env1 <- env1
    env2 <- -(env1*freq_env1/freq_env2)
    df$Env[df$genotype==0] <- c(rep(env1,MHF*freq_env1),rep(env2,MHF*freq_env2))
    df$Env[df$genotype==1] <- c(rep(env1,Heter*freq_env1),rep(env2,Heter*freq_env2))
    df$Env[df$genotype==2] <- c(rep(env1,MajorHomo*freq_env1),rep(env2,MajorHomo*freq_env2))
  }
  
  if(Exposure=="uniform"){
    df$Env <- round(runif(N, 20, 70))
  }
  
  if(Exposure=="normal"){
    df$Env <- rnorm(N, 25, 3)
  }
  
  return(df)
}

GenerateTrait <- function(df, g_mean, g_var, a_env, a_error, Dist){
  if(Dist=="Normal"){df$error <- rnorm(N, mean=0, sd=1)}
  if(Dist=="Chi-squared"){df$error <- rchisq(N, df = 1)}
  if(Dist=="Gamma"){df$error <- rgamma(N, shape = 2, scale = 0.5)}
  
  df$Gmean <- df$genotype
  df$Gvar <- df$genotype*df$Env

  df$Gmean <- ((df$Gmean - min(df$Gmean))/(max(df$Gmean)-min(df$Gmean)))*g_mean
  df$Gvar <- ((df$Gvar - min(df$Gvar)) / (max(df$Gvar) - min(df$Gvar))) *g_var
  df$Env <- ((df$Env - min(df$Env)) / (max(df$Env) - min(df$Env)))*a_env
  df$error <- ((df$error - min(df$error)) / (max(df$error) - min(df$error)))*a_error
  df$Trait <- df$Gmean + df$Gvar + df$Env + df$error
  return(df)
}

GenerateData <- function(N, MAF, freq_env, env1, g_mean, g_var, a_env, a_error, Dist){
  Variables <- GenerateVariable(N, MAF, freq_env, env1)
  SimulatedData <- GenerateTrait(Variables, g_mean, g_var, a_env, a_error, Dist)
  return(SimulatedData)
}

QUAIL <- function(x, y, quantiles = 100){
  x <- (x - mean(x))/sd(x)
  y <- lm(y~x)$residuals
  a_i_tau_list <- rep()
  k <- quantiles - 1 ### number of quantiles to fit the model
  d <- rnorm(length(x))
  for (i in 1:k){
    tau_curr <- i/(k+1)
    Qreg <- rq(y ~ d, tau = tau_curr, method = "fn")
    SE_d_cur <- summary(Qreg , se = "ker")$coefficients[2,2]
    balance <-   sqrt(-tau_curr^2 + tau_curr)/SE_d_cur    
    a_i_tau <- ifelse(residuals(Qreg) < 0 ,1,0)
    a_i_tau_list[[i]] <- (-a_i_tau + tau_curr)/balance
  }
  
  int_rank_score <- data.frame(score = rep(0, length(y)))
  for (i in 1:(k-1)){
    if ( i > quantiles/2){
      term2 <- 1
    }else {
      term2 <- -1
    }
    int_rank_score <- int_rank_score + term2*a_i_tau_list[[i]]
  }
  
  int_rank_score <- int_rank_score/(quantiles/2)
  Xstar <-  lm(x ~ 1)$residuals
  null_phi_1 <- int_rank_score$score*sqrt(sum(Xstar^2))
  res <- summary(lm(null_phi_1 ~ Xstar))
  return(res)
}


IndQTLmapping <- function(indData, type){
  indData$genotype <- as.factor(indData$genotype)
  Levene.P <- leveneTest(resi ~ genotype, data = indData, center = mean)$`Pr(>F)`[1]
  BF.P <- leveneTest(resi ~ genotype, data = indData, center = median)$`Pr(>F)`[1]
  Barlett.P <- bartlett.test(resi ~ genotype, data = indData)$p.value
  FK.P <- fligner.test(resi ~ genotype, data = indData)$p.value
  
  res <- summary(dglm(resi ~ as.numeric(genotype), ~as.numeric(genotype), data = indData, family = gaussian))
  DGLM.varP <- res$dispersion.summary$coefficients[2,4]
  DGLM.varEffect <- res$dispersion.summary$coefficients[2,1]
  
  Y.i <- tapply(indData$resi, indData$genotype, median)
  indData$DRM <- abs(indData$resi - Y.i[as.factor(indData$genotype)])
  DRM.P <- summary(lm(DRM ~ as.numeric(genotype), data = indData))$coef[2,4]
  DRM.effect <- summary(lm(DRM ~ as.numeric(genotype), data = indData))$coef[2,1]
  
  if(type=="raw"){
    indData$resid.CLS <- residuals(summary(lm(Trait ~ as.numeric(genotype) + Env, data=indData)))
    indData$resid.SVLM <- residuals(summary(lm(Trait ~ as.numeric(genotype) + Env, data = indData)))
    indData$InverseNorm <- qnorm((rank(indData$resi,na.last="keep") - 0.5)/sum(!is.na(indData$resi)))
    TSSR.P <- summary(lm(InverseNorm^2 ~ as.numeric(genotype), data=indData))$coef[2,4]
    TSSR.effect <- summary(lm(InverseNorm^2 ~ as.numeric(genotype), data=indData))$coef[2,1]
  }else{
    indData$resid.CLS <- qqnorm(residuals(summary(lm(Trait ~ as.numeric(genotype) + Env, data=indData))), plot.it = FALSE)$x
    indData$resid.SVLM <- qqnorm(residuals(summary(lm(Trait ~ as.numeric(genotype) + Env, data = indData))), plot.it = FALSE)$x
    TSSR.P <- NA
    TSSR.effect <- NA
  }
  CLS.P <- cor.test((indData$resid.CLS)^2, as.numeric(indData$genotype), method="spearman")$p.value
  SVLM.P <- summary(lm((resid.SVLM^2) ~ as.numeric(genotype), data = indData))$coef[2,4]
  SVLM.effect <- summary(lm((resid.SVLM^2) ~ as.numeric(genotype), data = indData))$coef[2,1]

  QUAILres <- QUAIL(x=as.numeric(indData$genotype), y=indData$resi)
  QUAIL.P <- QUAILres$coefficients[2,4]
  QUAIL.effect <- QUAILres$coefficients[2,1]
  
  out <- data.frame(t(data.frame(Pvalue=c(Levene.P, BF.P, Barlett.P, FK.P, DGLM.varP, DRM.P, CLS.P, SVLM.P, TSSR.P, QUAIL.P),
                                 effect=c(NA, NA, NA, NA, DGLM.varEffect, DRM.effect, NA, SVLM.effect, TSSR.effect, QUAIL.effect))))
  colnames(out) <- c("Levene", "BF", "Bartlett", "FK", "DGLM", "DRM", "CLS", "SVLM", "Zscore","QUAIL")
  return(out)
}  

RunSimulation <- function(N, MAF, freq_env, env1, g_mean, g_var, a_env, a_error, Dist){
  indData <- GenerateData(N, MAF, freq_env, env1, g_mean, g_var, a_env, a_error, Dist)
  
  ## unadjusted value
  # indData$resi <- residuals(lm(Trait ~ BMI + Age + Sex + Smoking, data = indData))
  indData$resi <- residuals(lm(Trait ~ Env, data = indData))
  output.raw <- IndQTLmapping(indData, type="raw")
  rownames(output.raw) <- c("unadjusted_p","unadjusted_effect")
  
  ## RINT value
  indData$resi <- qqnorm(residuals(lm(Trait ~ Env, data = indData)),plot.it = FALSE)$x
  output.rint <- IndQTLmapping(indData, type="RINT")
  rownames(output.rint) <- c("RINT_p","RINT_effect")
  
  output <- bind_rows(output.raw, output.rint)
  
  Indresult <- output %>% mutate(TotalSamples=N, MAF=MAF, mean_effect=g_mean, var_effect=g_var, env_effect=a_env, error_effect=a_error,
                                 ErrorDistribution = Dist, freq_env1=freq_env, env1=env1, freq_env2=1-freq_env, env2=-(env1*freq_env1/(1-freq_env))) %>%
    mutate(Gmean_mean = mean(indData$Gmean/indData$Trait,na.rm = T), Gmean_sd = sd(indData$Gmean/indData$Trait,na.rm = T),
           Gvar_mean = mean(indData$Gvar/indData$Trait, na.rm = T), Gvar_sd = sd(indData$Gvar/indData$Trait, na.rm = T),
           Env_mean = mean(indData$E/indData$Trait, na.rm=T), Env_sd = sd(indData$E/indData$Trait, na.rm=T),
           error_mean = mean(indData$error/indData$Trait, na.rm=T), error_sd = sd(indData$error/indData$Trait, na.rm=T))
  
  write.table(Indresult, paste0("/scratch/.../Xiaopu/SimRes.ideal.Sample",N,".MAF",MAF,".mean",g_mean,".var", g_var, ".", Dist, 
                                "fenv.", freq_env, ".E", Exposure, ".count",count), col.names=F ,row.names=T, sep="\t", quote=F, append=T)
}

for (i in c(1:1000)){
  RunSimulation(N, MAF, freq_env, env1, g_mean, g_var, a_env, a_error, Dist)
}
