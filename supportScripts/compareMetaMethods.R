 setwd("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/")




library(readr)
library(lme4)
library(RColorBrewer)
library(mvmeta)
library(progress)
library(parallel)
library(foreach)
library(doParallel)
library(vioplot)
library(moments)

readDoubleMatrix <- function(path){
  
  file <- ""
  if(endsWith(path, ".gz")){
    file <- gzfile(path)
  } else {
    file <- file(path)
  }
  
  firstColName <- scan(file, what = "character", n = 1, sep = '\t', quiet = T )
  if(firstColName == ""){
    # if there is no column name for the first column ...1is used by read_delim of the readr lib
    firstColName <- "...1"
  }
  
  l <- list()
  l[[firstColName]] <- col_character()
  colTypes <- structure(list(cols = l, default = col_double()), class = "col_spec")
  
  table_tmp <- read_delim(file, delim = "\t", quote = "", col_types = colTypes)
  table <- as.matrix(table_tmp[,-1])
  rownames(table) <- table_tmp[,1][[1]]
  return(table)
}

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10)


cohorts <- c("LL", "LLS_660Q", "LLS_OmniExpr", "NTR_Affy6", "NTR_GoNL", "RS")
palette(adjustcolor(brewer.pal(n = 3, name = "Accent"), alpha.f = 0.6))

stx3 <- "ENSG00000166900"
nod2 <- "ENSG00000167207"
nod2QtlSnp <- "rs1981760"





if(FALSE){
  genoList <- list()
  variantsTmp <- c()
  for(cohort in cohorts){
    genoList[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsPcCor/",cohort,"/merged2.dosages.txt"))
    variantsTmp <- c(variantsTmp,rownames(genoList[[cohort]]))
  }
  variantsCount <- table(variantsTmp)
  sharedVariants <- names(variantsCount)[variantsCount == length(cohorts)]
  
  for(cohort in cohorts){
    genoList[[cohort]] <- genoList[[cohort]][sharedVariants,]
  }
  
  expList <- list()
  sharedExpressionGenesTmp <- c()
  for(cohort in cohorts){
    expList[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsPcCor/",cohort,"/correctedExpression.txt"))
    sharedExpressionGenesTmp <- c(sharedExpressionGenesTmp, rownames(expList[[cohort]]))
  }
  sharedExpressionGenesCount <- table(sharedExpressionGenesTmp)
  sharedExpressionGenes <- names(sharedExpressionGenesCount)[sharedExpressionGenesCount == length(cohorts)]
  for(cohort in cohorts){
    expList[[cohort]] <- expList[[cohort]][sharedExpressionGenes,]
  }
  
  covList <- list()
  sharedCovariatsTmp <- c()
  for(cohort in cohorts){
    covList[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsPcCor/",cohort,"/covariates.combined.txt"))
    sharedCovariatsTmp <- c(sharedCovariatsTmp, colnames(covList[[cohort]]))
  }
  sharedCovariatesCount <- table(sharedCovariatsTmp)
  sharedCovariats <- names(sharedCovariatesCount)[sharedCovariatesCount == length(cohorts)]
  for(cohort in cohorts){
    covList[[cohort]] <- covList[[cohort]][,sharedCovariats]
  }
  str(covList)
  
  eqtls <- read.delim("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/eQtlgenLeadVariants.txt.gz")
  str(eqtls)
  
  eqtls <- eqtls[eqtls$feature_id %in% sharedExpressionGenes,]
  eqtls <- eqtls[eqtls$snp_id %in% sharedVariants,]
  
  for(cohort in cohorts){
    sampleTmp <- rownames(covList[[cohort]])
    sampleTmp <- c(sampleTmp, colnames(expList[[cohort]]))
    sampleTmp <- c(sampleTmp, colnames(genoList[[cohort]]))
    sampleCount <- table(sampleTmp)
    cohortSamples <- names(sampleCount)[sampleCount == 3]
    
    covList[[cohort]] <- covList[[cohort]][cohortSamples,]
    expList[[cohort]] <- expList[[cohort]][,cohortSamples]
    genoList[[cohort]] <- genoList[[cohort]][,cohortSamples]
    
  }
  
  cohortSampleCount <- sapply(cohorts, function(cohort){
    return(ncol(genoList[[cohort]]))
  })
  
  save(genoList, sharedVariants, sharedExpressionGenes, expList, sharedCovariats, covList, eqtls, cohortSampleCount, file = "biosInteractionsPcCor/compareMetaMethodsData.RData")
  
} else {
  load("biosInteractionsPcCor/compareMetaMethodsData.RData")
}

cov <- stx3
gene <- nod2
snp <- nod2QtlSnp


genotypes <- do.call(c, sapply(genoList, function(x){return(x[snp,])}))
expression <- do.call(c, sapply(expList, function(x){return(x[gene,])}))
covariate <- do.call(c, sapply(covList, function(x){return(x[,cov])}))
cohortLables <- factor(rep(names(cohortSampleCount), times  = cohortSampleCount))


#### mega anova
summary(lm(expression ~ genotypes * covariate + cohortLables))


#### mega anova with cohort interaction

summary(lm(expression ~ genotypes * covariate + genotypes * cohortLables))


#### mega LRT

fullModel <- lm(expression ~ genotypes * covariate + cohortLables)
baseModel <- lm(expression ~ genotypes + covariate + cohortLables)

logLik_full <- logLik(fullModel)
logLik_base <- logLik(baseModel)

LRT_stat <- -2 * (logLik_base - logLik_full)

df <- attr(logLik_full, "df") - attr(logLik_base, "df")
(p_value <- pchisq(LRT_stat, df, lower.tail = FALSE))
qnorm(p_value/2)

#### mega LRT with cohort * genotype interaction

fullModel <- lm(expression ~ genotypes * covariate + genotypes * cohortLables)
baseModel <- lm(expression ~ genotypes + covariate + genotypes * cohortLables)

logLik_full <- logLik(fullModel)
logLik_base <- logLik(baseModel)

LRT_stat <- -2 * (logLik_base - logLik_full)

df <- attr(logLik_full, "df") - attr(logLik_base, "df")
(p_value <- pchisq(LRT_stat, df, lower.tail = FALSE))
qnorm(p_value/2)


#### Meta LRT
cohort <- "LLS_OmniExpr"
cohort <- "LL"
cohortZ <- sapply(cohorts, function(cohort){
 
  fullModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
  baseModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] + covList[[cohort]][,cov])
 
  summary(fullModel)
  
  logLik_full <- logLik(fullModel)
  logLik_base <- logLik(baseModel)
  
  LRT_stat <- -2 * (logLik_base - logLik_full)
  
  # Compute p-value using Chi-square Distribution
  df <- attr(logLik_full, "df") - attr(logLik_base, "df")
  p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
  
  print(paste(coef(fullModel)[4], qnorm(p_value/2)))
  
  if(coef(fullModel)[4] < 0){
    return(qnorm(p_value/2))
  } else {
    return(-qnorm(p_value/2))
  }
  
})
cohortZ

(z <- sum(cohortZ * cohortSampleCount)/sqrt(sum(cohortSampleCount *cohortSampleCount)))
2*pnorm(-abs(z))

#### Mixed effect mega analysis

testModel <- lmer(expression ~ genotypes + covariate + (1 + genotypes + covariate | cohortLables), REML = F) 
fullModel <- lmer(expression ~ genotypes * covariate + (1 + genotypes + covariate | cohortLables), REML = F)
summary(fullModel)
anova(fullModel,testModel )


logLik_full <- logLik(fullModel)
logLik_base <- logLik(testModel)

LRT_stat <- -2 * (logLik_base - logLik_full)

df <- attr(logLik_full, "df") - attr(logLik_base, "df")
(p_value <- pchisq(LRT_stat, df, lower.tail = FALSE))
qnorm(p_value/2)

#### Multivariate meta analysis

betas <- matrix(NA, nrow = length(cohorts), ncol = 4)
betaCovariances <- list()

for(c in 1:length(cohorts)){
  
  cohort = cohorts[c]
  cohortFit <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
  
  betaCovariances[[cohort]] <- vcov(cohortFit)
  betas[c,] <- coef(cohortFit)
  
  
}
str(betas)
str(betaCovariances)

summary(mvmeta(betas, betaCovariances, method = "fixed"))

#### Multivariate meta analysis with LRT

betasFull <- matrix(NA, nrow = length(cohorts), ncol = 4)
betasBase <- matrix(NA, nrow = length(cohorts), ncol = 3)
betaCovariancesFull <- list()
betaCovariancesBase <- list()

for(c in 1:length(cohorts)){
  
  cohort = cohorts[c]
  cohortFitBase <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] + covList[[cohort]][,cov])
  cohortFitFull <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
  
  betaCovariancesBase[[cohort]] <- vcov(cohortFitBase)
  betaCovariancesFull[[cohort]] <- vcov(cohortFitFull)
  betasBase[c,] <- coef(cohortFitBase)
  betasFull[c,] <- coef(cohortFitFull)
  
}

metaFull <- mvmeta(betasFull, betaCovariancesFull, method = "ml")



metaBase <- mvmeta(betasBase, betaCovariancesBase, method = "ml")

summary(metaFull)
2*pnorm(-abs(13.6))
summary(metaBase)


logLik_full <- logLik(metaFull)
logLik_base <- logLik(metaBase)

LRT_stat <- -2 * (logLik_base - logLik_full)


df <- attr(logLik_full, "df") - attr(logLik_base, "df")
(p_value <- pchisq(LRT_stat, df, lower.tail = FALSE))
qnorm(p_value/2)



#### large compare

totalEqtls <- nrow(eqtls)
totalEqtls <- 100 #for testing

#zscores <- matrix(data =NA, nrow = totalEqtls, ncol = 4)


i =1

#pb <- progress_bar$new(
#  format = "  eQTLs [:bar] :percent (:eta remaining)",
#  total = totalEqtls, clear = FALSE, width = 100
#)

cl <- makeCluster(19)
clusterExport(cl, "eqtls")
clusterExport(cl, "genoList")
clusterExport(cl, "expList")
clusterExport(cl, "covList")
clusterExport(cl, "cohortSampleCount")
clusterExport(cl, "cohorts")

clusterEvalQ(cl, {library(mvmeta)})
i <- 62772 
#for(i in 1:totalEqtls){
#zscoresList <- parLapply(cl, 1:totalEqtls, function(i) {
zscoresList <- lapply( 1:totalEqtls, function(i) {
  
  zscores <- vector(mode = "numeric", length = 5)
  
  tryCatch(
    {
  
  gene <- eqtls[i,1]
  snp <- eqtls[i,2]
  
  
  genotypes <- do.call("c", sapply(genoList, function(x){return(x[snp,])}))
  expression <- do.call("c", sapply(expList, function(x){return(x[gene,])}))
  covariate <- do.call("c", sapply(covList, function(x){return(x[,cov])}))
  cohortLables <- factor(rep(names(cohortSampleCount), times  = cohortSampleCount))
  
  
  #### Mega anova
  
  # zscores[1] <- summary(lm(expression ~ genotypes * covariate + genotypes * cohortLables))$coefficients["genotypes:covariate", "t value"]
  
  #### Mega LRT
  
  fullModel <- lm(expression ~ genotypes * covariate + genotypes * cohortLables)
  baseModel <- lm(expression ~ genotypes + covariate + genotypes * cohortLables)
  
  logLik_full <- logLik(fullModel)
  logLik_base <- logLik(baseModel)
  
  LRT_stat <- -2 * (logLik_base - logLik_full)
  
  df <- attr(logLik_full, "df") - attr(logLik_base, "df")
  p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
  
  
  if(coef(fullModel)["genotypes:covariate"] < 0){
    zscores[2] <- qnorm(p_value/2)
  } else {
    zscores[2] <- -qnorm(p_value/2)
  }
  
  
  #### Meta LRT
  
  cohortZ <- sapply(cohorts, function(cohort){
    
    fullModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
    baseModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] + covList[[cohort]][,cov])
    
    
    logLik_full <- logLik(fullModel)
    logLik_base <- logLik(baseModel)
    
    LRT_stat <- -2 * (logLik_base - logLik_full)
    
    # Compute p-value using Chi-square Distribution
    df <- attr(logLik_full, "df") - attr(logLik_base, "df")
    p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
    
    if(coef(fullModel)[4] < 0){
      return(qnorm(p_value/2))
    } else {
      return(-qnorm(p_value/2))
    }
    
  })
  
  
  zscores[3] <- sum(cohortZ * cohortSampleCount)/sqrt(sum(cohortSampleCount * cohortSampleCount))
  
  zscores[4] <- sum(cohortZ * sqrt(cohortSampleCount))/sqrt(sum(cohortSampleCount))
  
  # cohortZ <- sapply(cohorts, function(cohort){
  #   
  #   fullModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
  #   baseModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] + covList[[cohort]][,cov])
  #   
  #   
  #   logLik_full <- logLik(fullModel)
  #   logLik_base <- logLik(baseModel)
  #   
  #   LRT_stat <- -2 * (logLik_base - logLik_full)
  #   
  #   # Compute p-value using Chi-square Distribution
  #   df <- attr(logLik_full, "df") - attr(logLik_base, "df")
  #   p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
  #   return(qnorm(p_value/2))
  #   
  #   if(coef(fullModel)[4] < 0){
  #     return(qnorm(p_value/2))
  #   } else {
  #     return(-qnorm(p_value/2))
  #   }
  #   
  # })
  # 
  # #### Meta multivariate anova
  # 
  # betas <- matrix(NA, nrow = length(cohorts), ncol = 4)
  # betaCovariances <- list()
  # 
  # for(c in 1:length(cohorts)){
  #   
  #   cohort = cohorts[c]
  #   cohortFit <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
  #   
  #   betaCovariances[[cohort]] <- vcov(cohortFit)
  #   betas[c,] <- coef(cohortFit)
  #   
  #   
  # }
  # 
  # zscores[4] <- summary(mvmeta(betas, betaCovariances, method = "fixed"))$coefficients[4, "z"]
  # 
  # zscores[5] <- summary(mvmeta(betas, betaCovariances, method = "reml"))$coefficients[4, "z"]
  # #pb$tick() 
  },
  error = function(cond) {
    message("error in ", i)
    })
  return(zscores)
})
length(zscoresList)

zscores <- do.call(rbind, zscoresList)

par(pty="s")
colnames(zscores) <- c("mega_anova", "mega_LRT", "meta_LRT", "meta_LRT_2", "meta_multivariate_anova_mixed")
pairs(zscores[,c(2,3,4)], upper.panel = NULL,  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.3), pch = 16)

cor.test(zscores[,2], zscores[,3])
cor.test(zscores[,2], zscores[,4])

colnames(zscores) <- c("mega_anova", "mega_LRT", "meta_LRT", "meta_multivariate_anova", "meta_multivariate_anova_mixed")

#saveRDS(zscores, "biosInteractionsPcCor/metaComparisonResults.rds")

zscores <- readRDS("biosInteractionsPcCor/metaComparisonResults.rds")

par(pty="s")
pairs(zscores, upper.panel = NULL,  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.3), pch = 16)

layout(1)
plot(zscores[,2], zscores[,3],  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, xlab = "LRT mega analysis (corrected for genotype * cohort)", ylab = "LRT meta analysis")
abline(h=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(v=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(0,1, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(h=c(-5,5), col = adjustcolor("firebrick", alpha.f = 0.4), lwd =2)
abline(v=c(-5,5), col = adjustcolor("firebrick", alpha.f = 0.4), lwd =2)

plot(zscores[,2], zscores[,5],  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, xlab = "LRT mega analysis (corrected for genotype * cohort)", ylab = "Mixed effect multivariate meta anova")
abline(h=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(h=-5, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(h=5, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(v=c(0,5,-5), col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(0,1, col = adjustcolor("black", alpha.f = 0.4), lwd =2)


bonfZ <- qnorm(0.05/nrow(eqtls)/2)
sum(abs(zscores[,3]) >=  -bonfZ)



BiosinteractionZ <- readRDS("biosInteractionsPcCor/metaZ.rds")

str(BiosinteractionZ[,stx3])

pipelineZ <- BiosinteractionZ[paste0(eqtls$snp_id, "-", eqtls$feature_id),stx3]

par(pty="s")
pairs(cbind(zscores,pipelineZ), upper.panel = NULL,  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.3), pch = 16)
range(pipelineZ)

layout(matrix(1:2,nrow =1))

plot(pipelineZ, zscores[,2],  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, xlab = "Pipeline Z-scores", ylab = "Mega LRT in R")
abline(h=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(v=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(0,1, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(h=c(-5,5), col = adjustcolor("firebrick", alpha.f = 0.4), lwd =2)
abline(v=c(-5,5), col = adjustcolor("firebrick", alpha.f = 0.4), lwd =2)

plot(pipelineZ, zscores[,"meta_LRT"],  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, xlab = "Pipeline Z-scores", ylab = "Meta LRT in R")
abline(h=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(v=0, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(0,1, col = adjustcolor("black", alpha.f = 0.4), lwd =2)
abline(h=c(-5,5), col = adjustcolor("firebrick", alpha.f = 0.4), lwd =2)
abline(v=c(-5,5), col = adjustcolor("firebrick", alpha.f = 0.4), lwd =2)


cor.test(pipelineZ, zscores[,"meta_LRT"])
cor.test(pipelineZ, zscores[,"meta_LRT_2"])

which(zscores[,"meta_LRT"] >= 5 & pipelineZ < 4)

pipelineZ[zscores[,"meta_LRT"] >= 5 & pipelineZ < 4]

zscores[zscores[,"meta_LRT"] >= 5 & pipelineZ < 4,]
#rs754281319-ENSG00000101347




str(pipelineZ)
par(pty="s")
pairs(cbind(zscores[,c(2,3)],pipelineZ), upper.panel = NULL,  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.3), pch = 16)







####  rs754281319-ENSG00000101347



cov <- stx3
gene <- "ENSG00000101347"
snp <- "rs754281319"



genotypes <- do.call("c", sapply(genoList, function(x){return(x[snp,])}))
expression <- do.call("c", sapply(expList, function(x){return(x[gene,])}))
covariate <- do.call("c", sapply(covList, function(x){return(x[,cov])}))
cohortLables <- factor(rep(names(cohortSampleCount), times  = cohortSampleCount))

plot(genotypes * covariate, col = round(genotypes) +1)

summary(lm(expression ~ genotypes * covariate + genotypes * cohortLables))

plot(covariate)
plot(expression)
kurtosis(expression)
kurtosis(`covariate` * genotypes)

cohort <- "LLS_OmniExpr"
cohort <- "LL"
#c("LLS_OmniExpr", "LL", "RS", "LLS_660Q")
layout(matrix(1:12,nrow =2))
cohortZ <- sapply(cohorts, function(cohort){
  
  
  d <- data.frame(exp = expList[[cohort]][gene,], cov = covList[[cohort]][,cov], geno = genoList[[cohort]][snp,])
 
  plot(d$geno * d$exp)
  plot(d$exp)
  
  k <- kurtosis(d$exp)
  
  n <- length(d$exp)
  SE_kurtosis <- sqrt((24 * n * (n - 1)^2) / ((n - 3) * (n - 2) * (n + 3) * (n + 5)))
  
  k / SE_kurtosis
   
  fullModel <- lm(exp ~ geno * cov, data = d)
  baseModel <- lm(exp ~ geno + cov, data = d)
  
  summary(fullModel)
  
  logLik_full <- logLik(fullModel)
  logLik_base <- logLik(baseModel)
  
  LRT_stat <- -2 * (logLik_base - logLik_full)
  
  # Compute p-value using Chi-square Distribution
  df <- attr(logLik_full, "df") - attr(logLik_base, "df")
  p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
  
  print(paste(coef(fullModel)[4], qnorm(p_value/2)))
  
  
  plot(genoList[[cohort]][snp,],expList[[cohort]][gene,], col = round(genoList[[cohort]][snp,]) + 1, main = cohort, xlab = "Genotype dosage", ylab = "eQTL gen expression")
  abline(lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,]))
  
  plot(covList[[cohort]][,cov], expList[[cohort]][gene,], col = round(genoList[[cohort]][snp,]) + 1, pch =16, main = cohort, xlab = "Covariate expression", ylab = "eQTL gen expression")

  
  dosage = 2
  for(dosage in 0:2){
    x <- data.frame("geno" = dosage, "cov" = range(covList[[cohort]][,cov]), check.names = F)
    y <- predict.lm(fullModel, x)  
    lines(x$cov,y, col = dosage + 1, lwd = 3)
  }

  if(coef(fullModel)[4] < 0){
    return(qnorm(p_value/2))
  } else {
    return(-qnorm(p_value/2))
  }
  


  
})
cohortZ
(z <- sum(cohortZ * cohortSampleCount)/sqrt(sum(cohortSampleCount *cohortSampleCount)))
2*pnorm(-abs(z))

sum(abs(cohortZ) * cohortSampleCount)/sqrt(sum(cohortSampleCount *cohortSampleCount))


cor.test(expList[[cohort]][gene,], covList[[cohort]][,cov])


listPipelinePerCorhort <- list()
for(cohort in cohorts){
  listPipelinePerCorhort[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsPcCor/",cohort,"/query2.txt"))[,stx3]
}






##### check distributions

totalEqtls <- nrow(eqtls)
totalEqtls <- 100 #for testing

#zscores <- matrix(data =NA, nrow = totalEqtls, ncol = 4)


i =1

#pb <- progress_bar$new(
#  format = "  eQTLs [:bar] :percent (:eta remaining)",
#  total = totalEqtls, clear = FALSE, width = 100
#)

#for(i in 1:totalEqtls){
distTestList <- lapply( 1:totalEqtls, function(i) {
  
  res <- vector(mode = "numeric", length = 4)
  
  tryCatch(
    {
      
      gene <- eqtls[i,1]
      snp <- eqtls[i,2]
      
      
      genotypes <- do.call("c", sapply(genoList, function(x){return(x[snp,])}))
      expression <- do.call("c", sapply(expList, function(x){return(x[gene,])}))
      covariate <- do.call("c", sapply(covList, function(x){return(x[,cov])}))
      cohortLables <- factor(rep(names(cohortSampleCount), times  = cohortSampleCount))
      
      #### Mega LRT
      
      fullModel <- lm(expression ~ genotypes * covariate + genotypes * cohortLables)
      baseModel <- lm(expression ~ genotypes + covariate + genotypes * cohortLables)
      
      logLik_full <- logLik(fullModel)
      logLik_base <- logLik(baseModel)
      
      LRT_stat <- -2 * (logLik_base - logLik_full)
      
      df <- attr(logLik_full, "df") - attr(logLik_base, "df")
      p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
      
      
      if(coef(fullModel)["genotypes:covariate"] < 0){
        res[1] <- qnorm(p_value/2)
      } else {
        res[1] <- -qnorm(p_value/2)
      }
      res[2] <- kurtosis(expression)
      
      #### Meta LRT
      
      cohortZ <- sapply(cohorts, function(cohort){
        
        fullModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
        baseModel <- lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,] + covList[[cohort]][,cov])
        
        
        logLik_full <- logLik(fullModel)
        logLik_base <- logLik(baseModel)
        
        LRT_stat <- -2 * (logLik_base - logLik_full)
        
        # Compute p-value using Chi-square Distribution
        df <- attr(logLik_full, "df") - attr(logLik_base, "df")
        p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
        
        if(coef(fullModel)[4] < 0){
          return(qnorm(p_value/2))
        } else {
          return(-qnorm(p_value/2))
        }
        
      })
      
      res[3] <- sum(cohortZ * cohortSampleCount)/sqrt(sum(cohortSampleCount * cohortSampleCount))
      
      kw <- sapply(cohorts, function(cohort){
        kurtosis(expList[[cohort]][gene,] ) * cohortSampleCount[cohort]
      })
      
      
      
      res[4] <- sum(kw)  / sum(cohortSampleCount)
      
    },
    error = function(cond) {
      message("error in ", i)
    })
  return(res)
})
str(distTestList)

distTest <- do.call(rbind, distTestList)



plot(distTest[,1],distTest[,2])
plot(distTest[,1],distTest[,2], ylim = c(0,10))

layout(matrix(1:2, nrow = 1))

plot(distTest[,3],distTest[,4], col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, ylab = "Mean Kurtosis per cohort\n(3 is expected for normal distribution)", xlab = "Interaction z-score")
abline(v=c(-5,5), col = "firebrick")
abline(h=3, col = "forestgreen")


plot(distTest[,3],distTest[,4], col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, ylim = c(0, 10), ylab = "Mean Kurtosis per cohort\n(3 is expected for normal distribution)", xlab = "Interaction z-score", main = "Zoom in of y-axis")
abline(v=c(-5,5), col = "firebrick")
abline(h=3, col = "forestgreen")

sum()
layout(1)
hist(distTest[abs(distTest[,3])>= 5,4], breaks = 100, main = "Kurtosis of eqtl gene expression of significant interactions", xlab = "Kurtosis")
abline(v=3, col = "forestgreen")










zScoreIntList <- lapply( 1:totalEqtls, function(i) {
  
  res <- vector(mode = "numeric", length =2)
   
  tryCatch(
    {
      
      gene <- eqtls[i,1]
      snp <- eqtls[i,2]
      
      
      genotypes <- do.call("c", sapply(genoList, function(x){return(x[snp,])}))
      expression <- do.call("c", sapply(expList, function(x){return(x[gene,])}))
      expression <- qnorm((rank(expression,na.last="keep")-0.5)/sum(!is.na(expression)))
      covariate <- do.call("c", sapply(covList, function(x){return(x[,cov])}))
      cohortLables <- factor(rep(names(cohortSampleCount), times  = cohortSampleCount))
      
      #### Mega LRT
      
      fullModel <- lm(expression ~ genotypes * covariate + genotypes * cohortLables)
      baseModel <- lm(expression ~ genotypes + covariate + genotypes * cohortLables)
      
      logLik_full <- logLik(fullModel)
      logLik_base <- logLik(baseModel)
      
      LRT_stat <- -2 * (logLik_base - logLik_full)
      
      df <- attr(logLik_full, "df") - attr(logLik_base, "df")
      p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
      
      
      if(coef(fullModel)["genotypes:covariate"] < 0){
        res[1] <- qnorm(p_value/2)
      } else {
        res[1] <- -qnorm(p_value/2)
      }

      #### Meta LRT
      
      cohortZ <- sapply(cohorts, function(cohort){
        
        exp <- expList[[cohort]][gene,]
        exp <- qnorm((rank(exp,na.last="keep")-0.5)/sum(!is.na(exp)))
        
        fullModel <- lm(exp ~ genoList[[cohort]][snp,] * covList[[cohort]][,cov])
        baseModel <- lm(exp ~ genoList[[cohort]][snp,] + covList[[cohort]][,cov])
        
        
        logLik_full <- logLik(fullModel)
        logLik_base <- logLik(baseModel)
        
        LRT_stat <- -2 * (logLik_base - logLik_full)
        
        # Compute p-value using Chi-square Distribution
        df <- attr(logLik_full, "df") - attr(logLik_base, "df")
        p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
        
        if(coef(fullModel)[4] < 0){
          return(qnorm(p_value/2))
        } else {
          return(-qnorm(p_value/2))
        }
        
      })
      
      res[2] <- sum(cohortZ * cohortSampleCount)/sqrt(sum(cohortSampleCount * cohortSampleCount))

      
    },
    error = function(cond) {
      message("error in ", i)
    })
  return(res)
})
str(zScoreIntList)

zScoreInt <- do.call(rbind, zScoreIntList)

plot(zScoreInt[,1], zScoreInt[,2])

par(pty="s")
plot(zScoreInt[,2], distTest[,3], xlab = "Interaction z-score with force normal", ylab = "Interaction z-score without force normal", col = adjustcolor(palette()[1], alpha.f = 0.6), pch = 16, asp = 1)
abline(0,1)


plot(distTest[,4], abs(zScoreInt[,2])-abs(distTest[,3]), xlab = "Kurtosis of eqtl gene expression", ylab = "Difference in interaction z-score", col = adjustcolor(palette()[1], alpha.f = 0.6), pch =16)
