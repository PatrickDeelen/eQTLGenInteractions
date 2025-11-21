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
    genoList[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsInt25pcIv/",cohort,"/merged2.dosages.txt"))
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
    expList[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsInt25pcIv/",cohort,"/correctedExpression.txt"))
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
    covList[[cohort]] <- readDoubleMatrix(paste0("biosInteractionsInt25pcIv/",cohort,"/covariates.combined.txt"))
    sharedCovariatsTmp <- c(sharedCovariatsTmp, colnames(covList[[cohort]]))
  }
  sharedCovariatesCount <- table(sharedCovariatsTmp)
  sharedCovariats <- names(sharedCovariatesCount)[sharedCovariatesCount == length(cohorts)]
  for(cohort in cohorts){
    covList[[cohort]] <- covList[[cohort]][,sharedCovariats]
  }
  str(covList)
  
  eqtls <- read.delim("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/independent_variants_filtered_lbf2_mlog10p5_annotated_20250325.b.txt.gz")
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
  
  save(genoList, sharedVariants, sharedExpressionGenes, expList, sharedCovariats, covList, eqtls, cohortSampleCount, file = "biosInteractionsInt25pcIv/compareMetaMethodsData.RData")
  
} else {
  load("biosInteractionsInt25pcIv/compareMetaMethodsData.RData")
}


BiosinteractionZ <- readRDS("biosInteractionsInt25pcIv/metaZTest.rds")

head(sort(abs(BiosinteractionZ[,stx3]),decreasing = T))

#top effect of stx3
cov <- stx3
gene <- "ENSG00000167528"
snp <- "rs2634667"

cov <- stx3
gene <- nod2
snp <- nod2QtlSnp

#Interaction much more signficant when do int
cov <- stx3
gene <- "ENSG00000229644"
snp <- "rs2302559"


layout(rbind(matrix(rep(1,length(cohorts)), nrow = 1),matrix(2:((length(cohorts)*2)+1), nrow = 2)), heights = c(1,5,5))

par(mar=c(0,0,0,0), pty = "m")
plot.new()
plot.window(0:1,0:1)
text(0.5,0.5,paste("eQTL: ", snp, "-", gene, " covariate: ", cov ),cex = 2)

par(mar=c(5, 4, 4, 2), pty = "s")

for(cohort in cohorts){
  d <- data.frame(exp = expList[[cohort]][gene,], cov = covList[[cohort]][,cov], geno = genoList[[cohort]][snp,])
  
  fullModel <- lm(exp ~ geno * cov, data = d)
  baseModel <- lm(exp ~ geno + cov, data = d)
  
  
  logLik_full <- logLik(fullModel)
  logLik_base <- logLik(baseModel)
  
  LRT_stat <- -2 * (logLik_base - logLik_full)
  
  # Compute p-value using Chi-square Distribution
  df <- attr(logLik_full, "df") - attr(logLik_base, "df")
  p_value <- pchisq(LRT_stat, df, lower.tail = FALSE)
  
  z<-qnorm(p_value/2)
  if(coef(fullModel)[4] < 0){
    z <- z*-1
  }
  print(paste0("cohort ", cohort, " Z-score ", z))
 
  
  #plot(genoList[[cohort]][snp,],expList[[cohort]][gene,], col = round(genoList[[cohort]][snp,]) + 1, main = cohort, xlab = "Genotype dosage", ylab = "eQTL gen expression")
  #abline(lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,]))
  
  vioplot(expList[[cohort]][gene,] ~ round(genoList[[cohort]][snp,]), main = cohort, xlab = "Genotype dosage", ylab = "eQTL gen expression", col = (min(round(genoList[[cohort]][snp,])):max(round(genoList[[cohort]][snp,])))+1, at = (min(round(genoList[[cohort]][snp,])):max(round(genoList[[cohort]][snp,])))  )
  abline(lm(expList[[cohort]][gene,] ~ genoList[[cohort]][snp,]))
  
  
  plot(covList[[cohort]][,cov], expList[[cohort]][gene,], col = round(genoList[[cohort]][snp,]) + 1, pch =16, main = cohort, xlab = "Covariate expression", ylab = "eQTL gen expression")
  
  for(dosage in (min(round(genoList[[cohort]][snp,])):max(round(genoList[[cohort]][snp,])))){
    x <- data.frame("geno" = dosage, "cov" = range(covList[[cohort]][,cov]), check.names = F)
    y <- predict.lm(fullModel, x)  
    lines(x$cov,y, col = dosage + 1, lwd = 3)
  }

}

#expListNoInt <- expList

plot(expList[["RS"]][gene,], expListNoInt[["RS"]][gene,], col = round(genoList[[cohort]][snp,]) + 1, pch =16, xlab = "Expression data with INT", ylab = "Expression data no INT")
abline(lm(expListNoInt[["RS"]][gene,]~expList[["RS"]][gene,]))


plot(expList[["RS"]][gene,], expListNoInt[["RS"]][gene,], col = round(genoList[["RS"]][snp,]) + 1, pch =16, xlab = "Expression data with INT", ylab = "Expression data no INT")


exp <- qnorm((rank(expListNoInt[["RS"]][gene,],na.last="keep")-0.5)/sum(!is.na(expListNoInt[["RS"]][gene,])))
plot(exp, expListNoInt[["RS"]][gene,])
plot(expList[["RS"]][gene,], qnorm((rank(expListNoInt[["RS"]][gene,],na.last="keep")-0.5)/sum(!is.na(expListNoInt[["RS"]][gene,]))), col = round(genoList[["RS"]][snp,]) + 1, pch =16, xlab = "Expression data INT before PC correction", ylab = "Expression data INT after PC correction")


