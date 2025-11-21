
library(readr)
#change X1 in case of specified header for row names. Use ...1 if first colun has no ID
colTypes <- cols( .default = col_double(),  `ID` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426//1_DataQC/out/LL/outputfolder_exp/exp_data_QCd/exp_data_preprocessed.txt", delim = "\t", quote = "", col_types = colTypes)
exp <- as.matrix(table_tmp[,-1])
rownames(exp) <- table_tmp[,1][[1]]
rm(table_tmp)


colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/testing/merged2.dosages.txt", delim = "\t", quote = "", col_types = colTypes)
geno <- as.matrix(table_tmp[,-1])
rownames(geno) <- table_tmp[,1][[1]]
rm(table_tmp)
geno <- geno[,rownames(exp)]


colTypes <- cols( .default = col_double(),  `...1` = col_character())
table_tmp <- read_delim("/groups/umcg-bios/tmp01/projects/BIOS_for_eQTLGenII/pipeline/20220426//1_DataQC/out/LL/outputfolder_gen/gen_PCs/GenotypePCs.txt", delim = "\t", quote = "", col_types = colTypes)
genoPcs <- as.matrix(table_tmp[,-1])
rownames(genoPcs) <- table_tmp[,1][[1]]
rm(table_tmp)
genoPcs <- genoPcs[rownames(exp),]

all(rownames(exp) == colnames(geno) )
all(rownames(exp) == rownames(genoPcs) )


id_nod2 <- "ENSG00000167207"
id_covGene <- "ENSG00000185219"  #ENSG00000166900
id_nod2Snp <- "rs1981760"

plot(exp[,id_covGene], exp[,id_nod2], col = geno[id_nod2Snp,] + 1)

expGeneMean <- apply(exp, 2, mean)
sampleQuality <- cor(t(exp), expGeneMean)


plot(exp[,id_covGene], sampleQuality)
cor.test(exp[,id_covGene], sampleQuality)

cor.test(exp[,id_covGene], genoPcs[,1])


eqtls <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline//data/sign_qtls_cistrans.txt.gz"))
str(eqtls)
eqtls <- eqtls[eqtls[,"feature_id"] %in% colnames(exp),]
eqtls <- eqtls[eqtls[,"snp_id"] %in% rownames(geno),]

eqtl <- eqtls[9966,]
eqtl <- eqtls[32325,]

sampleQualityInteractionZ <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + sampleQuality)
  model2 <- lm(exp[,gene] ~ geno[snp,] * sampleQuality)
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = F)
  
  z <- qnorm(logP/2, log.p = F)
  if(coefficients(model2)[4] > 0){
    z <- abs(z)
  } else {
    z <- -abs(z)
  }
  
  return(z)
  
})

hist(sampleQualityInteractionZ)

stx3InteractionZ <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + exp[,id_covGene])
  model2 <- lm(exp[,gene] ~ geno[snp,] * exp[,id_covGene])
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  pchisq(teststat, df = 1, lower.tail = FALSE, log.p = F)
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = F)
  
  z <- qnorm(logP/2, log.p = F)
  if(coefficients(model2)[4] > 0){
    z <- abs(z)
  } else {
    z <- -abs(z)
  }
  
  return(z)
  
})

hist(stx3InteractionZ)

genoPc1InteractionZ <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + genoPcs[,1])
  model2 <- lm(exp[,gene] ~ geno[snp,] * genoPcs[,1])
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = F)
  
  z <- qnorm(logP/2, log.p = F)
  
  if(coefficients(model2)[4] > 0){
    z <- abs(z)
  } else {
    z <- -abs(z)
  }
  
  summary(model2)$coefficients
  
  return(z)
  
})

hist(genoPc1InteractionZ)

plot(sampleQualityInteractionZ, stx3InteractionZ)
plot(sampleQualityInteractionZ, genoPc1InteractionZ)
plot(stx3InteractionZ, genoPc1InteractionZ)


stx3InteractionZfull <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + geno[snp,] * genoPcs[,1] + geno[snp,] * sampleQuality + exp[,id_covGene])
  model2 <- lm(exp[,gene] ~ geno[snp,] + geno[snp,] * genoPcs[,1] + geno[snp,] * sampleQuality + geno[snp,] * exp[,id_covGene])
  
  summary(model2)
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = F)
  
  z <- qnorm(logP/2, log.p = F)
  if(coefficients(model2)[8] > 0){
    z <- abs(z)
  } else {
    z <- -abs(z)
  }
  
  return(z)
  
})


hist(stx3InteractionZfull)

plot(stx3InteractionZ, stx3InteractionZfull)



correctionModel <- lm(stx3InteractionZ ~ sampleQualityInteractionZ + genoPc1InteractionZ)
summary(correctionModel)
stx3InteractionZReconstructed <- residuals(correctionModel)

plot(stx3InteractionZReconstructed, stx3InteractionZfull)
plot(stx3InteractionZ, stx3InteractionZReconstructed)


layout(matrix(1:3, nrow = 1))
par(mar = c(5, 5, 4, 2) + 0.1)
plot(stx3InteractionZ, stx3InteractionZfull, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, 
     col=adjustcolor("dodgerblue2", alpha.f = 0.5),
     xlab = "Interaction Z-scores STX3\nSimple model",
     ylab = "Interaction Z-scores STX3\nComplex model")

plot(stx3InteractionZ, stx3InteractionZReconstructed, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, 
     col=adjustcolor("dodgerblue2", alpha.f = 0.5),
     xlab = "Interaction Z-scores STX3\nSimple model",
     ylab = "Interaction Z-scores STX3\nReconstructed score")

plot(stx3InteractionZReconstructed, stx3InteractionZfull, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, 
     col=adjustcolor("dodgerblue2", alpha.f = 0.5),
     xlab = "Interaction Z-scores STX3\nReconstructed score",
     ylab = "Interaction Z-scores STX3\nComplex model")

cor.test(stx3InteractionZ, stx3InteractionZfull)

cor.test(stx3InteractionZReconstructed, stx3InteractionZfull)







sampleQualityInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + sampleQuality)
  model2 <- lm(exp[,gene] ~ geno[snp,] * sampleQuality)
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = T)
  
  z <- qnorm(logP, log.p = T)
  if(coefficients(model2)[4] > 0){
    z <- z * -1
  }
  
  #return(z)
  return(summary(model2)$coefficients[4,"t value"])
  
})

hist(sampleQualityInteractionZ)

stx3InteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + exp[,id_covGene])
  model2 <- lm(exp[,gene] ~ geno[snp,] * exp[,id_covGene])
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = F)
  
  
  z <- qnorm(logP/2, log.p = F)
  if(coefficients(model2)[4] > 0){
    z <- z * -1
  }
  
  #return(z)
  return(summary(model2)$coefficients[4,"t value"])
  
})

hist(stx3InteractionZ)

genoPc1InteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + genoPcs[,1])
  model2 <- lm(exp[,gene] ~ geno[snp,] * genoPcs[,1])
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = T)
  
  z <- qnorm(logP, log.p = T)
  if(coefficients(model2)[4] > 0){
    z <- z * -1
  }
  
  #return(z)
  
  return(summary(model2)$coefficients[4,"t value"])
  
  
  
})


stx3InteractionTfull <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  
  model1 <- lm(exp[,gene] ~ geno[snp,] + geno[snp,] * genoPcs[,1] + geno[snp,] * sampleQuality + exp[,id_covGene])
  model2 <- lm(exp[,gene] ~ geno[snp,] + geno[snp,] * genoPcs[,1] + geno[snp,] * sampleQuality + geno[snp,] * exp[,id_covGene])
  
  
  model1LogLik <- logLik(model1)
  model2LogLik <- logLik(model2)
  
  teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik))
  
  logP <- pchisq(teststat, df = 1, lower.tail = FALSE, log.p = T)
  
  z <- qnorm(logP, log.p = T)
  if(coefficients(model2)[8] > 0){
    z <- z * -1
  }
  
  #return(z)
  
  return(summary(model2)$coefficients[8,"t value"])
  
})

correctionModel <- lm(stx3InteractionZ ~ sampleQualityInteractionZ + genoPc1InteractionZ)
summary(correctionModel)
stx3InteractionZReconstructed <- residuals(correctionModel)

plot(stx3InteractionZReconstructed, stx3InteractionZfull)
plot(stx3InteractionZ, stx3InteractionZReconstructed)


layout(matrix(1:3, nrow = 1))
par(mar = c(5, 5, 4, 2) + 0.1)
plot(stx3InteractionZ, stx3InteractionZfull, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, 
     col=adjustcolor("dodgerblue2", alpha.f = 0.5),
     xlab = "Interaction Z-scores cov gene\nSimple model",
     ylab = "Interaction Z-scores cov gene\nComplex model")

plot(stx3InteractionZ, stx3InteractionZReconstructed, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, 
     col=adjustcolor("dodgerblue2", alpha.f = 0.5),
     xlab = "Interaction Z-scores cov gene\nSimple model",
     ylab = "Interaction Z-scores cov gene\nReconstructed score")



plot(stx3InteractionZReconstructed, stx3InteractionZfull, bg = adjustcolor("dodgerblue2", alpha.f = 0.3), pch = 21, 
     col=adjustcolor("dodgerblue2", alpha.f = 0.5),
     xlab = "Interaction Z-scores cov gene\nReconstructed score",
     ylab = "Interaction Z-scores cov gene\nComplex model")


cor.test(stx3InteractionZ, stx3InteractionZfull)

cor.test(stx3InteractionZReconstructed, stx3InteractionZfull)

         