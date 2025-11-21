
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
id_stx3 <- "ENSG00000166900"
id_nod2Snp <- "rs1981760"

plot(exp[,id_stx3], exp[,id_nod2], col = geno[id_nod2Snp,] + 1)

model1 <- lm(exp[,id_nod2] ~ exp[,id_stx3] + geno[id_nod2Snp,])
model2 <- lm(exp[,id_nod2] ~ exp[,id_stx3] * geno[id_nod2Snp,])

model1LogLik <- logLik(model1)
model2LogLik <- logLik(model2)


model1LogLik
model2LogLik


(teststat <- -2 * (as.numeric(model1LogLik)-as.numeric(model2LogLik)))

pchisq(teststat, df = 1, lower.tail = FALSE)


summary(lm(exp[,id_nod2] ~ exp[,id_stx3] * geno[id_nod2Snp,]))




expGeneMean <- apply(exp, 2, mean)
sampleQuality <- cor(t(exp), expGeneMean)


eqtls <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline//data/sign_qtls_cistrans.txt.gz"))
str(eqtls)
eqtls <- eqtls[eqtls[,"feature_id"] %in% colnames(exp),]
eqtls <- eqtls[eqtls[,"snp_id"] %in% rownames(geno),]

eqtl <- eqtls[9966,]
sampleQualityInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * sampleQuality))
  return(lmRes$coefficients["geno[snp, ]:sampleQuality","t value"])
  
})

hist(sampleQualityInteractionT, breaks = 100)
View(sampleQualityInteractionT)
head(order(abs(sampleQualityInteractionT), decreasing =T))

stx3InteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * exp[,id_stx3]))
  return(lmRes$coefficients["geno[snp, ]:exp[, id_stx3]","t value"])
  
})
hist(stx3InteractionT, breaks = 100)

stx3CovSampleQualityInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * exp[,id_stx3] + geno[snp,]  * sampleQuality))
  return(lmRes$coefficients["geno[snp, ]:exp[, id_stx3]","t value"])
  
})
hist(stx3CovSampleQualityInteractionT, breaks = 100)

plot(stx3InteractionT, sampleQualityInteractionT)
cor.test(stx3InteractionT, sampleQualityInteractionT)
plot(stx3InteractionT, stx3CovSampleQualityInteractionT)
cor.test(stx3InteractionT, stx3CovSampleQualityInteractionT)

correctionModel <- lm(stx3InteractionT ~ sampleQualityInteractionT)

stx3PostcorrectSampleQualityInteractionT <- residuals(correctionModel)

plot(stx3CovSampleQualityInteractionT, stx3PostcorrectSampleQualityInteractionT)
cor.test(stx3CovSampleQualityInteractionT, stx3PostcorrectSampleQualityInteractionT)


layout(matrix(1:2,nrow = 1))
plot(stx3InteractionT, stx3CovSampleQualityInteractionT, main = "normal model vs\n model with term for rna quality. r=0.95")
plot(stx3PostcorrectSampleQualityInteractionT, stx3CovSampleQualityInteractionT, main = "Corrected normal model vs\n model with term for rna quality. r=0.98")



genoPc1InteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * genoPcs[,1]))
  return(lmRes$coefficients["geno[snp, ]:genoPcs[, 1]","t value"])
  
})

hist(genoPc1InteractionT, breaks = 100)
plot(genoPc1InteractionT, )

plot(stx3InteractionT, genoPc1InteractionT)
cor.test(stx3InteractionT, genoPc1InteractionT)


stx3CovSampleQualityGenoPcsInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * exp[,id_stx3] + geno[snp,]  * sampleQuality + geno[snp,] * genoPcs[,1] + geno[snp,] * genoPcs[,2] + geno[snp,] * genoPcs[,3] + geno[snp,] * genoPcs[,4] ))
  return(lmRes$coefficients["geno[snp, ]:exp[, id_stx3]","t value"])
  
})
hist(stx3CovSampleQualityGenoPcsInteractionT, breaks = 100)
plot(stx3InteractionT, stx3CovSampleQualityGenoPcsInteractionT)
cor.test(stx3InteractionT, stx3CovSampleQualityGenoPcsInteractionT)

sampleQualityGenoPcsInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * sampleQuality + geno[snp,] * genoPcs[,1] + geno[snp,] * genoPcs[,2] + geno[snp,] * genoPcs[,3] + geno[snp,] * genoPcs[,4] ))
  return(lmRes$coefficients[8:12,"t value"])
  
})
str(sampleQualityGenoPcsInteractionT)
str(stx3InteractionT)

correctionModel2 <- lm(stx3InteractionT ~ t(sampleQualityGenoPcsInteractionT))
summary(correctionModel2)
stx3PostcorrectSampleQualityGenoPcsInteractionT <- residuals(correctionModel2)

plot(stx3CovSampleQualityGenoPcsInteractionT, stx3PostcorrectSampleQualityGenoPcsInteractionT)
cor.test(stx3CovSampleQualityInteractionT, stx3PostcorrectSampleQualityGenoPcsInteractionT)

plot(stx3InteractionT, stx3PostcorrectSampleQualityGenoPcsInteractionT)
cor.test(stx3InteractionT, stx3PostcorrectSampleQualityGenoPcsInteractionT)


#############


str(sampleQuality)
str(exp)
geneXsq <- cor(sampleQuality, exp)[1,]
hist()
which.max(geneXsq)

geneAsCov <- "ENSG00000185219"
geneXsq[geneAsCov]


geneInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * exp[,geneAsCov]))
  return(lmRes$coefficients["geno[snp, ]:exp[, geneAsCov]","t value"])
  
})
hist(geneInteractionT, breaks = 100)


geneCovSampleQualityGenoPcsInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * exp[,geneAsCov] + geno[snp,]  * sampleQuality + geno[snp,] * genoPcs[,1] + geno[snp,] * genoPcs[,2] + geno[snp,] * genoPcs[,3] + geno[snp,] * genoPcs[,4] ))
  return(lmRes$coefficients["geno[snp, ]:exp[, geneAsCov]","t value"])
  
})
hist(geneCovSampleQualityGenoPcsInteractionT, breaks = 100)


plot(geneInteractionT, geneCovSampleQualityGenoPcsInteractionT)
cor.test(geneInteractionT, geneCovSampleQualityGenoPcsInteractionT)
  

plot(geneInteractionT, t(sampleQualityGenoPcsInteractionT)[,1])
abline(lm(t(sampleQualityGenoPcsInteractionT)[,1] ~ geneInteractionT))
  cor(geneInteractionT, t(sampleQualityGenoPcsInteractionT)[,1])

summary(lm(geneInteractionT ~ t(sampleQualityGenoPcsInteractionT)[,1]))

geneCorrectionModel2 <- lm(geneInteractionT ~ t(sampleQualityGenoPcsInteractionT)[,1])
summary(geneCorrectionModel2)
genePostcorrectSampleQualityGenoPcsInteractionT <- residuals(geneCorrectionModel2)

plot(geneCovSampleQualityGenoPcsInteractionT, geneInteractionT)
cor.test(geneCovSampleQualityGenoPcsInteractionT, geneInteractionT)

plot(genePostcorrectSampleQualityGenoPcsInteractionT, geneCovSampleQualityGenoPcsInteractionT)
cor.test(genePostcorrectSampleQualityGenoPcsInteractionT, geneCovSampleQualityGenoPcsInteractionT)



layout(matrix(1:2,nrow = 1))
plot(geneInteractionT, geneCovSampleQualityGenoPcsInteractionT, xlab = "base model (1)", ylab = "model with tech covariats (3)", main = "normal model vs\n model with quality and geno PCs. r=0.79")
plot(genePostcorrectSampleQualityGenoPcsInteractionT, geneCovSampleQualityGenoPcsInteractionT, xlab = "base model (1) corrected for technical model (2)", ylab = "model with tech covariats (3)", main = "Corrected normal model vs\n model with quality and geno PCs. r=0.96")

genePreCorrectedCovInteractionT <- apply(eqtls, 1, function(eqtl){
  
  snp <- eqtl[2]
  gene <- eqtl[1]
  expCor <- residuals(lm(exp[,gene] ~ geno[snp,] * sampleQuality + geno[snp,] * genoPcs[,1] + geno[snp,] * genoPcs[,2] + geno[snp,] * genoPcs[,3] + geno[snp,] * genoPcs[,4] ))
      
  lmRes <- summary(lm(expCor ~ geno[snp,]  * exp[,geneAsCov]))
  return(lmRes$coefficients["geno[snp, ]:exp[, geneAsCov]","t value"])
  
})

plot(genePreCorrectedCovInteractionT, geneInteractionT)
cor.test(genePreCorrectedCovInteractionT, geneInteractionT)

plot(geneInteractionT, geneCovSampleQualityGenoPcsInteractionT)
cor.test(geneInteractionT, geneCovSampleQualityGenoPcsInteractionT)

layout(1)
plot(genePostcorrectSampleQualityGenoPcsInteractionT, genePreCorrectedCovInteractionT, xlab=  "Post corrected interaction", ylab ="Pre corrected interaction")
cor.test(genePostcorrectSampleQualityGenoPcsInteractionT, genePreCorrectedCovInteractionT)











library(parallel)

tableGeneCovSampleQualityGenoPcsInteractionT <- mclapply(colnames(exp)[1:10], function(geneAsCov){
  
  
  geneCovSampleQualityGenoPcsInteractionT <- apply(eqtls, 1, function(eqtl){
    
    snp <- eqtl[2]
    gene <- eqtl[1]
    lmRes <- summary(lm(exp[,gene] ~ geno[snp,]  * exp[,geneAsCov] + geno[snp,]  * sampleQuality + geno[snp,] * genoPcs[,1] + geno[snp,] * genoPcs[,2] + geno[snp,] * genoPcs[,3] + geno[snp,] * genoPcs[,4] ))
    return(lmRes$coefficients["geno[snp, ]:exp[, geneAsCov]","t value"])
    
  })
  
  return(geneCovSampleQualityGenoPcsInteractionT)
  
}, mc.cores = 10)

tableGeneCovSampleQualityGenoPcsInteractionT <- do.call(cbind, tableGeneCovSampleQualityGenoPcsInteractionT)

str(tableGeneCovSampleQualityGenoPcsInteractionT)



