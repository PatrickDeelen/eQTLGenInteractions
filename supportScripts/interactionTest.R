
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
summary(lm(exp[,id_nod2] ~ exp[,id_stx3] * geno[id_nod2Snp,]))


g1 <- "ENSG00000035115"
g2 <- "ENSG00000035115"
g2 <- "ENSG00000101842"
snp <- "rs7605824"

plot(exp[,g1], exp[,g2], col = geno[snp,] + 1)
summary(lm(exp[,g2] ~ exp[,g1] * geno[snp,]))


limix <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/e7/640c18222e79a4499aaaa08279a1f4/limix_out/tValueAssociations_ENSG00000035115.txt.gz", row.names = 1))

limixP <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/e7/640c18222e79a4499aaaa08279a1f4/limix_out/tValueAssociations_ENSG00000035115_perm.txt.gz", row.names = 1))

limix2 <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/1a/8d26668dbb9dde29fcbc801e28ed26/limix_out/tValueAssociations_ENSG00000035115.txt.gz", row.names = 1))

limix2P <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/1a/8d26668dbb9dde29fcbc801e28ed26/limix_out/tValueAssociations_ENSG00000035115_perm.txt.gz", row.names = 1))

limix3 <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/4e/21b98e5718b845c921ca4d128ae96d/limix_out/pValueAssociations_ENSG00000035115.txt.gz", row.names = 1))

limix3P <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/4e/21b98e5718b845c921ca4d128ae96d/limix_out/pValueAssociations_ENSG00000035115_perm.txt.gz", row.names = 1))

str(limix3)

layout(matrix(1:2,ncol =2))
hist(limix3, breaks = 100, xlab = "Log ratio test p-value", main = "Real run")
hist(limix3P, breaks = 100, xlab = "Log ratio test p-value", main = "Permutation run")

layout(matrix(1:2,ncol =2))
hist(limix3[snp,], breaks = 100, xlab = "Log ratio test p-value", main = "Real run")
hist(limix3P[paste0(snp,"_0"),], breaks = 100, xlab = "Log ratio test p-value", main = "Permutation run 1")

limix3Z <- -qnorm(limix3/2)
limix3PZ <- -qnorm(limix3P/2)


snpCheck <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/4e/21b98e5718b845c921ca4d128ae96d/test.raw", sep = " ", header = T))

dataForTest <- as.matrix(read.delim(paste0("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/4e/21b98e5718b845c921ca4d128ae96d/Cov.",g2,".ForTest.txt"), sep = ",", header = F))
phenoForTest <- as.matrix(read.delim("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/interactions/ieQTL_nextflow_pipeline/work/4e/21b98e5718b845c921ca4d128ae96d/phenoForTest.txt", sep = ",", header = T))

exp2 <- exp[match(phenoForTest[,"X"], rownames(exp)),]
geno2 <- geno[,match(phenoForTest[,"X"], colnames(geno))]

str(snpCheck)
snpCheck2 <- snpCheck[match(phenoForTest[,"X"], snpCheck[,"IID"]),]  
str(snpCheck2)

plot(phenoForTest[,"X0"], dataForTest[,"V2"])

str(dataForTest[,"V2"])
str(exp)


plot(dataForTest[,"V2"], exp2[,g2])


plot(dataForTest[,"V4"], exp2[,g2] * (((geno2[snp,])-1)*-1)+1)

plot(dataForTest[,"V4"], dataForTest[,"V2"] * dataForTest[,"V3"])

summary(lm(phenoForTest[,"X0"] ~ dataForTest[,"V2"] + dataForTest[,"V3"] + dataForTest[,"V4"]))

limix3[snp, g2]


table(dataForTest[,"V3"], (((geno2[snp,])-1)*-1)+1)

table(snpCheck2[,"rs7605824_A"], geno2[snp,])

table(snpCheck2[,"rs7605824_A"])


plot(abs(limix2[snp,]), limix3Z[snp,])
plot(limix2P[paste0(snp,"_0"),], limix3P[paste0(snp,"_0"),])

head(abs(limix2[snp,]))
head(-limix3Z[snp,])


str(limix[snp,])

hist(limix[snp,])

limix[snp,g2]
limixP[paste0(snp,"_0"),g2]


hist(limixP[paste0(snp,"_0"),])

tail(sort(abs(limix2P[paste0(snp,"_0"),])))

tail(sort(abs(limix2[snp,])))

plot(limix[snp,], limixP[paste0(snp,"_0"),], 	pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.5), cex = 0.5)

layout(1)
x <- cor.test(limix3[snp, ], limix3P[paste0(snp,"_0"), ])
plot(limix3Z[snp, ], limix3PZ[paste0(snp,"_0"), ], 	pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.8), cex = 0.3, xlab = "Log ratio Z-score of ENSG00000035115", ylab = "Log ratio Z-score of permutation round 1", main = paste0("r = ",round(x$estimate,3),", p-value = ",x$p.value,""))
abline(lm(limix3PZ[paste0(snp,"_0"), ] ~ limix3Z[snp,]))

layout(1)
x <- cor.test(limix3PZ[paste0(snp,"_0"), ], limix3P[paste0(snp,"_0"), ])
plot(limix3PZ[paste0(snp,"_1"), ], limix3PZ[paste0(snp,"_0"), ], 	pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.8), cex = 0.3, xlab = "Log ratio Z-score of permutation round 2", ylab = "Log ratio Z-score of permutation round 1", main = paste0("r = ",round(x$estimate,3),", p-value = ",x$p.value,""))
abline(lm(limix3PZ[paste0(snp,"_0"), ] ~ limix3PZ[paste0(snp,"_1"), ]))




plot(limix3[snp,-1586 ], limix3P[paste0(snp,"_0"),-1586 ], 	pch = 16, col=adjustcolor("dodgerblue2", alpha.f = 0.8), cex = 0.3, xlab = "Interaction T-stats of ENSG00000035115", ylab = "Interaction T-stats of permutations", main = "r = 0.03, p-value = 1.073e-05")
abline(lm(limix3P[paste0(snp,"_0"),-1586 ] ~ limix3[snp,-1586 ]))

cor.test(limix3[snp,-1586 ], limix3P[paste0(snp,"_0"),-1586 ])

str(limix3P)
cor.test(limix3P[paste0(snp,"_0"),-1586 ], limix3P[paste0(snp,"_1"),-1586 ])

which.max(limix[snp,])


str(exp)

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



