setwd("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/")




library(readr)

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

#BiosinteractionZ <- readDoubleMatrix("biosInteractionsInt25pcIv/metaZTest.txt")
#saveRDS(BiosinteractionZ, file = "biosInteractionsInt25pcIv/metaZTest.rds")
BiosinteractionZ <- readRDS("biosInteractionsInt25pcIv/metaZTest.rds")


BiosinteractionZNoInt <- readRDS("biosInteractionsPcCor/metaZ.rds")

BiosinteractionZNoInt[grep(nod2, rownames(BiosinteractionZNoInt)),stx3]

a <- intersect(rownames(BiosinteractionZNoInt), rownames(BiosinteractionZ))
layout(1)
par(pty="s")
plot(BiosinteractionZNoInt[a,stx3], BiosinteractionZ[a,stx3],  xlim = c(-17,17), ylim = c(-17,17),asp = 1, col = adjustcolor(palette()[1], alpha.f = 0.7), pch = 16, main = "STX3 interactions", xlab = "No INT", ylab = "Best normalization for replicating co-eqtls")
abline(0,1)
cor.test(BiosinteractionZNoInt[a,stx3], BiosinteractionZ[a,stx3])

a[BiosinteractionZNoInt[a,stx3] > -5 & BiosinteractionZ[a,stx3] < -10 ]

BiosinteractionZNoInt[ "rs2302559-ENSG00000229644",stx3]
BiosinteractionZ["rs2302559-ENSG00000229644",stx3]

BiosinteractionCount <- readDoubleMatrix("biosInteractionsPcCor/metaSampleCount.txt")


stx3 <- "ENSG00000166900"
nod2 <- "ENSG00000167207"
adm <- "ENSG00000148926"


hist(BiosinteractionZ[,stx3])

BiosinteractionZ[grep(nod2, rownames(BiosinteractionZ)),stx3]

head(sort(abs(BiosinteractionZ[,stx3]),decreasing = T))

interactionZ2017 <- readDoubleMatrix("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/interactions/results2017/topCovZ.txt.gz")

str(interactionZ2017)




interactionZ2017_2 <- data.frame(cisGenes = factor(sapply(strsplit(rownames(interactionZ2017), "_"), function(x){return(x[2])})), zscore = interactionZ2017[,1])
str(interactionZ2017_2)


interactionZ2017_3 <- aggregate(zscore ~ cisGenes,FUN = function(x){x[which.max(abs(x))]}, data = interactionZ2017_2)
rownames(interactionZ2017_3) <- interactionZ2017_3$cisGenes
str(interactionZ2017_3)

hist(interactionZ2017_3$zscore)

str(BiosinteractionZ[,stx3])

BiosinteractionZ_stx3_2 <- data.frame(cisGenes = factor(sapply(strsplit(names(BiosiznteractionZ[,stx3]), "-"), function(x){return(x[2])})), zscore = BiosinteractionZ[,stx3])

BiosinteractionZ_stx3_3 <- aggregate(zscore ~ cisGenes,FUN = function(x){x[which.max(abs(x))]}, data = BiosinteractionZ_stx3_2)
rownames(BiosinteractionZ_stx3_3) <- BiosinteractionZ_stx3_3$cisGenes

sharedGenes <- intersect(interactionZ2017_3$cisGenes, BiosinteractionZ_stx3_3$cisGenes)

length(sharedGenes)

library(RColorBrewer)
layout(matrix(1:2, nrow =1))
par(pty="s")
palette(adjustcolor(brewer.pal(n = 3, name = "Accent"), alpha.f = 0.6))
plot(interactionZ2017_3[sharedGenes, 2], BiosinteractionZ_stx3_3[sharedGenes, 2], pch = 16, col = 1, xlab = "2017 interaction z-scores of first module", ylab = "New bios interaction z-score", main = "STX3 covariate", asp = 1)
abline(lm(BiosinteractionZ_stx3_3[sharedGenes, 2] ~ interactionZ2017_3[sharedGenes, 2]))
cor.test(BiosinteractionZ_stx3_3[sharedGenes, 2], interactionZ2017_3[sharedGenes, 2])

palette(adjustcolor(brewer.pal(n = 3, name = "Accent"), alpha.f = 0.6))
plot(abs(interactionZ2017_3[sharedGenes, 2]), abs(BiosinteractionZ_stx3_3[sharedGenes, 2]), pch = 16, col = 1, xlab = "2017 interaction z-scores of first module", ylab = "New bios interaction z-score", main = "STX3 covariate abs Z-scores", asp = 1)
abline(lm( abs(BiosinteractionZ_stx3_3[sharedGenes, 2]) ~  abs(interactionZ2017_3[sharedGenes, 2])))
cor.test( abs(BiosinteractionZ_stx3_3[sharedGenes, 2]),  abs(interactionZ2017_3[sharedGenes, 2]))
 abline(0,1)

outliers <- abs(interactionZ2017_3[sharedGenes, 2]) > 10 & abs(BiosinteractionZ_stx3_3[sharedGenes, 2]) < 5 
plot(abs(interactionZ2017_3[sharedGenes[outliers], 2]), abs(BiosinteractionZ_stx3_3[sharedGenes[outliers], 2]), pch = 16, col = 1, xlab = "2017 interaction z-scores of first module", ylab = "New bios interaction z-score", main = "STX3 covariate abs Z-scores", asp = 1)

cat(sharedGenes[outliers], sep = "\n")

interactionZ2017_3["ENSG00000145730",2]
BiosinteractionZ_stx3_3["ENSG00000145730", 2]

str(BiosinteractionZ)

correctionModel <- lm(BiosinteractionZ[,stx3] ~ 
                        BiosinteractionZ[,"Comp1"] + 
                        BiosinteractionZ[,"Comp2"] + 
                        BiosinteractionZ[,"Comp3"] + 
                        BiosinteractionZ[,"Comp4"] + 
                        BiosinteractionZ[,"Comp5"] + 
                        BiosinteractionZ[,"Comp6"] + 
                        BiosinteractionZ[,"Comp7"] + 
                        BiosinteractionZ[,"Comp8"] + 
                        BiosinteractionZ[,"Comp9"] + 
                        BiosinteractionZ[,"Comp10"] + 
                        BiosinteractionZ[,"Comp11"] + 
                        BiosinteractionZ[,"Comp12"] + 
                        BiosinteractionZ[,"Comp13"] + 
                        BiosinteractionZ[,"Comp14"] + 
                        BiosinteractionZ[,"Comp15"] + 
                        BiosinteractionZ[,"Comp16"] + 
                        BiosinteractionZ[,"Comp17"] + 
                        BiosinteractionZ[,"Comp18"] + 
                        BiosinteractionZ[,"Comp19"] + 
                        BiosinteractionZ[,"Comp20"] + 
                        BiosinteractionZ[,"Comp21"] + 
                        BiosinteractionZ[,"Comp22"] + 
                        BiosinteractionZ[,"Comp23"] + 
                        BiosinteractionZ[,"Comp24"] + 
                        BiosinteractionZ[,"Comp25"]  )
plot(residuals(correctionModel))



correctionModel <- lm(BiosinteractionZ[,stx3] ~ 
                        BiosinteractionZ[,"Comp1"]
                        )
summary(correctionModel)

mean(BiosinteractionZ[,"Comp1"])

plot(residuals(correctionModel))

str(residuals(correctionModel))

summary(correctionModel)
residuals(correctionModel)[grep(nod2, names(residuals(correctionModel)))]

BiosinteractionZ[grep(nod2, rownames(BiosinteractionZ)),stx3]

str(BiosinteractionZ)
plot(BiosinteractionZ[,stx3],residuals(correctionModel))




ludeCisTransP <- read.delim("ludeCisTrans/Positives.txt")
ludeCisTransN <- read.delim("ludeCisTrans/Negatives.txt")

ludeCisTransN <- ludeCisTransN[,-7]


BiosinteractionZ_2 <- cbind(data.frame(cisGenes = factor(sapply(strsplit(rownames(BiosinteractionZ), "-"), function(x){return(x[2])}))), BiosinteractionZ)
str(BiosinteractionZ_2)


BiosinteractionZAgg <- aggregate(. ~ cisGenes, data = BiosinteractionZ_2, FUN = function(x){x[which.max(abs(x))]})

str(BiosinteractionZAgg)
rownames(BiosinteractionZAgg) <- BiosinteractionZAgg$cisGene

str(ludeCisTransP)
x <- ludeCisTransP[1,]
ludeCisTransP2 <- apply(ludeCisTransP, 1, function(x){
  a <- NA
  b <- NA
  if(x["CisGene"] %in% rownames(BiosinteractionZAgg) & x["TransGene"] %in% colnames(BiosinteractionZAgg)){
    a <- BiosinteractionZAgg[x["CisGene"],x["TransGene"]]
  } 
  if(x["TransGene"] %in% rownames(BiosinteractionZAgg) & x["CisGene"] %in% colnames(BiosinteractionZAgg)){
    b <- BiosinteractionZAgg[x["TransGene"],x["CisGene"]]
  } 
  c(x, a, b)
})

x <- ludeCisTransN2[1,]
ludeCisTransN2 <- apply(ludeCisTransN, 1, function(x){
  a <- NA
  b <- NA
  if(x["CisGene"] %in% rownames(BiosinteractionZAgg) & x["TransGene"] %in% colnames(BiosinteractionZAgg)){
    a <- BiosinteractionZAgg[x["CisGene"],x["TransGene"]]
  } 
  if(x["TransGene"] %in% rownames(BiosinteractionZAgg) & x["CisGene"] %in% colnames(BiosinteractionZAgg)){
    b <- BiosinteractionZAgg[x["TransGene"],x["CisGene"]]
  } 
  c(x, a, b)
})


ludeCisTransP2 <- data.frame(t(ludeCisTransP2))

colnames(ludeCisTransP2)[7] <- "transGeneCovariate"
colnames(ludeCisTransP2)[8] <- "cisGeneCovariate"

ludeCisTransN2 <- data.frame(t(ludeCisTransN2))
colnames(ludeCisTransN2)[7] <- "transGeneCovariate"
colnames(ludeCisTransN2)[8] <- "cisGeneCovariate"



ludeCisTransP2$transGeneCovariate <- as.numeric(ludeCisTransP2$transGeneCovariate)
ludeCisTransP2$cisGeneCovariate <- as.numeric(ludeCisTransP2$cisGeneCovariate)
ludeCisTransP2$Direction <- as.numeric(ludeCisTransP2$Direction       )

ludeCisTransN2$transGeneCovariate <- as.numeric(ludeCisTransN2$transGeneCovariate)
ludeCisTransN2$cisGeneCovariate <- as.numeric(ludeCisTransN2$cisGeneCovariate)
ludeCisTransN2$Direction <- as.numeric(ludeCisTransN2$Direction       )

plot(ludeCisTransP2[,"Direction"], ludeCisTransP2[,"cisGeneCovariate"])
plot(ludeCisTransP2[,"Direction"], ludeCisTransP2[,"transGeneCovariate"])




hist( ludeCisTransP2[,"cisGeneCovariate"], breaks = 100)
hist( ludeCisTransP2[,"transGeneCovariate"], breaks = 100)



library("vioplot")

vioplot( ludeCisTransP2[,"cisGeneCovariate"] ~ ludeCisTransP2[,"Direction"])

vioplot(abs(ludeCisTransP2[,"cisGeneCovariate"]), abs(ludeCisTransP2[,"transGeneCovariate"]), names = c("cis gene as covariate", "trans gene as covariate"))

abline(h= 1.595005)
abline(h= 1.966393)

t.test(abs(ludeCisTransP2[,"transGeneCovariate"]), abs(ludeCisTransP2[,"cisGeneCovariate"]))


vioplot(abs(ludeCisTransP2[,"cisGeneCovariate"]), abs(ludeCisTransN2[,"cisGeneCovariate"]), names = c("cis gene as covariate positive set", "cis gene as covariate negative set"))
t.test(abs(ludeCisTransP2[,"cisGeneCovariate"]), abs(ludeCisTransN2[,"cisGeneCovariate"]))



vioplot(abs(ludeCisTransN2[,"cisGeneCovariate"]), abs(ludeCisTransN2[,"transGeneCovariate"]), names = c("cis gene as covariate", "trans gene as covariate"), main = "negative set")


wilcox.test(abs(ludeCisTransP2[,"cisGeneCovariate"]), abs(ludeCisTransN2[,"cisGeneCovariate"]))

mean(abs(ludeCisTransP2[,"cisGeneCovariate"]), na.rm = T)
mean(abs(ludeCisTransN2[,"cisGeneCovariate"]), na.rm = T)


sum(abs(ludeCisTransP2[,"cisGeneCovariate"]) >=5, na.rm = T)
sum(abs(ludeCisTransN2[,"cisGeneCovariate"]) >=5, na.rm = T)



View(ludeCisTransP2[!is.na(ludeCisTransP2[,"cisGeneCovariate"]) & abs(ludeCisTransP2[,"cisGeneCovariate"]) >=5,])

write.table(ludeCisTransP2[!is.na(ludeCisTransP2[,"cisGeneCovariate"]) & abs(ludeCisTransP2[,"cisGeneCovariate"]) >=5,], file = "positivesWithStrongInteraction.txt", sep ="\t", quote = F, row.names = F)





