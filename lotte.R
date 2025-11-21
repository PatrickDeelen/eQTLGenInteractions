

library(data.table)
path <- "/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/input//L1/UT/B.Exp.txt.gz"
readDoubleMatrix <- function(path){
  
  file <- ""
  if(endsWith(path, ".gz")){
    file <- gzfile(path)
  } else {
    file <- file(path)
  }
  
  table_tmp <- fread(path, sep = "\t", quote = "", data.table = F)
  table <- as.matrix(table_tmp[,-1])
  rownames(table) <- table_tmp[,1]
  return(table)
}


#listPerCelltypeMedianList <- list()
listPerCelltypeMeanList <- list()

#cellDataPath <- "/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_oneK1k/input/"
cellDataPath <- "/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3/wg3_multiome/input/"

l <- "L1"
treatment <- "UT"
ct <- "B"

for(l in c("L1", "L2")){
  
  for(treatment in c("UT", "24hCA")){
    
    cellTypeFolder <- file.path(cellDataPath, l, treatment)
    
    if(dir.exists(cellTypeFolder)){
    
      cellTypes <- gsub("\\.Exp\\.txt.*", "", list.files(file.path(cellDataPath, l, treatment), pattern = ".*\\.Exp\\.txt.*"))
  
      for(ct in cellTypes){
        
        cat(paste0(l,"-",ct,"-",treatment,"\n"))
        cat(" ",paste0(cellTypeFolder, "/", ct, ".qtlInput.Pcs.txt"),"\n")
        if(!file.exists(paste0(cellTypeFolder, "/", ct, ".qtlInput.Pcs.txt"))){
          cat("  skip\n")
          
        } else {
          
          if(file.exists(paste0(cellTypeFolder, "/", ct, ".Exp.txt"))){
            exp <- readDoubleMatrix(paste0(cellTypeFolder, "/", ct, ".Exp.txt"))
          } else {
            exp <- readDoubleMatrix(paste0(cellTypeFolder, "/", ct, ".Exp.txt.gz"))
          }
        
          
          pcs <- readDoubleMatrix(paste0(cellTypeFolder, "/", ct, ".qtlInput.Pcs.txt"))
          
          rownames(pcs)[!rownames(pcs) %in% colnames(exp)]
          
          
          selectedSamples <- rownames(pcs)[rownames(pcs) %in% colnames(exp)]
          #selectedSamples <- rownames(pcs)
          exp <- exp[,selectedSamples]
          
          if(treatment == "UT"){
            combinedName <- paste0(l,"-",ct)
          } else {
            combinedName <- paste0(l,"-",ct,"-",treatment)
          }
          
          
          #listPerCelltypeMedianList[[combinedName]] <- apply(exp, 1, median)
          listPerCelltypeMeanList[[combinedName]] <- apply(exp, 1, mean)
          
        }
        
        
      }
    }
  }
}


names(listPerCelltypeMeanList)

listPerCelltypeMeanList[["L1-Platelet"]] <- NULL
listPerCelltypeMeanList[["L1-other_T"]] <- NULL
listPerCelltypeMeanList[["L1-unannotated"]] <- NULL
listPerCelltypeMeanList[["L1-plasmablast"]] <- NULL

geneNames <- sapply(listPerCelltypeMeanList, names)



#sharedGenes <-  Reduce(intersect, geneNames)
allGenes <- unique(do.call(c, geneNames))

any(is.na(allGenes))

str(allGenes)

#listPerCelltypeMedian <- lapply(listPerCelltypeMedianList, function(g){g[allGenes]})
listPerCelltypeMean <- lapply(listPerCelltypeMeanList, function(g){g[allGenes]})


perCelltypeMean <- do.call(cbind, listPerCelltypeMean)

rownames(perCelltypeMean) <- allGenes

sum(is.na(perCelltypeMean))
perCelltypeMean[is.na(perCelltypeMean)] <- 0

any(is.na(rownames(perCelltypeMean)))

geneMapping <- read.delim("/groups/umcg-fg/tmp04/projects/downstreamer/depict2_bundle/reference_datasets/human_b37/genes_Ensembl94.txt")

rownames(perCelltypeMean)[!(rownames(perCelltypeMean) %in% geneMapping$Gene.stable.ID | rownames(perCelltypeMean) %in% geneMapping$Gene.name)]

updatedRownames <- rownames(perCelltypeMean)

updatedRownames[updatedRownames %in% geneMapping$Gene.name] <- geneMapping$Gene.stable.ID[match(updatedRownames[updatedRownames %in% geneMapping$Gene.name], geneMapping$Gene.name)]

rownames(perCelltypeMean) <- updatedRownames
perCelltypeMean <- perCelltypeMean[(rownames(perCelltypeMean) %in% geneMapping$Gene.stable.ID),]



hist(apply(perCelltypeMean, 1, mean), breaks = 100)


hist(apply(perCelltypeMean, 1, function(a){sum(a>0.05)}), breaks = length(listPerCelltypeMeanList))

str(perCelltypeMedian)


#perCelltypeMedian2 <- perCelltypeMedian[apply(perCelltypeMedian, 1, mean) >= 0.5,]

#perCelltypeMean2 <- perCelltypeMean[apply(perCelltypeMean, 1, mean) >= 0.05,]

perCelltypeMean2 <- perCelltypeMean[apply(perCelltypeMean, 1, mean) >= 0.2,]
dim(perCelltypeMean2)

hist(perCelltypeMean2[,2], breaks = 100)
str(perCelltypeMean2)
colnames(perCelltypeMean2)
#subset to genes scored by pascalX
perCelltypeMean2 <- perCelltypeMean2[(rownames(perCelltypeMean2) %in% geneMapping$Gene.stable.ID) ,]


library(heatmap3)
library(RColorBrewer)

heatmap3(cor(perCelltypeMean2), scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024))


perCelltypeMean2Mean <- apply(perCelltypeMean2, 1 ,mean)

perCelltypeMean2Centered <- apply(perCelltypeMean2, 2, function(x){x - perCelltypeMean2Mean})

heatmap3(cor(perCelltypeMean2Centered), scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024))

str(perCelltypeMean2Centered)


write.table(perCelltypeMean2, file = gzfile("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/multiome_meanPerCelltype.txt.gz"), sep = "\t", quote = F, col.names = NA)
write.table(perCelltypeMean2Centered, file = gzfile("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/oneK1K_meanPerCelltypeCentered.txt.gz"), sep = "\t", quote = F, col.names = NA)

hist(perCelltypeMean2Centered[,6], breaks = 100)
colnames(perCelltypeMean2Centered)

recount3medianCenter <- readDoubleMatrix("/groups/umcg-fg/tmp04/projects/genenetwork/recount3/CombinedHealthyTissue/combinedHealthyTissue_medianPerTissueCentered.txt.gz")
recount3median <- readDoubleMatrix("/groups/umcg-fg/tmp04/projects/genenetwork/recount3/CombinedHealthyTissue/combinedHealthyTissue_medianPerTissue.txt.gz")

rownames(recount3medianCenter) <- gsub("\\..*", "", rownames(recount3medianCenter))
rownames(recount3median) <- gsub("\\..*", "", rownames(recount3median))

all(rownames(recount3medianCenter) == rownames(recount3median))


sharedGenes <- intersect(rownames(recount3medianCenter), rownames(perCelltypeMean2))

heatmap3(cor(perCelltypeMean2[sharedGenes,]), scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024))

sCvsRecount3medianCenter <- cor(perCelltypeMean2Centered[sharedGenes,], recount3medianCenter[sharedGenes,])
heatmap3(sCvsRecount3medianCenter, scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024), margins = c(10,5))

sCvsRecount3median <- cor(perCelltypeMean2[sharedGenes,], recount3median[sharedGenes,])
heatmap3(sCvsRecount3median, scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024), margins = c(10,5))
max(sCvsRecount3median)


ibdGwas <- read.delim("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/pascalXRes/IBD.txt.gz")
rownames(ibdGwas) <- ibdGwas$gene
ibdGwas$zcore <- -qnorm(ibdGwas$pvalue/2)


sharedGenes <- intersect(rownames(ibdGwas), rownames(perCelltypeMean2))
plot(ibdGwas[sharedGenes,"zcore"], perCelltypeMean2[sharedGenes,"L1-NK-24hCA"])
cor.test(ibdGwas[sharedGenes,"zcore"], perCelltypeMean2[sharedGenes,"L1-NK-24hCA"])


plot(perCelltypeMean2[sharedGenes,"L1-NK-24hCA"])

