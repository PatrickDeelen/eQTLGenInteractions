library(heatmap3)
library(RColorBrewer)
palette(brewer.pal(n = 3, name = "Accent"))

tmpFolder <- "/local/3688655/"
dataFolder <- "/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene"
#

datasets <- data.frame(
  name = c("COVID19", "AIDA", "OneK1K", "SLE", "Aging", "Tabula", "TabulaSmartSeq2"),
  file = c("40ee3b22-395b-451b-883d-fb7e34ee612f.h5ad", "3a8a77a6-f069-479c-a2c1-252feafd0106.h5ad", "078a26dc-0585-4b94-9252-ee2f2dda5742.h5ad", "d51627ad-0123-4eb9-82b0-75f017862307.h5ad", "e11a76f0-57ff-4358-8a91-008655475059.h5ad", "", ""),
  assayFilter = c("10x 5' v2", NA, NA, NA, "10x 5' v2", "10x 3' v3", "Smart-seq2"),
  celltypeCol = c("celltype", "Annotation_Level4", "predicted.celltype.l2", "cell_type", "", "", "")
)

rownames(datasets) <- datasets$name

datasets <- datasets[datasets$name != "SLE",]
datasets <- datasets[datasets$name != "Aging",]


if(FALSE){
  
  library(anndata)
  #for(d in 1:nrow(datasets)){
  for(d in 2:3){
    
    name <- datasets[d, "name"]
    assayFilter <- datasets[d, "assayFilter"]
    celltypeCol <- datasets[d, "celltypeCol"]
    
    print(name)
    
    # Load the file
    datasetH5ad <- read_h5ad(file.path(tmpFolder, datasets[d, "file"]))
    
    
    
    metaData <- datasetH5ad$obs
    
    
     exp <- datasetH5ad$chunk_X(select = 0:(datasetH5ad$n_obs-1))
    
    colnames(exp) <- datasetH5ad$var_names
    rownames(exp) <- datasetH5ad$obs_names
    
    
    
    
    if(!is.na(assayFilter)){
      print(paste0("Subsetting to: ", assayFilter))
      metaData <- metaData[metaData$assay == assayFilter,]
    } 
    
    if(name == "SLE"){
      metaData$tissueCell_type <- as.character(metaData$cell_type)
      metaData$tissueCell_type[!is.na(metaData$ct_cov)] <- as.character(metaData$ct_cov)[!is.na(metaData$ct_cov)]
      
      metaData$disease2 <- as.character(metaData$disease)
      metaData$disease2[metaData$disease2 == "systemic lupus erythematosus"] <- "-SLE-cases"
      metaData$disease2[metaData$disease2 == "normal"] <- ""
      
      metaData$tissueCell_type <- as.factor(paste0(metaData$tissueCell_type, metaData$disease2))
      
    } else if(name == "Aging"){
      metaData$tissueCell_type <- as.character(metaData$cell_type)
      
      metaData$cmv2 <- as.character(metaData$cmv)
      metaData$cmv2[metaData$cmv2 == "negative"] <- "-CMVneg"
      metaData$cmv2[metaData$cmv2 == "positive"] <- "-CMVpos"
      
      metaData$tissueCell_type <- as.factor(paste0(metaData$tissueCell_type, metaData$cmv2))
      
    } else {
      metaData$tissueCell_type <- paste0(metaData$tissue, "-", metaData[,celltypeCol])
    }
    
    
    tctsCounts <- table(metaData$tissueCell_type)
    tctsCounts <- tctsCounts[tctsCounts >= 100]
    
    
    if(name == "Aging"){
      
      #Exlucde cell types with corrupt data that could not be loaded
      tctsCounts <- tctsCounts[names(tctsCounts) != "memory B cell-CMVneg"]
      tctsCounts <- tctsCounts[names(tctsCounts) != "memory B cell-CMVpos"]
      
    }
    
    tcts <- names(tctsCounts)
    
    
    expSum <- sapply(tcts, function(tct){
      print(tct)
      
      
      tct_cells <- rownames(metaData)[metaData$tissueCell_type == tct ]
      #tct_exp <- exp[tct_cells,]
      
      cellIndex <- which(datasetH5ad$obs_names %in% tct_cells)
      
      tct_exp <- datasetH5ad$chunk_X(select = cellIndex)
      colnames(tct_exp) <- datasetH5ad$var_names
      
      tct_expSum <- apply(tct_exp, 2, sum)
      if(name == "SLE"){
        return(tct_expSum)
      } else{
        return(log2(tct_expSum+1))
      }
    })
    
    saveRDS(expSum, file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctSum.rds")))
    saveRDS(tctsCounts, file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctCounts.rds")))
    
    expPseudobulk <- lapply(tcts, function(tct){
      print(tct)
      
      metaDataTct <- metaData[metaData$tissueCell_type == tct,]
      
      donorCounts <- table(metaDataTct$donor_id)
      donorCounts <- donorCounts[donorCounts >= 100]
      tctDonors <- names(donorCounts)
      
      
      tct_cells <- rownames(metaData)[metaData$tissueCell_type == tct ]
      
      cellIndex <- which(datasetH5ad$obs_names %in% tct_cells)
      
     #tct_exp <- datasetH5ad$chunk_X(select = cellIndex)
      #colnames(tct_exp) <- datasetH5ad$var_names
     # rownames(tct_exp) <- datasetH5ad$obs_names[cellIndex]
      
      tctPseudobulk <- sapply(tctDonors, function(donor){
        
        donorCells <- rownames(metaDataTct)[metaDataTct$donor_id == donor]
        
        donorExp <- exp[donorCells,]
        
        donor_expSum <- apply(donorExp, 2, sum)
        
        return(log2(donor_expSum+1))
        
      })
      
      return(tctPseudobulk)
      
    })
    
    
  }
}

perDataSetTctsCounts <- lapply(datasets$name, function(name){
  return(readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctCounts.rds"))))
})

perDataSetExpSum <- lapply(datasets$name, function(name){
  return(readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctSum.rds"))))
})

names(perDataSetExpSum) <- datasets$name
names(perDataSetTctsCounts) <- datasets$name

if(FALSE){
  library(broman)
  
  
  # name <- datasets$name[5]
  
  sink <- lapply(datasets$name, function(name){
    
    cat(name, "\n")
    expSum <- perDataSetExpSum[[name]]
    tctsCounts <- perDataSetTctsCounts[[name]]
    
    
    zeroPerGene <- apply(expSum, 1,function(a){sum(a == 0) / length(a) * 100})
    meanPerGene <- apply(expSum, 1,mean)
    #hist(zeroPerGene, breaks = 100)
    #hist(meanPerGene, breaks = 100)
    #plot(zeroPerGene, meanPerGene)
    expSum2 <- expSum[zeroPerGene < 80 & meanPerGene > 1, ]
    dim(expSum2)
    
    
    zeroPerTct2 <- apply(expSum2, 2,function(a){sum(a == 0) / length(a) * 100})
    
    expSum3 <- expSum2[, zeroPerTct2 <= 25]
    
    dim(expSum3)
    
    meanPerGene3 <- apply(expSum3, 1,mean)
    zeroPerGene3 <- apply(expSum3, 1,function(a){sum(a == 0) / length(a) * 100})
    #hist(meanPerGene3, breaks = 100)
    #hist(zeroPerGene3, breaks = 100)
    
    expSum4 <- expSum3[meanPerGene3 >= 2 & zeroPerGene3 < 10, ]
    dim(expSum4)
    
    meanPerGene4 <- apply(expSum4, 1,mean)
    #hist(meanPerGene4, breaks = 100)
    
    meanPerTct4 <- apply(expSum4, 2,mean)
    sdPerTct4 <- apply(expSum4, 2,sd)
    zeroPerTct4 <- apply(expSum4, 2,function(a){sum(a == 0) / length(a) * 100})
    
    #hist(zeroPerTct4, breaks = 100)
    
    #plot(meanPerTct4, as.numeric(tctsCounts[names(meanPerTct4)]), ylab = "Number of cells per celltype", log =  "y", xlab = "Mean log2 exp")
    #plot(sdPerTct4, as.numeric(tctsCounts[names(meanPerTct4)]), ylab = "Number of cells per celltype", log =  "y", xlab = "SD log2 exp")
    #plot(zeroPerTct4, as.numeric(tctsCounts[names(meanPerTct4)]), ylab = "Number of cells per celltype", log =  "y", xlab = "Percentage of genes with zero expression")
    
    expSum5 <- expSum4[, zeroPerTct4 < 10]
    dim(expSum5)
    
    meanPerTct5 <- apply(expSum5, 2,mean)
    sdPerTct5 <- apply(expSum5, 2,sd)
    zeroPerTct5 <- apply(expSum5, 2,function(a){sum(a == 0) / length(a) * 100})
    
    #plot(meanPerTct5, as.numeric(tctsCounts[names(meanPerTct5)]), ylab = "Number of cells per celltype", log =  "y", xlab = "Mean log2 exp")
    #plot(sdPerTct5, as.numeric(tctsCounts[names(sdPerTct5)]), ylab = "Number of cells per celltype", log =  "y", xlab = "SD log2 exp")
    #plot(zeroPerTct5, as.numeric(tctsCounts[names(zeroPerTct5)]), ylab = "Number of cells per celltype", log =  "y", xlab = "Percentage of genes with zero expression")
    
    dim(expSum5)
    
    str(expSum5)
    
    expQcQq <- normalize(expSum5)
    dimnames(expQcQq) <- dimnames(expSum5)
    
    saveRDS(expQcQq, file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctSumQcQq.rds")))
    return()
    
    
    
  })
}

perDataSetExpSumQcQq <- lapply(datasets$name, function(name){
  return(readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctSumQcQq.rds"))))
})

names(perDataSetExpSumQcQq) <- datasets$name

layout(matrix(1:4, nrow  =2 ))
hist(perDataSetExpSumQcQq[["Tabula"]][,1])
hist(perDataSetExpSumQcQq[["Tabula"]][,2])
hist(perDataSetExpSumQcQq[["Tabula"]][,3])
hist(perDataSetExpSumQcQq[["Tabula"]][,4])

layout(1)
boxplot(perDataSetExpSumQcQq[["Tabula"]][,1:10])

sink <- lapply(datasets$name, function(name){
  pdf(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene/log2QcQq", name, ".pdf"), width = 90, height = 90)
  heatmap3(cor(perDataSetExpSumQcQq[[name]]), scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024), margins = c(25,25))
  dev.off()
  return()
})

sink <- lapply(datasets$name, function(name){
  
  
  geneMean <- apply(perDataSetExpSumQcQq[[name]], 1 ,mean)
  
  expSumQcQqCentered <- apply(perDataSetExpSumQcQq[[name]], 2, function(x){x - geneMean})
  
  
  saveRDS(expSumQcQqCentered, file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctSumQcQqMeanCentered.rds")))
  
  write.table(expSumQcQqCentered, file = gzfile(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene/processing/",name,"MeanCentered.txt.gz")), quote = F, sep ="\t", col.names = NA)
  
  tctDendro <- as.dendrogram(hclust(as.dist(1 - cor(expSumQcQqCentered)), method = "ward.D2"))
  
  saveRDS(tctDendro, file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctDendro.rds")))
  
  pdf(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctDendro.pdf")), width = 100)
  plot(tctDendro)
  dev.off()
  
})

datasetsToAlignToTabula <- c("COVID19", "AIDA", "OneK1K")
genesPerDataset <- lapply(c("Tabula", datasetsToAlignToTabula), function(name){
  return(rownames(perDataSetExpSumQcQq[[name]]))
})
genesPerDatasetCount <- table(do.call(c,genesPerDataset))
sharedGenes <- names(genesPerDatasetCount)[genesPerDatasetCount == length(datasetsToAlignToTabula)+1]
length(sharedGenes)


tabulaDist <- perDataSetExpSumQcQq[["Tabula"]][sharedGenes,1]
tabulaDistSort <- sort(tabulaDist)
tabulaGeneMean <- apply(perDataSetExpSumQcQq[["Tabula"]][sharedGenes,], 1 ,mean)


perDataSetExpSumQcQqAlignedToTabula <- lapply(datasetsToAlignToTabula, function(name){
  
  data <- perDataSetExpSumQcQq[[name]][sharedGenes,]
  
  
  data2 <- apply(data, 2, function(exp){
    return(tabulaDistSort[rank(exp)])
  })
  
  rownames(data2) <- rownames(data)
  return(data2)
  
  
})
names(perDataSetExpSumQcQqAlignedToTabula) <- datasetsToAlignToTabula

boxplot(cbind(perDataSetExpSumQcQq[["Tabula"]][sharedGenes,1:10], perDataSetExpSumQcQqAlignedToTabula[["COVID19"]][,1:10]))

sink <- lapply(datasetsToAlignToTabula, function(name){
  
  expSumQcQqCentered <- apply(perDataSetExpSumQcQqAlignedToTabula[[name]], 2, function(x){x - tabulaGeneMean})
  
  write.table(expSumQcQqCentered, file = gzfile(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene/processing/",name,"MeanCenteredTabula.txt.gz")), quote = F, sep ="\t", col.names = NA)
  
})



datasetsToMerge <- c("COVID19",         "AIDA",            "OneK1K"     ,           "Tabula" )

genesPerDataset <- lapply(datasetsToMerge, function(name){
  return(rownames(perDataSetExpSumQcQq[[name]]))
})
genesPerDatasetCount <- table(do.call(c,genesPerDataset))
sharedGenes <- names(genesPerDatasetCount)[genesPerDatasetCount == length(datasetsToMerge)]


subsetPerDataset <- lapply(datasetsToMerge, function(name){
  data <- perDataSetExpSumQcQq[[name]][sharedGenes,]
  
  colnames(data) <- paste0(colnames(data), "_(",name,")")
  
  return(data)
})

combinedExpSumQcQq <- do.call(cbind, subsetPerDataset)
library(broman)
combinedExpSumQcQq2 <- normalize(combinedExpSumQcQq)
rownames(combinedExpSumQcQq2) <- rownames(combinedExpSumQcQq)
colnames(combinedExpSumQcQq2) <- colnames(combinedExpSumQcQq)

geneMean <- apply(combinedExpSumQcQq2, 1 ,mean)

expSumQcQqCentered <- apply(combinedExpSumQcQq2, 2, function(x){x - geneMean})
study <- as.factor(gsub("\\)", "", gsub(".+\\(", "", colnames(expSumQcQqCentered))))
levels(study) <- brewer.pal(n = length(levels(study)), name = "Accent")
study <- as.character(study)
studyAnn <- rep("steelblue2", times = length(study))



pdf(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene/log2QcQq", "Combined", ".pdf"), width = 90, height = 90)
heatmap3(cor(expSumQcQqCentered), scale = "none", col = colorRampPalette(brewer.pal(n = 11, name = "Spectral"), bias = 2)(1024), margins = c(25,25), ColSideColors = study)
dev.off()


write.table(expSumQcQqCentered, file = gzfile(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene/processing/","Combined","MeanCentered.txt.gz")), quote = F, sep ="\t", col.names = NA)

tctDendro <- as.dendrogram(hclust(as.dist(1 - cor(expSumQcQqCentered)), method = "ward.D2"))
plot(tctDendro)

saveRDS(tctDendro, file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0("Combined", "_tctDendro.rds")))

