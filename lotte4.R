library(heatmap3)
library(RColorBrewer)
library(dendextend)
library(readxl)
library(ComplexHeatmap)
palette(brewer.pal(n = 8, name = "Accent"))


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





read.depict2 <- function(path, potential_traits=NULL) {
  if (is.null(potential_traits)) {
    potential_traits <- excel_sheets(path)
    potential_traits <- potential_traits[grep("Overview", potential_traits, invert=T)]
  }
  
  output <- list()
  for (sheet in potential_traits) {
    tmp <- tryCatch({data.frame(read_excel(path, sheet=sheet, col_types ="guess", trim_ws = T), stringsAsFactors=F)},
                    error=function(a){return(NA)},
                    warn=function(a){return(NA)})
    
    
    for (i in 1:ncol(tmp)) {
      if (class(tmp[,i]) == "character"){
        tmp[,i] <- type.convert(tmp[,i], as.is=T)
        
      }
    }
    
    rownames(tmp) <- tmp[,1]
    output[[sheet]] <- tmp
  }
  
  return(output)
}



traitDendro <- readRDS(file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/traitDendo.rds"))


traits <- read.delim("/groups/umcg-fg/tmp04/projects/downstreamer/depict2_bundle/scripts/trait_files/lotte.txt", header =F)$V1


perDataSetTctDendro <- lapply(datasets$name, function(name){
  return(readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0(name, "_tctDendro.rds"))))
})

names(perDataSetTctDendro) <- datasets$name

resultFolder <- "/groups/umcg-fg/tmp04/projects/downstreamer/lotte/results13/"

perTraitDsRes <- lapply(traits, function(trait){
  return(read.depict2(paste0(resultFolder, trait,"/",trait,"_TissueCelltypeEnrichment_enrichtments.xlsx")))
})
names(perTraitDsRes) <- traits

name <- "COVID19"
trait <- "SLE"

maxPerTrait <- 10
sigCol <- "Bonferroni.significant" #"FDR.5..significant"   "Bonferroni.significant"

library(circlize)

name <- "COVID19"
name <- "AIDA"
name <- "OneK1K"
name <- "Tabula"



for(name in c("COVID19", "AIDA", "OneK1K", "Tabula")){
  
  print(name)
  
  significantPerTrait <- lapply(traits, function(trait){
    datasetTraitRes <- perTraitDsRes[[trait]][[name]]
    
    datasetTraitResSig <- datasetTraitRes[datasetTraitRes$Enrichment.Z.score > 0 & datasetTraitRes[,sigCol],"Gene.set"]
    
    toReturn <- ifelse(length(datasetTraitResSig) > maxPerTrait, maxPerTrait, length(datasetTraitResSig))
    
    return(datasetTraitResSig[1:toReturn])
    
  })
  
  significantAnyTrait <- unique(do.call(c, significantPerTrait[!is.na(significantPerTrait)]))
  
  if(length(significantAnyTrait) <2){
    print(" - non signficant")
    return()
  }
  
  zscoresPerTrait <- sapply(traits, function(trait){
    datasetTraitRes <- perTraitDsRes[[trait]][[name]]
    
    return(datasetTraitRes[match(significantAnyTrait, datasetTraitRes$Gene.set),"Enrichment.Z.score"])
    
  })
  
  rownames(zscoresPerTrait) <- significantAnyTrait
  
  
  if(name == "Combined"){
    tctDendro <- readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0("Combined", "_tctDendro.rds")))
  } else {
    tctDendro <- perDataSetTctDendro[[name]]
  }
  
  
  toPrune <- labels(tctDendro)[!labels(tctDendro) %in% significantAnyTrait]
  tctDendroSig <- prune(tctDendro, toPrune)
  
  numberOfCells <- nrow(zscoresPerTrait)
  
  heigthRatio <-  nrow(zscoresPerTrait)  / ncol(zscoresPerTrait)
  col_fun = colorRamp2(c(0,1.5,3,3.5,4,4.5,5,5.5,6,8), c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#800026"))
  
  ht <-  Heatmap(zscoresPerTrait[labels(tctDendroSig),], width = unit(5, "cm"), height = unit(5 * heigthRatio, "cm"), col = col_fun, cluster_columns =  traitDendro, cluster_rows  = FALSE , rect_gp = gpar(col = "white", lwd = 1), row_names_max_width = max_text_width(
    labels(tctDendroSig), 
    gp = gpar(fontsize = 12)), row_names_side = "left")
  
  pdf(file.path(resultFolder, paste0(name, "_top2.pdf")), width = 20, height = 20)
  draw(ht)
  dev.off()
  
  
  # pdf(file.path(resultFolder, paste0(name, "_top.pdf")), width = 20, height = 20)
  #  heatmap3(zscoresPerTrait[labels(tctDendroSig),], scale = "none", balanceColor = F, margins = c(5,25), Rowv = NA, Colv = traitDendro, cexRow = 0.7, col = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"), bias = 1)(1024))
  # dev.off()
  
  
  
}




for(name in c("COVID19", "AIDA", "OneK1K", "Tabula")){
  
  print(name)

  
  celltypes <- perTraitDsRes[[1]][[name]][,"Gene.set"]
  
  zscoresPerTrait <- sapply(traits, function(trait){
    datasetTraitRes <- perTraitDsRes[[trait]][[name]]
    
    return(datasetTraitRes[match(celltypes, datasetTraitRes$Gene.set),"Enrichment.Z.score"])
    
  })
  
  rownames(zscoresPerTrait) <- celltypes
  
  
  if(name == "Combined"){
    tctDendro <- readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0("Combined", "_tctDendro.rds")))
  } else {
    tctDendro <- perDataSetTctDendro[[name]]
  }
  
  
  toPrune <- labels(tctDendro)[!labels(tctDendro) %in% celltypes]
  tctDendroSig <- prune(tctDendro, toPrune)
  
  numberOfCells <- nrow(zscoresPerTrait)
  
  heigthRatio <-  nrow(zscoresPerTrait)  / ncol(zscoresPerTrait)
  col_fun = colorRamp2(c(0,1.5,3,3.5,4,4.5,5,5.5,6,8), c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C", "#BD0026", "#800026", "#800026"))
  
  colnames(zscoresPerTrait)[colnames(zscoresPerTrait) == "SLE3"] <- "SLE"
  labels(traitDendro)[labels(traitDendro) == "SLE3"] <- "SLE"
  
  ht <-  Heatmap(zscoresPerTrait[labels(tctDendroSig),], width = unit(5, "cm"), height = unit(5 * heigthRatio, "cm"), col = col_fun, cluster_columns =  traitDendro, cluster_rows  = FALSE , rect_gp = gpar(col = "white", lwd = 1), row_names_max_width = max_text_width(
    labels(tctDendroSig), 
    gp = gpar(fontsize = 12)), row_names_side = "left",
  heatmap_legend_param = list(title = "Z-score"))
  
  pdf(file.path(resultFolder, paste0(name, "_all.pdf")), width = 10, height = (nrow(zscoresPerTrait) * 0.2) + 2)
  draw(ht)
  dev.off()
  
  write.table(zscoresPerTrait[labels(tctDendroSig), labels(traitDendro)], file = file.path(resultFolder, paste0(name, "_all.txt")), sep = "\t", quote = F, col.names = NA)
  
  print(-qnorm((0.05/nrow(zscoresPerTrait)) /2))
  
}


