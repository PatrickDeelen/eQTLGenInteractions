setwd("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/eqltVariantAnalysis")

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)

palette(brewer.pal(n = 8, name = "Accent"))

resultFolder <- "/groups/umcg-fg/tmp04/projects/downstreamer/lotte/results13/"

#  See Zhang et al, Nature Genetics volume 48, pages481–487(2016) for details
#Formula taken from Supplementary Information
zToBeta <- function(z, maf, n) {
  
  chi = z * z;
  a = 2 * maf * (1 - maf) * (n + chi);
  beta = z / Math.sqrt(a);
  se = 1 / Math.sqrt(a);
  return(c(beta, se));
}



source("Rb.R")
traitDendro <- readRDS(file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/traitDendo.rds"))

load("wip7.RData")
rm(c)
eqtls$gene_name2 <- eqtls$gene_name
eqtls$gene_name2[eqtls$gene_name == ""] <- eqtls$phenotype[eqtls$gene_name == ""]


celltypes <- c("B", "CD4_T", "CD8_T", "DC", "Mono", "NK")

celltype <- celltypes[1]

celltypesQtls <- lapply(celltypes, function(celltype){
  
  cellTypeQtls <- read.delim(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/rb/Meta_v5/UMCG_v2-UMCG_v3-oneK1k-Multiome-Wijst.",celltype , ".incPerm.tsv"))  
  cellTypeQtls <- cellTypeQtls[!is.na(cellTypeQtls$I_square) & cellTypeQtls$I_square <= 40, ]
  cellTypeQtls <- cellTypeQtls[!is.na(cellTypeQtls$n_datasets) & cellTypeQtls$n_datasets >= 3, ]
  
  rownames(cellTypeQtls) <- cellTypeQtls$alt_QTL
  
  cellTypeQtls <- cellTypeQtls[rownames(cellTypeQtls) %in% rownames(eqtls),]
  
  
  cellTypeQtls$assessed_alleleEqtlgen <- eqtls$eff_allele[match(rownames(cellTypeQtls), rownames(eqtls))]
  
  cellTypeQtls$betaCorrected <- cellTypeQtls$beta
  cellTypeQtls$betaCorrected[cellTypeQtls$assessed_alleleEqtlgen != cellTypeQtls$assessed_allele] <- cellTypeQtls$beta[cellTypeQtls$assessed_alleleEqtlgen != cellTypeQtls$assessed_allele] * -1
  
  cellTypeQtls$zscoreCorrected <- cellTypeQtls$recalibrated_z_score
  cellTypeQtls$zscoreCorrected[cellTypeQtls$assessed_alleleEqtlgen != cellTypeQtls$assessed_allele] <- cellTypeQtls$recalibrated_z_score[cellTypeQtls$assessed_alleleEqtlgen != cellTypeQtls$assessed_allele] * -1
  
  
  return(cellTypeQtls)
  
})
names(celltypesQtls) <- celltypes

sharedEqtls <- lapply(celltypes, function(celltype){
  rownames(celltypesQtls[[celltype]])
})
sharedEqtls <- do.call(c,sharedEqtls)
sharedEqtls <- table(sharedEqtls)
sharedEqtls <- names(sharedEqtls)[sharedEqtls == length(celltypes)]
length(sharedEqtls)

sharedEqtlsZscores <- sapply(celltypes, function(celltype){
  celltypesQtls[[celltype]][sharedEqtls,"zscoreCorrected"]
})

singleCellEqtlsDendro <- as.dendrogram(hclust(as.dist(1 - cor(sharedEqtlsZscores)), method = "ward.D2"))
plot(singleCellEqtlsDendro)

traits <- read.delim("/groups/umcg-fg/tmp04/projects/downstreamer/depict2_bundle/scripts/trait_files/lotte.txt", header =F)$V1

traitsMrFile <- c("asth", "cd", "ibd", "ms", "ra", "sle", "t1dm")
names(traitsMrFile) <- c("ASTH", "CD",  "IBD",  "MS",   "RA",   "SLE3", "T1DM")


traitsMr <- traits[traits %in% names(traitsMrFile)]

trait <- traitsMr[6]
trait <- "SLE3"
traitsPerCelltypeRb <- lapply(traitsMr, function(trait){
  print(trait)
  
  traitMrRes <- read.delim(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/rb/snp_gene_",traitsMrFile[trait] , ".tsv"), sep = "\t")  

  traitMrRes$qtl <- paste0(traitMrRes$variant,"-", traitMrRes$phenotype)
  
  traitMrRes <- traitMrRes[match(unique(traitMrRes$qtl),traitMrRes$qtl),]
  
  rownames(traitMrRes) <- traitMrRes$qtl
  
  layout(matrix(1:6, nrow = 1))
  
  celltype <- "NK"
  rbPerCelltype <- t(sapply(celltypes, function(celltype){
    
    sharedQtls <- intersect(rownames(traitMrRes), rownames(celltypesQtls[[celltype]]))
    
    celltypeQtls2 <- celltypesQtls[[celltype]][sharedQtls,]
    traitMrRes2 <- traitMrRes[sharedQtls,]
    
    rbAll <- calcu_cor_true(traitMrRes2$beta_eqtl, traitMrRes2$standard_error, celltypeQtls2$betaCorrected, celltypeQtls2$beta_se, 0)
    
    print(paste0(celltype, " - ", sum(traitMrRes2$mr_gene)))
    
    
    traitMrRes2$qtl[traitMrRes2$mr_gene][traitMrRes2$beta_eqtl[traitMrRes2$mr_gene] > 0]
    
    if(any(traitMrRes2$mr_gene)){
      varFactor <- as.factor(traitMrRes$variant[traitMrRes2$mr_gene])
      plot(traitMrRes2$beta_eqtl[traitMrRes2$mr_gene], celltypeQtls2$betaCorrected[traitMrRes2$mr_gene], xlab = "eQtlgen beta", ylab = paste0(celltype," beta"), main = trait, pch = 16, col = 1)
      if(sum(traitMrRes2$mr_gene) > 5){
        rbMr <- calcu_cor_true(traitMrRes2$beta_eqtl[traitMrRes2$mr_gene], traitMrRes2$standard_error[traitMrRes2$mr_gene], celltypeQtls2$betaCorrected[traitMrRes2$mr_gene], celltypeQtls2$beta_se[traitMrRes2$mr_gene], 0)
        mtext(paste0("Rb: ", rbMr[1]))
      } else {
        rbMr <- c(0,0,1)
        mtext(paste0("Rb: NA"))
      }
      
      
      
      
    } else {
      plot.new()
      rbMr <- c(0,0,1)
    }
    
    
    
    return(c(rbAll, rbMr))
    
  }))
  colnames(rbPerCelltype) <- c("rbAll", "rbSeAll", "rbPAll", "rbMr", "rbSeMr", "rbPMr")
  rbPerCelltype
  
  return(rbPerCelltype)
  
})
names(traitsPerCelltypeRb) <- traitsMr

traitsPerCelltypeRb

trait <- "SLE3"
trait <- "IBD"
traitReplicatingEqtls <- lapply(traitsMr, function(trait){
  print(trait)
  
  traitMrRes <- read.delim(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/rb/snp_gene_",traitsMrFile[trait] , ".tsv"), sep = "\t")  
  
  traitMrRes$qtl <- paste0(traitMrRes$variant,"-", traitMrRes$phenotype)
  
  traitMrRes <- traitMrRes[match(unique(traitMrRes$qtl),traitMrRes$qtl),]
  
  rownames(traitMrRes) <- traitMrRes$qtl
  
  layout(matrix(1:6, nrow = 1))
  
  celltype <- "CD8_T"
  traitReplicatingEqtlsPerCelltype <- lapply(celltypes, function(celltype){
    
    sharedQtls <- intersect(rownames(traitMrRes), rownames(celltypesQtls[[celltype]]))
    
    celltypeQtls2 <- celltypesQtls[[celltype]][sharedQtls,]
    traitMrRes2 <- traitMrRes[sharedQtls,]
    
    celltypeQtls2 <- celltypeQtls2[traitMrRes2$mr_gene,]
    traitMrRes2 <- traitMrRes2[traitMrRes2$mr_gene,]
    
    celltypeQtls2$SigReplication <- celltypeQtls2$recalibrated_p_value <= ( 0.05/nrow(celltypeQtls2))
    
    if(any(celltypeQtls2$SigReplication)){
      
      celltypeQtls3 <- celltypeQtls2[ celltypeQtls2$SigReplication ,]
      traitMrRes3 <- traitMrRes2[celltypeQtls2$SigReplication ,]
      
      diseaseDirectionCorrection <- ifelse(traitMrRes3$beta_gwas * traitMrRes3$r_tagging < 0, -1, 1) 
      
      celltypeQtls3$betaCorrectedDiseaseDirection <- celltypeQtls3$betaCorrected * diseaseDirectionCorrection 
      
   
      varFactor <- as.factor(celltypeQtls3$snp_id)
      
      plot(traitMrRes3$beta_eqtl, celltypeQtls3$betaCorrectedDiseaseDirection, xlab = "eQtlgen beta", ylab = paste0(celltype," beta"), main = trait, pch = 16, col = as.numeric(varFactor) )
      
      
      
      return(traitMrRes3$qtl)
      
    } else {
      
      plot.new()
      return()
    }
    
    zscoresCellTypes <- lapply(celltypes, function(celltype){
      
    })
    
  })
  names(traitReplicatingEqtlsPerCelltype) <- celltypes
  
  eqtlsReplicating <- do.call(c, traitReplicatingEqtlsPerCelltype)
  eqtlsReplicating <- unique(eqtlsReplicating)
  
  celltype <- "CD8_T"
  replicatingEqtlsZscores <-  sapply(celltypes, function(celltype){
    
    celltypeQtls2 <- celltypesQtls[[celltype]][eqtlsReplicating,]
    traitMrRes2 <- traitMrRes[eqtlsReplicating,]
    
    diseaseDirectionCorrection <- ifelse(traitMrRes2$beta_gwas * traitMrRes2$r_tagging < 0, -1, 1) 
    
    celltypeQtls2$zscoreCorrectedDiseaseDirection <- celltypeQtls2$zscoreCorrected * diseaseDirectionCorrection
  
    
    return(celltypeQtls2$zscoreCorrectedDiseaseDirection)
    
  })
  
  if(is.null(dim(replicatingEqtlsZscores))){
    replicatingEqtlsZscores <- t(as.matrix(replicatingEqtlsZscores))
  }
  
  
  eqtls2 <- eqtls[eqtlsReplicating,]
  traitMrRes2 <- traitMrRes[eqtlsReplicating,]
  
  diseaseDirectionCorrection <- ifelse(traitMrRes2$beta_gwas * traitMrRes2$r_tagging < 0, -1, 1) 

  eqtls2$tstatCorrectedDiseaseDirection <- eqtls2$tstat * diseaseDirectionCorrection
  
  replicatingEqtlsZscores <- cbind( eqtls2$tstatCorrectedDiseaseDirection, replicatingEqtlsZscores)
  
  rownames(replicatingEqtlsZscores) <- paste0(eqtls2$variant, " -> ", eqtls2$gene_name2)
  colnames(replicatingEqtlsZscores)[1] <- "Blood"

  
  
  maxAbsZ <- 5
  col_fun = colorRamp2(c(seq(-maxAbsZ,-1, length.out = 10), 0 , seq(1,maxAbsZ, length.out = 10)), c(rev(brewer.pal(n = 9, name = "Blues")), "#FFFFFF", "#FFFFFF", "#FFFFFF", brewer.pal(n = 9, name = "YlOrRd")))
  
  heigthRatio <-  nrow(replicatingEqtlsZscores)  / ncol(replicatingEqtlsZscores)
  
  ht <- Heatmap(replicatingEqtlsZscores[,c(labels(singleCellEqtlsDendro), "Blood"), drop = F], col = col_fun, na_col = "White", cluster_columns =  FALSE, cluster_rows = FALSE,
          width = unit(5, "cm"), height = unit(5 * heigthRatio, "cm"), rect_gp = gpar(col = "grey", lwd = 1), heatmap_legend_param = list(title = "Z-score"), column_title = trait)
  
  pdf(file.path(resultFolder, paste0(trait, "_MrCelltypeReplication.pdf")), width = 20, height = 20)
  draw(ht)
  dev.off()

  
  traitMrRes2[, c("beta_eqtl", "beta_gwas", "eff_allele_eqtl", "eff_allele_gwas", "variant", "variant_gwas", "r_tagging")]
  eqtls2[, c("beta", "tstat", "tstatCorrectedDiseaseDirection", "eff_allele", "variant")]
  
  dim(eqtls2)
  
  return(traitReplicatingEqtlsPerCelltype)
  
})

names(traitReplicatingEqtls) <- traitsMr



traitsPerCelltypeRb

traitsPerCelltypeRbZ <- sapply(traitsPerCelltypeRb, function(rbPerCelltype){
  z <- qnorm(rbPerCelltype[,3] /2)
  
  z[rbPerCelltype[,1] > 0] <- z[rbPerCelltype[,1] > 0] * -1
  return(z)
})

maxAbsZ <- max(abs(traitsPerCelltypeRbZ))
maxAbsZ <- 5
col_fun = colorRamp2(c(seq(-maxAbsZ,-1, length.out = 10), 0 , seq(1,maxAbsZ, length.out = 10)), c(rev(brewer.pal(n = 9, name = "Blues")), "#FFFFFF", "#FFFFFF", "#FFFFFF", brewer.pal(n = 9, name = "YlOrRd")))

Heatmap(traitsPerCelltypeRbZ, col= col_fun )


traitsPerCelltypeMrRbZ <- sapply(traitsPerCelltypeRb, function(rbPerCelltype){
  z <- qnorm(rbPerCelltype[,6] /2)
  
  z[rbPerCelltype[,4] > 0] <- z[rbPerCelltype[,4] > 0] * -1
  return(z)
})

maxAbsZ <- max(abs(traitsPerCelltypeMrRbZ))
maxAbsZ <- 6

toPrune <- labels(traitDendro)[!labels(traitDendro) %in% traitsMr]
traitDendro2 <- prune(traitDendro, toPrune)
plot(traitDendro2)

heigthRatio <-  nrow(traitsPerCelltypeMrRbZ)  / ncol(traitsPerCelltypeMrRbZ)



col_fun = colorRamp2(c( 0 , seq(1,maxAbsZ, length.out = 10)), c( "#FFFFFF", "#FFFFFF", brewer.pal(n = 9, name = "YlGn")))


ht <- Heatmap(traitsPerCelltypeMrRbZ[labels(singleCellEqtlsDendro),labels(traitDendro2)], col= col_fun, cluster_columns =  FALSE, cluster_rows = FALSE,
        width = unit(5, "cm"), height = unit(5 * heigthRatio, "cm"), rect_gp = gpar(col = "grey", lwd = 1))

pdf(file.path(resultFolder, paste0("celltypeRbReplication.pdf")), width = 20, height = 20)
draw(ht)
dev.off()


str(traitMrRes2)
celltypeQtls2

test <- eqtls[rownames(traitMrRes2),]


layout(1)
plot(test$scMeta_NK_beta, celltypeQtls2$beta)

plot(test$beta, traitMrRes2$beta_eqtl)
table(test$eff_allele, traitMrRes2$al)
plot(test$scMeta_NK_beta, celltypeQtls2$betaCorrected)




expSumQcQqCentered <- readRDS(file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0("COVID19", "_tctSumQcQqMeanCentered.rds")))

sdPerGene <- apply(expSumQcQqCentered, 1, sd)


expSumQcQqCenteredScaled <- expSumQcQqCentered / sdPerGene



traitsPerCelltypeRb <- lapply(traitsMr, function(trait){
  print(trait)
  
  traitMrRes <- read.delim(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/rb/snp_gene_",traitsMrFile[trait] , ".tsv"), sep = "\t")  
  
  
  mrGenes <- unique(traitMrRes[traitMrRes$mr_gene,"phenotype" ])
  return(mrGenes)
  
})

names(traitsPerCelltypeRb) <- traitsMr

str(expSumQcQqCentered[traitsPerCelltypeRb[["SLE3"]],])

traitDendro <- readRDS(file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/traitDendo.rds"))

tctDendro <- readRDS(file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/singleCellData/cellxgene", paste0("COVID19", "_tctDendro.rds")))

Heatmap(expSumQcQqCenteredScaled[traitsPerCelltypeRb[["SLE3"]],][,labels(tctDendro)], cluster_columns  = FALSE)


apply(expSumQcQqCentered[traitsPerCelltypeRb[["SLE3"]],][,labels(tctDendro)], 1, sd)

