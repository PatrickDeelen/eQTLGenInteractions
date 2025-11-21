setwd("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/eqltVariantAnalysis")
palette(brewer.pal(n = 3, name = "Accent"))


#load("wip6.RData")
ldPerChr <- readRDS("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/fix_20250509/ld_matrices_20250509.rds")
coexp <- readRDS("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/all_permuted_independentVariants_4GenPC20ExpPC_2025-02-05/exported_matrix/coexpression_matrix_16781_genes.rds")

str(eqtls)
eqtlsCisPerChr <- lapply(1:22, function(chr){
  as.character(eqtls[eqtls$type != "trans" & eqtls$chromosome == chr,"variant_index"])
})
names(eqtlsCisPerChr) <- 1:22


transWithHintOfCis <- read.delim("../transNoCis/transSnpsWithHintOfCisEffect.txt", header = F)$V1
str(transWithHintOfCis)

table(eqtls$type)
transSnp <- eqtls$variant[eqtls$type == "trans"]

table(transSnp %in% transWithHintOfCis)


transNoColoc <- transSnp[!transSnp %in% transWithHintOfCis]
str(transNoColoc)

eqtlsTransNoCisColoc <- eqtls[eqtls$type == "trans" & (eqtls$variant %in% transNoColoc),]

transNoLd <- sapply(1:nrow(eqtlsTransNoCisColoc), function(e){
  transVarIndex <- as.character(eqtlsTransNoCisColoc[e,"variant_index"])
  chr <- eqtlsTransNoCisColoc[e, "chromosome"]
  cisEqtls <- eqtlsCisPerChr[[chr]]
  
  ld <- ldPerChr[[chr]][transVarIndex,cisEqtls] ^ 2

  return(!any(ld >= 0.1))
  
})

table(transNoLd)

eqtlsTransNoCis <- eqtlsTransNoCisColoc[transNoLd,]

dim(eqtlsTransNoCis)
sort(table(eqtlsTransNoCis$consequence))

selection <- eqtlsTransNoCis[!is.na(eqtlsTransNoCis$consequence) & (eqtlsTransNoCis$consequence == "missense_variant" |  eqtlsTransNoCis$consequence == "inframe_insertion") , c("variant", "phenotype", "consequence", "consequenceGene")]


pickAnyGene2 <- apply(selection, 1, function(eqtl){
  
  gene <- eqtl["phenotype"]
  variant <- eqtl["variant"]
  
  a <- vipAnnotations[as.character(variant) == vipAnnotations$ID,]
  
  p <- which(!is.na(a$PICK))
  
  return(a[p,"Gene"])
  
})

selection$consequenceGene  <- sapply(pickAnyGene2, function(z ){
  if(length(z) == 0){
    return ("")
  } else {
    return (z)
  }
})
selection[order(selection$consequenceGene),]
cat(selection[selection$consequenceGene == "ENSG00000104976","phenotype"], sep = "\n")
sum(selection$consequenceGene == "ENSG00000187045")

vipAnnotations[vipAnnotations$ID == "rs2918299" & !is.na(vipAnnotations$PICK),]
vipAnnotations[vipAnnotations$ID == "rs117755721" & !is.na(vipAnnotations$PICK),]



dim(eqtlsTransNoCis)


sort(table(eqtlsTransNoCis$consequence))



sort(table(eqtlsTransNoCis$consequenceGene))


sum()
View(eqtls[eqtls$phenotype == "ENSG00000159216" & eqtls$type == "cis",])

View(eqtls[eqtls$phenotype == "ENSG00000103740",])

sum(eqtlsTransNoCis$consequenceGene == "ENSG00000159216", na.rm = T)
eqtlsTransNoCis[!is.na(eqtlsTransNoCis$consequenceGene) & eqtlsTransNoCis$consequenceGene == "ENSG00000159216" , c("variant", "phenotype", "consequence", "consequenceGene")]


transGenes <- eqtlsTransNoCis[!is.na(eqtlsTransNoCis$consequenceGene) & eqtlsTransNoCis$consequenceGene == "ENSG00000159216" , c("phenotype")]
transVariantPossibleCisDriven <- as.character(eqtls[eqtls$phenotype %in% transGenes & eqtls$chromosome == 21 & eqtls$bp >= 33787801 & eqtls$bp <= 37004667,"variant_index"])
transVariantPossibleCisDriven <- transVariantPossibleCisDriven[ ! transVariantPossibleCisDriven %in% transVariants]
sum(transVariantPossibleCisDriven %in% transVariants)
transVariants <- as.character(eqtlsTransNoCis[!is.na(eqtlsTransNoCis$consequenceGene) & eqtlsTransNoCis$consequenceGene == "ENSG00000159216" , c("variant_index")])
runx1CisVariants <- as.character(eqtls[eqtls$phenotype == "ENSG00000159216" & eqtls$type == "cis","variant_index"])

library(heatmap3)

ldBlock <- ldPerChr[[21]][c(transVariants, runx1CisVariants), c(transVariants, runx1CisVariants)]

#a <- data.frame(type = factor(c(rep("transSnp", length(transVariants)), rep("runx1CisSnp", length(runx1CisVariants)))))
cols <- brewer.pal(n = 3, name = "Accent")
c <- c(rep(cols[1], length(transVariants)), rep(cols[2], length(runx1CisVariants)))
heatmap3(ldBlock^2, RowSideColors = c, ColSideColors = c, balanceColor = T, legendfun= function() showLegend(c("Trans variants no cis", "runx1 cis variants"), col = cols[1:2]), scale = "none")

cols

ldBlock <- ldPerChr[[21]][c(transVariants, runx1CisVariants,transVariantPossibleCisDriven), c(transVariants, runx1CisVariants,transVariantPossibleCisDriven)]

#a <- data.frame(type = factor(c(rep("transSnp", length(transVariants)), rep("runx1CisSnp", length(runx1CisVariants)))))
cols <- brewer.pal(n = 3, name = "Accent")
c <- c(rep(cols[1], length(transVariants)), rep(cols[2], length(runx1CisVariants)), rep(cols[3], length(transVariantPossibleCisDriven)))
heatmap3(ldBlock^2, RowSideColors = c, ColSideColors = c, balanceColor = T, legendfun= function() showLegend(c("Trans variants no cis", "runx1 cis variants","Other trans variants"), col = cols), scale = "none")

cgeneAnnotations <- read.delim("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/public_data/Homo_sapiens.GRCh38.106.gtf.gz", comment.char = "#", header = F)
geneAnnotations <- geneAnnotations[geneAnnotations$V3 == "gene",]
matches <- regexec("gene_id (ENSG\\d*);",geneAnnotations$V9)
geneAnnotations$gene <- sapply(regmatches(geneAnnotations$V9, matches), function(a){a[2]})

str(geneAnnotations)
table(geneAnnotations$V1)

lonp1Trans <- read.delim("lonp1Trans.txt")
lonp1Trans <- as.matrix(lonp1Trans[lonp1Trans$variant == "47286831",-1])[1,]

tail(sort(abs(lonp1Trans)))

a <- as.character(eqtls$variant_index[eqtls$phenotype == "ENSG00000196365"])
b <- "47286831"
View(eqtls[eqtls$phenotype == "ENSG00000196365",])


CCDC26 <- "ENSG00000229140"
RUNX1 <- "ENSG00000159216"

noCistransGenes <- eqtlsTransNoCis[!is.na(eqtlsTransNoCis$consequenceGene) & eqtlsTransNoCis$consequenceGene == RUNX1 , c("phenotype")]

allTrans <- eqtls[eqtls$type == "trans" & eqtls$chromosome == 21 & eqtls$bp >= 33787801 & eqtls$bp <= 37004667,"phenotype"]

otherTrans <- allTrans[!allTrans %in% noCistransGenes]



transGenes %in% rownames(coexp)



coexpSubset <- coexp[c(RUNX1, noCistransGenes), c(RUNX1, noCistransGenes)]

cols <- brewer.pal(n = 3, name = "Accent")
c <- c(cols[1], rep(cols[2], length(noCistransGenes)))


heatmap3(coexpSubset^2, scale = "none", balanceColor = T, RowSideColors = c, ColSideColors = c, legendfun= function() showLegend(c("Cis gene", "no cis eqtl trans"), col = cols), distfun = function(x) as.dist(1 - coexpSubset))

coexpSubset <- coexp[c(RUNX1, noCistransGenes,otherTrans), c(RUNX1, noCistransGenes,otherTrans)]

cols <- brewer.pal(n = 3, name = "Accent")
c <- c(cols[1], rep(cols[2], length(noCistransGenes)), rep(cols[3], length(otherTrans)))


heatmap3(coexpSubset^2, scale = "none", balanceColor = T, RowSideColors = c, ColSideColors = c, legendfun= function() showLegend(c("Cis gene", "no cis eqtl trans","other trans"), col = cols), distfun = function(x) as.dist(1 - coexpSubset))

eqtlsTransNoCis$consequenceGene
geneCor <- apply(eqtlsTransNoCis, 1, function(eqtl){
  
  gene <- eqtl["phenotype"]
  originGene <- eqtl["consequenceGene"]
  
  if(is.na(originGene) | originGene == ""){
    return(NA)
  } else if(! originGene %in% rownames(coexp)){
    return(NA)
  } else {
    return(coexp[originGene, gene])
  }
  
})

hist(geneCor ^2, breaks = 100, xlab = "coexpression r2", main = "trans gene vs local gene")
which(geneCor ^2 > 0.15)
eqtlsTransNoCis["rs1657817-ENSG00000145335",]


singleCellColoc <- read.delim("trans_eQTLs_ColocClassesAdded_ToPatrick_2025-09-03.txt.gz")
table(singleCellColoc$variant_type)

rownames(singleCellColoc) <- paste0(singleCellColoc$variant,"-",singleCellColoc$phenotype)
dim(singleCellColoc)
length(row.names(singleCellColoc) %in% row.names(eqtlsTransNoCis))

singleCellColoc$transNoCis <- row.names(singleCellColoc) %in% row.names(eqtlsTransNoCis)

singleCellColocTrans <- singleCellColoc[singleCellColoc$type == "trans",]

table(singleCellColocTrans$transNoCis, singleCellColocTrans$variant_type)

for(type in unique(singleCellColocTrans$variant_type)){
  
  print(type)
  
  print(fisher.test(table(singleCellColocTrans$variant_type == type, singleCellColocTrans$transNoCis))$p.value)
  
}

table(singleCellColocTrans$transNoCis)

singleCellColocTrans2 <- singleCellColocTrans[match(unique(singleCellColocTrans$variant), singleCellColocTrans$variant),]

t <- table(singleCellColocTrans2$variant_type == "Some blood-cell" | singleCellColocTrans2$variant_type == "Exclusively blood-cell"| singleCellColocTrans2$variant_type == "Exclusively non-blood-cell", singleCellColocTrans2$transNoCis)
t
prop.table(t, margin = 2)
fisher.test(t)
table(singleCellColocTrans2$variant_type)
"Not GWAS colocalising"
selection$variant_type <- singleCellColoc[match(row.names(selection), row.names(singleCellColoc)),"variant_type"] 

length(unique(eqtlsTransNoCis$phenotype))


eqtls$cd4Sig <- !is.na(eqtls$scMeta_CD4_T_fdr) & eqtls$scMeta_CD4_T_fdr <= 0.05

singleCellColocTrans$cd4Sig <- eqtls[match(rownames(singleCellColocTrans), rownames(eqtls)),"cd4Sig"]

a <- table(singleCellColocTrans2$cd4Sig, singleCellColocTrans2$transNoCis)
fisher.test(a)



sum(!is.na(eqtls$scMeta_CD4_T_fdr) & eqtls$scMeta_CD4_T_fdr <= 0.05)

str(selection)
str(eqtls$sc)