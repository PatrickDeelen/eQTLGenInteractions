

ldPerChr <- readRDS("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/fix_20250509/ld_matrices_20250509.rds")
coexp <- readRDS("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/all_permuted_independentVariants_4GenPC20ExpPC_2025-02-05/exported_matrix/coexpression_matrix_16781_genes.rds")


cisTransColoc <- read.delim("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/cisTransCoexp/trans-eQTLs_annotated_with_ctc_class_and_cis-eGene_20251016.txt.gz")
cisTransColoc$qtl <- paste0(cisTransColoc$variant,"-", cisTransColoc$phenotype)
rownames(cisTransColoc) <- cisTransColoc$qtl
str(cisTransColoc)

load("/groups/umcg-fg/tmp04/projects/eqtlgen-phase2/freeze3/Interpretation/eqltVariantAnalysis/wip7.RData")


traitMrRes <- read.delim(paste0("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/rb/snp_gene_","sle" , ".tsv"), sep = "\t")  
traitMrRes$qtl <- paste0(traitMrRes$variant,"-", traitMrRes$phenotype)

traitMrRes <- traitMrRes[match(unique(traitMrRes$qtl),traitMrRes$qtl),]
traitMrRes <- traitMrRes[traitMrRes$mr_gene,]
rownames(traitMrRes) <- traitMrRes$qtl

View(eqtls[rownames(traitMrRes),])
View(traitMrRes)

all(rownames(traitMrRes) %in% rownames(cisTransColoc))

dim(traitMrRes)

qtl <- rownames(traitMrRes)[3]
sapply(rownames(traitMrRes), function(qtl){
  
  colocRes <- cisTransColoc[qtl,"cis_gene_list"]
  
})

cisTransColoc[rownames(traitMrRes),"cis_gene_list"]
