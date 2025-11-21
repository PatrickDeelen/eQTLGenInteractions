library(heatmap3)
library(RColorBrewer)
palette(brewer.pal(n = 3, name = "Accent"))


#

datasets <- data.frame(
  name = c("COVID19", "AIDA", "OneK1K", "SLE", "Aging", "Tabula", "TabulaSmartSeq2"),
  file = c("40ee3b22-395b-451b-883d-fb7e34ee612f.h5ad", "3a8a77a6-f069-479c-a2c1-252feafd0106.h5ad", "078a26dc-0585-4b94-9252-ee2f2dda5742.h5ad", "d51627ad-0123-4eb9-82b0-75f017862307.h5ad", "e11a76f0-57ff-4358-8a91-008655475059.h5ad", "", ""),
  assayFilter = c("10x 5' v2", NA, NA, NA, "10x 5' v2", "10x 3' v3", "Smart-seq2"),
  celltypeCol = c("celltype", "Annotation_Level4", "predicted.celltype.l2", "cell_type", "", "", "")
)

rownames(datasets) <- datasets$name


traitDendo <- readRDS(file = file.path("/groups/umcg-fg/tmp04/projects/downstreamer/lotte/traitDendo.rds"))




