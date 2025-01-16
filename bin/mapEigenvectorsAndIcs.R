library(optparse)
library(readr)
library(MASS)

option_list <- list(
  make_option(c("-e", "--expression"), type = "character",
              help = "Path to the eQTLgen normalized expression data"),
  make_option(c("-v", "--eigenvectors"), type = "character",
              help = "Path to the eigenvectors as calculated in the full eQTLgen data"),
  make_option(c("-i", "--ics"), type = "character",
              help = "Path to the independent components as calculated in the full eQTLgen data"),
  make_option(c("-o", "--out"), type = "character",
              help = "Output file name.")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser)


cat(paste0("args$",names(args), " = \"", unlist(args), "\"\n"))

readDoubleMatrix <- function(path){

  firstColName <- scan(path, what = "character", n = 1, sep = '\t', quiet = T )
  if(firstColName == ""){
    # if there is no column name for the first column ...1 is used by read_delim of the readr lib
    firstColName <- "...1"
  }

  l <- list()
  l[[firstColName]] <- col_character()
  colTypes <- structure(list(cols = l, default = col_double()), class = "col_spec")

  table_tmp <- read_delim(path, delim = "\t", quote = "", col_types = colTypes)
  table <- as.matrix(table_tmp[,-1])
  rownames(table) <- table_tmp[,1][[1]]
  return(table)
}

expression <- readDoubleMatrix(args$expression)
eigenvectors <- readDoubleMatrix(args$eigenvectors)
ics <- readDoubleMatrix(args$ics)

# Sometimes the expression matrix is transposed
if(length(intersect(rownames(eigenvectors), colnames(expression)) == 0){
    expression <- t(expression)
}

sharedGenes <- intersect(rownames(eigenvectors), colnames(expression))

length(sharedGenes)

if(length(sharedGenes) < 10000){
  stop ("Not enough genes matching")
}

mappedEigenvectors <- expression[,sharedGenes] %*% t(ginv(eigenvectors[sharedGenes,]))
colnames(mappedEigenvectors) <- colnames(eigenvectors)

mappedIcs <- expression[,sharedGenes] %*% t(ginv(ics[sharedGenes,]))
colnames(mappedIcs) <- paste0("IC", 1:ncol(ics))

mappedCombined <- cbind(mappedEigenvectors, mappedIcs)
write.table(mappedCombined, file = args$out, quote = F, sep = "\t", col.names = NA)
