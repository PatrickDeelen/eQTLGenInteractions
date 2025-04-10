
library(optparse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(caret)
library(readr)

option_list <- list(
  make_option(c("-s", "--general_covs"), type = "character",
              help = "Path to the covariate file with sex and age. Sample ids = genotype ids."),
  make_option(c("-c", "--cell_counts"), type = "character", default = NULL,
              help = "Path to the covariate file with (deconvoluted) cell counts. Sample ids - expression ids."),
  make_option(c("-g", "--genotype_pcs"), type = "character", default = NULL,
              help = "Path to the covariate file with genotype PCs. Sample ids = genotype ids."),
  make_option(c("-r", "--rna_qual"), type = "character", default = NULL,
              help = "Path to the file with RNA quality (for RNA-seq). Sample ids = expression ids."),
  make_option(c("-v", "--eigenAndIca"), type = "character", default = NULL,
              help = "Path to the file with eigenvectors and ICs"),
  make_option(c("-e", "--expression"), type = "character", default = NULL,
              help = "Path to the file with eigenvectors and ICs"),
  make_option(c("-i", "--gte"), type = "character",
              help = "Path to the genotype to expression id conversion file."),
  make_option(c("-o", "--out"), type = "character",
              help = "Output file name."),
  make_option(c("-n", "--rename"), type = "boolean",
              help = "Are raw expression sample ids different from genotype sample ids?.")
)

parser <- OptionParser(usage = "%prog [options] file", option_list = option_list)
args <- parse_args(parser)

getwd()
cat(paste0("args$",names(args), " = \"", unlist(args), "\"\n"))


rename_samples <- function(d, gte, key_col = 1, value_col = 2){
  if (length(intersect(row.names(d), gte[,key_col])) < length(intersect(row.names(d), gte[,value_col]))){
    cat("Sample ids are already converted to expression ids! Skipping id conversion.\n")
    return(d)
  }
  d_m <- d[match(gte[,key_col], row.names(d), nomatch = 0),, drop = F]
  gte_m <- gte[match(row.names(d_m), gte[,key_col], nomatch = 0),]
  row.names(d_m) <- gte_m[,value_col]
  cat("Coverted genotype sample ids to expression sample ids.\nN samples before: ", nrow(d), "\nN samples after: ", nrow(d_m), "\n")
  return(d_m)
}

readDoubleMatrix <- function(path){

  firstColName <- scan(path, what = "character", n = 1, sep = '\t', quiet = T )
  if(firstColName == ""){
    # if there is no column name for the first column ...1is used by read_delim of the readr lib
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

# Combine general covariates with cell counts if they are provided

covar_main <- read.delim(args$general_covs, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")


covar_main$age[is.na(as.numeric(covar_main$age)) & !is.na(covar_main$age)]

if(!all(c("gender", "age") %in% colnames(covar_main))){
  stop("Covariate file does not contain gender and sex column")
}

if(any(!unique(covar_main$gender) %in% c("M","F",NA))){
  stop("Gender column contains unexpected values, should be M F or NA")
}

#F1 M2
covar_main$gender[covar_main$gender=="F" ] <- 1
covar_main$gender[covar_main$gender=="M" ] <- 2
covar_main$gender <- as.numeric(covar_main$gender)


gte <- read.delim(args$gte, check.names = F, header = F, as.is = T, sep ="\t")


if (! is.null(args$cell_counts) && args$cell_counts != "NA"){
  covar_cell_counts <- read.delim(args$cell_counts, check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")
  covar_cell_counts <- rename_samples(covar_cell_counts, gte, key_col = 2, value_col = 1)
  covar <- merge(covar_main, covar_cell_counts, by = 0)
  row.names(covar) <- covar$Row.names
  covar$Row.names = NULL
  cat("Combined major covariates with cell counts\n")

} else {
  covar <- covar_main
  cat("No cell counts provided!\n")

}
str(covar)



# add genotype PCs
geno_pcs <- read.delim(args$genotype_pcs, check.names = F, header = T, row.names = 1, as.is = T)
colnames(geno_pcs) <- paste0("GenoPC", 1:ncol(geno_pcs))
#geno_pcs <- rename_samples(geno_pcs[,1:4], gte)

covar_merged <- merge(covar, geno_pcs, by = 0)
row.names(covar_merged) <- covar_merged$Row.names
covar_merged$Row.names = NULL

cat("Added genotype PCs\n")


if (!is.null(args$rna_qual)){
  rna_qual <- read.delim(args$rna_qual, check.names = F, header = T, row.names = 1, as.is = T)
  #rna_qual <- rename_samples(rna_qual, gte)
  covar_merged <- merge(covar_merged, rna_qual, by = 0)
  row.names(covar_merged) <- covar_merged$Row.names
  covar_merged$Row.names = NULL

  cat("Added RNA quality\n")
}

if (!is.null(args$eigenAndIca)){
  eigenAndICa <- read.delim(args$eigenAndIca, check.names = F, header = T, row.names = 1, as.is = T)
  covar_merged <- merge(covar_merged, eigenAndICa, by = 0)
  row.names(covar_merged) <- covar_merged$Row.names
  covar_merged$Row.names = NULL

  cat("Added mapped eigenvectors and ICs quality\n")

}

#Remove covariates with only NA values
fullMissingCovar <- apply(covar_merged, 2, function(x){all(is.na(x))})
if(any(fullMissingCovar)){
  write("Excluded covariates with full missingness", stderr())
  write(paste0(" - ", colnames(covar_merged)[fullMissingCovar]), stderr(), sep = "\n")
  covar_merged <- covar_merged[,colSums(is.na(covar_merged)) < nrow(covar_merged), drop = F]
}

# Remove covariates with near zero variance
near_zero_var <- nearZeroVar(covar_merged,  uniqueCut = 5)
if (length(near_zero_var) > 0){

  write("Excluded covariates with low variablility", stderr())
  write(paste0(" - ", colnames(covar_merged)[near_zero_var]), stderr(), sep = "\n")

  covar_merged <- covar_merged[, -near_zero_var]
}

# Remove samples with NA values
anyMissingSample <- apply(covar_merged, 1, function(x){any(is.na(x))})

cat("Test Number of covariates in the combined file:", ncol(covar_merged), "\n")


if(any(anyMissingSample)){
  write(paste0("Excluded samples with NA values in covariates: ",sum(anyMissingSample)), stderr())
    covar_merged <- covar_merged[!anyMissingSample, ]
}

cat("Number of covariates in the combined file:", ncol(covar_merged), "\n")
cat("Number of samples with covariate data available:", nrow(covar_merged), "\n")


# run INT on general covariates, cell counts and RNA quality
covar_names_for_INT <- names(which(apply(covar_merged, 2, function(x) length(unique(x)) > 3)))
covar_merged_int <- cbind(covar_merged[, !colnames(covar_merged) %in% covar_names_for_INT, drop =F], apply(covar_merged[,covar_names_for_INT, drop = F], 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) ))
#colnames(covar_merged_int) <- c(colnames(covar_merged[, !colnames(covar_merged) %in% covar_names_for_INT]), covar_names_for_INT)




#Expression data is already INT so we add that now
expression <- readDoubleMatrix(args$expression)
covar_merged_int <- merge(covar_merged_int, expression, by = 0)
row.names(covar_merged_int) <- covar_merged_int$Row.names
covar_merged_int$Row.names = NULL

cat("Added pre INT expression data\n")



# Plot covariate distributions prior to INT
# pdf(paste0(args$out, ".distributions.pdf"), height = 5, width = 10)
# #covar_long <- pivot_longer(covar_merged, cols = everything())
# for (covar_name in colnames(covar_merged)[1:2]){
#   p1 <- ggplot(covar_merged, aes(get(covar_name))) + geom_histogram(bins = 50, color = "black", alpha = 0.5, fill = "dodgerblue4") + theme_bw() + xlab(covar_name) + ggtitle(paste0("Raw ", covar_name))
#   p2 <- ggplot(covar_merged_int, aes(get(covar_name))) + geom_histogram(bins = 50, color = "black", alpha = 0.5, fill = "dodgerblue4") + theme_bw() + xlab(covar_name) + ggtitle(paste0("INT ", covar_name))
#   print(p1 + p2)
# }
# dev.off()

# Write INT covariates
write.table(covar_merged_int, file = args$out, sep = "\t", quote = F, col.names = NA)

write.table(colnames(covar_merged_int), file = "availableCovariates.txt", sep = "\t", quote = F, col.names = F, row.names = F)


