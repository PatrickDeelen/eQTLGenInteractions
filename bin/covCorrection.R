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

covariates <- readDoubleMatrix("covariates.combined.txt")
covariates <- covariates[,c(paste0("Comp", 1:25),"AvgExprCorrelation")]
covariates <- cbind(1,covariates)

geneExp <- readDoubleMatrix("expressionT.txt")
geneExp <- t(geneExp)
geneExp <- geneExp[rownames(covariates),]

correctedExp <- matrix(NA, nrow=nrow(geneExp), ncol=ncol(geneExp), dimnames = dimnames(geneExp))

for(i in 1:ncol(geneExp)){
  correctedExp[,i] <- residuals(lm.fit(covariates, geneExp[,i], ))
}

write.table(t(correctedExp), file = "correctedExpression.txt", sep = "\t", col.names = NA, quote =F )