


library(readr)


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

args <- commandArgs(trailingOnly = TRUE)

data <- readDoubleMatrix(args[1])

dataT <- t(data)

write.table(dataT, file = args[2], quote = F, sep = "\t", col.names = NA)