commandArgs()

for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
  } else {
    assign(ta[[1]][1],TRUE)
    cat("assigned ",ta[[1]][1]," the value of TRUE\n")
  }
}


N <- as.matrix(read.table("nSpets.txt"))
extensionLength <- as.numeric(extensionLengthStr)
genomeLength <- as.numeric(genomeLengthStr)
genomeCoverageRatio <- as.numeric(genomeCoverageRatioStr)

data <- as.matrix(read.table("data.txt"))
nRows <- length(data[,1])
lambdaGlobal <- array(N*2*extensionLength / (genomeLength*genomeCoverageRatio), c(nRows,1))
lambda10K <- (data[,2] - data[,1])*2*extensionLength/10000;
lambda20K <- (data[,3] - data[,1])*2*extensionLength/20000;
lambdaTemp <- array(0, c(nRows,3))
lambdaTemp[,1] <- lambdaGlobal
lambdaTemp[,2] <- lambda10K
lambdaTemp[,3] <- lambda20K
lambda <- apply(lambdaTemp, 1, max)
x <- ppois(data[,1]-1, lambda, lower.tail = FALSE)
fdr <- p.adjust(x, "BH")
x <- sapply(x, function(a){if(a<0.005){b = as.numeric(format(a, scientific = T, digits = 3))}else{b = round(a, 2)}})
fdr<- sapply(fdr, function(a){if(a<0.005){b = as.numeric(format(a, scientific = T, digits = 3))}else{b = round(a, 2)}})
#x <- as.numeric(format(x, scientific = F))
#fdr <- as.numeric(format(x, scientific = F))
a <- matrix(1, nRows, 2)
a[,1] = x
a[,2] = fdr
write.table(a, file="result.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
