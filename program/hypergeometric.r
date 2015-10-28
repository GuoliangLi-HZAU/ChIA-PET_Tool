N <- as.matrix(read.table("nTags.txt"))
###### input format for data
###### AB  A  B
###### AB: pet count between anchors A and B
###### A: tags in anchor A
###### B: tags in anchor B
data <- as.matrix(read.table("data.txt"))
nRows <- length(data[,1])
q <- data[,1] - 1
m <- data[,2]
n <- N - m
k <- data[,3]
x <- phyper(q, m, n, k, lower.tail = FALSE)

###### p-value adjustment with Benjamini-Hockberg method
fdr <- p.adjust(x, "BH")
a <- matrix(1, nRows, 4)
a[,1] <- x
a[,2] <- fdr
a[,3:4] <- -log10(a[,1:2]) ###### -log10(p-value)
a[is.infinite(a[,3]),3] <- 1000 ###### replace the extreme values with 1000
a[is.infinite(a[,4]),4] <- 1000 ###### replace the extreme values with 1000
a <- as.data.frame(a)
for (i in 1:2) {
  a[[i]] <- as.numeric(format(a[[i]], digits = 3))
}
for (i in 3:4) {
  a[[i]] <- round(a[[i]], 2)
}
write.table(a, file="result.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

