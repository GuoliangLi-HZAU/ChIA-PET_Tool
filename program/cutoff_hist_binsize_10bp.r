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

#### plus_plus
data_file <- sprintf("%s.bedpe.selected.intra-chrom.distance.plusplus.txt", filePrefix)
plusplus <- as.matrix(read.table(data_file))
#### plus_minus
data_file <- sprintf("%s.bedpe.selected.intra-chrom.distance.plusminus.txt", filePrefix)
plusminus <- as.matrix(read.table(data_file))
#### minus_plus
data_file <- sprintf("%s.bedpe.selected.intra-chrom.distance.minusplus.txt", filePrefix)
minusplus <- as.matrix(read.table(data_file))
#### minus_minus
data_file <- sprintf("%s.bedpe.selected.intra-chrom.distance.minusminus.txt", filePrefix)
minusminus <- as.matrix(read.table(data_file))

breaks <- c(-3E9, seq(-100000,100000,10), 3E9)
mids <- seq(-100000,100000-10,10) + 10/2

h.plusplus   <- hist(plusplus,   breaks=breaks)$counts
h.plusminus  <- hist(plusminus,  breaks=breaks)$counts
h.minusplus  <- hist(minusplus,  breaks=breaks)$counts
h.minusminus <- hist(minusminus, breaks=breaks)$counts
max_count <- max(c(h.plusplus, h.plusminus, h.minusplus, h.minusminus))

write.table(h.plusplus,   file="h.plusplus.binsize_10bp.txt")
write.table(h.plusminus,  file="h.plusminus.binsize_10bp.txt")
write.table(h.minusplus,  file="h.minusplus.binsize_10bp.txt")
write.table(h.minusminus, file="h.minusminus.binsize_10bp.txt")

distribution_diff <- h.minusplus - (h.plusplus + h.plusminus + h.minusminus)/3
write.table(distribution_diff, file="h.diff.binsize_10bp.txt")


hist_plot_file <- sprintf("%s.hist.binsize_10bp.png", filePrefix)
png(hist_plot_file, width = 960, height = 480)
plot(  mids, h.plusplus[2:20001], xlim=c(-100000,100000), ylim=c(1,max_count), xlab="pet span (bp)", ylab="frequency of pets", col="red", log="y")
points(mids, h.plusminus[2:20001],  col="green")
points(mids, h.minusplus[2:20001],  col="blue")
points(mids, h.minusminus[2:20001], col="black")
legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), col=c("red", "green", "blue", "black"))
dev.off()

indexes_1 <- 10001:20001
logMids <- log10(mids[indexes_1]+10/2)

hist_plot_file.loglog <- sprintf("%s.hist.loglog.binsize_10bp.png", filePrefix)
png(hist_plot_file.loglog, unit = "in", width = 7, height = 5, res = 300, pointsize = 4)
indexes <- (h.plusplus[indexes_1] > 0)
plot(  logMids[indexes], (h.plusplus[indexes_1])[indexes],  col="red", xlim=c(1.5,5), ylim=c(1,max_count), xlab="Pet span (bp)", ylab="frequency of pets", xaxt="n", log="y")
indexes <- (h.plusminus[indexes_1] > 0)
points(logMids[indexes], (h.plusminus[indexes_1])[indexes], col="green")
indexes <- (h.minusplus[indexes_1] > 0)
points(logMids[indexes], (h.minusplus[indexes_1])[indexes], col="blue")
indexes <- (h.minusminus[indexes_1] > 0)
points(logMids[indexes], (h.minusminus[indexes_1])[indexes], col="black")
axis(1, at=seq(-5,5,1), lab=c(-100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000))
legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), text.col=c("red", "green", "blue", "black"))
for (i in -4:4) {
    abline(v=i, lty="dashed")
}
dev.off()

indexes_1 <- 10001:20001
logMids <- log10(mids[indexes_1]+10/2)

hist_diff_plot_file.loglog <- sprintf("%s.hist_diff.loglog.binsize_10bp.png", filePrefix)
png(hist_diff_plot_file.loglog, unit = "in", width = 7, height = 5, res = 300, pointsize = 4)
indexes <- (distribution_diff[indexes_1] > 0)
plot(  logMids[indexes], (distribution_diff[indexes_1])[indexes],  col="red", xlim=c(1.5,5), ylim=c(1,max_count), xlab="Pet span (bp)", ylab="frequency of pets", xaxt="n", log="y")
axis(1, at=seq(-5,5,1), lab=c(-100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000))
legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), text.col=c("red", "green", "blue", "black"))
for (i in -4:4) {
    abline(v=i, lty="dashed")
}
#abline(v=log10(8000), lty="dashed")
dev.off()

