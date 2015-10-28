## =================================================================== ##
## This file generates all the images needed for the report
## =================================================================== ##

# General settings and parameters passing
options(warn = -1)
args <- commandArgs(trailingOnly = TRUE)

program_dir <- args[1]
setwd(args[2])
input_prefix <- args[3]
output_dir <- paste(args[2], "/ChIA-PET_Tool_Report/images/", sep = "")
cytoband_file_path <- paste(args[1], "/Rscript_and_genome_data/", args[4], sep = "")
species <- as.numeric(args[5])
#program_dir <- "/home/lycai/Project/test1/program"
#output_dir <- "/home/lycai/Project/test1/Output_dir_2/ChIA-PET_Tool_Report/images/"
#input_prefix <- "Bcell"
#setwd("/home/lycai/Project/test1/Output_dir_2")
#cytoband_file_path <- "/home/lycai/Project/test1/program/Rscript_and_genome_data/Mouse_cytoBandIdeo.txt"


# sink(paste(args[2], "/ChIA-PET_Tool_Report/log.txt", sep = ""), type = c("output", "message"), append = T)

# Check if necessary libraries are installed 
check_pkg <- function(pkg) {
    if(require(pkg, character.only = TRUE)){
        print(paste("Package", pkg, "is loaded correctly", sep = " "))
    } else {
        print(paste("Trying to install package", pkg, sep = " "))
        install.packages(pkg, repos="http://cran.us.r-project.org", dep = TRUE)
        if(require(pkg, character.only = TRUE)){
            print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
        } else{
            install.packages(pkg, repos="http://cran.rstudio.com/", dep = TRUE)
            if(require(pkg, character.only = TRUE)){
               print(paste("Package", pkg, "is installed and loaded correctly", sep = ""))
               } else{
                stop(paste("Couldn't install package", pkg, sep = " "));
              }
        }
    }
}
check_pkg("grid")
check_pkg("xtable")
check_pkg("RCircos")


# source plotting fucntions
source(paste(program_dir, "/Rscript_and_genome_data/Plotting_functions.R", sep = ""))

## =================================================================== ##
## Set color palettes for later use
## =================================================================== ##
fill_palette_1 <- c("#4D4D4D", "#C3C3C3")
fill_palette_2 <- c("#4D4D4D", "#969696", "#C3C3C3", "#E6E6E6", "#F5F5F5")
# fill_palette_3 <- c("#2A286B", "#4E69B2", "#EDA1AD", "#B9539F")
#fill_palette_3 <- c("#2A286B", "#4E69B2", "#72A7E5", "#E0ECF9")
fill_palette_3 <- c("#2F7ED8", "#72A7E5", "#95BDEB", "#DCE9F8")
fill_palette_4 <- c("#4D4D4D", "#969696", "#C3C3C3", "#E6E6E6")



format_axis_labels <- function(a) {
  sapply(a, function(x) {if (x >= 1000000) {x = paste(x/1000000, "M", sep = "")}
                         else if (x>=1000) {x = paste(x/1000, "K", sep = "")}
                         else {x = x}
})
}


## =================================================================== ##
## .summary.txt
## =================================================================== ##
summary_info <- read.table(paste(input_prefix, ".summary.txt", sep = ""), sep = "\t", header = F)
colnames(summary_info) <- c("Title", "Percentage")
summary_info$Percentage <- round(summary_info$Percentage, 4)
summary_info$Percentage <- paste(summary_info$Percentage * 100, "%", sep = "")
#name.width <- max(sapply(summary_info$Percentage, nchar))
#summary_info$Percentage <- format(summary_info$Percentage, justify = "right")
summary_info_tb <- xtable(summary_info, align = c("c", "l", "r"))
print(summary_info_tb, type='html', file=paste(output_dir, "Rplot25.html", sep = ""), include.rownames = F, append = T,
  html.table.attributes = "class = 'mytable'")



## =================================================================== ##
## .runningInformation.LinkerFiltering_FastQ_PET.txt
## =================================================================== ##
running_info <- read.table(paste(input_prefix, ".runningInformation.LinkerFiltering_FastQ_PET.txt", sep = ""), sep = "\t", header = F)
colnames(running_info) <- c("Title", "Information")
# running_info[, 1] <- gsub("_", " ", running_info[, 1])

# png graph
# png(paste(output_dir, "Rplot1.png", sep = ""), width = 650, height = 450, res = 80)
# grid.table(running_info, gpar.corefill = gpar(fill = "#E6E6E6", col = NA), gpar.colfill = gpar(fill = NA, col = NA),
#            gpar.coltext = gpar(col = "#585858", cex = 1.4), gpar.coretext = gpar(col = "black"),
#            h.odd.alpha = 0, h.even.alpha = 0.7, v.odd.alpha = 1, v.even.alpha = 1,
#            show.rownames = FALSE, padding.v = unit(6, "mm"))
# dev.off()

running_info_tb <- xtable(running_info, align = rep("c", 3))
print(running_info_tb, type='html', file=paste(output_dir, "Rplot1.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")



## =================================================================== ##
## .runningInformation.Mapping.txt
## =================================================================== ##
mapping_info <- read.table(paste(input_prefix, ".runningInformation.Mapping.txt", sep = ""), sep = "\t", header = F)
colnames(mapping_info) <- c("Title", "Information")
mapping_info_tb <- xtable(mapping_info, align = rep("c", 3))
print(mapping_info_tb, type='html', file=paste(output_dir, "Rplot18.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")



## =================================================================== ##
## .runningInformation.RemovingRedundancy.txt
## =================================================================== ##
remove_redundancy_info <- read.table(paste(input_prefix, ".runningInformation.RemovingRedundancy.txt", sep = ""), sep = "\t", header = F)
colnames(remove_redundancy_info) <- c("Title", "Information")
remove_redundancy_info_tb <- xtable(remove_redundancy_info, align = rep("c", 3))
print(remove_redundancy_info_tb, type='html', file=paste(output_dir, "Rplot19.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")
	  
	  
	  
## =================================================================== ##
## .runningInformation.Categories.txt
## =================================================================== ##
category_info <- read.table(paste(input_prefix, ".runningInformation.Categories.txt", sep = ""), sep = "\t", header = F)
colnames(category_info) <- c("Title", "Information")
category_info_tb <- xtable(category_info, align = rep("c", 3))
print(category_info_tb, type='html', file=paste(output_dir, "Rplot20.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")	  
	  
	  
## =================================================================== ##
## .runningInformation.Peakcalling.txt
## =================================================================== ##
peakcalling_info <- read.table(paste(input_prefix, ".runningInformation.Peakcalling.txt", sep = ""), sep = "\t", header = F)
colnames(peakcalling_info) <- c("Title", "Information")
peakcalling_info$Information <- as.character(peakcalling_info$Information)
peakcalling_info_tb <- xtable(peakcalling_info, align = rep("c", 3))
print(peakcalling_info_tb, type='html', file=paste(output_dir, "Rplot21.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")	


	  
## =================================================================== ##
## .runningInformation.clustering.txt
## =================================================================== ##
clustering_info <- read.table(paste(input_prefix, ".runningInformation.clustering.txt", sep = ""), sep = "\t", header = F)
colnames(clustering_info) <- c("Title", "Information")
clustering_info$Information <- as.character(clustering_info$Information)
clustering_info_tb <- xtable(clustering_info, align = rep("c", 3))
print(clustering_info_tb, type='html', file=paste(output_dir, "Rplot22.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")	
	  
	  
	  
## =================================================================== ##
## .runningInformation.Report.txt
## =================================================================== ##
report_info <- read.table(paste(input_prefix, ".runningInformation.Report.txt", sep = ""), sep = "\t", header = F)
colnames(report_info) <- c("Title", "Information")
report_info_tb <- xtable(report_info, align = rep("c", 3))
print(report_info_tb, type='html', file=paste(output_dir, "Rplot23.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")	
	  
	  
	  
## =================================================================== ##
## .basic_statistics.txt
## =================================================================== ##
basic_stat <- read.table(paste(input_prefix, ".basic_statistics.txt", sep = ""), sep = "\t", header = F)
colnames(basic_stat) <- c("PET category", "Number")
#basic_stat$Number <- format(basic_stat$Number, big.mark = ",")
#class(Number)
basic_stat$Number <- as.numeric(basic_stat$Number)
#class(Number)
#print(Number[2]/Number[1])
Percentage <- round(c(basic_stat$Number[2]/basic_stat$Number[1], basic_stat$Number[3]/basic_stat$Number[2], basic_stat$Number[4]/basic_stat$Number[3], basic_stat$Number[5]/basic_stat$Number[4], basic_stat$Number[6]/basic_stat$Number[5], basic_stat$Number[7]/basic_stat$Number[5], basic_stat$Number[8]/basic_stat$Number[5]) * 100, 2)
#Percentage <- round(c(Number[2]/Number[1], Number[3]/Number[2], Number[4]/Number[3], Number[5]/Number[4], Number[6]/Number[4], Number[7]/Number[4]) * 100, 2)
#print(Percentage)
Percentage <- paste(Percentage, "%", sep = "")
#print(Percentage)
#Percentage <- c("\", paste(Percentage[1], "of", "1.", sep = " "), paste(Percentage[2], "of", "2.", sep = " "), paste(Percentage[3], "of", "3.", sep = " "), paste(Percentage[4], "of", "4.", sep = " "), paste(Percentage[5], "of", "4.", sep = " "), paste(Percentage[6], "of", "4.", sep = " "), "\", "\")
Percentage <- c("NA", paste(Percentage[1], "of (1)", sep = " "), paste(Percentage[2], "of (2)", sep = " "), paste(Percentage[3], "of (3)", sep = " "), paste(Percentage[4], "of (4)", sep = " "), paste(Percentage[5], "of (5)", sep = " "), paste(Percentage[6], "of (5)", sep = " "), paste(Percentage[7], "of (5)", sep = " "), "NA", "NA")
#print(Percentage)

Percentage_2 <- round(c(basic_stat$Number[2]/basic_stat$Number[1], basic_stat$Number[3]/basic_stat$Number[1], basic_stat$Number[4]/basic_stat$Number[1], basic_stat$Number[5]/basic_stat$Number[1], basic_stat$Number[6]/basic_stat$Number[1], basic_stat$Number[7]/basic_stat$Number[1], basic_stat$Number[8]/basic_stat$Number[1]) * 100, 2)
#Percentage <- round(c(Number[2]/Number[1], Number[3]/Number[2], Number[4]/Number[3], Number[5]/Number[4], Number[6]/Number[4], Number[7]/Number[4]) * 100, 2)
#print(Percentage)
Percentage_2 <- paste(Percentage_2, "%", sep = "")
#print(Percentage)
#Percentage <- c("\", paste(Percentage[1], "of", "1.", sep = " "), paste(Percentage[2], "of", "2.", sep = " "), paste(Percentage[3], "of", "3.", sep = " "), paste(Percentage[4], "of", "4.", sep = " "), paste(Percentage[5], "of", "4.", sep = " "), paste(Percentage[6], "of", "4.", sep = " "), "\", "\")
Percentage_2 <- c("NA", paste(Percentage_2[1], "of (1)", sep = " "), paste(Percentage_2[2], "of (1)", sep = " "), paste(Percentage_2[3], "of (1)", sep = " "), paste(Percentage_2[4], "of (1)", sep = " "), paste(Percentage_2[5], "of (1)", sep = " "), paste(Percentage_2[6], "of (1)", sep = " "), paste(Percentage_2[7], "of (1)", sep = " "), "NA" , "NA")
Order <- paste("(", c(1:10), ")", sep = "")
basic_stat <- cbind(basic_stat, Percentage, Percentage_2, Order)
basic_stat$Number <- format(basic_stat$Number, big.mark = ",")
colnames(basic_stat)[4] <- "Percentage of total PETs"
#print(basic_stat)
basic_stat_tb <- xtable(basic_stat, align = c("c", "l", rep("r", 4)))
print(basic_stat_tb, type='html', file=paste(output_dir, "Rplot2.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")

# png graph
# png(paste(output_dir, "Rplot2.png", sep = ""), width = 650, height = 250, res = 80)
# grid.table(basic_stat, gpar.corefill = gpar(fill = "#E6E6E6", col = NA), gpar.colfill = gpar(fill = NA, col = NA),
#            gpar.coltext = gpar(col = "#585858", cex = 1.4), gpar.coretext = gpar(col = "#585858"),
#            h.odd.alpha = 0, h.even.alpha = 0.7, v.odd.alpha = 1, v.even.alpha = 1,
#            show.rownames = FALSE, padding.v = unit(6, "mm"), padding.h = unit(5.8, "cm"))
# dev.off()



## =================================================================== ##
## .running_time.txt
## =================================================================== ##
report_info <- read.table(paste(input_prefix, ".running_time.txt", sep = ""), sep = "\t", header = F)
colnames(report_info) <- c("Module", "Running time (min)")
report_info_tb <- xtable(report_info, align = rep("c", 3))
print(report_info_tb, type='html', file=paste(output_dir, "Rplot24.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")  



## =================================================================== ##
## .mapping_statistics.txt
## =================================================================== ##
files_order <- c("1_1", "2_2", "1_2", "2_1")
file_names <- paste(input_prefix, files_order, "mapping_statistics.txt", sep = ".")
data_list <- lapply(file_names, read.table, header = T, row.names = 1, sep = "\t")
#print(data_list[[1]])
for (i in 1:4) {
        colnames(data_list[[i]]) <- rownames(data_list[[i]])
        for (j in 1:3) {
          data_list[[i]][[j]] <- format(data_list[[i]][[j]], big.mark = ",")
        }
        #data_list[[i]] <- data_list[[i]][c(1, 3, 2), c(1, 3, 2)]
        #data_list[[i]] <- as.data.frame(sapply(data_list[[i]], function(q){format(q, big.mark = ",")}))
        print(xtable(data_list[[i]], align = rep("c", 4)), type='html', file=paste(output_dir, "Rplot3_", i, ".html", sep = ""), include.rownames = T, append = T,
              html.table.attributes = "class = 'mytable'")
}
# draw mapping statistics of 1_1, 2_2, 1_2, 2_1 on the same page
# png(paste(output_dir, "Rplot3.png", sep = ""), width = 650, height = 600, res = 80)
# grid.newpage()
# plot_layout(5, 3, heights = unit(c(0.05, rep(0.95 / 5, 5)), "npc"), widths = unit(c(0.05, 0.85, 0.1), "npc"))
# for (i in 1:4) {
#     colnames(data_list[[i]]) <- rownames(data_list[[i]])
#     data_list[[i]] <- data_list[[i]][c(1, 3, 2), c(1, 3, 2)]
#     if (i == 1) {colname_col = "#585858"} else {colname_col = "white"}
#     grid.table(data_list[[i]], gpar.corefill = gpar(fill = "#E6E6E6", col = NA), gpar.colfill = gpar(fill = NA, col = NA),
#                gpar.rowfill = gpar(fill = NA, col = NA), gpar.rowtext = gpar(col = "#585858", cex = 1.1),
#                gpar.coltext = gpar(col = colname_col, cex = 1.1), gpar.coretext = gpar(col = "#585858"),
#                h.odd.alpha = 0, h.even.alpha = 0.7, v.odd.alpha = 1, v.even.alpha = 1,
#                padding.v = unit(6, "mm"), padding.h = unit(1, "cm"), vp = viewport(layout.pos.row = i + 1, layout.pos.col = 2))
#     grid.text(files_order[i], x = 0.5, y = 0.4, vp = viewport(layout.pos.row = i + 1, layout.pos.col = 3), gp = gpar(cex = 1.5))
# 
# }
# grid.text("Read Two", x = 0.65, y = 1, vp = viewport(layout.pos.row = 1, layout.pos.col = 2), gp = gpar(col = "#585858", cex = 1.8))
# grid.text("Read One", x = 0.5, y = 0, vp = viewport(layout.pos.row = 3, layout.pos.col = 1), rot = 90, gp = gpar(col = "#585858", cex = 1.8))
# dev.off()



## =================================================================== ##
## .linker_alignment_score_distribution.txt
## =================================================================== ##
linker_align_score_dis <- read.table(paste(input_prefix, ".linker_alignment_score_distribution.txt", sep = ""), header = F)
temp <- t(as.matrix(linker_align_score_dis))
row.names(temp) <- c("Type", "Read1", "Read2")
png(paste(output_dir, "Rplot4.png", sep = ""), width = 680, height  = 350)
#pdf(paste(output_dir, "Rplot4.pdf", sep = ""), width = 9, height  = 6)
# draw_back_to_back(linker_align_score_dis, "#C77965", "#7E827A", "Score")   
#barplot(temp[-1, ], legend.txt = T, beside = T, names.arg = temp[1, ], col = c("#2F7ED8", "#7E827A"), yaxt="n", 
#        args.legend = list(x = "topright", inset=c(0.03, -0.2), cex = 1, bty = "n", border = NA))
par(xpd = T)
barplot(temp[-1, ], names.arg = temp[1, ], beside = T, col = c("#2F7ED8", "#7E827A"), yaxt="n")
legend("topright", inset = c(0.03, -0.2), legend = row.names(temp)[2:3], fill = c("#2F7ED8", "#7E827A"), cex = 1, bty = "n", border = NA) 
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)
dev.off()


## =================================================================== ##
## .linker_alignment_score_difference_distribution.txt
## =================================================================== ##
linker_align_score_diff <- read.table(paste(input_prefix, ".linker_alignment_score_difference_distribution.txt", sep = ""), header = F)
temp <- t(as.matrix(linker_align_score_diff))
row.names(temp) <- c("Type", "Read1", "Read2")
png(paste(output_dir, "Rplot5.png", sep = ""), width = 680, height = 350)
#pdf(paste(output_dir, "Rplot5.pdf", sep = ""), width = 9, height  = 6)
# draw_back_to_back(linker_align_score_diff, "#C77965", "#7E827A", "Score")  
#barplot(temp[-1, ], legend.txt = T, beside = T, names.arg = temp[1, ], col = c("#2F7ED8", "#7E827A"), yaxt="n", 
#        args.legend = list(x = "topright", inset=c(0.03, -0.2), cex = 1, bty = "n", border = NA))
par(xpd = T) 
barplot(temp[-1, ], names.arg = temp[1, ], beside = T, col = c("#2F7ED8", "#7E827A"), yaxt="n")
legend("topright", inset = c(0.03, -0.2), legend = row.names(temp)[2:3], fill = c("#2F7ED8", "#7E827A"), cex = 1, bty = "n", border = NA) 
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)
dev.off()


## =================================================================== ##
## .tag_length_distribution.txt
## =================================================================== ##
tag_len_dis <- read.table(paste(input_prefix, ".tag_length_distribution.txt", sep = ""))
temp <- t(as.matrix(tag_len_dis))
row.names(temp) <- c("Type", "Read1", "Read2")
png(paste(output_dir, "Rplot6.png", sep = ""), width = 680, height = 350)
#pdf(paste(output_dir, "Rplot6.pdf", sep = ""), width = 9, height  = 6)
# draw_back_to_back(tag_len_dis, "#C77965", "#7E827A", "Length (bp)")
# barplot(temp[-1, ], legend.txt = T, beside = T, names.arg = temp[1, ], col = c("#2F7ED8", "#7E827A"), yaxt="n", 
#        args.legend = list(x = "topright", inset=c(0.03, -0.2), cex = 1, bty = "n", border = NA))
par(xpd = T)
barplot(temp[-1, ], names.arg = temp[1, ], beside = T, col = c("#2F7ED8", "#7E827A"), yaxt="n")
legend("topright", inset = c(0.03, -0.2), legend = row.names(temp)[2:3], fill = c("#2F7ED8", "#7E827A"), cex = 1, bty = "n", border = NA)
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)   
dev.off()



## =================================================================== ##
## .linker_composition_distribution.txt
## =================================================================== ##
linker_compo <- read.table(paste(input_prefix, ".linker_composition_distribution.txt", sep = ""), check.names= F, header = T, stringsAsFactors = F)
#linker_no <- NULL
#for (i in 1:2) {
#    linker_no <- c(linker_no, strsplit(linker_compo[i], "\t")[[1]][1], strsplit(linker_compo[i], "\t")[[1]][2])
#}
#linker_no <- as.numeric(c(linker_no, strsplit(linker_compo[6], "\t")[[1]][[2]]))
#linker_no[6] <- sum(linker_no[1:5])
#names(linker_no) <- c("1_1", "1_2", "2_1", "2_2", "Ambigous", "Total")
#linker_no <- as.data.frame(t(linker_no))
#linker_no <- linker_no[c(1, 4, 2, 3, 5, 6)]
#print(linker_compo)
for (i in 1:ncol(linker_compo)) {
  linker_compo[[i]][1] <- format(as.numeric(linker_compo[[i]][1]), big.mark = ",")
}
linker_no_tb <- xtable(linker_compo, align = rep("c", 7))
print(linker_no_tb, type='html', file=paste(output_dir, "Rplot7.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")

# png graph
# png(paste(output_dir, "Rplot7.png", sep = ""), width = 650, height = 80, res = 80)
# grid.table(linker_no, gpar.corefill = gpar(fill = "#E6E6E6", col = NA), gpar.colfill = gpar(fill = NA, col = NA),
#            gpar.coltext = gpar(col = "#585858", cex = 1.3), gpar.coretext = gpar(col = "#585858"),
#            h.odd.alpha = 0, h.even.alpha = 0.7, v.odd.alpha = 1, v.even.alpha = 1,
#            show.rownames = FALSE, padding.v = unit(8, "mm"))
# dev.off()


## =================================================================== ##
## .bedpe.qc.dist.txt
## =================================================================== ##
#png(paste(output_dir, "Rplot8.png", sep = ""), width = 350, height = 350)
png(paste(output_dir, "Rplot8.png", sep = ""), unit = "in", width = 2.5, height = 1.85, res = 300, pointsize = 4)
#png(paste(output_dir, "Rplot8.png", sep = ""), width = 450, height = 350)
draw_barplot("bedpe.qc.dist", "V2", col = fill_palette_3, border = NA, yaxt="n", 
             xlab = "Quality Score")
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)  
dev.off()


## =================================================================== ##
## .bedpe.strand.dist.txt
## =================================================================== ##
# png(paste(output_dir, "Rplot9.png", sep = ""), width = 650, res = 100)
png(paste(output_dir, "Rplot9.png", sep = ""), unit = "in", width = 2.5, height = 1.85, res = 300, pointsize = 4)
# png(paste(output_dir, "Rplot9.png", sep = ""), width = 450, height = 350)
par(mar = c(5, 4, 4, 4))
draw_barplot("bedpe.strand.dist", c("V2", "V3"), col = fill_palette_3, border = NA, yaxt="n",
             xlab = "Strands information of two corresponding reads")
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)
dev.off()


## =================================================================== ##
## .bedpe.selected.dist.txt
## =================================================================== ##
# png(paste(output_dir, "Rplot10.png", sep = ""), width = 650, height = 470)
png(paste(output_dir, "Rplot10.png", sep = ""), unit = "in", width = 2.5, height = 1.85, res = 300, pointsize = 4)
# png(paste(output_dir, "Rplot10.png", sep = ""), width = 450, height = 350)
par(mar = c(5, 4, 4, 4))
draw_barplot("bedpe.selected.dist", "V2", col = fill_palette_3, border = NA, yaxt="n",  
             xlab = "Number of occurrences")
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)
dev.off()


## =================================================================== ##
## .bedpe.selected.unique.intra-chrom.strand.dist.txt
## =================================================================== ##
# png(paste(output_dir, "Rplot11.png", sep = ""), width = 650, height = 470)
png(paste(output_dir, "Rplot11.png", sep = ""), unit = "in", width = 2.5, height = 1.85, res = 300, pointsize = 4)
par(mar = c(5, 4, 4, 4))
draw_barplot("bedpe.selected.unique.intra-chrom.strand.dist", c("V2", "V3"), col = fill_palette_3, border = NA, yaxt="n", 
             xlab = "Number of occurrences")
axis(2, axTicks(4), unlist(format_axis_labels(axTicks(4))), las = 1)
dev.off()


## =================================================================== ##
## .bedpe.selected.intra-chrom.distance.minusplus.txt
## =================================================================== ##
fragment_len_stat <- scan(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.minusplus.txt", sep = ""))
fragment_len_stat <- as.numeric(fragment_len_stat)
fragment_len_stat <- fragment_len_stat[fragment_len_stat < 8000]
fragment_len_stat <- data.frame(mean(fragment_len_stat), sd(fragment_len_stat), median(fragment_len_stat), mad(fragment_len_stat, constant = 1))
colnames(fragment_len_stat) <- c("Mean", "Standard deviation", "Median", "Median absolute deviation")
fragment_len_stat <- xtable(fragment_len_stat, align = rep("c", 5))
print(fragment_len_stat, type='html', file=paste(output_dir, "Rplot17.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")





## =================================================================== ##
## .bedpe.selected.intra-chrom.distance.minusplus.txt
## =================================================================== ##
#### plus_plus
data_file <- read.table(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.plusplus.txt", sep = ""), header = F)
plusplus <- as.matrix(data_file)
#### plus_minus
data_file <- read.table(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.plusminus.txt", sep = ""), header = F)
plusminus <- as.matrix(data_file)
#### minus_plus
data_file <- read.table(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.minusplus.txt", sep = ""), header = F)
minusplus <- as.matrix(data_file)
#### minus_minus
data_file <- read.table(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.minusminus.txt", sep = ""), header = F)
minusminus <- as.matrix(data_file)

breaks <- c(-3E9, seq(-100000,100000,10), 3E9)
mids <- seq(-100000,100000-10,10) + 10/2

h.plusplus <- hist(plusplus,   breaks=breaks)$counts
h.plusplus <- h.plusplus[-length(h.plusplus)]
h.plusminus <- hist(plusminus,  breaks=breaks)$counts
h.plusminus <- h.plusminus[-length(h.plusminus)]
h.minusplus <- hist(minusplus,  breaks=breaks)$counts
h.minusplus <- h.minusplus[-length(h.minusplus)]
h.minusminus <- hist(minusminus, breaks=breaks)$counts
h.minusminus <- h.minusminus[-length(h.minusminus)]
max_count <- max(c(h.plusplus, h.plusminus, h.minusplus, h.minusminus))

write.table(h.plusplus,   file="h.plusplus.binsize_10bp.txt")
write.table(h.plusminus,  file="h.plusminus.binsize_10bp.txt")
write.table(h.minusplus,  file="h.minusplus.binsize_10bp.txt")
write.table(h.minusminus, file="h.minusminus.binsize_10bp.txt")

distribution_diff <- h.minusplus - (h.plusplus + h.plusminus + h.minusminus)/3
write.table(distribution_diff, file="h.diff.binsize_10bp.txt")


# hist_plot_file <- sprintf("%s.hist.binsize_10bp.png", filePrefix)
#png(paste(output_dir, "Rplot11.png", sep = ""), width = 960, height = 480)
#plot(  mids, h.plusplus[2:20001], xlim=c(-100000,100000), ylim=c(1,max_count), xlab="pet span (bp)", ylab="frequency of pets", col="red", log="y")
#points(mids, h.plusminus[2:20001],  col="green")
#points(mids, h.minusplus[2:20001],  col="blue")
#points(mids, h.minusminus[2:20001], col="black")
#legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), col=c("red", "green", "blue", "black"))
#dev.off()

indexes_1 <- 10001:20001
logMids <- log10(mids[indexes_1]+10/2)

# hist_plot_file.loglog <- sprintf("%s.hist.loglog.binsize_10bp.png", filePrefix)
png(paste(output_dir, "Rplot12.png", sep = ""), width = 680, height = 420)
indexes <- (h.plusplus[indexes_1] > 0)
plot(  logMids[indexes], (h.plusplus[indexes_1])[indexes],  col="red", xlim=c(1.5,5), ylim=c(1,max_count), xlab="PET span (bp)", ylab="Frequency of PETs", xaxt="n", cex = 1)
indexes <- (h.plusminus[indexes_1] > 0)
points(logMids[indexes], (h.plusminus[indexes_1])[indexes], col="green")
indexes <- (h.minusplus[indexes_1] > 0)
points(logMids[indexes], (h.minusplus[indexes_1])[indexes], col="blue")
indexes <- (h.minusminus[indexes_1] > 0)
points(logMids[indexes], (h.minusminus[indexes_1])[indexes], col="black")
axis(1, at=seq(-5,5,1), lab=c(-100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000))
legend("topright", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), text.col=c("red", "green", "blue", "black"), border = NA, cex = 1)
for (i in -4:4) {
    abline(v=i, lty="dashed")
}
dev.off()

indexes_1 <- 10001:20001
logMids <- log10(mids[indexes_1]+10/2)




# png(paste(output_dir, "Rplot24.png", sep = ""), unit = "in", width = 7, height = 5, res = 300, pointsize = 4)
# indexes <- (distribution_diff[indexes_1] > 0)
# plot(  logMids[indexes], (distribution_diff[indexes_1])[indexes],  col="red", xlim=c(1.5,5), ylim=c(1,max_count), xlab="Pet span (bp)", ylab="frequency of pets", xaxt="n", log="y")
# axis(1, at=seq(-5,5,1), lab=c(-100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000))
# legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), text.col=c("red", "green", "blue", "black"))
# for (i in -4:4) {
#     abline(v=i, lty="dashed")
# }
# #abline(v=log10(8000), lty="dashed")
# dev.off()



png(paste(output_dir, "Rplot16.png", sep = ""), width = 680, height = 420)
indexes <- (distribution_diff[indexes_1] > 0)
plot(  logMids[indexes], (distribution_diff[indexes_1])[indexes],  col="#7E827A", xlim=c(1.5,5), xlab="PET span (bp)", ylab="Difference of PET frequencies", xaxt="n", cex = 1)
axis(1, at=seq(1,5,1), lab=c(10, 100, 1000, 10000, 100000))
#axis(1, at=seq(-5,5,1), lab=c(-100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000))
#legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), text.col=c("red", "green", "blue", "black"))
for (i in -4:4) {
    abline(v=i, lty="dashed")
}
abline(h = 0, lty = "dashed")
#abline(v=log10(8000), lty="dashed")
dev.off()




# hist_diff_plot_file.loglog <- sprintf("%s.hist_diff.loglog.binsize_10bp.png", filePrefix)
png(paste(output_dir, "Rplot24.png", sep = ""), width = 680, height = 420)
indexes <- (distribution_diff[indexes_1] > 0)
plot(  logMids[indexes], (distribution_diff[indexes_1])[indexes],  col="#7E827A", xlim=c(1.5,5), ylim=c(1,max_count), xlab="PETs span (bp)", ylab="Frequency of PETs", xaxt="n", log="y", cex = 1)
axis(1, at=seq(-5,5,1), lab=c(-100000, -10000, -1000, -100, -10, 0, 10, 100, 1000, 10000, 100000))
#legend("topleft", c("plus,plus", "plus,minus", "minus,plus", "minus,minus"), text.col=c("red", "green", "blue", "black"))
for (i in -4:4) {
    abline(v=i, lty="dashed")
}
#abline(v=log10(8000), lty="dashed")
dev.off()





## =================================================================== ##
## .bedpe.selected.intra-chrom.distance.gap.txt
## =================================================================== ##
# png(paste(output_dir, "Rplot16.png", sep = ""), unit = "in", width = 3.5, height = 2.5, res = 300, pointsize = 4)
# dis_gap <- read.table(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.gap.txt", sep = ""), header = F)
# dis_gap <- dis_gap[dis_gap$V1 > 0, ]
# #dis_gap <- dis_gap[, dis_gap$V2 > 0]
# plot(log10(dis_gap$V1), dis_gap$V2, xlab = "log10(Span)", ylab = "Counts", cex = 0.5, lwd = 0.5, xlim = c(0, 6))
# abline(h = 0, lwd = 0.5, col = "red")
# dev.off()



## =================================================================== ##
## .bedpe.selected.intra-chrom.distance.gap.txt
## =================================================================== ##
# png(paste(output_dir, "Rplot24.png", sep = ""), unit = "in", width = 3.5, height = 2.5, res = 300, pointsize = 4)
# dis_gap <- read.table(paste(input_prefix, ".bedpe.selected.intra-chrom.distance.gap.txt", sep = ""), header = F)
# dis_gap <- dis_gap[dis_gap$V1 > 0 & dis_gap$V2 >0, ]
# #print(dis_gap)
# #dis_gap <- dis_gap[dis_gap$V2 > 0, ]
# plot(log10(dis_gap$V1), log10(dis_gap$V2), xlab = "log10(Span)", ylab = "log10(Counts)", cex = 0.5, lwd = 0.5)
# dev.off()



## =================================================================== ##
## .PET_count_distribution.txt
## =================================================================== ##
pet_count <- read.table(paste(input_prefix, ".PET_count_distribution.txt", sep = ""), sep = "\t", header = T)
#colnames(mapping_info) <- c("Title", "Information")
pet_count[[5]] <- round(pet_count[[5]], 4)
pet_count[[5]] <- paste(pet_count[[5]] * 100, "%", sep = "")

# get the PET conuts that have number of clusters more than 10000
pet_count[[2]] <- as.numeric(pet_count[[2]])
pet_10000 <- as.character(pet_count[-10, ][which(pet_count[[2]] > 20000), ][[1]])
pet_count[[2]] <- as.character(pet_count[[2]])
colnames(pet_count) <- c("PET.counts", "No. of clusters", "No intra chrom clusters", "No inter chrom clusters", "Percent of intra chrom clusters")
pet_count_tb <- xtable(pet_count, align = rep("c", 6))
print(pet_count_tb, type='html', file=paste(output_dir, "Rplot30.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")



## =================================================================== ##
## Basic statistics of interactions
## =================================================================== ##
interactions <- read.table(paste(input_prefix, ".cluster.FDRfiltered.txt", sep = ""))
interactions_stat <- interactions[, 1:7]
interactions_stat[, c(1, 4)] <- sapply(interactions_stat[, c(1, 4)], as.character)
interactions_stat[, c(2, 3, 5, 6, 7)] <- sapply(interactions_stat[, c(2, 3, 5, 6, 7)], as.numeric)
interactions_same_chr <- interactions_stat[which(interactions_stat[[1]] == interactions_stat[[4]]), ]
interactions_same_dis <- apply(interactions_same_chr, 1, function(x) {x = as.numeric(x[c(2, 3, 5, 6)]) 
                                                                      a = abs((x[1] + x[2] - x[3] - x[4]) / 2)})
interactions_same_dis_type <- sapply(interactions_same_dis / 1000, function(x) {x = as.numeric(x)
                                                                                if (x < 100) {a = "< 100Kb"}
                                                                                else if (x >= 100 & x < 1000) {a = "[100Kb, 1Mb)"}
                                                                                else if (x >= 1000 & x < 10000) {a = "[1Mb, 10Mb)"}
                                                                                else {a = "> 10Mb"}})

interactions_same_dis_type <- as.data.frame(table(interactions_same_dis_type))
colnames(interactions_same_dis_type) <- c("Distance", "Frequency")
interactions_same_dis_type[["Interaction_type"]] <- "Intra-chromosomal"
# interactions_same_dis_type <- interactions_same_dis_type[c(2, 1, 3, 4), c(3, 1, 2)]
# interactions_same_dis_type[5, ] <- c("Inter-chromosomal", "NOT AVAILABLE", nrow(interactions_stat) - nrow(interactions_same_chr))
#interactions_same_dis_type <- rbind(interactions_same_dis_type, c("Inter-chromosomal", "NOT AVAILABLE", nrow(interactions_stat) - nrow(interactions_same_chr)))
#for (i in seq_along(c("< 1Kb", "[1Kb, 100Kb)", "[100Kb, 1Mb)", "[1Mb, 10Mb)"))) {
#    if (!(i %in% interactions_same_dis_type$Distance)) {
#         interactions_same_dis_type <- rbind(interactions_same_dis_type, c("Inter-chromosomal", i, 0))              
#         }
#}
#print(c("Not available", nrow(interactions_stat) - nrow(interactions_same_chr), "Inter-chromosomal"))
intra_chrom_interactions <- data.frame(Distance = "Different chromosomes", Frequency = nrow(interactions_stat) - nrow(interactions_same_chr), Interaction_type = "Inter-chromosomal")
#names(intra_chrom_interactions) <- c("Distance", "Frequency", "Interaction type")
#intra_chrom_interactions <- as.data.frame(intra_chrom_interactions)
#print(intra_chrom_interactions)
interactions_same_dis_type <- rbind(interactions_same_dis_type, intra_chrom_interactions)
#print(interactions_same_dis_type)
#print(ncol(interactions_same_dis_type))
for (j in c("< 100Kb", "[100Kb, 1Mb)", "[1Mb, 10Mb)", "> 10Mb")) {
    if (!(j %in% interactions_same_dis_type$Distance)) {
         temp <- data.frame(Distance = j, Frequency = 0, Interaction_type = "Intra-chromosomal")
         interactions_same_dis_type <- rbind(interactions_same_dis_type, temp)
         }
}
#print(interactions_same_dis_type)
# mv ${DATA_FOLDER}/${OUTPUT_PREFIX}.hist.loglog.png ${DATA_FOLDER}/ChIA-PET_Tool_Report/images/Rplot12.png
order <- c("< 100Kb", "[100Kb, 1Mb)", "[1Mb, 10Mb)", "> 10Mb", "Different chromosomes")
#interactions_same_dis_type$sort <- factor(interactions_same_dis_type$Distance, levels = order)
interactions_same_dis_type <- interactions_same_dis_type[match(order, interactions_same_dis_type$Distance), ]
#print(interactions_same_dis_type)
#interactions_same_dis_type$Frequency <- sapply(interactions_same_dis_type$Frequency, function(q){round(as.numeric(q), 0)}) 
interactions_same_dis_type$Frequency <- as.character(interactions_same_dis_type$Frequency)
#interactions_same_dis_type$Frequency <- sapply(interactions_same_dis_type$Frequency, function(q){round(as.numeric(q), 0)}) 
colnames(interactions_same_dis_type) <- c("Distance", "Frequency", "Interaction type")
interactions_same_dis_type_tb <- xtable(interactions_same_dis_type, align = rep("c", 4))
print(interactions_same_dis_type_tb, type='html', file=paste(output_dir, "Rplot13.html", sep = ""), include.rownames = F, append = T,
      html.table.attributes = "class = 'mytable'")

# png(paste(output_dir, "Rplot13.png", sep = ""), width = 650, height = 250, res = 90)
# grid.table(interactions_same_dis_type, gpar.corefill = gpar(fill = "#E6E6E6", col = NA), gpar.colfill = gpar(fill = NA, col = NA),
#            gpar.coltext = gpar(col = "#585858", cex = 1.3), gpar.coretext = gpar(col = "#585858"),
#            h.odd.alpha = 0, h.even.alpha = 0.7, v.odd.alpha = 1, v.even.alpha = 1,
#            show.rownames = FALSE, padding.v = unit(8, "mm"))
# dev.off()


## =================================================================== ##
## Plot intra chromosomal interaction
## =================================================================== ##
# h19 cytoband data downlaod
# http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
# mm10 cytoband data download
# http://hgdownload.cse.ucsc.edu/goldenpath/mm10/database/
interactions_same_chr <- interactions_same_chr[which(!(interactions_same_chr[[7]] %in% pet_10000)), ]
clusters <- interactions_same_chr
peaks <- read.table(paste(input_prefix, ".peak.FDRfiltered.txt", sep = ""))
cyto <- read.table(cytoband_file_path, fill = T, header = T)
peaks[[4]] <- (peaks[[4]] - min(peaks[[4]]))/(max(peaks[[4]]) - min(peaks[[4]]))

# filter cytoband data
cyto <- cyto[nchar(as.character(cyto[[1]])) <= 5, ]
layout_nrow <- length(unique(cyto[[1]])) - 1
png(paste(output_dir, "Rplot14.png", sep = ""), width = 800, height = 1000)
# pdf(paste(output_dir, "Rplot14.pdf", sep = ""), width = 10, height = 15)
plot_layout(layout_nrow * 3, 2, heights = unit(rep(c(1.3, 0.4, 1.3)/(layout_nrow * 3), layout_nrow), "npc"), widths = unit(c(0.1, 0.9), "npc"))
plot_intra_chr_interaction(clusters, cyto, peaks)
dev.off()


## =================================================================== ##
## Plot circular view
## =================================================================== ##
# Set plotting parameters
interactions_stat <- interactions_stat[which(!(interactions_stat[[7]] %in% pet_10000)), ]
interactions_stat <- plot_intra_chr_interaction(interactions_stat, cyto, peaks, return = TRUE)
interactions_diff <- interactions_stat[which(as.character(interactions_stat[[1]]) != as.character(interactions_stat[[4]])), ]
# interactions  <- interactions[1:1000, ]

# data(UCSC.HG19.Human.CytoBandIdeogram)
# data(UCSC.Mouse.GRCm38.CytoBandIdeogram)


if (species == 1) {
    data(UCSC.HG19.Human.CytoBandIdeogram)
    cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
} else if (species == 2) {
    data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
    cyto.info <- UCSC.Mouse.GRCm38.CytoBandIdeogram
} else if (species == 3) {
    cyto.info <- cyto
}
# cyto.info <- UCSC.Mouse.GRCm38.CytoBandIdeogram
# cyto.info <- UCSC.HG19.Human.CytoBandIdeogram

cyto_color <- data.frame(Stain = c("gpos100", "gpos75", "gpos66", "gpos50", "gpos33", "gpos25", "gneg", "acen", "gvar", "stalk"),
                         BandColor = c("#4D4D4D", "#969696", "#AAAAAA", "#C3C3C3", "#C2C2C2", "#E6E6E6", "white", "white", "#F3E6C4", "#5A81AA"))


chr.exclude <- NULL
tracks.inside <- 2
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)

params <- RCircos.Get.Plot.Parameters()
params$radius.len <- 4
params$chr.ideog.pos <- 4.3
params$chr.name.po <- 4.8
params$plot.radius <- 5
params$chrom.width <- 0.2
params$track.in.start <- 3
params$text.size <- 20
params$track.padding <- 0.3
params$track.height <- 0.6

#params$hist.width <- 2000
params$track.background <- NA
params$grid.line.color <- NA
params$hist.color <- "#77C4D2"
#params$sub.tracks <- 1
RCircos.Reset.Plot.Parameters(params)


png(paste(output_dir, "Rplot15.png", sep = ""), unit = "in", width = 7, height = 7, res = 350, pointsize = 4)
# RCircos.Set.Plot.Area()
#pdf("Rplot16.pdf", width = 8, height = 8)
par(mai=c(0.01, 0.01, 0.01, 0.01), cex = 2)
plot.new()
plot.window(c(-5,5), c(-5, 5))
RCircos.Chromosome.Ideogram.Plot()
# data(RCircos.Ribbon.Data);
# RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data,
#                       track.num=11, by.chromosome=FALSE, twist=FALSE);

RCircos.Histogram.Plot(peaks, 4, 1, "in")
if (nrow(interactions_diff) == 0) {
  print ("There is no inter-chromosomal interaction.")
  dev.off()
  } else {
    RCircos.Link.Plot(interactions_diff, 2, FALSE)
# data(RCircos.Ribbon.Data)
# RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data,
#                      + track.num=11, by.chromosome=FALSE, twist=FALSE)
  dev.off()
}

