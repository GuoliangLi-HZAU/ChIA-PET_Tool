# =================================================================== ##
## Preprocess data and plot bar graph 
## =================================================================== ##
# import multiple files into a list orderly
# merge 1_1 1_2 2_1 2_2 objects and draw barplot
draw_barplot <- function(pattern, merge_col, ...) {
    
    # Import file in right order
    files_order <- c("1_1", "2_2", "1_2", "2_1")
    file_names <- paste(input_prefix, files_order, pattern, "txt", sep = ".")
    data_list <- lapply(file_names, read.table, header = F)
    data_merged <- Reduce(function(x, y) merge(x, y, by = merge_col, all = TRUE), data_list)
    data_merged[is.na(data_merged)] <- 0
    
    # Merge the first two columns of "Strand distribution" file
    if(length(merge_col) == 2) {
        data_merged[, merge_col[1]] <- paste(data_merged[, merge_col[1]], data_merged[, merge_col[2]], sep = " ")
        data_merged <- data_merged[!(names(data_merged) %in% merge_col[2])]
    }
    colnames(data_merged) <- c("value", "A_A", "B_B", "A_B", "B_A")
    
    # Order rows
    if (is.numeric(data_merged[[1]])) {
        data_merged <- data_merged[order(data_merged$value), ]
    } else {
        value_order <- c("+ +", "+ -", "- -", "- +", "+ .", "- .", ". +", ". -", ". .")
        data_merged[[1]] <- factor(data_merged[[1]], levels = value_order)
        data_merged <- data_merged[order(data_merged[[1]]), ]
    }
    
    # Set cex according to number of rows
    if (nrow(data_merged) >= 20) {
        cex.names = 0.6
    } else if (nrow(data_merged) < 10) {
        cex.names = 1.5
    } else { cex.names = 1 }
    
    # Plot bar graph
    barplot(t(data_merged[, -1]), legend.text = T, args.legend = list(x = "topright", inset=c(-0.08, 0), cex = 1, bty = "n", border = NA),
            names.arg = data_merged[, 1], cex.names = cex.names, ...)
}


## =================================================================== ##
## Plot back to back barplot 
## =================================================================== ##
draw_back_to_back <- function(data, left_col, right_col, y_label) {
    
    # Convert variable to the proper type
    data <- as.data.frame(sapply(data, as.numeric))
    score_max<- max(data[[2]], data[[3]])
    n <- nrow(data)
    
    # Set layout for plotting elements
    grid.newpage()
    vp_layout <- viewport(height = unit(0.9, "npc"), width = unit(0.9, "npc"), 
                          layout = grid.layout(nrow = 1, ncol = 2, width = unit(c(0.1, 0.9), "npc")))
    pushViewport(vp_layout)
    
    # Plot back to back horizental bar plot 
    vp_barplot <- viewport(layout.pos.col = 2, xscale = c(-score_max, score_max), yscale = c(-n / 10, n))
    pushViewport(vp_barplot)
    grid.lines(x = unit(rep(0, 2), "native"), y = unit(c(-n / 10, n), "native"), gp = gpar(col = "grey"))
    grid.segments(x0 = rep(-score_max, n), y0 = data[[1]] + 0.4, x1 = rep(score_max, n), y1 = data[[1]] + 0.4, default.units = "native",
                  gp = gpar(lty = "dotted"))  
    grid.rect(x = -data[[2]], y = data[[1]], width = data[[2]], height = rep(0.8, n), default.units = "native", 
              just = c("left", "bottom"), gp = gpar(fill = left_col, col = NA))
    grid.rect(x = rep(0, n), y = data[[1]], width = data[[3]], height = rep(0.8, n), default.units = "native", 
              just = c("left", "bottom"), gp = gpar(fill = right_col, col = NA))
    
    # Plot text boxes in the bottom
    box_width = 1.2*stringWidth("Right end")
    box_height = 1.5*stringHeight("Right end")
    grid.rect(x = 0, y = unit(-n / 20, "native") - box_height * 0.25, width = box_width, height = box_height, default.units = "native",
              just = c("right", "bottom"), gp = gpar(fill = left_col, col = NA))
    grid.rect(x = 0, y = unit(-n / 20, "native") - box_height * 0.25, width = box_width, height = box_height, default.units = "native",
              just = c("left", "bottom"), gp = gpar(fill = right_col, col = NA))
    grid.text(label = "Read one", x = 0, y = -n / 20, default.units = "native", 
              just = c("right", "bottom"), gp = gpar(col = "white", cex = 1))
    grid.text(label = "Read two", x = 0, y = -n / 20, default.units = "native", 
              just = c("left", "bottom"), gp = gpar(col = "white", cex = 1))
    
    # Plot xaxis
    ticks <- round(seq(-score_max, score_max, length.out = 11), digits = 0)
    grid.xaxis(at = ticks, label = abs(ticks), gp = gpar(lineheight = 0.8, cex = 0.6))
    upViewport()
    
    # Plot score on the left 
    vp_text <- viewport(layout.pos.col = 1, yscale = c(-n / 10, n))
    pushViewport(vp_text)
    draw_text <- function(x) {
        grid.text(label = x[1], x = 0.5, y = unit(x[1] + 0.25, "native"), just = c("right", "bottom"))
    }
    apply(data, 1, draw_text)
    grid.text(label = y_label, x = 0.5, y = unit(-n / 20, "native"), just = c("center", "bottom"))
    upViewport()
}


## =================================================================== ##
## Set layout for plotting
## =================================================================== ##
plot_layout <- function(nrow, ncol, heights = NULL, widths = NULL) {
    
    # Function to test if sth is a interger
    is.wholenumber <-
        function(x, tol = .Machine$double.eps ^ 0.5)  abs(x - round(x)) < tol
    
    # Check if nrow or ncol are "null unit"
    if (is.null(heights)) {
        heights = unit(rep_len(1, nrow), "null")
    } 
    if (is.null(widths)) {
        widths = unit(rep_len(1, ncol), "null")
    }
    
    # Set layput for later plotting
    grid.newpage()
    layout <- grid.layout(nrow = nrow, ncol = ncol, 
                          heights = heights,
                          widths = widths)
    vp_overall <- viewport(name = "vp_overall", width = unit(0.95, "npc"), height = unit(0.95, "npc"), layout = layout)
    pushViewport(vp_overall)
}


## =================================================================== ##
## function to draw intra chromosomal interactions
## =================================================================== ##

plot_intra_chr_interaction <- function(interaction, cyto, peaks, return = FALSE) {
    
    # Remove interactions on chrM
    interaction <- interaction[as.character(interaction[[1]]) != "chrM" & as.character(interaction[[4]]) != "chrM", ]
    peaks <- peaks[peaks[[1]] != "chrM", ]
    cyto <- cyto[cyto[[1]] != "chrM", ]
    
    # Set color scheme for each "stain" in cytoband data
    colnames(cyto) <- c("chrom", "start", "end", "name", "stain")
    cyto_color <- data.frame(stain = c("gpos100", "gpos75", "gpos66", "gpos50", "gpos33", "gpos25", "gneg", "acen", "gvar", "stalk"),
                             col = c("#4D4D4D", "#969696", "#AAAAAA", "#C3C3C3", "#C2C2C2", "#E6E6E6", "white", "white", "#F3E6C4", "#5A81AA"))
    cyto <- merge(cyto, cyto_color, by = "stain", all.x = T)
    
    # Order chromosome from chr1 to chrY
    idx <- grepl("[0-9]+", cyto$chrom)
    chrom_no <- paste("chr", unique(c(sort(as.numeric(gsub("chr", "", cyto$chrom[idx]))), gsub("chr", "", sort(cyto$chrom[!idx])))), sep = "")
    
    # Whether return color columns for RCiscos plot
    # interaction <- interaction[, 1:7]  
    interaction$PlotColor <- apply(interaction, 1, function(x){x[7] <- as.numeric(x[7]); if(x[7] == 2){a = "#EF667A"}
                                                               else if(x[7] == 3){a = "#E92F4A"}
                                                               else if(x[7] == 4 | x[7] == 5){a = "#E92F1E"}
                                                               else if(x[7] >=6 & x[7] <= 10){a = "#AE2F4A"}
                                                               else if(x[7] >= 11 & x[7] <= 20){a = "#7C2F4A"}
                                                               else(a = "#5A2F63")})
    
    if (return == TRUE) {
        #col_div <- split(chrom_no, ceiling(seq_along(chrom_no)/4))
#         names(col_div) <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
#                             "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")[1:length(col_div)]
#        names(col_div) <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
#                            "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")[1:length(col_div)]
        #names(col_div) <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
        #                     "#66A61E", "#E6AB02", "#A6761D", "#666666")[1:length(col_div)]
        ##freq <- unique(as.numeric(interaction[[7]]))[order(unique(as.numeric(interaction[[7]])))]
        ##color = colorRampPalette(c("#E92F4A", "#F499A6"))(length(freq)) 
        #names(col_1) <- color
        ##freq_color <- data.frame(PlotColor = color, freq = freq)

#         interaction$PlotColor <- apply(interaction, 1, function(x){if((x[1] %in% col_div[[1]]) | (x[4] %in% col_div[[1]])) {a = names(col_div)[1]}
#                                                                    else if((x[1] %in% col_div[[2]]) | (x[4] %in% col_div[[2]])) {a = names(col_div)[2]}
#                                                                    else if((x[1] %in% col_div[[3]]) | (x[4] %in% col_div[[3]])) {a = names(col_div)[3]}
#                                                                    else if((x[1] %in% col_div[[4]]) | (x[4] %in% col_div[[4]])) {a = names(col_div)[4]}
#                                                                    else if((x[1] %in% col_div[[5]]) | (x[4] %in% col_div[[5]])) {a = names(col_div)[5]}
#                                                                    else {a = names(col_div)[6]}
#                                                                   }
#                                        )
        

        #interaction$PlotColor <- apply(interaction, 1, function(x){a = col_1[which(names(col_1) == x[7])]})
        ## interaction <- merge(interaction, freq_color, by.x = "V7", by.y = "freq", all.x = T)[, c("V1", "V2", "V3", "V4", "V5", "V6", "PlotColor")]
        #interaction$PlotColor <- apply(interaction, 1, function(x){freq_color[which(freq_color$freq == as.character(x[7])), ][["color"]]})
#interaction <- interaction[, c(1:6, 8)]
        return(interaction)
    } else {
      #  interaction <- interaction[, c(1:6, 8)]
        # Order cytoband data
        cyto$chrom_no_f <- factor(cyto$chrom, levels = chrom_no)
        cyto <- cyto[with(cyto, order(chrom_no_f, start)), ]
        # Seperate centromere data from cytoband data
        centro <- cyto[cyto$stain == "acen", ]
        centro <- centro[, c("chrom", "start", "end")]
        cyto <- cyto[cyto$stain != "acen", ]
        cyto <- cyto[, c("chrom", "start", "end", "col")]
        
        # Convert variables to correct type
       # interaction <- interaction[which(as.character(interaction[[1]]) == as.character(interaction[[4]])), ]
        cyto[, c(1, 4)] <- sapply(cyto[, c(1, 4)], as.character)
        cyto[, 2:3] <- sapply(cyto[, 2:3], as.numeric)
        centro[, 1] <- as.character(centro[, 1])
        centro[, 2:3] <- sapply(centro[, 2:3], as.numeric)
        
        peaks[, 1] <- as.character(peaks[, 1])
        peaks[, 2:4] <- sapply(peaks[, 2:4], as.numeric)
        
        # Calculate max coordinate for each chromosome and the max peak signal
        cyto_max <- max(cyto[["end"]])
        max_score <- max(peaks[[4]])
        
        # Draw ideogram, interaction curves and peaks signals for each chromosome
        for (i in seq.int(length(chrom_no))) { 
            # Get chromosome name
            if (grepl("X", chrom_no[i])) {
                chrom <- paste("chr", "X", sep = "")
            } else if (grepl("Y", chrom_no[i])) {
                chrom <- paste("chr", "Y", sep = "") 
            }  else { chrom <- paste("chr", i, sep = "") }
            
            # Filter data for ideogram
            chrom_cyto <- cyto[which(cyto[[1]] == chrom), ]
            chrom_cyto <- chrom_cyto[order(chrom_cyto[[2]]), ]
            chrom_max <- max(chrom_cyto[["end"]])
            chrom_centro <- centro[which(centro[["chrom"]] == chrom), ]
            
            # Filter data for peaks
            chrom_peaks <- peaks[which(peaks[[1]] == chrom), ]
            chrom_peaks <- chrom_peaks[order(chrom_peaks[[2]]), ]
            chrom_peaks[[4]] <- chrom_peaks[[4]] / max_score * 10
            chrom_peaks[[4]][chrom_peaks[[4]] > 1] <- 1
            # chrom_peaks[[4]][chrom_peaks[[4]] > quantile(chrom_peaks[[4]], 0.95)] <- max_score
            
            
            # Filter data for inetraction curves
            chrom_interaction <- interaction[which(interaction[[1]] == chrom), ]
            chrom_interaction <- chrom_interaction[, c(2, 3, 5, 6, 8)]
            #print(head(chrom_interaction))
            
            # Plot chromosome name
            vp_chr = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 1)
            grid.text(chrom, x = unit(0.5, "npc"), y = unit(0, "npc"),
                      just = c("right", "bottom"), gp = gpar(col = "#323745", fontsize = 15) , vp = vp_chr) 
            
            # Plot chromosome ideogram
            vp = viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp", xscale = c(0, cyto_max))

            
            # Cytoband information
         #   chrom_cyto <- rectGrob(x = unit(chrom_cyto[["start"]], "native"), y = unit(0.07, "npc"), width = unit((chrom_cyto[["end"]] - chrom_cyto[["start"]]) * 0.95, "native"), height = unit(0.86, "npc"),
         #                          just = c("left", "bottom"), gp = gpar(fill = chrom_cyto[["col"]], col = NA), vp = vp)            
         #   grid.draw(chrom_cyto)
         #   
         #    
         #   if (nrow(chrom_centro) == 0) {
         #       chrom_rect <- roundrectGrob(x = unit(-0.0015, "npc"), y = unit(0, "npc"), width = unit(chrom_max, "native") + unit(0.0015, "npc"),
         #                                   r = unit(1.3, "mm"),
         #                                   just = c("left", "bottom"), gp = gpar(fill = NA, lex = 2.5, col = "#9F9F9F"), vp = vp)   
         #       grid.draw(chrom_rect)
         #   } else {
         #       # Left part
         #       chrom_left <- chrom_centro[1, ]
         #       chrom_left$start <- 0
         #       chrom_rect_left <- roundrectGrob(x = unit(-0.0015, "npc"), y = unit(0, "npc"), width = unit(chrom_left[["end"]] - chrom_left[["start"]], "native") + unit(0.0015, "npc"),
         #                                        r = unit(1.3, "mm"),
         #                                        just = c("left", "bottom"), gp = gpar(fill = NA, lex = 2.5, col = "#9F9F9F"), vp = vp)   
         #       grid.draw(chrom_rect_left)
         #       
         #       # Right part
         #       chrom_right <- chrom_centro[2, ]
         #       chrom_right$end <- chrom_max
         #       chrom_rect_right <- roundrectGrob(x = unit(chrom_right[["start"]], "native"), y = unit(0, "npc"), width = unit((chrom_right[["end"]] - chrom_right[["start"]]) * 1.05, "native"),
         #                                         r = unit(1.3, "mm"),
         #                                         just = c("left", "bottom"), gp = gpar(fill = NA, lex = 2.5, col = "#9F9F9F"), vp = vp)            
         #       grid.draw(chrom_rect_right)
         #   }           


        grid.rect(x = unit(chrom_cyto[["start"]], "native"), y = unit(1, "npc"), width = unit(chrom_cyto[["end"]] - chrom_cyto[["start"]], "native"), height = unit(1, "npc"),
                               just = c("left", "top"), gp = gpar(fill = chrom_cyto[["col"]], col = NA, col = "#9F9F9F"), vp = vp)    
               

        if (nrow(chrom_centro) == 0) {
            grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_max, "native"), height = unit(1, "npc"),
                                        just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   

        } else {
            # Left part
            chrom_left <- chrom_centro[1, ]
            chrom_left$start <- 0
            grid.rect(x = unit(0, "npc"), y = unit(1, "npc"), width = unit(chrom_left[["end"]] - chrom_left[["start"]], "native"), height = unit(1, "npc"),
                                             just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)   

            
            # Right part
            chrom_right <- chrom_centro[2, ]
            chrom_right$end <- chrom_max
            grid.rect(x = unit(chrom_right[["start"]], "native"), y = unit(1, "npc"), width = unit(chrom_right[["end"]] - chrom_right[["start"]], "native"), height = unit(1, "npc"),
                                              just = c("left", "top"), gp = gpar(fill = NA, col = "#9F9F9F"), vp = vp)            

        }


            # Plot interaction curves
#             draw_curves <- function(vector) {
#                 x_start = (vector[1] + vector[2]) / 2
#                 x_end = (vector[3] + vector[4]) / 2
#                 angle <- seq(0, pi, length=50)
#                 polygonGrob(x = unit(seq(x_start, x_end, length=50), "native"),
#                             y = unit(-sin(angle) , "npc"),
#                             gp = gpar(fill = "#E92F4A", col = NA, alpha = 0.3, lwd = 0.001)
#                 )
#             } 
    
            draw_Bezier_Curve <- function(vector) {
              x_start = (as.numeric(vector[1]) + as.numeric(vector[2])) / 2
              x_end = (as.numeric(vector[3]) + as.numeric(vector[4])) / 2
              x1 <- x_start; y1 <- -0.2;
              x4 <- x_end; y4 <- -0.2;
              if ((x4 - x1) <= 1000000) {
                  x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 50; y2 <- -3;
                  x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 50; y3 <- -3;
              } else if((x4 - x1) <= 10000000) {
                x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 5; y2 <- -3;
                x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 5; y3 <- -3;  
              } else {
                x2 <- (x_start + x_end)/2 + (x_start - x_end)/2 * 0.5; y2 <- -3;
                x3 <- (x_start + x_end)/2 - (x_start - x_end)/2 * 0.5; y3 <- -3;
              }
              bezierGrob(x = unit(c(x1, x2, x3, x4), "native"), y = c(y1, y2, y3, y4),
                         gp = gpar(col = as.character(vector[5]), fill = as.character(vector[5]), alpha = 0.5, lwd = 0.5))
            }
            
            curve <- apply(chrom_interaction, 1, draw_Bezier_Curve)
            curves <- gList()
            for (k in seq.int(curve)) {
                curves[[k]] <- gList(curve[[k]])
            }
            vp_interactions <- viewport(layout.pos.row = 3 * i - 1, layout.pos.col = 2, name = "vp_interactions", xscale = c(0, cyto_max), just = c("left", "top"))
            pushViewport(vp_interactions)
            grid.draw(curves)
            upViewport()
            
            #  Plot peaks
            vp_peaks <- viewport(layout.pos.row = 3 * i - 2, layout.pos.col = 2, name = "vp_peaks", xscale = c(0, cyto_max), yscale = c(0, max_score), just = c("left", "bottom"))
            if (nrow(chrom_peaks) == 0) {
                next
            } else {
                grid.rect(x = unit(chrom_peaks[[2]], "native"), y = unit(0.09, "npc"), width = unit(chrom_peaks[[3]] - chrom_peaks[[2]], "native"), height = unit(chrom_peaks[[4]], "npc"),
                          just = c("left", "bottom"), gp = gpar(fill = "#77C4D2", col = "#77C4D2"), vp = vp_peaks)
            }


        }
    }
}




## =================================================================== ##
## Overide RCircos.Chromosome.Ideogram.Plot function in RCircos package
## =================================================================== ##
RCircos.Chromosome.Ideogram.Plot <-function()
{
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
  RCircos.Pos <- RCircos.Get.Plot.Positions()
  RCircos.Par <- RCircos.Get.Plot.Parameters()
  
  
  #  Plot chromosome outlines, chromosome names, and 
  #	chromosome highlights
  #	_________________________________________________________
  #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  right.side <- nrow(RCircos.Pos)/2;
  outer.location <- RCircos.Par$chr.ideog.pos + RCircos.Par$chrom.width;
  inner.location <-RCircos.Par$chr.ideog.pos
  
  chroms <- unique(RCircos.Cyto$Chromosome);
  
  
  
  #	Add chromosome bands (Giema stain positive only)
  #	_________________________________________________________
  #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  RCircos.Cyto <- RCircos.Cyto[, -6]
  RCircos.Cyto <- merge(RCircos.Cyto, cyto_color, by = "Stain")
  RCircos.Cyto <- RCircos.Cyto[with(RCircos.Cyto, order(Chromosome, ChromStart)), ]
  RCircos.Cyto$BandColor <- as.character(RCircos.Cyto$BandColor)
  
  for(a.band in 1:nrow(RCircos.Cyto))
  {
    a.color <- RCircos.Cyto$BandColor[a.band];
    if(a.color=="white") { next; }
    
    start <- RCircos.Cyto$Location[a.band]-RCircos.Cyto$Unit[a.band]+1;
    end   <- RCircos.Cyto$Location[a.band];
    
    pos.x<- c(RCircos.Pos[start:end,1]*outer.location, 
              RCircos.Pos[end:start,1]*inner.location);
    pos.y<- c(RCircos.Pos[start:end,2]*outer.location, 
              RCircos.Pos[end:start,2]*inner.location);
    polygon(pos.x, pos.y, col=a.color, border=NA);
    
  }
  
  
  
  for(a.chr in 1:length(chroms))
  {
    the.chr  <- RCircos.Cyto[RCircos.Cyto$Chromosome==chroms[a.chr],];
    start <- the.chr$Location[1]- the.chr$Unit[1] + 1;
    end   <- the.chr$Location[nrow(the.chr)];
    mid <- round((end-start+1)/2, digits=0)+start;
    
    chr.color <- the.chr$ChrColor[nrow(the.chr)];
    
    #  Draw chromosome outlines
    #  _________________________________________________________
    #  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    pos.x<- c(RCircos.Pos[start:end,1]*outer.location, 
              RCircos.Pos[end:start,1]*inner.location);
    pos.y<- c(RCircos.Pos[start:end,2]*outer.location, 
              RCircos.Pos[end:start,2]*inner.location);
    polygon(pos.x, pos.y, lwd = 1.5, border = "#9F9F9F");
    
    
    #	Add chromosome names
    #	_________________________________________________________
    #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    chr.name <- sub(pattern="chr", replacement="", chroms[a.chr]);
    text(RCircos.Pos[mid,1]*RCircos.Par$chr.name.pos,
         RCircos.Pos[mid,2]*RCircos.Par$chr.name.pos,
         label=chr.name, srt=RCircos.Pos$degree[mid]);
    
    
    #	Add chromosome highlights
    #	_________________________________________________________
    #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    #    lines(RCircos.Pos[start:end,]*RCircos.Par$highlight.pos, 
    #          col=chr.color, lwd=RCircos.Par$highlight.width);
  }
}




## =================================================================== ##
## Overide RCircos.Histogram.Plot function in RCircos package
## =================================================================== ##
RCircos.Histogram.Plot<-function(hist.data, data.col, track.num, side)
{
  RCircos.Pos <- RCircos.Get.Plot.Positions();
  RCircos.Par <- RCircos.Get.Plot.Parameters();
  RCircos.Cyto <- RCircos.Get.Plot.Ideogram();
  
  
  #  Convert raw data to plot data. The raw data will be validated
  #	first during the convertion
  #	____________________________________________________________
  #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  hist.data <- RCircos.Get.Plot.Data(hist.data, "plot");
  
  
  #	Each heatmap cell is centered on data point location. Make
  #	sure each one will be in the range of it chromosome
  #	____________________________________________________________
  #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  hist.locations <- as.numeric(hist.data[, ncol(hist.data)]);
  start <- hist.locations - RCircos.Par$hist.width;
  end   <- hist.locations + RCircos.Par$hist.width;
  
  data.chroms <- as.character(hist.data[,1]);
  chromosomes <- unique(data.chroms);
  cyto.chroms <- as.character(RCircos.Cyto$Chromosome);
  
  for(a.chr in 1:length(chromosomes))
  {
    cyto.rows <- which(cyto.chroms==chromosomes[a.chr]);
    locations <- as.numeric(RCircos.Cyto$Location[cyto.rows]);
    chr.start <- min(locations) - RCircos.Cyto$Unit[cyto.rows[1]];
    chr.end   <- max(locations);
    
    data.rows <- which(data.chroms==chromosomes[a.chr]);
    start[data.rows[start[data.rows]<chr.start]] <- chr.start;
    end[data.rows[end[data.rows]>chr.end]] <- chr.end;
  }
  
  
  #	Plot position for current track. 
  #	___________________________________________________________
  #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  locations <- RCircos.Track.Positions(side, track.num);

  out.pos <- locations[1] - RCircos.Par$track.padding  /3;
  in.pos <- locations[2] - RCircos.Par$track.padding  /3;
  
  
  #       Label colors
  #       _________________________________________________________
  #       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  hist.colors <- RCircos.Get.Plot.Colors(hist.data, RCircos.Par$hist.color); 
  
  
  #	Draw histogram
  #	___________________________________________________________
  #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  num.subtrack <- RCircos.Par$sub.tracks;
  # RCircos.Track.Outline(out.pos, in.pos, RCircos.Par$sub.tracks);
  
  for(a.point in 1:nrow(hist.data))
  {
    hist.height <- hist.data[a.point, data.col];
    
    the.start <- start[a.point];
    the.end <- end[a.point];
    
    #	Plot rectangle with specific height for each 
    #	data point
    #	_______________________________________________
    #	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    
    height <- in.pos + RCircos.Par$track.height*hist.height;
    polygon.x<- c(RCircos.Pos[the.start:the.end,1]*height, 
                  RCircos.Pos[the.end:the.start,1]*in.pos);
    polygon.y<- c(RCircos.Pos[the.start:the.end,2]*height, 
                  RCircos.Pos[the.end:the.start,2]*in.pos);
    polygon(polygon.x, polygon.y, col=hist.colors[a.point], 
            border=NA);
  }
}