## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function(coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}


draw_matrix2 <- function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color) {
    
    n = nrow(matrix)
    m = ncol(matrix)
    
    coord_x = pheatmap:::find_coordinates(m, gaps_cols)
    coord_y = pheatmap:::find_coordinates(n, gaps_rows)
    
    x = coord_x$coord - 0.5 * coord_x$size
    y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)
    
    coord = expand.grid(y = y, x = x)
    
    res = gList()
    
    res[["rect"]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, col = border_color))
    
    if (attr(fmat, "draw")) {
        res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
    }
    
    if (!is.null(gaps_cols)) {
        for (i in 1:(length(gaps_cols) - 1)) {
            res[[paste0("line", as.character(i))]] <- linesGrob(x = rep(x[gaps_cols[i] + 1], 2) - unit(1.5, "bigpts"), y = unit(c(0, 1), "npc"), 
                gp = gpar(col = "black", lwd = 1, lty = 2))
        }
    }
    res = gTree(children = res)
    
    return(res)
}


plotHeatmap <- function(object, fileName, cnv.scale = 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T) {
    
    assignInNamespace(x = "draw_matrix", value = draw_matrix2, ns = asNamespace("pheatmap"))
    assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
    
    breaks <- seq(-2, 2, by = 0.2)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    
    idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:22)])
    xlabel <- rep("", length(rownames(object@data)))
    half <- round(table(object@annotation.filt$Chr)[as.character(1:22)]/2)[-1]
    xpos <- c(half[1], (idx[-22] + half))
    xlabel[xpos] <- 1:22
    
    data <- object@control.normalized[[cnv.scale]]
    
    if (only_soi) 
        data <- data[, !(colnames(data) %in% object@control.sample.ids)]
    
    pheatmap(t(data), cluster_cols = cluster_cols, cluster_rows = cluster_rows, gaps_col = idx, color = color, 
        breaks = breaks, labels_col = xlabel, show_rownames = show_rownames, filename = fileName)
    
}

## practical for around 20 samples
plotGEAndGT <- function(chrMat, genoMat, fileName) {
    rownames(genoMat) <- paste0(rownames(genoMat), "_GT")
    rownames(chrMat) <- paste0(rownames(chrMat), "_RNASeq")
    
    chrMatAll <- rbind(chrMat, genoMat)
    plot.data <- melt(chrMatAll)
    
    plot.data$value2 <- "neutral"
    plot.data$value2[plot.data$value > 0] <- "amplification"
    plot.data$value2[plot.data$value < 0] <- "deletion"
    plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
    plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
    plot.data$sampleType <- gsub("_RNASeq", "", gsub("_GT", "", plot.data$X1))
    
    p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
        deletion = muted("blue"), neutral = "white")) + facet_wrap(~plot.data$sampleType, scales = "free_y", ncol = 1) + theme_grey(base_size = 6) + 
        theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank(), strip.text.x = element_blank(), legend.text = element_text(colour = "black", 
            size = 7, face = "bold"), legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 6, 
            colour = "black", angle = -45, hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
            plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
    ggsave(fileName, plot = p)
    
}

plotLargeScaleEvent2 <- function(chrMat, fileName) {
    

    plot.data <- melt(chrMat)
    
    plot.data$value2 <- "neutral"
    plot.data$value2[plot.data$value > 0] <- "amplification"
    plot.data$value2[plot.data$value < 0] <- "deletion"
    plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
    plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
    
    p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
        deletion = muted("blue"), neutral = "white")) + theme_grey(base_size = 6) + theme(legend.position = "right", legend.direction = "vertical", 
        legend.title = element_blank(), strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, face = "bold"), 
        legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 5, colour = "black", 
            angle = -45, hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
    ggsave(fileName, plot = p)
    
}



plotLargeScaleEvent <- function(object, fileName) {
    
    samps <- object@large.scale.cnv.events
    chrs <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
    chrMat <- matrix(0, ncol = 44, nrow = length(rownames(samps)))
    colnames(chrMat) <- chrs
    rownames(chrMat) <- rownames(samps)
    
    
    for (x in 1:dim(samps)[1]) {
        
        chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleAmp[x], collapse = " ")), split = " ")))
        chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- 1
        
        chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(samps$LargeScaleDel[x], collapse = " ")), split = " ")))
        chrMat[x, match(intersect(chrWithEvents, chrs), colnames(chrMat))] <- (-1)
        
    }
    
    plot.data <- melt(chrMat)
    
    plot.data$value2 <- "neutral"
    plot.data$value2[plot.data$value > 0] <- "amplification"
    plot.data$value2[plot.data$value < 0] <- "deletion"
    plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
    plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
    
    p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", y = "") + scale_fill_manual(values = c(amplification = muted("red"), 
        deletion = muted("blue"), neutral = "white")) + theme_grey(base_size = 6) + theme(legend.position = "right", legend.direction = "vertical", 
        legend.title = element_blank(), strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, face = "bold"), 
        legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 5, colour = "black", 
            angle = -45, hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
    ggsave(fileName, plot = p)
    
}


plotBAFInSeperatePages <- function(loh, folderName) {
    
    dir.create(folderName)
    for (i in 1:length(names(loh))) {
        sample_name <- names(loh)[i]
        maf <- loh[[i]]
        
        all_Chr <- NULL
        all_BAF_D <- NULL
        
        overlays1 <- list()
        textOverlays <- list()
        max_num <- 1
        
        for (m in 1:22) {
            
            chrBAF <- maf[which(as.character(maf$chr) == as.character(m)), ]
            chrBAF <- chrBAF[order(as.numeric(as.character(chrBAF$pos))), ]
            
            all_Chr <- c(all_Chr, as.numeric(as.character(chrBAF$chr)))
            all_BAF_D <- c(all_BAF_D, (as.numeric(as.character(chrBAF$dev))))
            
            max_num <- (length(chrBAF$chr) + max_num)
            hr2 <- makeRectangleOverlay(start = max_num, end = max_num, dp = DisplayPars(lty = "dashed", lwd = 1, color = "grey"))
            overlays1 = c(overlays1, list(hr2))
            
            if ((m%%2) == 0) {
                hr1 <- makeTextOverlay(as.character(m), xpos = max_num - (length(chrBAF$chr)/2), ypos = 0.8, coords = "genomic", dp = DisplayPars(pointSize = 1, 
                  cex = 0.7))
                textOverlays = c(textOverlays, list(hr1))
            }
            if ((m%%2) == 1) {
                hr1 <- makeTextOverlay(as.character(m), xpos = max_num - (length(chrBAF$chr)/2), ypos = 0.75, coords = "genomic", dp = DisplayPars(pointSize = 1, 
                  cex = 0.7))
                textOverlays = c(textOverlays, list(hr1))
                
            }
            
        }
        
        
        pList = list(makeTitle(paste(sample_name, sep = ""), cex = 1), makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", 
            sep = "")), makeTitle(paste(" ", sep = "")), `Tumor BAF` = makeGenericArray(probeStart = 1:length(all_BAF_D), intensity = as.matrix(all_BAF_D), 
            dp = DisplayPars(color = "darkred", ylim = c(0, 0.5), type = "dot", pointSize = 0.1)), makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", 
            sep = "")), makeTitle(paste(" ", sep = "")))
        
        pdf(paste(folderName, "/", sample_name, ".pdf", sep = ""))
        gdPlot(pList, overlays = c(textOverlays, overlays1))
        dev.off()
    }
}


plotBAFAllSamples <- function(loh, fileName) {
    
    to.plot <- NULL
    for (i in 1:length(names(loh))) {
        to.plot <- rbind(to.plot, data.frame(loh[[i]], sampleId = names(loh)[i]))
        
    }
    
    to.plot <- to.plot[order(to.plot$position), ]
    idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$chr) == x)))))))
    to.plot <- to.plot[idx, ]
    to.plot$index <- 1:nrow(to.plot)
    
    last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$chr == x]))
    
    p <- ggplot(to.plot, aes(x = index, y = dev)) + geom_point(col = "darkred", size = 0.05) + scale_x_continuous(breaks = last_bp_vec, minor_breaks = NULL, 
        labels = seq(22)) + theme_bw() + facet_wrap(~sampleId, ncol = 1, strip.position = "left") + ylab("BAF deviation") + theme(legend.position = "none", 
        axis.text.x = element_text(size = 6, colour = "black", angle = -45, hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), 
        axis.ticks = element_line(size = 0.4), plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"), strip.text.x = element_text(size = 6), 
        strip.background = element_rect(fill = "white"))
    
    ggsave(fileName, plot = p)
}

plotGEAndBAFOneSample <- function(object, cnv.scale, loh.scale, sample, n = 50, scale.iteration = 50) {
    
    row_median <- apply(object@centered.data, 2, median)
    centered.control.normalized <- t(apply(object@centered.data, 1, "-", row_median))
    centered.control.normalized <- AverageReference(data = centered.control.normalized, ref_ids = object@control.sample.ids)
    centered.control.normalized[centered.control.normalized < (-2)] <- (-2)
    centered.control.normalized[centered.control.normalized > 2] <- 2
    centered.control.normalized[abs(centered.control.normalized) < 0.3] <- 0
    
    samples <- names(object@loh.median.filtered.data)
    plot_list = list()
    
    to.plot <- rbind(data.frame(variable = "expr.original", value=centered.control.normalized[, sample], object@annotation.filt[, c("Chr", "Position", "cytoband")]), 
       data.frame(variable = paste0("expr.scale", cnv.scale), value=object@control.normalized.visbound.noiseRemoved[[cnv.scale]][, sample], object@annotation.filt[, c("Chr", "Position", "cytoband")]) )
   # to.plot <- melt(data, measure.vars = colnames(data)[1:2])
    
    to.plot <- to.plot[order(to.plot$Position), ]
    idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$Chr) == x)))))))
    to.plot <- to.plot[idx, ]
    
    to.plot$index <- 1:nrow(to.plot)
    
    last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$Chr == x]))
    cnv.plot <- ggplot(to.plot, aes(x = index, y = value)) + geom_area(colour = "darkblue", fill = "darkblue") + scale_x_continuous(breaks = last_bp_vec, 
        minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~variable, ncol = 1, strip.position = "left") + ylab("Expression") + 
        ggtitle(paste("Expression", sample))
    
    maf <- object@loh[[sample]]
    data_smoothed <- maf$dev
    maf_temp <- maf
    for (i in 1:loh.scale) {
        
        maf_2 <- NULL
        for (m in 1:22) {
            chrBAF <- maf_temp[which(as.character(maf_temp$chr) == as.character(m)), ]
            chrBAF <- chrBAF[order(as.numeric(as.character(chrBAF$pos))), ]
            data_smoothed <- chrBAF$dev
            data_smoothed <- round(signal::filter(MedianFilter(n + 1), data_smoothed), digits = 2)
            chrBAF$dev <- data_smoothed
            maf_2 <- rbind(maf_2, chrBAF)
            
        }
        
        maf_temp <- maf_2
        n <- n + scale.iteration
    }
    scale<- maf_2

    to.plot <- rbind(original = data.frame(sample = "baf.original", object@loh[[sample]]), scale = data.frame(sample = paste0("baf.scale", loh.scale), 
        scale))
    
    to.plot <- to.plot[order(to.plot$position), ]
    idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$chr) == x)))))))
    to.plot <- to.plot[idx, ]
    to.plot$index <- 1:nrow(to.plot)
    
    last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$chr == x]))
    
    loh.plot <- ggplot(to.plot, aes(x = index, y = dev)) + geom_point(colour = "darkred", size = 0.05) + scale_x_continuous(breaks = last_bp_vec, 
        minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~sample, ncol = 1, strip.position = "left") + ylab("BAF") + ggtitle(paste("BAF", 
        sample))
    
    p <- ggarrange(loh.plot, cnv.plot, labels = c("A", "B"), ncol = 1, nrow = 2)
    ggsave(filename=paste0(sample, "_GE_BAF.png"), p)
}

plotGEAllSamples <- function(object, fileName = fileName, cnv.scale) {
    
    row_median <- apply(object@centered.data, 2, median)
    centered.control.normalized <- t(apply(object@centered.data, 1, "-", row_median))
    centered.control.normalized <- AverageReference(data = centered.control.normalized, ref_ids = object@control.sample.ids)
    centered.control.normalized[centered.control.normalized < (-2)] <- (-2)
    centered.control.normalized[centered.control.normalized > 2] <- 2
    centered.control.normalized[abs(centered.control.normalized) < 0.3] <- 0
    
    samples <- names(object@loh.median.filtered.data)
    plot_list = list()
    for (i in 1:length(samples)) {
        data <- data.frame(original = centered.control.normalized[, samples[i]], scale = object@control.normalized.visbound.noiseRemoved[[cnv.scale]][, 
            samples[i]], object@annotation.filt[, c("Chr", "Position", "cytoband")])
        
        
        to.plot <- melt(data, measure.vars = colnames(data)[1:2])
        
        
        to.plot <- to.plot[order(to.plot$Position), ]
        idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$Chr) == x)))))))
        to.plot <- to.plot[idx, ]
        
        
        to.plot$index <- 1:nrow(to.plot)
        
        to.plot$index <- 1:nrow(to.plot)
        
        
        
        last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$Chr == x]))
        plot_list[[i]] <- ggplot(to.plot, aes(x = index, y = value)) + geom_area(colour = "darkblue", fill = "darkblue") + scale_x_continuous(breaks = last_bp_vec, 
            minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~variable, ncol = 1, strip.position = "left") + ylab("Expression") + 
            ggtitle(samples[i])
    }
    
    pdf(fileName)
    for (i in 1:length(samples)) {
        print(plot_list[[i]])
    }
    dev.off()
}


plotBAFOneSample <- function(object, fileName) {
    
    object <- lohCallMedianFilterByChr(object, loh.scale = 1)
    scale1 <- object@loh.median.filtered.data
    object <- lohCallMedianFilterByChr(object, loh.scale = 2)
    scale2 <- object@loh.median.filtered.data
    object <- lohCallMedianFilterByChr(object, loh.scale = 3)
    scale3 <- object@loh.median.filtered.data
    
    samples <- names(object@loh.median.filtered.data)
    plot_list = list()
    for (i in 1:length(samples)) {
        to.plot <- rbind(original = data.frame(sample = "baf.original", object@loh[[samples[i]]]), scale1 = data.frame(sample = "baf.scale1", 
            scale1[[samples[i]]]), scale2 = data.frame(sample = "baf.scale2", scale2[[samples[i]]]), scale3 = data.frame(sample = "baf.scale3", 
            scale3[[samples[i]]]))
        
        to.plot <- to.plot[order(to.plot$position), ]
        idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$chr) == x)))))))
        to.plot <- to.plot[idx, ]
        to.plot$index <- 1:nrow(to.plot)
        
        last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$chr == x]))
        
        plot_list[[i]] <- ggplot(to.plot, aes(x = index, y = dev)) + geom_point(colour = "darkred", size = 0.05) + scale_x_continuous(breaks = last_bp_vec, 
            minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~sample, ncol = 1, switch = "y") + ylab("BAF") + ggtitle(samples[i]) + 
            theme(strip.text.y = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
    }
    
    pdf(fileName)
    for (i in 1:length(samples)) {
        print(plot_list[[i]])
    }
    dev.off()
    
}

plotSingleCellLargeScaleEventHeatmap <- function(finalChrMat, sampleName, chrs)
{

    title <- paste0(sampleName, "_", paste0(chrs, collapse="_"), sep="")

    finalChrMat<- finalChrMat[grep(sampleName, rownames(finalChrMat)), chrs]
    ord <- hclust( dist(finalChrMat, method = "euclidean"))$order

    plot.data <- melt(t(finalChrMat))

    plot.data$value2 <- "neutral"
    plot.data$value2[plot.data$value>0] <- "amplification"
    plot.data$value2[plot.data$value<0 ] <- "deletion"

    plot.data$value2 <- factor(plot.data$value2, levels=c("amplification", "deletion", "neutral"))
    plot.data$X1 <- factor(plot.data$X1, levels=colnames(finalChrMat))
    plot.data$X2 <- factor(plot.data$X2, levels=rev(rownames(finalChrMat)[ord]))

    p <- ggplot(aes(x=X1, y=X2, fill = value2), data=plot.data)
    p <- p + geom_tile(size=0.01) + labs(x="",y="")+
          scale_fill_manual(values = c("amplification"= muted("red"), "deletion"=muted("blue"), "neutral"="white"))  + #set a base size for all fonts
        theme_grey(base_size = 6)+ggtitle(title) + 
        theme(legend.position = "none",legend.direction = "vertical",
              legend.title = element_blank(),
               strip.text.x = element_blank(),
              legend.text = element_text(colour = "black",size = 7,face = "bold"),
              legend.key.height = grid::unit(0.8,"cm"),
              legend.key.width = grid::unit(0.5,"cm"),
              axis.text.x = element_text(size =6,colour = "black", angle = -45, hjust = 0),
              axis.text.y = element_text(size=2, vjust  =  0.2,colour = "black"),
              axis.ticks = element_line(size = 0.4),
              plot.title = element_text(colour = "black",hjust = 0,size = 6,face = "bold"))


    ggsave(paste0(title, ".png"), plot = p, width=1, height=3)


}

plotMUAndCooccurence <- function(results)
{
    for (j in 1:length(results))
    {
        mat <- results[[j]]       
        label <- unique(c(as.character(mat$Node1), as.character(mat$Node2)))
        inc.matrix <- matrix(0, nrow=length(label), ncol=length(label))
        rownames(inc.matrix) <- label
        colnames(inc.matrix) <- label
        for ( k in 1:nrow(mat))
        {   
            if(mat[k,5]== "occurence") inc.matrix[label==mat[k,1], label==mat[k,2]] <- log2(mat[k,3])
            if(mat[k,5]== "MU") {
                temp <- inc.matrix[label==mat[k,1], label==mat[k,2]] 
                if(is.na(temp)) {
                 inc.matrix[label==mat[k,1], label==mat[k,2]] <- (-log2(mat[k,3]))
                }
                else if(abs(temp)<(-log2(mat[k,3])))
                {
                    inc.matrix[label==mat[k,1],label==mat[k,2]] <- (-log2(mat[k,3]))
                }
            }
        }
 
        temp <- inc.matrix
        inc.matrix[abs(inc.matrix)<(-log2(0.01))] <- 0
        inc.matrix[which(is.na(inc.matrix))] <- 0
        #inc.matrix[inc.matrix<0] <- 0
        filt1 <-  which(apply(inc.matrix, 1, function(x) length(which(x==0))/length(x))==1)
        filt2 <-  which(apply(inc.matrix, 2, function(x) length(which(x==0))/length(x))==1)
        filter <- intersect(names(filt1), names(filt2))
        inc.matrix <- inc.matrix[!(rownames(inc.matrix) %in% filter), !(colnames(inc.matrix) %in% filter)]

        net=graph.adjacency(inc.matrix,mode="undirected",weighted=TRUE,diag=T) 
        V(net)$color=ifelse(grepl("amp", names(V(net))), "amp", "del")
          V(net)$vertex.names <- gsub("amp", "", gsub("del", "", V(net)$vertex.names))
        V(net)$size <- 30
        E(net)$type <- ifelse(E(net)$weight >0, "mu", "co")
        E(net)$weight <- abs(E(net)$weight)
        net2 <- asNetwork(net)
        
        p <- ggplot(net2, aes(x = x, y = y, xend = xend, yend = yend)) +
        geom_edges(aes(linetype=type, size=weight), color = "grey50") +
        geom_nodes(aes(color = color), size=15, alpha = 0.5) + 
        geom_nodetext(aes(color = "black", label = vertex.names),color = "black",
                      fontface = "bold", size=7)+ 
        scale_color_brewer(palette = "Set1") + ggtitle(names(results)[j])+
        theme_blank()

        ggsave(paste0(names(results)[j], "network.pdf"), plot = p, width=6, height=6)

        write.table(mat, paste0(names(results)[j], "_mat_edges.txt", sep=""), sep="\t", quote=F)
    }

}

plotSCellCNVTree <- function(finalChrMat, sampleName, path="C:\\Users\\aharmanci\\Downloads\\phylip-3.695\\phylip-3.695\\exe", fileName)
 {
#     #library(philentropy)
#     library(Rphylip)
#     library(phyloseq)

#     finalChrMat<- finalChrMat[grep(sampleName, rownames(finalChrMat)),]
#     chr <- names(which(apply(finalChrMat, 2, function(x) length(which(x!=0)))/nrow(finalChrMat)>0.1))
#     finalChrMat <- finalChrMat[, match(chr, colnames(finalChrMat))]
    
#     m<- distance(finalChrMat, method = "jaccard")
#     colnames(m) <- rownames(finalChrMat)
#     rownames(m) <- rownames(finalChrMat)
#     tree<-Rfitch(na.omit(m), global=F, path=path, outgroup=T)


#     pdf(fileName)
#     X <- as.matrix(finalChrMat)[,which(apply(finalChrMat, 2, sum)!=0)]
#     plot(tree, x.lim = 40, align.tip = TRUE, adj = 1, cex=0.2)
#     f <- function(n) c(muted("blue"), "white",muted("red"))
#     phydataplot(X, tree, "m", 2, width = 2, border = "white", lwd = 2,
#                 legend = "side", funcol = f)
#     dev.off()

    
#     del <- data.frame(apply(finalChrMat, 1, function(x) paste(names(which(x==(-1))), collapse="del,")))

#     amp <- data.frame(apply(finalChrMat, 1, function(x) paste(names(which(x==(1))), collapse="amp,")))

#     all <- paste0(del[,1], "del,", amp[,1], "amp")
#     all[which(all=="del,amp")] <- ""
#     all <- gsub(",amp", "", all)
#     all <- gsub("^del,", "", all)

#     data(GlobalPatterns)
#      physeq = prune_taxa(taxa_names(GlobalPatterns)[1:50], GlobalPatterns)
#     finalChrMat[finalChrMat==(-1)] <- 2
#      physeq@phy_tree <- tree

#     physeq@otu_table@.Data <- finalChrMat
#     physeq@tax_table <- NULL
#       physeq@sam_data <-  physeq@sam_data[1:length(colnames(finalChrMat)), ]
#      physeq@sam_data$X.SampleID <- colnames(finalChrMat)
#      physeq@sam_data$SampleType <- colnames(finalChrMat)
#      #physeq@sam_data$SampleType[!physeq@sam_data$SampleType %in% c("5q", "14q")] <- "others"
#       physeq = prune_taxa(taxa_names(physeq), physeq)
#       h <- plot_tree(physeq, ladderize="left", color="SampleType") +  scale_fill_brewer(palette="Spectral")
#      rownames(physeq@sam_data)  <- colnames(finalChrMat)

#     plot_tree(physeq)


}