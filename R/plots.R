## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function(coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3, "bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)
}

#@description helper function for plotHeatmap
# 
draw_matrix2 <- function(matrix, border_color, gaps_rows, gaps_cols, fmat, fontsize_number, number_color) {
    
    n = nrow(matrix)
    m = ncol(matrix)
    
    coord_x = pheatmap:::find_coordinates(m, gaps_cols)
    coord_y = pheatmap:::find_coordinates(n, gaps_rows)
    
    x = coord_x$coord - 0.5 * coord_x$size
    y = unit(1, "npc") - (coord_y$coord - 0.5 * coord_y$size)
    
    coord = expand.grid(y = y, x = x)
    
    res = gList()
    
    res[["rect"]] = rectGrob(x = coord$x, y = coord$y, width = coord_x$size, height = coord_y$size, gp = gpar(fill = matrix, 
        col = border_color))
    
    if (attr(fmat, "draw")) {
        res[["text"]] = textGrob(x = coord$x, y = coord$y, label = fmat, gp = gpar(col = number_color, fontsize = fontsize_number))
    }
    
    if (!is.null(gaps_cols)) {
        for (i in 1:(length(gaps_cols) - 1)) {
            res[[paste0("line", as.character(i))]] <- linesGrob(x = rep(x[gaps_cols[i] + 1], 2) - unit(1.5, "bigpts"), y = unit(c(0, 
                1), "npc"), gp = gpar(col = "black", lwd = 1, lty = 2))
        }
    }
    res = gTree(children = res)
    
    return(res)
}

#' @title plotHeatmap()
#'
#' @description  Visualization of the genomewide gene expression signal plot at different smoothing scales 
#'
#' @param object casper object
#' 
#' @param fileName fileName of the putput image
#'
#' @param cnv.scale expression.scale for the expression signal 
#'
#' @param cluster_cols boolean values determining if columns should be clustered 
#'
#' @param cluster_rows boolean values determining if rows should be clustered 
#'
#' @param show_rownames boolean values determining if rownames should be plotted
#'
#' @param only_soi boolean values determining if only samples of interest without control samples should be plotted
#'
#' @export
#'
#'
plotHeatmap10x <- function(object, fileName, cnv.scale = 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T, ...) {
    
    assignInNamespace(x = "draw_matrix", value = draw_matrix2, ns = asNamespace("pheatmap"))
    assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
    
     
    data <- object@control.normalized.noiseRemoved[[cnv.scale]]

    x.center <- mean(data)
    quantiles = quantile(data[data != x.center], c(0.01, 0.99))

    delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
    low_threshold = x.center - delta
    high_threshold = x.center + delta
    x.range = c(low_threshold, high_threshold)

    data[data < low_threshold] <- low_threshold
    data[data > high_threshold] <- high_threshold

    breaks <- seq(x.range[1], x.range[2], length=16)
      #  breaks <- seq(0, 2, length=16)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    
    idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:22)])
    xlabel <- rep("", length(rownames(object@data)))
    half <- round(table(object@annotation.filt$Chr)[as.character(1:22)]/2)[-1]
    xpos <- c(half[1], (idx[-22] + half))
    xlabel[xpos] <- 1:22
   
    if (only_soi) 
        data <- data[, !(colnames(data) %in% object@control.sample.ids)]
    
    pheatmap(t(data), cluster_cols = cluster_cols, cluster_rows = cluster_rows, gaps_col = idx, color = color, breaks = breaks, 
        labels_col = xlabel, show_rownames = show_rownames, filename = fileName, ...)
    
}

plotHeatmap <- function(object, fileName, cnv.scale = 3, cluster_cols = F, cluster_rows = T, show_rownames = T, only_soi = T,...) {
    
    assignInNamespace(x = "draw_matrix", value = draw_matrix2, ns = asNamespace("pheatmap"))
    assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", ns = asNamespace("pheatmap"))
    
    breaks <- seq(-2, 2, by = 0.2)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    
    idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:22)])
    xlabel <- rep("", length(rownames(object@data)))
    half <- round(table(object@annotation.filt$Chr)[as.character(1:22)]/2)[-1]
    xpos <- c(half[1], (idx[-22] + half))
    xlabel[xpos] <- 1:22
    
    data <- log2(object@control.normalized.noiseRemoved[[cnv.scale]])
    
    if (only_soi) 
        data <- data[, !(colnames(data) %in% object@control.sample.ids)]
    
    pheatmap(t(data), cluster_cols = cluster_cols, cluster_rows = cluster_rows, gaps_col = idx, color = color, breaks = breaks, 
        labels_col = xlabel, show_rownames = show_rownames, filename = fileName, ...)
    
}

#' @title plotGEAndGT()
#'
#' @description  Heatmap plot for large scale event calls identified by CaSpER and genotyping array.
#'
#' @param chrMat large scale events identified from CaSpER represented as matrix. Rows indicates samples (cells) whereas columns indicates chromosome arms
#' 
#' @param genoMat large scale events identified from genotyping array represented as matrix. Rows indicates samples (cells) whereas columns indicates chromosome arms
#'
#' @param fileName fileName of the putput image
#'
#' @export
#'
#'
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
    
    p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", 
        y = "") + scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), neutral = "white")) + 
        facet_wrap(~plot.data$sampleType, scales = "free_y", ncol = 1) + theme_grey(base_size = 6) + theme(legend.position = "right", 
        legend.direction = "vertical", legend.title = element_blank(), strip.text.x = element_blank(), legend.text = element_text(colour = "black", 
            size = 7, face = "bold"), legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"), 
        axis.text.x = element_text(size = 6, colour = "black", angle = -45, hjust = 0), axis.text.y = element_text(size = 6, 
            vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), plot.title = element_text(colour = "black", 
            hjust = 0, size = 6, face = "bold"))
    ggsave(fileName, plot = p)
    
}

#' @title plotLargeScaleEvent2()
#'
#' @description Visualization of the  large-scale CNV events among all the samples/cells
#'
#' @param chrMat large scale events identified from CaSpER represented as matrix. Rows indicates samples (cells) whereas columns indicates chromosome arms
#' 
#' @param fileName fileName of the output image
#'
#' @export
#'
#'
plotLargeScaleEvent2 <- function(chrMat, fileName) {
    
    
    plot.data <- melt(chrMat)
    
    plot.data$value2 <- "neutral"
    plot.data$value2[plot.data$value > 0] <- "amplification"
    plot.data$value2[plot.data$value < 0] <- "deletion"
    plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
    plot.data$X2 <- factor(plot.data$X2, levels = colnames(chrMat))
    
    p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", 
        y = "") + scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), neutral = "white")) + 
        theme_grey(base_size = 6) + theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank(), 
        strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, face = "bold"), legend.key.height = grid::unit(0.8, 
            "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 5, colour = "black", angle = -45, 
            hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
    ggsave(fileName, plot = p)
    
}

#' @title plotLargeScaleEvent()
#'
#' @description  Visualization of the  large-scale CNV events among all the samples/cells
#'
#' @param object casper object
#'
#' @param fileName fileName of the output image
#'
#' @export
#'
#'
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
    
    p <- ggplot(aes(x = X2, y = X1, fill = value2), data = plot.data) + geom_tile(colour = "gray50", size = 0.01) + labs(x = "", 
        y = "") + scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), neutral = "white")) + 
        theme_grey(base_size = 6) + theme(legend.position = "right", legend.direction = "vertical", legend.title = element_blank(), 
        strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, face = "bold"), legend.key.height = grid::unit(0.8, 
            "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 5, colour = "black", angle = -45, 
            hjust = 0), axis.text.y = element_text(size = 6, vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), 
        plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
    ggsave(fileName, plot = p)
    
}

#' @title plotBAFInSeperatePages()
#'
#' @description  Visualization of BAF deviation for each sample in separate pages 
#'
#' @param loh baf signal, user can either give smoothed baf signal or original baf signal as an input. 
#' 
#' @param folderName folder name for the output images  
#'
#' @return object
#'
#' @export
#'
#'
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
                hr1 <- makeTextOverlay(as.character(m), xpos = max_num - (length(chrBAF$chr)/2), ypos = 0.8, coords = "genomic", 
                  dp = DisplayPars(pointSize = 1, cex = 0.7))
                textOverlays = c(textOverlays, list(hr1))
            }
            if ((m%%2) == 1) {
                hr1 <- makeTextOverlay(as.character(m), xpos = max_num - (length(chrBAF$chr)/2), ypos = 0.75, coords = "genomic", 
                  dp = DisplayPars(pointSize = 1, cex = 0.7))
                textOverlays = c(textOverlays, list(hr1))
                
            }
            
        }
        
        
        pList = list(makeTitle(paste(sample_name, sep = ""), cex = 1), makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", 
            sep = "")), makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", sep = "")), `Tumor BAF` = makeGenericArray(probeStart = 1:length(all_BAF_D), 
            intensity = as.matrix(all_BAF_D), dp = DisplayPars(color = "darkred", ylim = c(0, 0.5), type = "dot", pointSize = 0.1)), 
            makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", sep = "")), makeTitle(paste(" ", sep = "")))
        
        pdf(paste(folderName, "/", sample_name, ".pdf", sep = ""))
        gdPlot(pList, overlays = c(textOverlays, overlays1))
        dev.off()
    }
}

#' @title plotBAFAllSamples()
#'
#' @description  Visualization of BAF shift signal for all samples together 
#'
#' @param loh baf signal, user can either give smoothed baf signal or original baf signal as an input. 
#' 
#' @param fileName fileName of the putput image
#'
#' @export
#'
#'
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
    
    p <- ggplot(to.plot, aes(x = index, y = dev)) + geom_point(col = "darkred", size = 0.05) + scale_x_continuous(breaks = last_bp_vec, 
        minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~sampleId, ncol = 1, strip.position = "left") + ylab("BAF deviation") + 
        theme(legend.position = "none", axis.text.x = element_text(size = 6, colour = "black", angle = -45, hjust = 0), axis.text.y = element_text(size = 6, 
            vjust = 0.2, colour = "black"), axis.ticks = element_line(size = 0.4), plot.title = element_text(colour = "black", 
            hjust = 0, size = 6, face = "bold"), strip.text.x = element_text(size = 6), strip.background = element_rect(fill = "white"))
    
    ggsave(fileName, plot = p)
}

#' @title plotGEAndBAFOneSample()
#'
#' @description  Gene expression and BAF signal for one sample in one plot
#'
#' @param object casper object
#' 
#' @param cnv.scale expression.scale for the expression signal 
#'
#' @param sample sample name 
#'
#' @param n window length used for median filtering
#' 
#' @param length.iterations increase in window length at each scale iteration 
#' 
#' @export
#'
#'
plotGEAndBAFOneSample <- function(object, cnv.scale, loh.scale, sample, n = 50, scale.iteration = 50) {
    
    row_median <- apply(object@centered.data, 2, median)
    centered.control.normalized <- t(apply(object@centered.data, 1, "-", row_median))
    centered.control.normalized <- AverageReference(data = centered.control.normalized, ref_ids = object@control.sample.ids)
    centered.control.normalized[centered.control.normalized < (-2)] <- (-2)
    centered.control.normalized[centered.control.normalized > 2] <- 2
    centered.control.normalized[abs(centered.control.normalized) < 0.3] <- 0
    
    samples <- names(object@loh.median.filtered.data)
    plot_list = list()
    
    to.plot <- rbind(data.frame(variable = "expr.original", value = (centered.control.normalized[, sample]), object@annotation.filt[, 
        c("Chr", "Position", "cytoband")]), data.frame(variable = paste0("expr.scale", cnv.scale), value = log2(object@control.normalized.noiseRemoved[[cnv.scale]][, 
        sample]), object@annotation.filt[, c("Chr", "Position", "cytoband")]))
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
    scale <- maf_2
    
    to.plot <- rbind(original = data.frame(sample = "baf.original", object@loh[[sample]]), scale = data.frame(sample = paste0("baf.scale", 
        loh.scale), scale))
    
    to.plot <- to.plot[order(to.plot$position), ]
    idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$chr) == x)))))))
    to.plot <- to.plot[idx, ]
    to.plot$index <- 1:nrow(to.plot)
    
    last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$chr == x]))
    
    loh.plot <- ggplot(to.plot, aes(x = index, y = dev)) + geom_point(colour = "darkred", size = 0.05) + scale_x_continuous(breaks = last_bp_vec, 
        minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~sample, ncol = 1, strip.position = "left") + ylab("BAF") + 
        ggtitle(paste("BAF", sample))
    
    p <- ggarrange(loh.plot, cnv.plot, labels = c("A", "B"), ncol = 1, nrow = 2)
    ggsave(filename = paste0(sample, "_GE_BAF.png"), p)
}

#' @title plotGEAllSamples()
#'
#' @description  plot gene expression signal for each sample seperately
#'
#' @param object casper object 
#'
#' @param fileName fileName of the putput image
#'
#' @param cnv.scale  expression.scale for the expression signal  
#'
#' @export
#'
#'
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
        data <- data.frame(original = (centered.control.normalized[, samples[i]]), scale = log2(object@control.normalized.noiseRemoved[[cnv.scale]][, 
            samples[i]]), object@annotation.filt[, c("Chr", "Position", "cytoband")])
        
        
        to.plot <- melt(data, measure.vars = colnames(data)[1:2])
        
        
        to.plot <- to.plot[order(to.plot$Position), ]
        idx = unlist(as.vector((sapply(1:22, function(x) as.vector(unlist(which(as.character(to.plot$Chr) == x)))))))
        to.plot <- to.plot[idx, ]
        
        
        to.plot$index <- 1:nrow(to.plot)
        
        to.plot$index <- 1:nrow(to.plot)
        
        
        
        last_bp_vec <- sapply(1:22, function(x) min(to.plot$index[to.plot$Chr == x]))
        plot_list[[i]] <- ggplot(to.plot, aes(x = index, y = value)) + geom_area(colour = "darkblue", fill = "darkblue") + 
            scale_x_continuous(breaks = last_bp_vec, minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~variable, 
            ncol = 1, strip.position = "left") + ylab("Expression") + ggtitle(samples[i])
    }
    
    pdf(fileName)
    for (i in 1:length(samples)) {
        print(plot_list[[i]])
    }
    dev.off()
}

#' @title plotBAFOneSample()
#'
#' @description  Visualization of BAF shift signal in different scales for one sample
#'
#' @param object  casper object
#' 
#' @param fileName fileName of the output image
#'
#' @export
#'
#'

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
            minor_breaks = NULL, labels = seq(22)) + theme_bw() + facet_wrap(~sample, ncol = 1, switch = "y") + ylab("BAF") + 
            ggtitle(samples[i]) + theme(strip.text.y = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
    }
    
    pdf(fileName)
    for (i in 1:length(samples)) {
        print(plot_list[[i]])
    }
    dev.off()
    
}

#' @title plotSingleCellLargeScaleEventHeatmap()
#'
#' @description  Visualization of large scale event summary for selected samples and chromosomes
#'
#' @param finalChrMat large scale events identified from CaSpER represented as matrix. Rows indicates samples (cells) whereas columns indicates chromosome arms
#' 
#' @param sampleName sample name
#'
#' @param chrs chromosome names
#'
#' @return object
#'
#' @export
#'
#'
plotSingleCellLargeScaleEventHeatmap <- function(finalChrMat, sampleName, chrs) {
    
    title <- paste0(sampleName, "_", paste0(chrs, collapse = "_"), sep = "")
    
    finalChrMat <- finalChrMat[grep(sampleName, rownames(finalChrMat)), chrs]
    ord <- hclust(dist(finalChrMat, method = "euclidean"))$order
    
    plot.data <- melt(t(finalChrMat))
    
    plot.data$value2 <- "neutral"
    plot.data$value2[plot.data$value > 0] <- "amplification"
    plot.data$value2[plot.data$value < 0] <- "deletion"
    
    plot.data$value2 <- factor(plot.data$value2, levels = c("amplification", "deletion", "neutral"))
    plot.data$X1 <- factor(plot.data$X1, levels = colnames(finalChrMat))
    plot.data$X2 <- factor(plot.data$X2, levels = rev(rownames(finalChrMat)[ord]))
    
    p <- ggplot(aes(x = X1, y = X2, fill = value2), data = plot.data)
    p <- p + geom_tile(size = 0.01) + labs(x = "", y = "") + scale_fill_manual(values = c(amplification = muted("red"), deletion = muted("blue"), 
        neutral = "white")) + theme_grey(base_size = 6) + ggtitle(title) + theme(legend.position = "none", legend.direction = "vertical", 
        legend.title = element_blank(), strip.text.x = element_blank(), legend.text = element_text(colour = "black", size = 7, 
            face = "bold"), legend.key.height = grid::unit(0.8, "cm"), legend.key.width = grid::unit(0.5, "cm"), axis.text.x = element_text(size = 6, 
            colour = "black", angle = -45, hjust = 0), axis.text.y = element_text(size = 2, vjust = 0.2, colour = "black"), 
        axis.ticks = element_line(size = 0.4), plot.title = element_text(colour = "black", hjust = 0, size = 6, face = "bold"))
    
    
    ggsave(paste0(title, ".pdf"), plot = p, width = 1, height = 3)
    
    
}

#' @title plotMUAndCooccurence()
#'
#' @description  Visualization of mutually exclusive and co-occuring events
#'
#' @param results output of extractMUAndCooccurence() function
#'
#' @export
#'
#'
plotMUAndCooccurence <- function(results) {
    for (j in 1:length(results)) {
        mat <- results[[j]]
        label <- unique(c(as.character(mat$Node1), as.character(mat$Node2)))
        inc.matrix <- matrix(0, nrow = length(label), ncol = length(label))
        rownames(inc.matrix) <- label
        colnames(inc.matrix) <- label
        for (k in 1:nrow(mat)) {
            if (mat[k, 5] == "occurence") 
                inc.matrix[label == mat[k, 1], label == mat[k, 2]] <- log2(mat[k, 3])
            if (mat[k, 5] == "MU") {
                temp <- inc.matrix[label == mat[k, 1], label == mat[k, 2]]
                if (is.na(temp)) {
                  inc.matrix[label == mat[k, 1], label == mat[k, 2]] <- (-log2(mat[k, 3]))
                } else if (abs(temp) < (-log2(mat[k, 3]))) {
                  inc.matrix[label == mat[k, 1], label == mat[k, 2]] <- (-log2(mat[k, 3]))
                }
            }
        }
        
        temp <- inc.matrix
        inc.matrix[abs(inc.matrix) < (-log2(0.01))] <- 0
        inc.matrix[which(is.na(inc.matrix))] <- 0
        # inc.matrix[inc.matrix<0] <- 0
        filt1 <- which(apply(inc.matrix, 1, function(x) length(which(x == 0))/length(x)) == 1)
        filt2 <- which(apply(inc.matrix, 2, function(x) length(which(x == 0))/length(x)) == 1)
        filter <- intersect(names(filt1), names(filt2))
        inc.matrix <- inc.matrix[!(rownames(inc.matrix) %in% filter), !(colnames(inc.matrix) %in% filter)]
        
        net = graph.adjacency(inc.matrix, mode = "undirected", weighted = TRUE, diag = T)
        V(net)$color = ifelse(grepl("amp", names(V(net))), "amp", "del")
        V(net)$vertex.names <- gsub("amp", "", gsub("del", "", V(net)$vertex.names))
        V(net)$size <- 30
        E(net)$type <- ifelse(E(net)$weight > 0, "mu", "co")
        E(net)$weight <- abs(E(net)$weight)
        net2 <- asNetwork(net)
        
        p <- ggplot(net2, aes(x = x, y = y, xend = xend, yend = yend)) + geom_edges(aes(linetype = type, size = weight), color = "grey50") + 
            geom_nodes(aes(color = color), size = 15, alpha = 0.5) + geom_nodetext(aes(color = "black", label = vertex.names), 
            color = "black", fontface = "bold", size = 7) + scale_color_brewer(palette = "Set1") + ggtitle(names(results)[j]) + 
            theme_blank()
        
        ggsave(paste0(names(results)[j], "network.pdf"), plot = p, width = 6, height = 6)
        
        write.table(mat, paste0(names(results)[j], "_mat_edges.txt", sep = ""), sep = "\t", quote = F)
    }
    
}

#' @title plotSCellCNVTree()
#'
#' @description  Pyhlogenetic tree-based clustering and visualization of the cells based on the CNV events from single cell RNA-seq Data.
#'
#' @param finalChrMat large scale events identified from CaSpER represented as matrix. Rows indicates samples (cells) whereas columns indicates chromosome arms
#' 
#' @param sampleName sample name
#'
#' @param path path to the executable containing fitch. If path = NULL, the R will search several commonly used directories for the correct executable file. More information about installing PHYLIP can be found on the PHYLIP webpage: http://evolution.genetics.washington.edu/phylip.html.
#'
#' @export
#'
#'
plotSCellCNVTree <- function(finalChrMat, sampleName, path = "C:\\Users\\aharmanci\\Downloads\\phylip-3.695\\phylip-3.695\\exe", 
    fileName) {
    finalChrMat <- finalChrMat[grep(sampleName, rownames(finalChrMat)), ]
    chr <- names(which(apply(finalChrMat, 2, function(x) length(which(x != 0)))/nrow(finalChrMat) > 0.1))
    finalChrMat <- finalChrMat[, match(chr, colnames(finalChrMat))]
    
    m <- distance(finalChrMat, method = "jaccard")
    colnames(m) <- rownames(finalChrMat)
    rownames(m) <- rownames(finalChrMat)
    tree <- Rfitch(na.omit(m), global = T, path = path, outgroup = T)
  
    pdf(fileName)
    X <- as.matrix(finalChrMat)[, which(apply(finalChrMat, 2, sum) != 0)]
    plot(tree, x.lim = 40, align.tip = TRUE, adj = 1, cex = 0.2)
    f <- function(n) c(muted("blue"), "white", muted("red"))
    phydataplot(X, tree, "m", 2, width = 2, border = "white", lwd = 2, legend = "side", funcol = f)
    dev.off()
    
    
}
