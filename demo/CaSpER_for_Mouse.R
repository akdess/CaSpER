#Load libraries
library(CaSpER)
library(GenomicRanges)
library(biomaRt)
library(gprofiler2)
library(plyr)


plotHeatmapMM <- function (object, fileName, cnv.scale = 3, cluster_cols = F, 
    cluster_rows = T, show_rownames = T, only_soi = F) 
{
    assignInNamespace(x = "draw_matrix", value = draw_matrix2, 
        ns = asNamespace("pheatmap"))
    assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", 
        ns = asNamespace("pheatmap"))
    breaks <- seq(-2, 2, by = 0.2)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:19)])
    xlabel <- rep("", length(rownames(object@data)))
    half <- round(table(object@annotation.filt$Chr)[as.character(1:19)]/2)[-1]
    xpos <- c(half[1], (idx[-19] + half))
    xlabel[xpos] <- 1:19
    data <- log2(object@control.normalized.noiseRemoved[[cnv.scale]])
    p<- pheatmap(t(data), cluster_cols = F, cluster_rows = T, gaps_col = idx, 
        color = color, breaks = breaks, labels_col = xlabel, 
        show_rownames = T, filename = fileName)
    return(p)
}

plotHeatmap10x <- function (object, fileName, cnv.scale = 3, cluster_cols = F, 
    cluster_rows = T, show_rownames = T, only_soi = T) 
{
    assignInNamespace(x = "draw_matrix", value = draw_matrix2, 
        ns = asNamespace("pheatmap"))
    assignInNamespace(x = "draw_colnames", value = "draw_colnames_45", 
        ns = asNamespace("pheatmap"))
    data <- object@control.normalized.noiseRemoved[[cnv.scale]]
    x.center <- mean(data)
    quantiles = quantile(data[data != x.center], c(0.01, 0.99))
    delta = max(abs(c(x.center - quantiles[1], quantiles[2] - 
        x.center)))
    low_threshold = x.center - delta
    high_threshold = x.center + delta
    x.range = c(low_threshold, high_threshold)
    data[data < low_threshold] <- low_threshold
    data[data > high_threshold] <- high_threshold
    breaks <- seq(x.range[1], x.range[2], length = 16)
    color <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(length(breaks))
    idx <- cumsum(table(object@annotation.filt$Chr)[as.character(1:19)])
    xlabel <- rep("", length(rownames(object@data)))
    half <- round(table(object@annotation.filt$Chr)[as.character(1:19)]/2)[-1]
    xpos <- c(half[1], (idx[-19] + half))
    xlabel[xpos] <- 1:19
    if (only_soi) 
        data <- data[, !(colnames(data) %in% object@control.sample.ids)]
    pheatmap(t(data), cluster_cols = F, cluster_rows = T, gaps_col = idx, 
        color = color, breaks = breaks, labels_col = xlabel, 
        show_rownames = T, filename = "heatmap.png")
}
#Generate Annotations
generateAnnotation <- function (id_type = "ensembl_gene_id", genes)
{
  id_type = "ensembl_gene_id"
  mart <- useDataset("mmusculus_gene_ensembl", useEnsembl(biomart = "ensembl"))
  
  G_list <- getBM(filters = "ensembl_gene_id", attributes = c(id_type,
                                                              "mgi_symbol", "ensembl_gene_id", "chromosome_name", "start_position",
                                                              "end_position", "band"), values = genes, mart = mart)
  common <- intersect(genes, G_list[, id_type])
  ord <- match(common, G_list[, id_type])
  annotation <- G_list[ord, ]
  annotation <- annotation[order(annotation$start_position),
                           ]
  annotation$cytoband <- paste0(annotation$chromosome_name,
                                substr((annotation$band), 0, 1))
  chr.list <- as.character(1:19)
  idx <- unlist(unique(as.vector(sapply(chr.list, function(x) as.vector(unlist(which(as.character(annotation$chromosome_name) ==
                                                                                       x)))))))
  annotation <- as.data.frame(annotation[idx, ])
  colnames(annotation)[c(1:6, 8)] <- c("Gene", "GeneSymbol",
                                       "EntrezID", "Chr", "start", "end", "cytoband")
  annotation$isCentromer <- rep("no", nrow(annotation))
  
  
  annotation$Position <- (as.numeric(annotation$start) + as.numeric(annotation$end))/2
  annotation$new_positions <- as.vector(unlist(lapply(lapply(split(annotation$cytoband,
                                                                   annotation$cytoband), length)[unique(annotation$cytoband)],
                                                      function(x) 1:x)))
  return(annotation)
}

lohCallMedianFilterByChr <- function(object, loh.scale, n = 50, scale.iteration = 50) {
        
    object@loh.median.filtered.data <- list()
    
    for (j in 1:length(object@loh)) {
        window <- n
        maf <- object@loh[[j]]
        data_smoothed <- maf$dev
        maf_temp <- maf
        for (i in 1:loh.scale) {
            
            maf_2 <- NULL
            for (m in 1:19) {
                chrBAF <- maf_temp[which(as.character(maf_temp$chr) == as.character(m)), ]
                chrBAF <- chrBAF[order(as.numeric(as.character(chrBAF$pos))), ]
                data_smoothed <- chrBAF$dev
                data_smoothed <- round(signal::filter(MedianFilter(window + 1), data_smoothed), digits = 2)
                chrBAF$dev <- data_smoothed
                maf_2 <- rbind(maf_2, chrBAF)
                
            }
            
            maf_temp <- maf_2
            window <- window + scale.iteration
        }
        object@loh.median.filtered.data[[j]] <- maf_2
    }
    names(object@loh.median.filtered.data) <- names(object@loh)
    
    return(object)
}

#'
PerformSegmentationWithHMM <- function(object, cnv.scale) {
    
    object <- generateParam(object, cnv.scale = cnv.scale)
    data <- object@control.normalized[[cnv.scale]]
    annotation <- object@annotation.filt[]
    annotation$cytoband <- annotation$Chr

    indices <- annotation$Chr %in% names(which(table(annotation$Chr) > 1))
    annotation <- annotation[indices, , drop = FALSE]
    data <- data[indices, , drop = FALSE]
    
    segments <- NULL
    for (i in 1:dim(data)[2]) {
        rdata <- GRanges(ranges=IRanges(start = annotation$start,
            end = annotation$end), seqnames = annotation$Chr,
            copy = data[, i], chr=annotation$Chr ,  space=annotation$Chr)
        rdata <- data.frame(rdata)
        rdata$chr <- as.factor(rdata$chr)

        hmm.segments <- HMMsegment(correctOut=rdata, param = object@hmmparam, verbose = F)
        segments <- rbind(segments, data.frame(ID = colnames(data)[i], hmm.segments$segs))
    }
    
    cytoband <- object@cytoband
    arms <- cytoband$V1
    arm_sizes <- cytoband$V3
    
    object@segments <- segments
    object@segments$event_scale <- rep("", nrow(object@segments))
    
    for (i in 1:dim(object@segments)[1]) {
        ind <- which(annotation$Position >= object@segments$start[i] & annotation$cytoband == object@segments[i, "chr"] & annotation$Position <= 
            object@segments$end[i])
        object@segments$size[i] <- object@segments$end[i] - object@segments$start[i]
        object@segments$num.marks[i] <- length(ind)
        pair_arm_sizes <- arm_sizes[match(as.character(object@segments[i, "chr"]), arms)]
        object@segments$arm.size.perc[i] <- object@segments$size[i]/pair_arm_sizes
        if (object@segments$size[i] > pair_arm_sizes * (1/3)) 
            object@segments$event_scale[i] <- "large_scale"
        if ((object@segments$size[i] > pair_arm_sizes * (1/10)) & (object@segments$size[i] < pair_arm_sizes * (1/3))) 
            object@segments$event_scale[i] <- "focal"
        
    }
    
    return(object)
}

generateLargeScaleEvents <- function (object) 
{
    amp <- extractEvents(segments = object@segments, cytoband = object@cytoband, 
        type = "amp")
    del <- extractEvents(segments = object@segments, cytoband = object@cytoband, 
        type = "del")
    final <- data.frame(amp = amp, del = del)
    colnames(final) <- c("LargeScaleAmp", "FocalAmp", 
        "LargeScaleAmpNum", "FocalAmpNum", "LargeScaleDel", 
        "FocalDel", "LargeScaleDelNum", "FocalDelNum")
    object@large.scale.cnv.events <- final
    return(object)
}
extractEvents <- function(segments, cytoband, type) {
    
    sample_ids <- unique(segments$ID)
    mat_focal_paired <- matrix(0, nrow = length(sample_ids), ncol = 19)
    mat_large_paired <- matrix(0, nrow = length(sample_ids), ncol = 19)
    
    colnames(mat_focal_paired) <- 1:19
    colnames(mat_large_paired) <- 1:19
    rownames(mat_focal_paired) <- sample_ids
    rownames(mat_large_paired) <- sample_ids
    
    arms <- cytoband$V1
    arm_sizes <- cytoband$V3
    
    for (i in 1:length(sample_ids)) {
        
        k <- which(segments$ID == as.character(sample_ids[i]))
        
        sub_segments <- segments[k, ]
        sub_segments$chr <- as.character(sub_segments$chr)
        uniq_arms <- unique(as.character(sub_segments$chr))
        
        for (t in 1:length(uniq_arms)) {
            
            t1 <- which(as.character(sub_segments$chr) == uniq_arms[t] & as.character(sub_segments$states2) == type)
            armSize <- sum(sub_segments$size[t1], na.rm = T)
            pair_arm_sizes <- arm_sizes[uniq_arms[t] == arms]
            isLargeScale <- armSize >= (pair_arm_sizes * (1/3))
            isFocalScale <- (armSize < pair_arm_sizes * (1/3)) & armSize >= (pair_arm_sizes * (1/10))
            
            if (isFocalScale) {
                mat_focal_paired[i, uniq_arms[t]] <- 1
            }
            if (isLargeScale) {
                mat_large_paired[i, uniq_arms[t]] <- 1
            }
            
        }
    }
    
    
    result_focal <- apply(mat_focal_paired, 1, function(x) paste(colnames(mat_focal_paired)[which(x == 1)], collapse = " "))
    result_large <- apply(mat_large_paired, 1, function(x) paste(colnames(mat_large_paired)[which(x == 1)], collapse = " "))
    
    result_focal_number <- apply(mat_focal_paired, 1, function(x) length(which(x == 1)))
    result_large_number <- apply(mat_large_paired, 1, function(x) length(which(x == 1)))
    
    
    summary_events <- cbind(result_large, result_focal, result_large_number, result_focal_number)
    return(summary_events)
}

runCaSpER <- function(object, removeCentromere = T, cytoband = object@cytoband, method = "iterative") {
    final.objects <- list()

    loh.list <- list()
    cnv.list <- list()
    
    message("Performing recursive median filtering...")
    
    for (i in 1:object@loh.scale) {
        loh.list[[i]] <- lohCallMedianFilterByChr(object, loh.scale = i)
    }
       
    message("Performing HMM segmentation...")
    
    for (i in 1:object@cnv.scale) {
        cnv.list[[i]] <- PerformSegmentationWithHMM(object, cnv.scale = i)
    }
    
    combin <- expand.grid(1:object@cnv.scale, 1:object@loh.scale)
    list.names <- apply(combin, 1, function(x) paste(x[1], x[2], sep = "_vs_"))
    
    for (i in 1:nrow(combin)) {
        loh.scale <- combin[i, 2]
        cnv.scale <- combin[i, 1]
        message("Processing cnv.scale:", cnv.scale, " loh.scale:", loh.scale, "...")
        object <- cnv.list[[cnv.scale]]
        object@loh.median.filtered.data <- loh.list[[loh.scale]]@loh.median.filtered.data
        object <- calculateLOHShiftsForEachSegment(object)
        object <- assignStates(object)
        final.objects[[i]] <- generateLargeScaleEvents(object)
    }
    names(final.objects) <- list.names

   return(final.objects)
}



#' @title extractLargeScaleEventsMouse()
#'
#' @description generates coherent set of large scale CNV events using the pairwise comparison of all scales from BAF and expression signals
#'
#' @param final.objects casper object
#' 
#' @param thr gamma threshold determining the least number of scales required to support 
#'
#' @return final large scale event summary reported as a matrix 
#'
#' @export
#'
#'
extractLargeScaleEventsMouse <- function(final.objects, thr = 0.5) {
    
    mergeScales <- mergeScalesAndGenerateFinalEventSummary(final.objects)
    mergeScalesAmp <- mergeScales$mergeScalesAmp
    mergeScalesDel <- mergeScales$mergeScalesDel

    finalChrMat <- matrix(0, ncol =19, nrow = length(rownames(mergeScales$mergeScalesAmp)))
    colnames(finalChrMat) <- 1:19
    rownames(finalChrMat) <- rownames(mergeScales$mergeScalesAmp)
    
    finalChrMat[(mergeScalesAmp/length(final.objects)) >= thr] <- 1
    finalChrMat[(mergeScalesDel/length(final.objects)) >= thr] <- (-1)
    
    return(finalChrMat)
}
#' @title mergeScalesAndGenerateFinalEventSummary()
#'
#' @description  helper function for extractLargeScaleEvents()
#'
#' @param final.objects list of casper objects
#'
#' @return list of objects
#'
#' @export
#'
#'
mergeScalesAndGenerateFinalEventSummary <- function(final.objects) {
    sampleNames <- rownames(final.objects[[1]]@large.scale.cnv.events)
    mergeScalesAmp <- matrix(0, ncol = 19, nrow = length(sampleNames))
    colnames(mergeScalesAmp) <- 1:19
    rownames(mergeScalesAmp) <- sampleNames
    
    mergeScalesDel <- matrix(0, ncol = 19, nrow = length(sampleNames))
    colnames(mergeScalesDel) <- 1:19
    rownames(mergeScalesDel) <- sampleNames
    
    for (i in 1:length(final.objects)) {
        
        object <- final.objects[[i]]
        
        for (x in 1:nrow((mergeScalesDel))) {
            
            chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(object@large.scale.cnv.events$LargeScaleAmp[sampleNames[x] == 
                rownames(object@large.scale.cnv.events)], collapse = " ")), split = " ")))
            if (length(chrWithEvents) > 0) 
                mergeScalesAmp[x, match(intersect(chrWithEvents, colnames(mergeScalesAmp)), colnames(mergeScalesAmp))] <- mergeScalesAmp[x, 
                  match(intersect(chrWithEvents, colnames(mergeScalesAmp)), colnames(mergeScalesAmp))] + 1
            
            chrWithEvents <- as.vector(unlist(strsplit(as.vector(paste(object@large.scale.cnv.events$LargeScaleDel[sampleNames[x] == 
                rownames(object@large.scale.cnv.events)], collapse = " ")), split = " ")))
            if (length(chrWithEvents) > 0) 
                mergeScalesDel[x, match(intersect(chrWithEvents, colnames(mergeScalesDel)), colnames(mergeScalesDel))] <- mergeScalesDel[x, 
                  match(intersect(chrWithEvents, colnames(mergeScalesDel)), colnames(mergeScalesDel))] + 1
            
        }
        
    }
    
    return(list(mergeScalesAmp = mergeScalesAmp, mergeScalesDel = mergeScalesDel))
}

# load CreateCasperObject object 
load("object.rda")
library(GenomicRanges)

p<- plotHeatmapMM(object, fileName="heatmap.png")
## runCaSpER
object@annotation$cytoband <- object@annotation$Chr
object@annotation.filt$cytoband <- object@annotation.filt$Chr

final.objects <- runCaSpER(object)

finalChrMat <- extractLargeScaleEventsMouse  (final.objects, thr=0.75) 

k<- pheatmap(finalChrMat[p$tree_row$labels[p$tree_row$order], ], 
  color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
  show_rownames=T, cluster_rows=F, cluster_cols=F)
