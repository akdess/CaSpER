
generateParam <- function(object, scale = 4) {
    param <- data.frame(strength = 1e+07, e = 0.9999999, mu = quantile(object@control.normalized[[scale]], na.rm = TRUE, prob = c(0.01, 
        0.05, 0.5, 0.95, 0.99)), lambda = 20, nu = 2.1, kappa = c(0.05, 0.05, 0.8, 0.05, 0.05) * 1000, m = 0, eta = c(5, 5, 
        50, 5, 5) * 10000, gamma = 3, S = 0)
    param$m <- param$mu
    param$S <- ((sd(2^object@control.normalized[[scale]], na.rm = TRUE)/sqrt(nrow(param)))^2)
    rownames(param) <- seq(1, 5)
    object@hmmparam <- param
    return(object)
}

runCaSPER <- function(object, removeCentromere = T, cytoband = object@cytoband, method = "iterative") {
    final.objects <- list()
    
    if (method == "iterative") {
        loh.list <- list()
        cnv.list <- list()
        
        message("Performing recursive median filtering...")
        
        for (i in 1:object@loh.scale) {
            loh.list[[i]] <- lohCallMedianFilterByChr(object, loh.scale = i)
        }
        
        
        message("Performing HMM segmentation...")
        
        for (i in 1:object@cnv.scale) {
            cnv.list[[i]] <- PerformSegmentationWithHMM(object, cnv.scale = i, removeCentromere = T, cytoband = cytoband)
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
    } else if (method == "fixed") {
        object <- PerformSegmentationWithHMM(object, cnv.scale = object@cnv.scale, removeCentromere = T, cytoband = cytoband)
        object <- lohCallMedianFilterByChr(object, loh.scale = object@loh.scale)
        object <- calculateLOHShiftsForEachSegment(object)
        object <- assignStates(object)
        final.objects[[1]] <- generateLargeScaleEvents(object)
    }
    return(final.objects)
}

extractLargeScaleEvents <- function(final.objects, thr = 0.5) {
    
    mergeScales <- mergeScalesAndGenerateFinalEventSummary(final.objects)
    mergeScalesAmp <- mergeScales$mergeScalesAmp
    mergeScalesDel <- mergeScales$mergeScalesDel
    
    chrs <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
    finalChrMat <- matrix(0, ncol = 44, nrow = length(rownames(mergeScales$mergeScalesAmp)))
    colnames(finalChrMat) <- chrs
    rownames(finalChrMat) <- rownames(mergeScales$mergeScalesAmp)
    
    finalChrMat[(mergeScalesAmp/length(final.objects)) >= thr] <- 1
    finalChrMat[(mergeScalesDel/length(final.objects)) >= thr] <- (-1)
    
    return(finalChrMat)
}


mergeScalesAndGenerateFinalEventSummary <- function(final.objects) {
    sampleNames <- rownames(final.objects[[1]]@large.scale.cnv.events)
    mergeScalesAmp <- matrix(0, ncol = 44, nrow = length(sampleNames))
    colnames(mergeScalesAmp) <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
    rownames(mergeScalesAmp) <- sampleNames
    
    mergeScalesDel <- matrix(0, ncol = 44, nrow = length(sampleNames))
    colnames(mergeScalesDel) <- as.vector(sapply(1:22, function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
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



PerformSegmentationWithHMM <- function(object, cnv.scale, removeCentromere = T, cytoband) {
    
    object <- generateParam(object, scale = cnv.scale)
    data <- object@control.normalized[[cnv.scale]]
    annotation <- object@annotation.filt
    
    if (removeCentromere) {
        isCentromer <- annotation$isCentromer == "no"
        data <- data[isCentromer, ]
        annotation <- annotation[isCentromer, ]
    }
    
    filt <- annotation$cytoband %in% names(which(table(annotation$cytoband) > 1))
    annotation <- annotation[filt, ]
    data <- data[filt, ]
    
    segments <- NULL
    for (i in 1:dim(data)[2]) {
        rdata <- RangedData(IRanges(start = annotation$start, end = annotation$end), space = annotation$cytoband, copy = data[, 
            i])
        hmm.segments <- HMMsegment(rdata, param = object@hmmparam, verbose = F)
        segments <- rbind(segments, data.frame(ID = colnames(data)[i], hmm.segments$segs))
    }
    
    arms <- paste(cytoband$V1, cytoband$V4, sep = "")
    arm_sizes <- cytoband$V3 - cytoband$V2
    
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


calculateLOHShiftsForEachSegment <- function(object) {
    segments <- object@segments
    loh <- object@loh.median.filtered.data
    mapping <- object@loh.name.mapping
    segments$MeanDev <- rep(NA, nrow(segments))
    segments$medianDev <- rep(NA, nrow(segments))
    segments$nVarSite <- rep(NA, nrow(segments))
    
    for (i in 1:length(names(loh))) {
        sampleName <- as.character(mapping$sample.name[names(loh)[i] == mapping$loh.name])
        
        sub.segments <- segments[segments$ID %in% sampleName, ]
        sub.segments$chr <- as.character(sub.segments$chr)
        loh.i <- loh[[i]]
        
        m <- apply(sub.segments, 1, function(x, loh.i) {
            
            chrom <- gsub("p", "", gsub("q", "", as.character(x["chr"])))
            start <- as.numeric(x["start"])
            end <- as.numeric(x["end"])
            ind <- which(loh.i$position >= start & loh.i$chr == chrom & loh.i$position <= end)
            
            if (length(ind) > 0) {
                x["MeanDev"] <- mean(loh.i$dev[ind])
                x["medianDev"] <- median(loh.i$dev[ind])
                x["nVarSite"] <- length(ind)
            }
            x[c("MeanDev", "medianDev", "nVarSite")]
        }, loh.i = loh[[i]])
        sub.segments[c("MeanDev", "medianDev", "nVarSite")] <- t(m)
        segments[segments$ID %in% sampleName, ] <- sub.segments
        
    }
    object@segments <- segments
    return(object)
    
}


assignStates <- function(object) {
    object@loh.shift.thr <- 0.15
    X <- na.omit(as.numeric(object@segments$medianDev[object@segments$event_scale == "large_scale"]))
    X <- X + rnorm(length(X), 0, 1e-11)
    mod4 <- densityMclust(na.omit(X), modelName = "V", verbose = F)
    
    class2 <- as.numeric(table(mod4$classification)["2"])
    if (!is.na(class2) & class2 > 0) {
        if (object@sequencing.type == "single-cell") {
            
            object@loh.shift.thr <- min(mod4$data[mod4$classification == 2])
        }
        if (object@sequencing.type == "bulk") {
            object@loh.shift.thr <- median(mod4$data[mod4$classification == 2])
        }
    }
    
    object@segments$states2 <- rep("neut", length(object@segments$state))
    
    object@segments$states2[as.numeric(as.character(object@segments$state)) == 1] <- "del"
    object@segments$states2[as.numeric(as.character(object@segments$state)) == 5] <- "amp"
    
    object@segments$states2[as.numeric(as.character(object@segments$state)) == 2 & object@segments$medianDev > object@loh.shift.thr] <- "del"
    object@segments$states2[as.numeric(as.character(object@segments$state)) == 4 & object@segments$medianDev > object@loh.shift.thr] <- "amp"
    
    return(object)
}

generateLargeScaleEvents <- function(object) {
    amp <- extractEvents(segments = object@segments, cytoband = object@cytoband, type = "amp")
    del <- extractEvents(segments = object@segments, cytoband = object@cytoband, type = "del")
    final <- data.frame(amp = amp, del = del)
    colnames(final) <- c("LargeScaleAmp", "FocalAmp", "LargeScaleAmpNum", "FocalAmpNum", "LargeScaleDel", "FocalDel", "LargeScaleDelNum", 
        "FocalDelNum")
    object@large.scale.cnv.events <- final
    
    return(object)
}

extractEvents <- function(segments, cytoband, type) {
    
    sample_ids <- unique(segments$ID)
    mat_focal_paired <- matrix(0, nrow = length(sample_ids), ncol = 46)
    mat_large_paired <- matrix(0, nrow = length(sample_ids), ncol = 46)
    
    colnames(mat_focal_paired) <- as.vector(sapply(c(1:22, "X"), function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
    colnames(mat_large_paired) <- as.vector(sapply(c(1:22, "X"), function(x) c(paste(x, "p", sep = ""), paste(x, "q", sep = ""))))
    
    rownames(mat_focal_paired) <- sample_ids
    rownames(mat_large_paired) <- sample_ids
    
    arms <- paste(cytoband$V1, cytoband$V4, sep = "")
    arm_sizes <- cytoband$V3 - cytoband$V2
    
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
