
CreateCasperObject <- function(raw.data, annotation, control.sample.ids, cytoband, loh.name.mapping, cnv.scale, loh.scale, 
    method, loh, project = "casperProject", sequencing.type, expr.cutoff = 4.5, display.progress = TRUE, log.transformed = TRUE, 
    centered.threshold = 3, window.length = 50, length.iterations = 50, vis.bound = 2, noise.thr = 0.3, genomeVersion = "hg19", 
    ...) {
    
    object <- new(Class = "casper", raw.data = raw.data, loh = loh, annotation = annotation, sequencing.type = sequencing.type, 
        control.sample.ids = control.sample.ids, project.name = project, version = casper.version, cytoband = cytoband, loh.name.mapping = loh.name.mapping, 
        cnv.scale = cnv.scale, loh.scale = loh.scale, method = method, window.length = window.length, length.iterations = length.iterations, 
        vis.bound = vis.bound, noise.thr = noise.thr, genomeVersion = genomeVersion)
    # filter cells on number of genes detected modifies the raw.data slot as well now
    if (!log.transformed) {
        object@raw.data <- log2(object@raw.data + 1)
    }
    object.raw.data <- object@raw.data
    
    gene.average <- log2(rowMeans((((2^data) - 1) * 10), na.rm = TRUE) + 1)
    genes.use <- gene.average > expr.cutoff
    object@data <- object.raw.data[genes.use, ]
    object.data <- object@data
    object@annotation.filt <- object@annotation[genes.use, ]
    
    centered.data <- object.data - rowMeans(object.data, na.rm = TRUE)
    centered.data[centered.data > centered.threshold] <- centered.threshold
    centered.data[centered.data < (-centered.threshold)] <- (-centered.threshold)
    centered.data <- centered.data - rowMeans(centered.data, na.rm = TRUE)
    object@centered.data <- centered.data
    
    object <- processMedianFiltering(object)
    
    return(object)
}

processMedianFiltering <- function(object) {
    
    object <- PerformMedianFilterByChr(object, window.length = object@window.length, length.iterations = object@length.iterations)
    
    object <- CenterSmooth(object)
    
    object <- ControlNormalize(object, vis.bound = object@vis.bound, noise.thr = object@noise.thr)
    # parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals('CreateCasperObject'))]
    # parameters.to.store$raw.data <- NULL
    return(object)
}


PerformMedianFilter <- function(object, window.length = 50, length.iterations = 50) {
    
    median.filtered.data <- list()
    
    for (i in 1:object@cnv.scale) {
        if (i == 1) {
            median.filtered.data[[i]] <- apply(object@centered.data, 2, function(x) round(filter(MedianFilter(window.length + 
                1), x), digits = 2))
            rownames(median.filtered.data[[i]]) <- rownames(object@centered.data)
            window.length <- window.length + length.iterations
        } else {
            median.filtered.data[[i]] <- apply(median.filtered.data[[i - 1]], 2, function(x) round(filter(MedianFilter(window.length + 
                1), x), digits = 2))
            rownames(median.filtered.data[[i]]) <- rownames(object@centered.data)
            window.length <- window.length + length.iterations
        }
    }
    object@median.filtered.data <- median.filtered.data
    return(object)
}

PerformMedianFilterByChr <- function(object, window.length = 50, length.iterations = 50) {
    
    median.filtered.data <- list()
    
    for (i in 1:object@cnv.scale) {
        if (i == 1) {
            median.filtered.data[[i]] <- apply(object@centered.data, 2, function(x) {
                dataByChr <- split(x, object@annotation.filt[, "Chr"])[unique(object@annotation.filt[, "Chr"])]
                dataAll <- c()
                sapply(dataByChr, function(y) {
                  r <- as.vector(y)
                  f <- round(filter(MedianFilter(window.length + 1), r), digits = 2)
                  dataAll <<- c(dataAll, f)
                  # return(dataAll)
                })
                return(dataAll)
            })
            
            rownames(median.filtered.data[[i]]) <- rownames(object@centered.data)
            window.length <- window.length + length.iterations
        } else {
            median.filtered.data[[i]] <- apply(median.filtered.data[[i - 1]], 2, function(x) {
                dataByChr <- split(x, object@annotation.filt[, "Chr"])[unique(object@annotation.filt[, "Chr"])]
                dataAll <- c()
                sapply(dataByChr, function(y) {
                  r <- as.vector(y)
                  f <- round(filter(MedianFilter(window.length + 1), r), digits = 2)
                  dataAll <<- c(dataAll, f)
                  # return(dataAll)
                })
                return(dataAll)
            })
            rownames(median.filtered.data[[i]]) <- rownames(object@centered.data)
            window.length <- window.length + length.iterations
        }
    }
    object@median.filtered.data <- median.filtered.data
    return(object)
}


CenterSmooth <- function(object) {
    object@center.smoothed.data <- list()
    for (i in 1:object@cnv.scale) {
        row_median <- apply(object@median.filtered.data[[i]], 2, median)
        object@center.smoothed.data[[i]] <- t(apply(object@median.filtered.data[[i]], 1, "-", row_median))
    }
    return(object)
    
}


ControlNormalize <- function(object, vis.bound, noise.thr) {
    object@control.normalized <- list()
    object@control.normalized.visbound.noiseRemoved <- list()
    object@control.normalized.visbound <- list()
    
    for (i in 1:object@cnv.scale) {
        object@control.normalized[[i]] <- AverageReference(data = object@center.smoothed.data[[i]], ref_ids = object@control.sample.ids)
        data <- object@control.normalized[[i]]
        data[data < (-vis.bound)] <- (-vis.bound)
        data[data > vis.bound] <- vis.bound
        object@control.normalized.visbound[[i]] <- data
        data[abs(data) < noise.thr] <- 0
        object@control.normalized.visbound.noiseRemoved[[i]] <- data
        
    }
    
    return(object)
    
}

AverageReference <- function(data, ref_ids) {
    average_reference_obs <- data[, ref_ids, drop = FALSE]
    # Reference gene within reference groups for (ref_group in ref_groups){
    grp_average <- rowMeans(average_reference_obs, na.rm = TRUE)
    
    data <- data - grp_average
    return(data)
}
