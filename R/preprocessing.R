
#' @title CreateCasperObject
#'
#' @param raw.data  the matrix of genes (rows) vs. cells (columns) containing the raw counts
#' 
#' @param annotation data.frame containing positions of each gene along each chromosome in the genome
#'
#' @param control.sample.ids  vector containing the  reference (normal) cell (sample) names 
#'
#' @param cytoband cytoband information downloaded from UCSC hg19: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz hg38:http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
#' 
#' @param loh.name.mapping  contains the cell (sample) name and the matching baf signal sample name
#' 
#' @param cnv.scale  maximum expression scale
#' 
#' @param loh.scale  maximum baf scale
#' 
#' @param method analysis type: itereative  or fixed (default: iterative)
#' 
#' @param loh The original baf signal
#'
#' @param sequencing.type sequencing.type sequencing type: bulk or single-cell
#'
#' @param expr.cutoff expression cutoff for lowly expressed genes
#'
#' @param log.transformed indicates if the data log2 transformed or not. (default:TRUE)
#' 
#' @param centered.threshold
#' 
#' @param window.length window length used for median filtering (default: 50)
#' 
#' @param length.iterations increase in window length at each scale iteration (default: 50)
#' 
#' @param vis.bound threshold for control normalized data for better visualization (default: 2)
#' 
#' @param genomeVersion genomeVersion: hg19 or hg38 (default: hg19)
#' 
#' @description Creation of a casper object.
#'
#' @return casper
#'
#' @export
#'
CreateCasperObject <- function(raw.data, annotation, control.sample.ids, cytoband, loh.name.mapping, cnv.scale, loh.scale, 
    method, loh, project = "casperProject", sequencing.type, expr.cutoff = 4.5, display.progress = TRUE, log.transformed = TRUE, 
    centered.threshold = 3, window.length = 50, length.iterations = 50, vis.bound = 2, noise.thr = 0.3, genomeVersion = "hg19", 
    ...) {
    
    object <- new(Class = "casper", raw.data = raw.data, loh = loh, annotation = annotation, sequencing.type = sequencing.type, 
        control.sample.ids = control.sample.ids, project.name = project,  cytoband = cytoband, loh.name.mapping = loh.name.mapping, 
        cnv.scale = cnv.scale, loh.scale = loh.scale, method = method, window.length = window.length, length.iterations = length.iterations, 
        vis.bound = vis.bound, noise.thr = noise.thr, genomeVersion = genomeVersion)
    # filter cells on number of genes detected modifies the raw.data slot as well now
    if (!log.transformed) {
        object@raw.data <- log2(object@raw.data + 1)
    }
    object.raw.data <- object@raw.data
    
    gene.average <- log2(rowMeans((((2^object.raw.data) - 1) * 10), na.rm = TRUE) + 1)
    genes.use <- gene.average > expr.cutoff
    object@data <- object.raw.data[genes.use, ]
    object.data <- object@data
    object@annotation.filt <- object@annotation[genes.use, ]
    
    centered.data <- object.data - rowMeans(object.data, na.rm = TRUE)
    centered.data[centered.data > centered.threshold] <- centered.threshold
    centered.data[centered.data < (-centered.threshold)] <- (-centered.threshold)
    centered.data <- centered.data - rowMeans(centered.data, na.rm = TRUE)
    object@centered.data <- centered.data
    
    object <- ProcessData(object)
    
    return(object)
}

#' @title ProcessData()
#'
#' @description Processing expression signal. Step 1. Recursively iterative median filtering  Step 2. Center Normalization Step 3. Control Normalization
#'
#' @param object casper object
#' 
#' @return object
#'
#' @export
#'
#' 
ProcessData <- function(object) {
    
    object <- PerformMedianFilterByChr(object, window.length = object@window.length, length.iterations = object@length.iterations)
    
    object <- CenterSmooth(object)
    
    object <- ControlNormalize(object, vis.bound = object@vis.bound, noise.thr = object@noise.thr)
    # parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals('CreateCasperObject'))]
    # parameters.to.store$raw.data <- NULL
    return(object)
}

#' @title PerformMedianFilter()
#'
#' @description Recusive iterative median filtering is applied to whole genome 
#'
#' @param object casper object
#'
#' @param window.length window length used for median filtering
#' 
#' @param length.iterations increase in window length at each scale iteration 
#'
#' @return object
#'
#' @export

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

#' @title PerformMedianFilterByChr()
#'
#' @description Recusive iterative median filtering is applied for each chromosome
#'
#' @param object casper object
#'
#' @param window.length window length used for median filtering
#' 
#' @param length.iterations increase in window length at each scale iteration 
#'
#' @return object
#'
#' @export
#'
 
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


#' @title CenterSmooth()
#'
#' @description Cell centric expression centering is performed. For each cell (or sample), we compute the mid-point of the expression level then we subtract the mid-point expression from the expression levels of all the genes for the corresponding cell
#'
#' @param object casper object
#' 
#' @return object
#'
#' @export
#'
#' 
CenterSmooth <- function(object) {
    object@center.smoothed.data <- list()#' 
    for (i in 1:object@cnv.scale) {
        row_median <- apply(object@median.filtered.data[[i]], 2, median)
        object@center.smoothed.data[[i]] <- t(apply(object@median.filtered.data[[i]], 1, "-", row_median))
    }
    return(object)
    
}

#' @title ControlNormalize()
#'
#' @description  The control normalization is performed by subtracting reference expression values from the tumor expression values.
#'
#' @param object casper object
#' 
#' @return object
#'
#' @export
#'
#'
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


#' @title AverageReference()
#'
#' @description the mean the expression level for each gene across all the reference cells (samples) are computed.
#'
#' @param object casper object
#' 
#' @return object
#'
#' @export
#'
#'
AverageReference <- function(data, ref_ids) {
    average_reference_obs <- data[, ref_ids, drop = FALSE]
    grp_average <- rowMeans(average_reference_obs, na.rm = TRUE)   
    data <- data - grp_average
    return(data)
}
