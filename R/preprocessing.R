# setwd('C:\\Users\\aharmanci\\Google Drive\\uthealth\\codebase_yale\\codebase_uthealth') install('CaSpER')
#' Initialize and setup the casper object
#'
#' Initializes the casper object and some optional filtering
#' @param raw.data Raw input data
#' @param project Project name (string)
#' @param min.cells Include genes with detected expression in at least this
#' many cells. Will subset the raw.data matrix as well. To reintroduce excluded
#' genes, create a new object with a lower cutoff.
#' @param min.genes Include cells where at least this many genes are detected.
#' @param is.expr Expression threshold for 'detected' gene. For most datasets, particularly UMI
#' datasets, will be set to 0 (default). If not, when initializing, this should be set to a level
#' based on pre-normalized counts (i.e. require at least 5 counts to be treated as expresesd) All
#' values less than this will be set to 0 (though maintained in object@raw.data).
#' @param normalization.method Method for cell normalization. Default is no normalization.
#' In this case, run NormalizeData later in the workflow. As a shortcut, you can specify a
#' normalization method (i.e. LogNormalize) here directly.
#' @param scale.factor If normalizing on the cell level, this sets the scale factor.
#' @param do.scale In object@@scale.data, perform row-scaling (gene-based
#' z-score). FALSE by default. In this case, run ScaleData later in the workflow. As a shortcut, you
#' can specify do.scale = TRUE (and do.center = TRUE) here.
#' @param do.center In object@@scale.data, perform row-centering (gene-based centering)
#' @param names.field For the initial identity class for each cell, choose this field from the
#' cell's column name
#' @param names.delim For the initial identity class for each cell, choose this delimiter from the
#' cell's column name
#' @param meta.data Additional metadata to add to the casper object. Should be a data frame where
#' the rows are cell names, and the columns are additional metadata fields
#' @param display.progress display progress bar for normalization and/or scaling procedure.
#' @param ... Ignored
#'
#' @return Returns a casper object with the raw data stored in object@@raw.data.
#' object@@data, object@@meta.data, object@@ident, also initialized.
#'
#' @import stringr
#' @import pbapply
#' @importFrom utils packageVersion
#' @importFrom Matrix colSums rowSums
#'
#' @export
#'
#' @examples
#' pbmc_raw <- read.table(
#'   file = system.file('extdata', 'pbmc_raw.txt', package = 'casper'),
#'   as.is = TRUE
#' )
#' pbmc_small <- CreateCasperObject(raw.data = pbmc_raw)
#' pbmc_small
#'
CreateCasperObject <- function(raw.data, annotation, control.sample.ids, cytoband, loh.name.mapping, cnv.scale, loh.scale, 
    method, loh, project = "casperProject", sequencing.type, expr.cutoff = 4.5, display.progress = TRUE, log.transformed = TRUE, 
    centered.threshold = 3, window.length = 50, length.iterations = 50, vis.bound = 2, noise.thr = 0.3, genomeVersion = "hg19", 
    ...) {
    
    casper.version <- packageVersion("casper")
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

#' Old R based implementation of ScaleData. Scales and centers the data
#'
#' @param object casper object
#' @param genes.use Vector of gene names to scale/center. Default is all genes in object@@data.
#' @param data.use Can optionally pass a matrix of data to scale, default is object@data[genes.use,]
#' @param do.scale Whether to scale the data.
#' @param do.center Whether to center the data.
#' @param scale.max Max value to accept for scaled data. The default is 10. Setting this can help
#' reduce the effects of genes that are only expressed in a very small number of cells.
#'
#' @return Returns a casper object with object@@scale.data updated with scaled and/or centered data.
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' \dontrun{
#' pbmc_small <- ScaleDataR(object = pbmc_small)
#' }
#'
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
