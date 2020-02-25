
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
    method, loh, project = "casperProject", matrix.type="raw", sequencing.type, expr.cutoff = 0.1, log.transformed = TRUE, 
    window.length = 50, length.iterations = 50,  genomeVersion = "hg19", filter= "median",
    ...) {
    
    object <- new(Class = "casper", raw.data = raw.data, loh = loh, annotation = annotation, sequencing.type = sequencing.type, 
        control.sample.ids = control.sample.ids, project.name = project,  cytoband = cytoband, loh.name.mapping = loh.name.mapping, 
        cnv.scale = cnv.scale, loh.scale = loh.scale, matrix.type = matrix.type, method = method, window.length = window.length, length.iterations = length.iterations, 
         genomeVersion = genomeVersion, filter = filter)

    object.raw.data <- object@raw.data
    if(matrix.type == "normalized") {
         gene.average <- rowMeans(((2^object.raw.data) - 1))
    }
    if(matrix.type == "raw") {
         gene.average <- rowMeans(object.raw.data)
    }

    indices <- which(gene.average < expr.cutoff)
    
    if(length(indices) >0) {
        object@data <- object.raw.data[-1 * indices, , drop=FALSE]
        object.data <- object@data
        object@annotation.filt <- object@annotation[-1 * indices, , drop=FALSE]
    } else {
        object@data <- object.raw.data
        object.data <- object@data
        object@annotation.filt <- object@annotation

    }
    
   if(is.null( object.data) || nrow( object.data) < 1 || ncol( object.data) < 1){  stop("Error, data has no rows or columns") }

    # processing steps adapted from infercnv R package
    if(matrix.type == "raw") {
        cs <- colSums(object.data)
        object.data <- sweep(object.data, STATS=cs, MARGIN=2, FUN="/")
        normalize_factor <-  median(cs)
        object.data  <- object.data * normalize_factor
        object.data <- log2(object.data + 1)
        #centered.data <-  t(scale(t(object.data)))
    }

    centered.data <-  t(apply(t(object.data), 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y))))

    object@centered.data <- AverageReference(data = centered.data, ref_ids = object@control.sample.ids)
    lower_bound <- mean(apply(object@centered.data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(object@centered.data, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))

    threshold = mean(abs(c(lower_bound, upper_bound)))

    object@centered.data[object@centered.data > threshold] <- threshold
    object@centered.data[object@centered.data < (-1 * threshold)] <- -1 * threshold

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
    
    if (object@filter=="median")   object <- PerformMedianFilterByChr(object)
    if (object@filter=="mean")   object <- PerformMeanFilterByChr(object)
    
    object <- CenterSmooth(object)
    
    object <- ControlNormalize(object)

    return(object)
}

# adapted from infercnv R package
#' @title PerformMeanFilter()
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
    
    filtered.data <- list()
    
    for (i in 1:object@cnv.scale) {
        if (i == 1) {
            filtered.data[[i]] <- apply(object@centered.data, 2, function(x) round(filter(MedianFilter(window.length + 
                1), x), digits = 2))
            rownames(filtered.data[[i]]) <- rownames(object@centered.data)
            window.length <- window.length + length.iterations
        } else {
            filtered.data[[i]] <- apply(filtered.data[[i - 1]], 2, function(x) round(filter(MedianFilter(window.length + 
                1), x), digits = 2))
            rownames(filtered.data[[i]]) <- rownames(object@centered.data)
            window.length <- window.length + length.iterations
        }
    }
    object@filtered.data <- filtered.data
    return(object)
}


#' @title PerformMeanFilterByChr()
#'
#' @description Recusive iterative median filtering is applied for each chromosome
#'
#' @param object casper object
#'
#' @return object
#'
#' @export
#'
 
PerformMeanFilterByChr <- function(object) {
    
    window.length = object@window.length 
    length.iterations = object@length.iterations

    filtered.data <- list()
    gene_chr_listing = object@annotation.filt[, "Chr"]
    chrs = unlist(unique(gene_chr_listing))

           
    for (i in 1:object@cnv.scale) {
          message("Performing Mean filtering...")
          message(paste0("Scale:", i, "..."))
          if (i == 1) { 
            filtered.data[[i]] <- object@centered.data
            for (chr in chrs) {
                chr_genes_indices <-  which(gene_chr_listing == chr)
             
                chr_data <-   filtered.data[[i]] [chr_genes_indices, , drop=FALSE]
                if (nrow(chr_data) > 1) {
                    chr_data <-  apply(chr_data, 2, caTools::runmean, k=window.length + 1)
                    filtered.data[[i]][chr_genes_indices, ] <- chr_data
                }
            }
            window.length <- window.length + length.iterations
        } else {
        
        filtered.data[[i]] <- filtered.data[[i-1]] 
            for (chr in chrs) {
                chr_genes_indices <-  which(gene_chr_listing == chr)
                chr_data <-  filtered.data[[i]][chr_genes_indices, , drop=FALSE]
                    if (nrow(chr_data) > 1) {
                        chr_data <- apply(chr_data, 2, caTools::runmean, k=window.length + 1)
                        filtered.data[[i]][chr_genes_indices, ] <- chr_data
                    }
            }
            window.length <- window.length + length.iterations
        } 
    }

    object@filtered.data <- filtered.data
    return(object)
}


#' @title PerformMedianFilterByChr()
#'
#' @description Recusive iterative median filtering is applied for each chromosome
#'
#' @param object casper object
#'
#' @return object
#'
#' @export
#'
 
PerformMedianFilterByChr <- function(object) {
    
    window.length = object@window.length 
    length.iterations = object@length.iterations

    filtered.data <- list()
    gene_chr_listing = object@annotation.filt[, "Chr"]
    chrs = unlist(unique(gene_chr_listing))

           
    for (i in 1:object@cnv.scale) {
          message("Performing Median Filtering...")
          message(paste0("Scale:", i, "..."))
          if (i == 1) { 
            filtered.data[[i]] <- object@centered.data
            for (chr in chrs) {
                chr_genes_indices <-  which(gene_chr_listing == chr)
             
                chr_data <-   filtered.data[[i]] [chr_genes_indices, , drop=FALSE]
                if (nrow(chr_data) > 1) {
                    chr_data <- apply(chr_data, 2, function(x) round(filter(MedianFilter(window.length + 1), x), digits = 2))
                        
                    filtered.data[[i]][chr_genes_indices, ] <- chr_data
                }
            }
            window.length <- window.length + length.iterations
        } else {
        
        filtered.data[[i]] <- filtered.data[[i-1]] 
            for (chr in chrs) {
                chr_genes_indices <-  which(gene_chr_listing == chr)
                chr_data <-  filtered.data[[i]][chr_genes_indices, , drop=FALSE]
                    if (nrow(chr_data) > 1) {
                        chr_data <- apply(chr_data, 2, function(x) round(filter(MedianFilter(window.length + 1), x), digits = 2))
                        filtered.data[[i]][chr_genes_indices, ] <- chr_data
                    }
            }
            window.length <- window.length + length.iterations
        } 
    }
    object@filtered.data <- filtered.data
    return(object)
}

 # processing steps adapted from infercnv R package
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
        row_median <- apply(object@filtered.data[[i]], 2, function(x) { median(x, na.rm=TRUE) } )
        object@center.smoothed.data[[i]] <- t(apply(object@filtered.data[[i]], 1, "-", row_median))
    }
    return(object)
    
}

 # processing steps adapted from infercnv R package
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
ControlNormalize <- function(object) {
    object@control.normalized <- list()
    object@control.normalized.noiseRemoved <- list()
    
    for (i in 1:object@cnv.scale) {
        object@control.normalized[[i]] <- AverageReference(data = object@center.smoothed.data[[i]], ref_ids = object@control.sample.ids)
        
        data <- 2^object@control.normalized[[i]]

        # adapted from infercnv R package
        lower_bound <- mean(apply(data, 2,
                                  function(x) quantile(x, na.rm=TRUE)[[1]]))
        upper_bound <- mean(apply(data, 2,
                                  function(x) quantile(x, na.rm=TRUE)[[5]]))

                                            # apply bounds
        data[data < lower_bound] <- lower_bound
        data[data > upper_bound] <- upper_bound

        
        vals = data[, colnames(data) %in% object@control.sample.ids]
        mean_ref_vals = mean(as.matrix(vals), na.rm=T)
        mean_ref_sd <- mean(apply(as.matrix(vals), 2, function(x) sd(x, na.rm = T))) * 1.5
        upper_bound = mean_ref_vals + mean_ref_sd
        lower_bound = mean_ref_vals - mean_ref_sd
        data[data > lower_bound & data < 
          data] = mean_ref_vals


        object@control.normalized.noiseRemoved[[i]] <- data
        
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
    average_reference_obs <- data[, colnames(data) %in% ref_ids, drop = FALSE]
    grp_average <- rowMeans(average_reference_obs, na.rm = TRUE)   
    data <- data - grp_average
    return(data)
}
