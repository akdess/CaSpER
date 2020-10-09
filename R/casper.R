

#' @details
#' The main functions you will need to use are CreateCasperObject() and runCaSpER(casper_object).
#' For additional details on running the analysis step by step, please refer to the example vignette.
#' @aliases CaSpER-package
"_PACKAGE"


#' The CaSpER Class
#'
#'
#' The CaSpER Class
#' The casper object is required for performing CNV analysis on single-cell and bulk RNA-Seq. It stores all information
#' associated with the dataset, including data, smoothed data, baf values, annotations, scale specific segments, scale specific large scale events etc. 
#'
#' @slot raw.data  raw project data
#' 
#' @slot data lowly expressed genes are filtered from the data
#' 
#' @slot loh  original baf signal 
#' 
#' @slot filtered.data  median filtered expression signal 
#' 
#' @slot loh.median.filtered.data  median filtered baf signal 
#' 
#' @slot centered.data  gene expression levels are centered around the mid-point. For each gene, the mid-point of expression level is computed among all the cells (or samples in bulk RNA-seq), then the mid-point expression level is subtracted from the expression levels
#' 
#' @slot center.smoothed.data  cell centric expression centering is performed. For each cell (or sample), we compute the mid-point of the expression level then we subtract the mid-point expression from the expression levels of all the genes for the corresponding cel
#' 
#' @slot control.normalized  control normalization is performed by subtracting reference expression values from the tumor expression values. 
#' 
#' @slot control.normalized.noiseRemoved  noise is removed from control normalized and thresholded data.
#' 
#' @slot large.scale.cnv.events large scale CNV events identified by CaSpER
#' 
#' @slot segments CNV segments identified  by CaSpER
#' 
#' @slot cytoband cytoband information downloaded from UCSC hg19: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz hg38:http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
#' 
#' @slot annotation positions of each gene along each chromosome in the genome
#' 
#' @slot annotation.filt  lowly expressed genes are filtered from gene annotation data.frame
#' 
#' @slot control.sample.ids vector containing the  reference (normal) cell (sample) names 
#' 
#' @slot project.name project name
#' 
#' @slot genomeVersion genomeVersion: hg19 or hg38
#' 
#' @slot hmmparam initial hmm parameters estimated from data
#' 
#' @slot plotorder cell (sample) ordering for heatmap plots 
#' 
#' @slot vis.bound threshold for control normalized data for better visualization
#' 
#' @slot noise.thr noise threshold for better visualization
#' 
#' @slot loh.name.mapping containing the cell (sample) name and the matching baf signal sample name
#' 
#' @slot sequencing.type sequencing type: bulk or single-cell
#' 
#' @slot cnv.scale maximum expression scale
#' 
#' @slot loh.scale maximum baf scale
#' 
#' @slot loh.shift.thr baf shift threshold estimated from baf signal using gaussian mixture models
#' 
#' @slot window.length window length used for median filtering
#' 
#' @slot length.iterations increase in window length at each scale iteration 
#' 
#' 
#' @name casper
#' @rdname casper
#' @aliases casper-class
#' @exportClass casper

casper <- methods::setClass("casper", slots = c(raw.data = "ANY", data = "ANY", filtered.data = "ANY", loh.median.filtered.data = "ANY", control.normalized.visbound.noiseRemoved="ANY",
    centered.data = "ANY", matrix.type="ANY", center.smoothed.data = "ANY", control.normalized = "ANY", control.normalized.noiseRemoved = "ANY", 
    large.scale.cnv.events = "ANY", segments = "ANY", cytoband = "ANY", annotation = "ANY", annotation.filt = "ANY", control.sample.ids = "character", 
    project.name = "character", genomeVersion = "ANY", hmmparam = "ANY", plotorder = "ANY",
	 loh.name.mapping = "ANY", sequencing.type = "ANY", cnv.scale = "ANY", loh.scale = "ANY", method = "ANY", 
    loh.shift.thr = "ANY", window.length = "ANY", length.iterations = "ANY", loh = "ANY", filter= "ANY"))

setMethod(f = "show", signature = "casper", definition = function(object) {
    cat("An object of class", class(object), "in project", object@project.name, "\n", nrow(x = object@raw.data), "genes across", 
        ncol(x = object@raw.data), "samples.\n")
    invisible(x = NULL)
})
