################################################################################ casper

#' The casper Class
#'
#' The casper object is required for single-cell CNV analysis. It stores all information
#' associated with the dataset, including data, annotations, analyes, etc. 
#'

casper <- methods::setClass("casper", slots = c(raw.data = "ANY", data = "ANY", median.filtered.data = "ANY", loh.median.filtered.data = "ANY", 
    centered.data = "ANY", center.smoothed.data = "ANY", control.normalized = "ANY", control.normalized.visbound = "ANY", control.normalized.visbound.noiseRemoved = "ANY", 
    large.scale.cnv.events = "ANY", segments = "ANY", cytoband = "ANY", annotation = "ANY", annotation.filt = "ANY", control.sample.ids = "character", 
    project.name = "character", genomeVersion = "ANY", version = "ANY", hmmparam = "ANY", plotorder = "ANY", vis.bound = "ANY", noise.thr = "ANY", 
    loh.name.mapping = "ANY", sequencing.type = "ANY", cnv.scale = "ANY", loh.scale = "ANY", method = "ANY", loh.shift.thr = "ANY", window.length = "ANY", 
    length.iterations = "ANY", loh = "ANY"))

setMethod(f = "show", signature = "casper", definition = function(object) {
    cat("An object of class", class(object), "in project", object@project.name, "\n", nrow(x = object@raw.data), "genes across", ncol(x = object@raw.data), 
        "samples.\n")
    invisible(x = NULL)
})
