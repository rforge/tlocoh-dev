#' Convert isopleths to rasters
#'
#' Adds rasterized version of isopleth(s) to a LoCoH-hullset object 
#'
#' @param lhs A \link{LoCoH-hullset} object
#' @param id The id(s) of the hullsets to create isopleths for 
#' @param k The k value of hullsets to create isopleths for
#' @param r The r value of hullsets to create isopleths for
#' @param a The a value of hullsets to create isopleths for
#' @param s The s value of hullsets to create isopleths for
#' @param hs.names The name(s) of saved hullsets to create isopleths for
#' @param sort.metric The name(s) of hull metric(s) used to form isopleths that rasters should be created for
#' @param iso.method The method(s) used to define isopleths that will be converted to raster 
#' @param raster A RasterLayer object that will be used to set the extent and cell size of the rasterized isopleth
#' @param dimSize The number of cells along the largest dimension of the track. The according raster will be calculated internally. Default is 100. Ignored if \code{raster} is passed. 
#' @param cell.size The size of each square cell in map units. Ignored if \code{raster} is passed. 
#' @param sf.cell.size The number of significant figures to use if the cell size has be computed based on \code{dimSize}. Default=2. Ignored if \code{raster} is passed. 
#' @param ll.round Whether to round the lower left coordinate to the lowest multiple of \code{cell.size}. Ignored if \code{raster} is passed. (T/F)
#' @param status Show status messages. T/F
#'  
#' @return A LoCoH-hullset object
#'
#' @details
#' This will take exising isopleths and create raster versions of them.
#'
#' @examples
#' \dontrun{
#' lhs <- lhs.iso.add(lhs)
#' lhs <- lhs.iso.rast(lhs)
#' }
#'
#' @export
#' @seealso \code{\link{lhs.iso.add}}


lhs.iso.rast <- function(lhs, id=NULL, k=NULL, r=NULL, a=NULL, s=NULL, hs.names = NULL, 
                        sort.metric=NULL, iso.method="pt.quantiles", raster=NULL, dimSize=100,  
                        cell.size=NULL, sf.cell.size=2, ll.round=TRUE, status=TRUE) {
         
    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (!requireNamespace("raster")) stop("package raster required")
    if (!is.null(sort.metric)) {
        if (FALSE %in% (sort.metric %in% hm.expr(names.only=T, desc=F, print=F))) stop("Unknown value(s) for sort.metric")  
    }
    
    start.time <- Sys.time()
    if (is.null(id) && is.null(r) && is.null(k) && is.null(a) && is.null(s) && is.null(hs.names)) {
        hs.matching.idx <- 1:length(lhs)
    } else {    
        hs.matching.idx <- lhs.select.which(lhs, id=id, r=r, k=k, a=a, s=s, hs.names=hs.names)
    }
    if (length(hs.matching.idx)==0) stop("No sets of hulls found matching those criteria")
    
    isos.converted <- 0
    if (status) {cat("Converting isopleths to rasters\n"); flush.console()}
    
    for (hs.idx in hs.matching.idx) {
        if (status) {cat(names(lhs)[hs.idx] , "\n", sep=""); flush.console()}
        if (!is.null(lhs[[hs.idx]][["isos"]])) {
            for (iso.idx in 1:length(lhs[[hs.idx]][["isos"]])) {
                blnCont <- TRUE
                if (!is.null(sort.metric)) {
                    blnCont <- blnCont && lhs[[hs.idx]][["isos"]][[iso.idx]][["sort.metric"]] %in% sort.metric
                }
                if (!is.null(iso.method)) {
                  blnCont <- blnCont && lhs[[hs.idx]][["isos"]][[iso.idx]][["iso.method"]] %in% iso.method
                }
                
                if (blnCont) {
                    iso.rast <- iso2raster(lhs[[hs.idx]][["isos"]][[iso.idx]][["polys"]], raster=raster, dimSize=dimSize, cell.size=cell.size, sf.cell.size=sf.cell.size, ll.round=ll.round, status=status)
                    lhs[[hs.idx]][["isos"]][[iso.idx]][["rast"]] <- iso.rast
                    isos.converted <- isos.converted + 1
                }
            }
        
        }

            
    }
    
    if (status) {
        time.taken = difftime(Sys.time(), start.time, units="auto")
        cat("Total time:", round(time.taken,1), units(time.taken), "\n", sep = " ")    
    }
    
    return(lhs)  

}
