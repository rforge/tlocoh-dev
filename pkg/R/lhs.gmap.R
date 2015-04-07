#' Download background images for plots
#'
#' @param lhs A \link{LoCoH-hullset} object
#' @param id The name(s) of individuals to extract
#' @param k The k value of hullsets to extract
#' @param r The r value of hullsets to extract
#' @param a The a value of hullsets to extract
#' @param s The s value of hullsets to extract
#' @param hs.names The name(s) of saved hullsets to extract
#' @param gmap The type of image to download: \code{"roadmap"}, \code{"satellite"}, \code{"hybrid"}, or \code{"terrain"}
#' @param status Show status messages, T/F
#'
#' @return A list containing one element for each id in \code{lhs}. The list is of class \code{locoh.gmap}.
#' Each element is another named list containing i) a raster (in the same projection as \code{lhs}, ii) a vector
#' of color values, and iii) the type of image (e.g., hybrid, terrain, etc.).
#'
#' @details This function can be used to download images for plotting hullsets. By downloading them once and
#' storing them as a variable, plotting is faster. One image is obtained for each id (individual).
#'
#' @examples
#' \dontrun{
#' require(raster)
#' require(dismo)
#' toni.bg <- lhs.gmap(toni.lhs, type="terrain")
#' plot(toni.lhs, iso=T, gmap=toni.bg)
#' }
#' @export

lhs.gmap <- function(lhs, id=NULL, k=NULL, r=NULL, a=NULL, s=NULL, hs.names = NULL, gmap=c("roadmap", "satellite", "hybrid", "terrain")[3], status=TRUE) {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (!requireNamespace("raster", quietly=TRUE)) stop("package raster required, please install")
    if (!requireNamespace("dismo", quietly=TRUE)) stop("package dismo required, please install")
    if (!gmap %in% c("roadmap", "satellite", "hybrid", "terrain")) stop("Unknown value for gmap")

    if (is.null(id) && is.null(r) && is.null(k) && is.null(a) && is.null(s) && is.null(hs.names)) {
        hs.matching.idx <- 1:length(lhs)
    } else {
        hs.matching.idx <- tlocoh:::lhs.select.which(lhs, id=id, r=r, k=k, a=a, s=s, hs.names=hs.names)
    }
    if (length(hs.matching.idx)==0) stop("No sets of hulls found matching those criteria")
    
    ids.all <- unique(sapply(lhs, function(hs) hs$id))
    if (!is.null(id)) ids.all <- intersect(id, ids.all)
    if (length(ids.all)==0) stop("No matching ids found")

    res <- list()
    
    for (idVal in ids.all) {
      hs.idx <- which(sapply(lhs, function(hs) hs$id)==idVal)

      extLatLong <- raster::projectExtent(lhs[[hs.idx[1]]]$pts, CRS("+proj=longlat +datum=WGS84"))

      if (status) cat("Downloading a ", gmap, " background image for ", idVal, "...", sep="")
      ## Download a basemap from Google
      base.map.merc <- dismo::gmap(extLatLong, type=gmap)
      col.merc <- base.map.merc@legend@colortable

      ## Project the downloaded basemap
      base.map.rast <- raster::projectRaster(base.map.merc, crs=lhs[[1]]$pts@proj4string, method="ngb")
      if (status) cat("Done \n")

      res[[idVal]] <- list(bg.rast=base.map.rast, bg.col=col.merc, type=gmap)
    }
    
    class(res) <- c("locoh.gmap", "list")
    return(invisible(res))

}