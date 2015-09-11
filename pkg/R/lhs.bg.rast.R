lhs.bg.rast <- function(lhs, type=c("roadmap", "satellite", "hybrid", "terrain")[2]) {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (!requireNamespace("raster", quietly=TRUE)) stop("package raster required, please install")
    if (!requireNamespace("dismo", quietly=TRUE)) stop("package dismo required, please install")
    if (!type %in% c("roadmap", "satellite", "hybrid", "terrain")) stop("Unknown value for type")
    
    bg.rast <- dismo::gmap(x, exp=1.05, type=type, filename='', style=NULL, scale=2, zoom=NULL, size=c(640, 640), rgb=FALSE, lonlat=TRUE)


}