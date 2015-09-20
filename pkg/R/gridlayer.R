#' Returns a square grid covering a bounding box
#' If both mindim and cellsize are passed, cellsize trumps mindim

gridlayer <- function(bbox, proj4string=CRS(as.character(NA)), mindim=50, 
                      cellsize=NULL, plotme=FALSE, status=TRUE) {
    
    ## Return a square grid as a SpatialPolygons object encompassing bbox
    ## bbox is a matrix with columns min and max, and rows x and y
    
    require(sp)
    if (!require(pbapply)) error("package pbapply required")
    
    if (is.null(mindim) && is.null(cellsize)) stop("Either mindim or cellsize is required")

    ## Compute the cellsize
    if (is.null(cellsize)) cellsize <- min(apply(bbox,1, diff))/ mindim
    
    ## Compute the lower left
    lower_left <- c(bbox[1,1], bbox[2,1]) # - (cellsize / 2)
    
    ## Compute the number of rows and colums needd
    ncol <- 1 + ceiling(diff(bbox[1,]) / cellsize)
    nrow <- 1 + ceiling(diff(bbox[2,]) / cellsize)
    if (status) cat("Constructing ", ncol * nrow, " cells \n", sep="")
    
    ## Compute the x and y offsets, making a two-column matrix
    xoffsets <- 0:(ncol-1) * cellsize
    yoffsets <- 0:(nrow-1) * cellsize
    cellctr <- expand.grid(xoffsets, yoffsets)

    ## Make a list of Polygons objects
    Polys_lst <- pblapply(1:nrow(cellctr), function(i) sp::Polygons(list(Polygon(
        cbind(lower_left[1] + cellctr[i,1] + (cellsize / 2) * c(-1,1,1,-1), 
              lower_left[2] + cellctr[i,2] + (cellsize / 2) * c(1,1,-1,-1)))),
              ID=as.character(i)))
    
    ## Package Polys_lst into a SpatialPolygons object
    Sr <- sp::SpatialPolygons(Polys_lst, proj4string=proj4string)
    
    ## Plot if needed
    if (plotme) {
        plot(Sr, axes=T, asp=1)
        points(x=bbox[1,c(1,2,2,1)], y=bbox[2,c(2,2,1,1)], pch=16, cex=2, col="blue")
        points(x=lower_left[1] + cellctr[,1], y=lower_left[2] + cellctr[,2], pch=16, col="red")
    }
    
    return(Sr)
    
}

