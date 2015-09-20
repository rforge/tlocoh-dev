#' Returns a hexagonal grid covering a bounding box
#' If both mindim and cellsize are passed, cellsize trumps mindim

hexlayer <- function(bbox, proj4string=CRS(as.character(NA)), mindim=50, 
                     cellsize=NULL, plotme=FALSE, status=TRUE) {
    
    ## Return a hexagonal grid as a SpatialPolygons object encompassing bbox
    ## bbox is a matrix with columns min and max, and rows x and y
    
    if (!require(pbapply)) error("package pbapply required")
    
    if (is.null(mindim) && is.null(cellsize)) stop("Either mindim or cellsize is required")
    
    ## Compute the hex radius
    if (!is.null(cellsize)) {
        hexrad <- cellsize / 2
    } else {
        hexrad <- min(apply(bbox,1, diff))/ mindim   
    }
    
    ## Compute the lower left, reducing the y dimension by the hexradius
    ## so that all points along the lower boundary are included
    lower_left <- c(bbox[1,1], bbox[2,1] - hexrad / 2)
    
    ## Define a vector of angles that will be used to compute the nodes of each hexagon
    theta <- c(0:5,0) * pi / 3

    ## Compute the number of rows and colums needd
    ncol <- 1 + ceiling(diff(bbox[1,]) / (hexrad * (1 + cos(pi/3))))
    nrow <- 1 + ceiling(diff(bbox[2,]) / (hexrad * 2 * sin(pi/3)))
    if (status) cat("Constructing ", ncol * nrow, " hexagons \n", sep="")
    
    ## Compute the x and y offsets, making a two-column matrix
    xoffsets <- 0:(ncol-1) * hexrad * (cos(pi/3) + 1)
    yoffsets <- 0:(nrow-1) * hexrad * 2 * sin(pi/3)
    hexctr <- expand.grid(xoffsets, yoffsets)

    ## Increate the y-value of the even columns by rad * sin(60)
    even_cols <- which(hexctr[,1] %in% xoffsets[seq(from=2, to=ncol, by=2)])
    hexctr[even_cols, 2] <- hexctr[even_cols, 2] + hexrad * sin(pi/3)

    ## Make a list of Polygons objects
    Polys_lst <- pblapply(1:nrow(hexctr), function(i) sp::Polygons(list(Polygon(cbind(lower_left[1] + hexctr[i,1] + hexrad * cos(theta), lower_left[2] + hexctr[i,2] + hexrad * sin(theta)))), ID=as.character(i)))
    
    ## Package Polys_lst into a SpatialPolygons object
    Sr <- sp::SpatialPolygons(Polys_lst, proj4string=proj4string)
    
    ## Plot if needed
    if (plotme) {
        plot(Sr, axes=T, asp=1)
        points(x=bbox[1,c(1,2,2,1)], y=bbox[2,c(2,2,1,1)], pch=16, cex=2, col="blue")
        points(x=lower_left[1] + hexctr[,1], y=lower_left[2] + hexctr[,2], pch=16, col="red")
    }
    
    return(Sr)
    
}

