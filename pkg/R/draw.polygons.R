#' Draw polygons on the plot window
#'
#' @description Returns a SpatialPolygonsDataFrame object containing polygons drawn with the mouse on a plot
#'
#' @param n Number of polygons to draw.
#' @param draw.reg Whether to draw the polygons on the plot. TRUE/FALSE.
#' @param col Color values for the polygons to be drawn. Ignored if \code{draw.reg=FALSE}.
#' @param alpha Transparency of the colors (0..1). Ignored if \code{draw.reg=FALSE}.
#' @param prompt.labels Whether to prompt for a label after each polygon is drawn (TRUE/FALSE). See Details.
#' @param ID A character vector of length \code{n} of the ID values for each polygon
#' @param proj4string An object of class \code{CRS} containing the coordinate system of the drawn polygons
#'
#' @details The number of polygons to be drawn must be specified by the argument \code{n}. If
#' \code{prompt.labels=FALSE}, default IDs will be constructed for each polygon.
#' 
#' \code{proj4string} can be used if the current plot window is displaying geographic data.
#'
#' @return An object of class \code{SpatialPolygonsDataFrame}
#'
#' @export

draw.polygons <- function(n, draw.reg=TRUE, col=NULL, alpha=0.5, prompt.labels=TRUE, ID=NULL, proj4string = CRS(as.character(NA))) {
    ## Returns a list with n items where each item is a list containing two elements:
    ##     1) a data frame of points of a closed polygon drawn by the user on the current plot window
    ##     2) a color object (string)
    ## The user will be given the chance to picks points for the polygon from the plot window using the mouse
    ## This requires a plot window to be active
    ## The purpose of this function is to help build a "regions" object which is one piece of a scatter plot 
    
    
    if (is.null(n)) stop("n (number of polygons) is a required argument")
    if (dev.cur()==1) stop("To use this function, a plot window must be active")
    blnRStudio <- names(dev.cur())=="RStudioGD"
    
    lst <- list()
    if (is.null(col)) {
        col <- rainbow(n, end=5/6)
        col.str <- rep("", n)
    } else {
        if (length(col) != n) stop("The number of colors must equal the number of regions to be drawn")
        col.str <- paste(" (", as.character(col), ")", sep="")
    }
    
    if (is.null(ID)) {
        ID <- 1:n
    } else {
        if (length(ID) != n) stop("ID should be the same length as a number of polygons")
    }
    
    if (blnRStudio) {
        msg_end <- "press escape"
    } else {
        msg_end <- "right-click and select 'stop'"
    }
    
    poly.labels <- NULL
    Polys.lst <- list()
    
    cat(cw(paste("To define a polygon region, click on the active plot window using a mouse.\nWhen finished, ", msg_end, ". The polygon will be automatically 'closed'\n\n", sep="")))
    for (i in 1:n) {
        cat("Please draw polygon #", i, " of ", n, col.str[i], "\n", sep="")
      	flush.console()
        pts <- tlocoh:::poly.from.plot(draw.poly=FALSE, status=FALSE)
        
        if (draw.reg) {
          col.use <- col2rgb(col[i])
          col.use <- rgb(red=col.use[1], green=col.use[2], blue=col.use[3], maxColorValue = 255, alpha=alpha * 255)
          polygon(pts, col=col.use)
        }
        
        if (prompt.labels) {
            if (!blnRStudio) bringToTop(-1)
            label <- readline(prompt = paste("Label for this region [", i, "]: ", sep=""))
            if (label=="") label <- as.character(i)
        } else {
            label <- as.character(i)
        }       
        poly.labels <- c(poly.labels, label)
        Polys.lst[[i]] <- Polygons(list(Polygon(coords=pts)), ID = ID[i])
    }
    
    # Create Spatial Polygons Data Frame
    sites.sp <- SpatialPolygons(Polys.lst, proj4string=proj4string)
    sites.spdf <- SpatialPolygonsDataFrame(sites.sp, data=data.frame(name=poly.labels, id=ID, col=col))
    
    return(sites.spdf)
}
