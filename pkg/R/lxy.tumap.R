#' Computes the time use statistics for a LoCoH-xy object using a gridded surface
#'
#' @param lxy A LoCoH-xy object
#' @param id The name(s) of individuals to analyze
#' @param ivg Inter-visit gap
#' @param grid Grid type ('hex' or 'square')
#' @param cellsize Cell size in map units
#' @param mindim The minimum number of rows or columns
#' @param gridtype No longer used  Deprecated
#' @param stats Whether to display messages
#' @export

lxy.tumap <- function(lxy, id=NULL, ivg=NULL, grid=c("hex", "square")[1], cellsize=NULL, mindim=20, gridtype=NULL, status=TRUE) {

    if (!missing(gridtype)) {
      warning("argument 'gridtype' is deprecated; please use 'grid' instead.", call. = FALSE)
      grid <- gridtype
    }
    
    if (!require(pbapply)) error("package pbapply required")
    if (!inherits(lxy, "locoh.lxy")) stop("lxy should be of class \"locoh.lxy\"")
    if (is.character(grid)) {
        if (!grid %in% c("hex", "square")) stop("unknown value for grid")
    } else if (!inherits(grid, "SpatialPolygons")) {
        stop("unknown value for grid")
    }
    
    if (is.null(ivg)) stop("Need a value for the inter-visit gap (ivg)")
    if (is.null(lxy[["pts"]][["dt"]])) stop("timestamps needed for this function")
    
    if (is.null(id)) {
        id <- levels(lxy[["pts"]][["id"]])
    } else {
        if (FALSE %in% (id %in% levels(lxy[["pts"]][["id"]]))) stop(paste("Can't find the following ID(s) in lxy: ", paste(id[!(id %in% levels(lxy[["pts"]][["id"]]))], collapse=", ", sep=""), sep=""))
    }

    res <- list()
    prj <- lxy$pts@proj4string
    
    for (idVal in id) {
        
        idVal.idx <- which(lxy[["pts"]][["id"]] == idVal)
        idVal.num.pts <- length(idVal.idx)
        xys.idVal <- coordinates(lxy[["pts"]])[idVal.idx,,drop=FALSE]
        
        ## Compute the bounding box for this individual
        xrange <- range(xys.idVal[,1])
        yrange <- range(xys.idVal[,2])
        bbox <- rbind(xrange, yrange)

        if (identical(grid,"hex")) {
            thisgrid <- hexlayer(bbox, proj4string=prj, mindim=mindim, cellsize=cellsize, plotme=FALSE)
        } else if (identical(grid,"square")) {
            thisgrid <- gridlayer(bbox, proj4string=prj, mindim=mindim, cellsize=cellsize, plotme=FALSE)
        } else {
            thisgrid <- grid
        }
        
        tusdata <- NULL
        
        for (ivg.idx in 1:length(ivg)) {
            ivgVal <- ivg[ivg.idx]
            if (status) cat(ivg.idx, " of ", length(ivg), ". Computing the number of visits in each hull for ivg=", ivgVal, " (", secs.fmt(ivgVal), ")\n", sep="")
            
            ## First, we identify the points which fall in each grid cell
            if (status) cat("  Identifying enclosed points...", sep=""); flush.console()
            enc.pts.idx.lst <- over(thisgrid, as(lxy[["pts"]][idVal.idx,], "SpatialPoints"), returnList=TRUE)
            if (status) cat("Done.\n"); flush.console()
            
            ## Identify the cells that have at least one point enclosed
            enc.pts.idx.lst.nonzero <- (sapply(enc.pts.idx.lst, length) != 0)
            
            ## Create a blank list that will contain the visit number of points in each cell
            ivg.tab.lst <- vector("list", length(enc.pts.idx.lst))

            ## Next, we update the non-zero elements of ivg.tab.lst with a vector of the 'visit number' of each point in there
            ivg.tab.lst[enc.pts.idx.lst.nonzero] <- pblapply(enc.pts.idx.lst[enc.pts.idx.lst.nonzero], function(x) as.numeric(table(cumsum(c(1, diff(as.numeric(lxy[["pts"]][["dt"]][idVal.idx] [x])) >= ivgVal)))))
            
            ## Calculate the number of separate visits and the mean number of locations per visit
            nsv <- rep(0, length(enc.pts.idx.lst))
            nsv[enc.pts.idx.lst.nonzero] <- sapply(ivg.tab.lst[enc.pts.idx.lst.nonzero], length)
            mnlv <- rep(0, length(enc.pts.idx.lst))
            mnlv[enc.pts.idx.lst.nonzero] <- sapply(ivg.tab.lst[enc.pts.idx.lst.nonzero], mean)
            
            ## Save these columns to the tusdata data frame
            tus_thisivg <- data.frame(nsv, mnlv)
            names(tus_thisivg) <- paste0(c("nsv","mnlv"), ".", ivg)
            if (is.null(tusdata)) {
                tusdata <- tus_thisivg
            } else {
                tusdata <- cbind(tus_thisivg, tus_thisivg)
            }
        }
        
        ## Attach the dataframe to the SpatialPolygons object
        thisgrid <- SpatialPolygonsDataFrame(thisgrid, data=tusdata, match.ID = TRUE)
        
        ## Save this SpatialPolygonsDataFrame to the result list
        res[[idVal]] <- thisgrid

    }
    
    class(res) <- c("locoh.tumap", "list")
    
    return(res)
    
}