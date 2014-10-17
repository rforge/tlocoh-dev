#' Convert isopleths to raster
#'
#' Converts isopleths in a SpatialPolygonsDataFrame object to a RasterLayer
#' 
#' @param polys A SpatialPolygonsDataFrame containing isopleths sorted by isopleth level (smallest to largest)
#' @param raster A RasterStack object to be used to set the extent and resolution for the output raster.
#' @param ext An extent object or NULL. Ignored if \code{raster} is passed. 
#' @param dimSize Numeric value used as the number of cells along the largest dimension of the data is numeric. Ignored if \code{raster} is passed.
#' @param cell.size . Ignored if \code{raster} is passed.
#' @param sf.cell.size . Ignored if \code{raster} is passed.
#' @param ll.round Anchor the lower left coordinate of the raster extent to a multiple of the cell size. Ignored if \code{raster} is passed.
#' @param status Show status messages and progress bar
#' @param debug Show debugging info
#'
#' @details
#' This presumes the SPDF contains isopleths ordered from lowest level to highest. The cell values of the resulting raster will
#' sum up to the largest isopleth leve. In order for the resulting raster to sum to 1, the 100% isopleth must be part of the input.
#
#' @import sp

iso2raster <- function(polys, raster=NULL, ext=NULL, dimSize=100, cell.size=NULL, sf.cell.size=2, ll.round=TRUE, status=TRUE, debug=FALSE) {

    if (!requireNamespace("raster")) stop("raster package required")

    if (!is(polys, "SpatialPolygonsDataFrame")) stop("polys must be class SpatialPolygonsDataFrame")
    if (TRUE %in% diff(polys@data[["ptp"]]) < 0 ) stop("Isopleth levels must be increasing order to rasterize ")
      
    if (is.null(raster)) {    
    
      if (is.null(ext)) {
          ext <- extent(polys)
      } else {
          if (!is(ext, "Extent")) stop('ext must be of class extent')
          if (ll.round) {
            if (status) cat('Setting ll.round to false because ext is passed \n')
            ll.round <- FALSE
          }
      }
      crs <- polys@proj4string
      
      # Define the cell size
      if (is.null(cell.size)) {
        if (is.null(dimSize)) stop("Must pass cell.size or dimSize")
          cell.size = signif(min((ext@xmax - ext@xmin) / dimSize, (ext@ymax - ext@ymin) / dimSize), digits = sf.cell.size)
      } 
      
      # Next, in order to get nice neat coordinates for our cells, we'll round xmin and ymin
      # *down* to the nearest multiple of cell.size, and round xmax and ymax *up* to the nearest
      # multiple of cell.size.
      if (ll.round) {
          ext@xmin <- floor(ext@xmin / cell.size) * cell.size
          ext@ymin <- floor(ext@ymin / cell.size) * cell.size
          ext@xmax <- ceiling(ext@xmax / cell.size) * cell.size
          ext@ymax <- ceiling(ext@ymax / cell.size) * cell.size
      }
      
      # Now we can compute the number of rows and columns needed.
      nrow <- (ext@ymax - ext@ymin) / cell.size
      ncol <- (ext@xmax - ext@xmin) / cell.size
    
    } else {
      if (!is(raster, "RasterLayer")) stop("If used, raster must be class RasterLayer")
      if (!is.null(cell.size) || !is.null(ext)) stop(cw("If you pass a raster object, cell.size and extent must be left null", final.cr=FALSE))
      if (!identical(raster@crs, polys@proj4string)) stop(cw("Isopleth coordinate system is different than the raster coordinate system", final.cr=FALSE))
      ext <- raster@extent
      ncols <- raster@ncols
      nrows <- raster@nrows
      crs <- raster@crs
      
      if (!identical(polys@proj4string, raster@crs)) warning("Input raster CRS is different than isopleths")
    } 
    
    # Create a blank raster, the cell values will be NA (missing data).
    rast.iso <- raster(ext, ncol=ncol, nrow=nrow, crs=crs)
    rast.mask <- rast.iso
    rast.blank <- rast.iso

    ## Assign starting values
    rast.iso[] <- 0
    rast.mask[] <- 1

    mymask <- list()
    
    ## Keep track of proportion of total points that have been added to rast.iso
    ptp.total <- 0
    if (status) pb <- txtProgressBar(min = 0, max = nrow(polys), style = 3)
    
    for (i in 1:nrow(polys)) {
    
      if (status) setTxtProgressBar(pb, i)  
      
      ## Compute the total of proportion of total points that should be distributed within the remaining cells
      this.iso.level <- polys@data[i, "ptp"]
      ptp.this.layer <- this.iso.level - ptp.total

      ## If there are no new points enclosed by this isopleth level, skip it
      if (ptp.this.layer > 0) {
      
        ## Rasterize the polygons in this isopleth level, creating a boolean raster
        this.iso.rast <- rasterize(polys[i,], rast.blank, field=1, fun='count', background=0, silent=TRUE)
        
        ## Set to 0 any cells that have already been assigned values from previous
        ## iterations
        this.iso.rast <- this.iso.rast * rast.mask
        
        ## Get the number of cells not masked out by previous isopleths
        this.iso.sum <- cellStats(this.iso.rast, stat="sum")
        
        if (this.iso.sum==0) {
          ## No new area.
          ## so what we need to do is to put these points in the prevous iso points.
                    
          ## Spread the additional point among the last isopleth added
          this.iso.rast <- last.iso.rast.bln * (ptp.this.layer / last.iso.sum)
          
          ## Add values to rast.iso
          rast.iso <- rast.iso + this.iso.rast
  
        } else {
          
          ## Update mask for the next iteration of loop
          rast.mask <- rast.mask - this.iso.rast
          
          if (debug) {
            mymask[[as.character(i)]] <- rast.mask
            rast.mask.min <- cellStats(rast.mask, stat="min")
            cat(i, "rast.mask.min=", rast.mask.min, "\n")
            plot(this.iso.rast)
            browser()
          }
          
          ## Save a copy of this boolean raster  in case the next one doesn't include 
          ## any new additional cells
          last.iso.rast.bln <- this.iso.rast
          last.iso.sum <- this.iso.sum
          
          ## Assign values to the non-zero cells
          this.iso.rast <- this.iso.rast * (ptp.this.layer / this.iso.sum)
          
          ## Add values to rast.iso
          rast.iso <- rast.iso + this.iso.rast
        
        }
        
        ## Update ptp for next pass
        ptp.total <- this.iso.level
      }
      
    }
    if (status) close(pb)
    
    if (debug) {
      print("done with that, lets look at mymask");browser()      
      par(mfrow = n2mfrow(4))
      for (i in 1:9) plot(mymask[[i]])
    }
    
    return(rast.iso)

}