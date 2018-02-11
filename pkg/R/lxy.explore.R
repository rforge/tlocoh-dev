#' Explore a LoCoH-xy object in a interactive map
#'
#' Explore a LoCoH-xy object in a zoomable clickable interactive map
#'
#' @param lxy A \link[tlocoh]{LoCoH-xy} object
#' @param id The name(s) of individuals to export
#' @param bg The background layer (i.e., tile) to display, or \code{'none'}
#' @param connect.dots Whether to draw a line between segments, T/F
#' @param line.popup Enable popup balloons on the line segments, T/F
#' @param width Width of the leaflet map including 'px' for units
#' @param height Height of the leaflet map including 'px' for units
#' @param status Display status messages, T/F
#'
#' @details This function displays a Locoh-xy object in an interactive window with a satellite image in the background.
#' You may zoom in and out, pan, and click on locations.
#'
#' To use this function, you must be using RStudio with the leaflet package installed. To display a satellite image in
#' the background, you must also be connected to the internet.
#'
#' @seealso \code{\link[tlocoh]{lxy.exp.kml}}
#'
#' @examples
#' \dontrun{
#' if (!require('devtools')) install.packages('devtools')
#' devtools::install_github('rstudio/leaflet')
#' mycon <- url("http://tlocoh.r-forge.r-project.org/toni.n5775.2005-08-22.2006-04-23.lxy.RData")
#' load(mycon)
#' close(mycon)
#' lxy.explore(toni.lxy, connect.dots=TRUE, line.popup=TRUE)
#' }
#'
#' @export


lxy.explore <- function(lxy, id=NULL, bg=c("esri_world_imagery","none")[1], connect.dots=FALSE, line.popup=FALSE, width="300px", height=width, status=TRUE) {

    if (!inherits(lxy, "locoh.lxy")) stop("lxy should be of class \"locoh.lxy\"")
    if (!requireNamespace("rgdal")) stop("package rgdal required")
    if (!require("leaflet")) stop("package leaflet required")
    if (!bg %in% c("esri_world_imagery","none")) stop("unknown option for bg")
    if (is.na(proj4string(lxy$pts))) stop("lxy doesn't have a coordinate reference system. See lxy.proj.add()")
    if (connect.dots && is.null(lxy$pts$dt)) stop("Time stamps not found, so can not connect the dots")

    ## Initialize a new leaflet object
    m <- leaflet::leaflet(width=width, height=height)

    ## Add tiles if needed
    if (bg != "none") {
        tiles.prp <- list()

        tiles.prp[["esri_world_imagery"]] <- list(
          tilesURL = "http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
          tilesAttr = "Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community"
          )

        tiles.prp[["mapquest_open"]] <- list(
            tilesURL = "http://otile{s}.mqcdn.com/tiles/1.0.0/sat/${z}/${x}/${y}.jpg",
            tilesAttr = "Portions Courtesy NASA/JPL-Caltech and U.S. Depart. of Agriculture, Farm Service Agency"
            )

        if (status) {
            cat("Downloading background map, please wait...\n")
            cat("(if your patience runs out, rerun with bg='none')\n")
            flush.console()
        }
        ## For other tile providers, see also http://leaflet-extras.github.io/leaflet-providers/preview/
        
        ## Add Tiles
        m <- m %>% leaflet::addTiles(tiles.prp[[bg]]$tilesURL, attribution = tiles.prp[[bg]]$tilesAttr)

    }

    id <- tlocoh::vectorize.parameter(id, type="character", sort=FALSE)
    if (is.null(id)) {
        id <- levels(lxy[["pts"]][["id"]])
        idx.all <- 1:nrow(lxy[["pts"]])
    } else {
        if (FALSE %in% (id %in% levels(lxy[["pts"]][["id"]]))) stop("id value(s) not found")
        idx.all <- which(lxy[["pts"]][["id"]] %in% id)
    }

    ## Get pts and project to lat-long
    crs.latlong <- CRS("+proj=longlat +datum=WGS84")
    pts.latlong <- spTransform(as(lxy$pts[idx.all,], "SpatialPoints"), crs.latlong)

    ## Set the map bounds
    pts.bbox <- bbox(pts.latlong)
    m <- m %>% leaflet::fitBounds(pts.bbox[1,1], pts.bbox[2,1], pts.bbox[1,2], pts.bbox[2,2])

    ## Create colors that will be used for points
    if (length(id)==1) {
        pts.cols <- topo.colors(length(pts.latlong))
    } else {

        ## Create a vector that will map the id level of each point to a color
        id.levels <- levels(lxy$pts@data[["id"]])
        idlev.col <- numeric(length(id.levels))
        for (k in 1:length(id)) {
            idlev.col[which(id.levels ==id[k])] <- k
        }
        
        ## Create the vector of colors
        pts.cols <- rainbow(length(id), end=5/6)[ idlev.col[as.numeric(lxy$pts@data[idx.all, "id"])]]
        
    }

    ## Remove final two characters of color code which leaflet doesn't recognize
    pts.cols <- substr(pts.cols, 1, 7)

    ## Create the line segments if needed
    if (connect.dots) {

        if (status) {
            cat("Creating SpatialLines...")
            flush.console()
            start.time <- Sys.time()
        }

        ## Create an empty variable to hold the line segements as a SpatialLines object
        segmentsall.sl <- NULL
        segmentsall.Lines <- NULL
        segmentsall.popup <- NULL
        
        ## Create an empty vector to hold the segment colors (used only if plotting multiple ids)
        if (length(id) > 1) {
            segmentsall.col <- NULL
            col.id <- rainbow(length(id), end=5/6)
            col.id <- substr(col.id, 1, 7)
        }

        seg.idx.base <- 0

        for (j in 1:length(id)) {
            idVal <- id[j]

            ## Get the coordinates for this individual
            idVal.idx <- which(lxy$pts@data[idx.all,"id"]==idVal)
            pts.coords.thisid <- coordinates(pts.latlong)[idVal.idx,]
            
            ## Make list of Lines objects
            seg.start.idx <- 1:(nrow(pts.coords.thisid)-1)
            segments.Lines <- lapply(seg.start.idx, function(i) Lines(list(Line(pts.coords.thisid[c(i,i+1),])), ID=as.character(seg.idx.base + i)) )
            
            ## Append to master list of Lines objects
            segmentsall.Lines = append(segmentsall.Lines, segments.Lines)

            ## Create the segment popup
            if (connect.dots && line.popup) {
                #segmentsall.popup <- c(segmentsall.popup, paste(pts.dt[idVal.idx][seg.start.idx], " - ", pts.dt[idVal.idx][seg.start.idx+1], sep=""))
                #segmentsall.popup <- c(segmentsall.popup, paste(lxy$pts$dt[idx.all][idVal.idx][seg.start.idx], " - ", lxy$pts$dt[idx.all][idVal.idx][seg.start.idx+1], sep=""))
                #segmentsall.popup <- c(segmentsall.popup, paste(pts.ptid[idVal.idx][seg.start.idx], " - ", pts.ptid[idVal.idx][seg.start.idx+1], sep=""))
                
                #WORKS: segmentsall.popup <- c(segmentsall.popup, paste("<table cellspacing=5><tr align='center'><td>", format(lxy$pts@data[idx.all[idVal.idx][seg.start.idx], "dt"], "%d-%b-%Y<br>%H:%M:%S"), "</td><td>&ndash;</td><td>", format(lxy$pts@data[idx.all[idVal.idx][seg.start.idx+1], "dt"], "%d-%b-%Y<br>%H:%M:%S"), "</td></tr></table>", sep=""))
                segmentsall.popup <- c(segmentsall.popup, paste("<table><tr align='center'><td colspan=3><strong>", idVal, "</strong></td></tr><tr align='center'><td>", format(lxy$pts@data[idx.all[idVal.idx][seg.start.idx], "dt"], "%d-%b-%y<br>%H:%M:%S"), "</td><td>&mdash;</td><td>", format(lxy$pts@data[idx.all[idVal.idx][seg.start.idx+1], "dt"], "%d-%b-%y<br>%H:%M:%S"), "</td></tr></table>", sep=""))
            }

            ## Add color values to the segment color vector
            if (length(id) > 1) segmentsall.col <- c(segmentsall.col, rep(col.id[j], nrow(pts.coords.thisid)-1))
            
            seg.idx.base <- seg.idx.base + length(seg.start.idx)
        }
        
        ## Create a SpatialLines object
        segmentsall.sl = SpatialLines(segmentsall.Lines, proj4string=crs.latlong)

        if (status) {
            time.taken <- difftime(Sys.time(), start.time, units="auto")
            cat("Done. Time:", round(time.taken,1), units(time.taken), "\n", sep = " ")

            cat("Adding polylines to map..."); flush.console()
            start.time <- Sys.time()

        }

        ## Define color vector for the segments
        if (length(id)==1) {
            seg.col <- pts.cols[1:(length(pts.cols)-1)]
        } else {
            seg.col <- segmentsall.col
        }

        ## Set the popup options
        if (line.popup) {
            line.options <- leaflet::pathOptions()
        } else {
            line.options <- list(clickable=FALSE)
        }

        ## Add polylines to map
        m <- m %>% leaflet::addPolylines(data=segmentsall.sl, weight=1, color=seg.col, popup=segmentsall.popup, options=line.options)
        
        if (status) {
            time.taken <- difftime(Sys.time(), start.time, units="auto")
            cat("Done. Time:", round(time.taken,1), units(time.taken), "\n", sep = " "); flush.console()
        }
        
    }

    ## Add circle markers to map
    if (status) cat("Adding circle markers..."); flush.console()
    
    ## Define balloon text
    balloon_text <- paste(lxy$pts@data[idx.all, "id"], lxy$pts@data[idx.all, "ptid"], sep=":")
    if (!is.null(lxy$pts[["dt"]])) balloon_text <- paste("<center><strong>", balloon_text, "</strong><br>", format(lxy$pts@data[idx.all, "dt"], "%d-%b-%Y<br>%H:%M:%S"), "</center>", sep="")

    ## Add circle markers
    #m <- m %>% addCircleMarkers(data=pts.latlong, radius=2, stroke=FALSE, fillColor=pts.cols, fillOpacity=1, popup=htmltools::htmlEscape(as.character(balloon_text)))
    m <- m %>% addCircleMarkers(data=pts.latlong, radius=2, stroke=FALSE, fillColor=pts.cols, fillOpacity=1, popup=as.character(balloon_text))

    if (status) {
        cat("Done\n")
        cat("Preparing map for display \n")
    }

    ## Display the map
    ## print(m)
    
    ## Return the leaflet object
    return(m)


}