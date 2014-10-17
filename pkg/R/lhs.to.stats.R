#' Computes statistics on the centroid distances for time-overlapped hulls
#'
#' Computes statistics on the hull-to-hull centroid distance for time-overlapped hulls, with an option to plot and overlay the distribution of centroid distances for random pairs of hulls 
#'
#' @param lhs A LoCoH-hullset object
#' @param id1 Hullset 1 id value(s). Can also be \code{'all'}.
#' @param id2 Hullset 2 id value(s). Can also be \code{'all'}.
#' @param n The number of randomly selected paired hulls to use as the NULL model of no association. Can also be "all"
#' @param iso.oz Apply isopleth filter to hs1 using the isopleths from hs2. T/F.
#' @param iso.lower The lower level isopleth for the isopleth filter
#' @param iso.upper The upper level isopleth for the isopleth filter
#' @param iso.sort.metric The name of a hull metric that was used to sort hulls in the construction of the
#' isopleths to be used as filters. If \code{auto} (default) it will pick the default sort metric used for density
#' isopleths (i.e., area for the k-method, and number of enclosed points for the a and r methods)
#' @param to.comp.hist Draw a histogram of the centroid distances of temporally overlapping hulls
#' @param breaks The number of breaks in the histgram (or another valid value for breaks, see \code{\link{hist}}).
#' @param to.mcd.outline.only Show the outline only of the histogram of the mean centroid distance for temporally overlapping hulls. T/F.
#' @param lwd.outline The line width of the histogram outline (ignored if \code{to.mcd.outline.only=F}).
#' @param hist.type The type of histogram to plot: 'density' or 'counts'.
#' @param col.to.mcd The color of the outline of the histogram of the distribution of the centroid distances of temporally overlapping hulls.
#' @param col.h2h.cd The color of the outline of the histogram of the distribution of the centroid distances for randomly paired hulls.
#' @param title The title to be displayed. Character. If NULL a title will be constructed.
#' @param title.show Whether to show the title. T/F.
#' @param title.id.only Whether to construct the title from the id values only. T/F. Ignored if \code{title} is passed or \code{title.show=FALSE}
#' @param title.sub.iso.enc Whether to include the isopleth filter information as the second line of the title, T/F.
#' @param mar The plot margins. A four item numeric vector
#' @param mgp The distance away from the edge of the plot for the 1) label, 2) tick marks, and 3) axis line. A three-item numeric vector
#' @param figs.per.page The number of plots per page.
#' @param panel.num A number or letter to display in the upper left hand corner of the plot when the plot will be used as part of a
#' multi-frame graphic (as in publications). Character
#' @param panel.num.inside.plot Whether to display panel.num inside the plot area itself, as opposed to the title area. T/F
#' @param png.dir The directory for a PNG file (filename will be constructed automatically). Ignored if png.fn is passed
#' @param png.dir.make Whether to create png.dir if it doesn't exist. T/F
#' @param png.width The width of the PNG image
#' @param png.height The height of the PNG image
#' @param png.overwrite Whether to overwrite an existing PNG file if it exists. T/F
#' @param png.pointsize The pointsize (in pixels) for the PNG image, equivalent to the height or width of a character in 
#' pixels (increase to make labels appear larger)
#' @param status Display status messages. T/F
#' @param ... Additional parameters that will be passed to the \code{\link{plot}} function
#'
#' @return A list object
#' @details
#' This returns a list object containing the centroid-to-centroid distances of a random selection of hulls from two individuals, which serves
#' as a NULL model of no interaction. It can also plot the histogram of mean-centroid-distance of time-overlapped hulls with the outline of 
#' the centroid distance of random-pairs of hulls overlain on top, to visually see how close the distributions match.
#'
#' It will also compute the Welch Two Sample t-test to see if the distribution for time-overlapped and randomly paired hulls have statistically significant means, 
#' and the Two-sample Kolmogorov-Smirnov test which tells you how likely the two distributions are the same.
#'
#' Note that before you can use this function, the mean-centroid-distance for time-overlapped hulls must be computed using \code{\link{lhs.to.add}}
#'
#' You can apply an isopleth filter by passing values for iso.lower and iso.upper. These should be the isopleth level (normally between 0 and 1) that the hull parent point must fall in to be included in the analysis. If, for example, you wanted to see whether association in the core area is significantly different from nuetral interaction, you would pass \code{iso.upper=0.5} and leave \code{iso.lower} NULL. Note that the isopleths with matching isopleth levels must already be present.
#'
#' @seealso \code{\link{lhs.to.add}}, \code{\link{lhs.iso.add}}
#' @export


lhs.to.stats <- function(lhs, id1="all", id2="all", n="all", iso.lower=NULL, iso.upper=NULL, iso.oz=TRUE, iso.sort.metric="auto", 
                       to.comp.hist=TRUE, breaks=20, to.mcd.outline.only=FALSE, lwd.outline=3, hist.type=c("density","counts")[1],
                       col.to.mcd="blue", col.h2h.cd="red", 
                       title=NULL, title.show=TRUE, title.id.only=FALSE, title.sub.iso.enc=TRUE, 
                       mar=c(3,3, if (title.show) 1.5 + (if(title.sub.iso.enc) 1.3 else 0) else 0.5, 0.5), mgp=c(1.8, 0.5, 0), figs.per.page=1,
                       panel.num=NULL, panel.num.inside.plot=!title.show, 
                       png.dir=NULL, png.dir.make=TRUE, png.width=800, png.height=png.width, png.overwrite=TRUE, png.pointsize=12+(png.width-480)/80, 
                       status=TRUE, ...) {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (!hist.type %in% c("density","counts")) stop("Unknown value for hist.type")
    
    ## Make a data frame of all of the hullsets and their properties
    lhs.info <- do.call(rbind, lapply(lhs, function(x) data.frame(id=x$id, s=x$s, k=n2z(x$k), a=n2z(x$a), r=n2z(x$r), stringsAsFactors=FALSE)))
    id.all <- unique(lhs.info[["id"]])
    if (length(id.all)==1) stop("This hull metric requires a hullset collection with more than one ID")

    if (identical(id1, "all")) id1 <- id.all
    if (identical(id2, "all")) id2 <- id.all
    res <- list()
    pngs.made <- NULL
        
    for (idVal1 in id1) {
        ## Save the indices of the hullset(s) which have this id
        hsi.idVal1 <- which(lhs.info[["id"]] == idVal1)
        
        for (h1 in hsi.idVal1) {
            if (status) cat(names(lhs)[h1], "\n");flush.console()

            ## Apply isopleth filter
            if (!is.null(iso.lower) || !is.null(iso.upper)) {
                if (!is.null(iso.lower) && !is.null(iso.upper) && iso.lower >= iso.upper) stop("iso.lower must be less than iso.upper")
                
                ## Grab the correct iso
                if (iso.sort.metric == "auto") {
                    iso.sort.metric.use <- if (lhs[[h1]][["mode"]] == "k") "area" else "nep"
                } else {
                    iso.sort.metric.use <- iso.sort.metric
                
                }
                
                iso.idx <- which(sapply(lhs[[h1]][["isos"]], function(x) x$sort.metric == iso.sort.metric.use))[1]
                if (is.na(iso.idx)) stop("Cant find isopleth")
                
                ## Create a SpatialPoints object for the hull parent points, will need it for the overlay
                h1.hpp.sp <- as(lhs[[h1]][["pts"]][lhs[[h1]][["hulls"]][["pts.idx"]],], "SpatialPoints")
                
                if (is.null(iso.upper)) {
                    hpp1.in.iso1.upper <- rep(TRUE, nrow(lhs[[h1]][["hulls"]]))
                } else {
                    iso.poly.idx <- which(sapply(lhs[[h1]][["isos"]][[iso.idx]][["polys"]][["iso.level"]], function(x) isTRUE(all.equal(x, iso.upper))))
                    if (length(iso.poly.idx)==0) stop("Can't find an isopleth matching iso.upper")
                    h1.iso.up.sp <- as(lhs[[h1]][["isos"]][[iso.idx]][["polys"]][iso.poly.idx,], "SpatialPolygons")
                    hpp1.in.iso1.upper <- !is.na(over(h1.hpp.sp, h1.iso.up.sp))
                }
                
                if (is.null(iso.lower)) {
                    hpp1.in.iso1.lower <- rep(FALSE, nrow(lhs[[h1]][["hulls"]]))
                } else {
                    iso.poly.idx <- which(sapply(lhs[[h1]][["isos"]][[iso.idx]][["polys"]][["iso.level"]], function(x) isTRUE(all.equal(x, iso.lower))))
                    if (length(iso.poly.idx)==0) stop("Can't find an isopleth matching iso.lower")
                    h1.iso.low.sp <- as(lhs[[h1]][["isos"]][[iso.idx]][["polys"]][iso.poly.idx,], "SpatialPolygons")
                    hpp1.in.iso1.lower <- !is.na(over(h1.hpp.sp, h1.iso.low.sp))
                }
                
                h1.hulls.to.include <- hpp1.in.iso1.upper
                h1.hulls.to.include[hpp1.in.iso1.lower] <- FALSE
                
                #rm(iso.sp, hpp1.in.upper, hpp1.in.iso1.lower)
                
            
            } else {
                ## No iso filter
                h1.hulls.to.include <- rep(TRUE, nrow(lhs[[h1]][["hulls"]]))
            }

            hs1.name <- paste(idVal1, ".s", lhs[[h1]][["s"]], ".", lhs[[h1]][["mode"]], lhs[[h1]][[  lhs[[h1]][["mode"]]  ]], sep="")
            
            ## Get centroid of each hull saved as a xy matrix
            h1.ctr <- t(sapply(lhs[[h1]][["hulls"]]@polygons, function(x) apply(x@Polygons[[1]]@coords,2,mean)))

            for (idVal2 in id2) {
                if (idVal2 != idVal1) {
                    hsi.idVal2 <- which(lhs.info[["id"]] == idVal2)
                    for (h2 in hsi.idVal2) {                    
                        if (status) cat(" - ", names(lhs)[h2], "\n");flush.console()
                        
                        ## Apply isopleth filter to hpp1 using the isopleths from h2 if iso.oz = T
                        if (iso.oz && (!is.null(iso.lower) || !is.null(iso.upper))) {
                            
                            ## Grab the correct iso
                            if (iso.sort.metric == "auto") {
                                iso.sort.metric.use <- if (lhs[[h2]][["mode"]] == "k") "area" else "nep"
                            } else {
                                iso.sort.metric.use <- iso.sort.metric
                            }
                            
                            iso.idx <- which(sapply(lhs[[h2]][["isos"]], function(x) x$sort.metric == iso.sort.metric.use))[1]
                            if (is.na(iso.idx)) stop("Cant find isopleth")
                            
                            ## Create a SpatialPoints object for the hull parent points, will need it for the overlay
                            ### hpp.sp <- as(lhs[[h1]][["pts"]][lhs[[h1]][["hulls"]][["pts.idx"]],], "SpatialPoints")
                            
                            if (is.null(iso.upper)) {
                                hpp1.in.iso2.upper <- rep(TRUE, nrow(lhs[[h1]][["hulls"]]))
                            } else {
                                iso.poly.idx <- which(sapply(lhs[[h2]][["isos"]][[iso.idx]][["polys"]][["iso.level"]], function(x) isTRUE(all.equal(x, iso.upper))))
                                if (length(iso.poly.idx)==0) stop("Can't find an isopleth matching iso.upper")
                                h2.iso.up.sp <- as(lhs[[h2]][["isos"]][[iso.idx]][["polys"]][iso.poly.idx,], "SpatialPolygons")
                                hpp1.in.iso2.upper <- !is.na(over(h1.hpp.sp, h2.iso.up.sp))
                            }
                            
                            if (is.null(iso.lower)) {
                                hpp1.in.iso2.lower <- rep(FALSE, nrow(lhs[[h1]][["hulls"]]))
                            } else {
                                iso.poly.idx <- which(sapply(lhs[[h2]][["isos"]][[iso.idx]][["polys"]][["iso.level"]], function(x) isTRUE(all.equal(x, iso.lower))))
                                if (length(iso.poly.idx)==0) stop("Can't find an isopleth matching iso.lower")
                                h2.iso.low.sp <- as(lhs[[h2]][["isos"]][[iso.idx]][["polys"]][iso.poly.idx,], "SpatialPolygons")
                                hpp1.in.iso2.lower <- !is.na(over(h1.hpp.sp, h2.iso.low.sp))
                            }
                            
                            h1.hulls.in.iso2.ul <- hpp1.in.iso2.upper
                            h1.hulls.in.iso2.ul[hpp1.in.iso2.lower] <- FALSE
                            
                            #print("here is the gist");browser()
                            
                            #rm(h1.iso.sp, hpp.sp, hpp1.in.iso2.upper, hpp1.in.iso1.lower)
                            h1.hulls.to.include <- h1.hulls.to.include & h1.hulls.in.iso2.ul 
                        
                        } 


                        ## Apply isopleth filter to hpp2 if iso.oz = T
                        if (iso.oz && (!is.null(iso.lower) || !is.null(iso.upper))) {

                            ## Create a SpatialPoints object for the hull parent points, will need it for the overlay
                            h2.hpp.sp <- as(lhs[[h2]][["pts"]][lhs[[h2]][["hulls"]][["pts.idx"]],], "SpatialPoints")

                            if (is.null(iso.upper)) {
                                hpp2.in.iso1.upper <- hpp2.in.iso2.upper <- rep(TRUE, nrow(lhs[[h2]][["hulls"]]))
                            } else {
                                hpp2.in.iso1.upper <- !is.na(over(h2.hpp.sp, h1.iso.up.sp))
                                hpp2.in.iso2.upper <- !is.na(over(h2.hpp.sp, h2.iso.up.sp))
                            }
                            
                            if (is.null(iso.lower)) {
                                hpp2.in.iso1.lower <- hpp2.in.iso2.lower <- rep(FALSE, nrow(lhs[[h2]][["hulls"]]))
                            } else {
                                hpp2.in.iso1.lower <- !is.na(over(h2.hpp.sp, h1.iso.low.sp))
                                hpp2.in.iso2.lower <- !is.na(over(h2.hpp.sp, h2.iso.low.sp))
                            }
                            
                            h2.hulls.to.include <- hpp2.in.iso1.upper
                            h2.hulls.to.include[hpp2.in.iso1.lower] <- FALSE

                            h2.hulls.in.iso2.ul <- hpp2.in.iso2.upper
                            h2.hulls.in.iso2.ul[hpp2.in.iso2.lower] <- FALSE
                            
                            #print("here is the gist");browser()
                            #rm(h1.iso.sp, hpp.sp, hpp1.in.iso2.upper, hpp1.in.iso1.lower)
                            
                            h2.hulls.to.include <- h2.hulls.to.include & h2.hulls.in.iso2.ul 
                        } else {
                            h2.hulls.to.include <- rep(TRUE, nrow(lhs[[h2]][["hulls"]]))
                        }

                        if (sum(h1.hulls.to.include)==0) {
                            cat("No hull parent points from h1 meet the criteria \n"); flush.console()
                        } else if (sum(h2.hulls.to.include)==0) {
                            cat("No hull parent points from h2 meet the criteria \n"); flush.console()
                        } else {
    
                            ## Pick a random sample of *all* hulls for use as the NULL model of no interaction.
                            ## For the sample size, either use the fixed sample size (if n is numeric)
                            ## or the number of temporally overlapping hulls that will be examined
                            #print("lets draw a sample from h1");browser()
                            
                            if (identical(n, "all")) {
                                #h1.hidx <- (1:nrow(lhs[[h1]][["hulls"]]))[order(runif(nrow(lhs[[h1]][["hulls"]])))][1:sum(h1.hulls.to.include)]
                                h1.hidx <- which(h1.hulls.to.include)[order(runif(sum(h1.hulls.to.include)))]
                            } else {
                                #h1.hidx <- sample(1:nrow(lhs[[h1]][["hulls"]]), min(n, nrow(lhs[[h1]][["hulls"]])), replace=FALSE)                        
                                h1.hidx <- sample(which(h1.hulls.to.include), min(n, sum(h1.hulls.to.include)), replace=FALSE)                        
                            }

                            ## Pick a sample of hulls for h2
                            if (identical(n, "all")) {
                                #h2.hidx <- (1:nrow(lhs[[h2]][["hulls"]]))[order(runif(nrow(lhs[[h2]][["hulls"]])))]
                                h2.hidx <- which(h2.hulls.to.include)[order(runif(sum(h2.hulls.to.include)))]
                            } else {
                                #h2.hidx <- sample(1:nrow(lhs[[h2]][["hulls"]]), min(n, nrow(lhs[[h2]][["hulls"]])), replace=FALSE)
                                h2.hidx <- sample(which(h2.hulls.to.include), min(n, sum(h2.hulls.to.include)), replace=FALSE)                        
                            }
                            
                            # Make sure we are pairing up the same number of hulls
                            maxn <- min(length(h1.hidx), length(h2.hidx))
                            h1.hidx.use <- h1.hidx[1:maxn]
                            h2.hidx.use <- h2.hidx[1:maxn]
                            
                            ## Create a xy matrix of the centroids of h2
                            h2.ctr <- t(sapply(lhs[[h2]][["hulls"]]@polygons, function(x) apply(x@Polygons[[1]]@coords,2,mean)))

                            ## Compute the centroid distances
                            h2h.cd <- sqrt((h1.ctr[h1.hidx.use,1] - h2.ctr[h2.hidx.use,1])^2 + (h1.ctr[h1.hidx.use,2] - h2.ctr[h2.hidx.use,2])^2)
                            
                            ## Construct a name for hullset 2 that will be used as part of the name for hull metric
                            hs2.name <- paste(idVal2, ".s", lhs[[h2]][["s"]], ".", lhs[[h2]][["mode"]], lhs[[h2]][[  lhs[[h2]][["mode"]]  ]], sep="")
    
                            ## Pull out the mean centroid distance for overlapping hulls
                            to.mcd <- lhs[[h1]][["hulls"]][[paste("to.mcd.", hs2.name, sep="")]][h1.hulls.to.include]
                            if (is.null(to.mcd)) stop("Can't find the mean centroid distance for temporally overlapping hulls")
                            
                            ## Draw histogram
                            if (to.comp.hist) {
                            
                                ## Create separate histogram objects for the two distributions
                                h2h.cd.hist <- hist(h2h.cd, breaks=breaks, plot=FALSE) 
                                to.mcd.hist <- hist(to.mcd, breaks=breaks, plot=FALSE) 
                                
                                ## Make a blank plot for the combined histograms
                                xlim <- c(min(min(h2h.cd.hist[["breaks"]]), min(to.mcd.hist[["breaks"]])), max(max(h2h.cd.hist[["breaks"]]), max(to.mcd.hist[["breaks"]])))
                                ylim <- c(0, max(max(h2h.cd.hist[[hist.type]]), max(to.mcd.hist[[hist.type]])))
    
                                if (title.show) {
                                    if (is.null(title)) {
                                        if (title.id.only) {
                                            title.str <- paste(idVal1, " and ", idVal2, sep="")
                                        } else {                                    
                                            title.str <- paste(hs1.name, " and ", hs2.name, sep="")
                                        }
                                        
                                        if (title.sub.iso.enc && (!is.null(iso.lower) || !is.null(iso.upper))) {
                                            title.str <- paste(title.str, "\n within isopleth ", if (is.null(iso.lower)) "0" else iso.lower, " and ", if (is.null(iso.upper)) "1" else iso.upper, sep="")
                                        }
                                    } else {
                                        title.str <- title
                                    }
                                } else {
                                    title.str <- NULL
                                }
    
    
                                if (is.null(png.dir)) {
                                    par(mfrow = n2mfrow(figs.per.page), mar=mar, mgp=mgp)
                                    fn <- NULL
                                } else {
                                    ## Create png folder if needed
                                    if (!file.exists(png.dir)) {
                                        if (png.dir.make) {
                                            dir.made <- dir.create(png.dir)
                                            if (!dir.made) stop(paste("Couldn't make output folder ", png.dir, sep=""))
                                        } else {
                                            stop(paste(png.dir, " doesn't exist and png.dir.make is False, so can not continue."))
                                        }
                                    }
                            
                                    if (figs.per.page != 1) stop("Can only output one page per plot")
                                    iso.lower.fn <- if (is.null(iso.lower)) NULL else paste(".ilow-", iso.lower, sep="")
                                    iso.upper.fn <- if (is.null(iso.upper)) NULL else paste(".iup-", iso.upper, sep="")
                                    iso.oz.fn <- if (iso.oz) "-oz" else NULL
                                    hist.type.fn.lst <- list(density="dens", counts="freq")
                                    
                                    fn <- file.path(png.dir, paste(hs1.name, ".", hs2.name, iso.lower.fn, iso.upper.fn, iso.oz.fn, ".to-mcd.hist-", hist.type.fn.lst[[hist.type]], ".png", sep=""))
                                    if (file.exists(fn) && !png.overwrite) stop(paste(fn, "exists"))
                                    png(filename=fn, height=png.height, width=png.width, bg="white", pointsize=png.pointsize)
                                    opar <- par(mfrow = n2mfrow(figs.per.page), mar=mar, mgp=mgp)
                                    pngs.made <- c(pngs.made, fn)
                                    
                                }
    
                                plot(NULL, xlim=xlim, ylim=ylim, xlab="hull-to-hull centroid distance", ylab=hist.type, main=title.str)
                                
                                if (to.mcd.outline.only) {
                                    ## Plot the outline of the histogram of the mean centroid distance for temporally overlapping hulls
                                    to.mcd.mat <- do.call(rbind, lapply(1:length(to.mcd.hist[[hist.type]]), function(i) matrix(c(to.mcd.hist[["breaks"]][i:(i+1)], rep(to.mcd.hist[[hist.type]][i],2)), ncol=2 )))  
                                    lines(to.mcd.mat, lwd=2, col=col.to.mcd)
                                    lines(matrix(c(rep(to.mcd.hist[["breaks"]][1],2), 0, to.mcd.hist[[hist.type]][1]), ncol=2 ), lwd=2, col=col.to.mcd)
                                    lines(matrix(c(rep(to.mcd.hist[["breaks"]][length(to.mcd.hist[["breaks"]])],2), 0, to.mcd.hist[[hist.type]][length(to.mcd.hist[[hist.type]])]), ncol=2 ), lwd=2, col=col.to.mcd)
                                } else {
                                    plot(to.mcd.hist, col="gray80", add=TRUE, freq=(hist.type=="counts"))
                                }
    
                                ## Plot the outline of the histogram of the hull2hull centroid distance for randomly paired hulls
                                h2h.cd.mat <- do.call(rbind, lapply(1:length(h2h.cd.hist[[hist.type]]), function(i) matrix(c(h2h.cd.hist[["breaks"]][i:(i+1)], rep(h2h.cd.hist[[hist.type]][i],2)), ncol=2 )))  
                                lines(h2h.cd.mat, lwd=lwd.outline, col=col.h2h.cd)
                                lines(matrix(c(rep(h2h.cd.hist[["breaks"]][1],2), 0, h2h.cd.hist[[hist.type]][1]), ncol=2 ), lwd=lwd.outline, col=col.h2h.cd)
                                lines(matrix(c(rep(h2h.cd.hist[["breaks"]][length(h2h.cd.hist[["breaks"]])],2), 0, h2h.cd.hist[[hist.type]][length(h2h.cd.hist[[hist.type]])]), ncol=2 ), lwd=lwd.outline, col=col.h2h.cd)
    
                                ## Add panel.num
                                if (!is.null(panel.num)) {
                                    if (panel.num.inside.plot) {
                                        text(x=par("usr")[1], y=par("usr")[4], labels=panel.num, cex=2, adj=c(-0.3,1.2), font=2)
                                    } else {
                                        mar.old <- par("mar")
                                        par(mar=c(0, 0.3, 0.2, 0))
                                        title(main=panel.num, adj=0, cex.main=2, line=-1.5)
                                        par(mar=mar.old)
                                    }
                                }
                                if (!is.null(png.dir)) dev.off()
    
                            }
                            
                            ttest <- t.test(x=h2h.cd, y=to.mcd, alternative='two.sided', conf.level=.95, var.equal=FALSE)
                            kstest <- ks.test (x=h2h.cd, y=to.mcd)
    
                            res[[paste(hs1.name, ".", hs2.name, sep="") ]] <- list(id1=idVal1, id2=idVal2, hs1=names(lhs)[h1], hs2=names(lhs)[h2],
                                                                                   id1.idx=h1.hidx.use, id2.idx=h2.hidx.use, h2h.cd=h2h.cd, to.mcd=to.mcd,
                                                                                   ttest=ttest, kstest=kstest, iso.lower=iso.lower, iso.upper=iso.upper, fn=fn) 
                        }
                    }
                }
            }
        }
    }
    
    if (!is.null(png.dir) && status) {
        cat("png file(s) made: \n")
        print(pngs.made)
    }

    class(res) <- "locoh.tostats"
    return(invisible(res))

}

summary.locoh.tostats <- function(tostats) {

    p <- do.call(rbind, lapply(tostats, function(x) data.frame(id1=x$id1, id2=x$id2, iso.lower=max(x$iso.lower,0), iso.upper=min(x$iso.upper,1), null.model.n=length(x$id1.idx),
            to.mcd.n=length(x$to.mcd), ttest.p=x$ttest$p.value, kstest.p=x$kstest$p.value )))
    rownames(p) <- NULL
    print(p)
}




