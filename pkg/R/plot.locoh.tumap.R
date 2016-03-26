#' Plot time use maps

#' @param tumap A locoh-time use map 
#' @param breaks The number of breaks for the symbology
#' @param nsv Plot the revisitation values (number of separate visits), T/F
#' @param mnlv Plot the proxy for visit duration (mean number of locations per visit), T/F
#' @param ivg The intervisit-gap interval to plot (in seconds)
#' @param legend Where to place the legend (set to FALSE to omit a legend)
#' @param ... Other arguments to pass to the plot function

#' @details
#' If \code{cex} is passed, it will change the text size for the plot and legend

#' @seealso \code{\link{lxy.tumap}}

#' @export

plot.locoh.tumap <- function(tumap, breaks=10, nsv=TRUE, mnlv=TRUE, ivg=NULL, legend=c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")[7], ...) {

    if (!inherits(tumap, "locoh.tumap")) stop("tumap should be of class \"locoh.tumap\"")

    ## See if the user passed a value for cex, in which case we'll use it for the legend also
    params <- list(...)
    if ("cex" %in% names(params)) {
        cex <- params$cex    
    } else {
        cex <- 1
    }

    ## Make a white to red color ramp for NSV
    if (nsv) col.nsv <- hsv(h=360/360,v=1, s=seq(0,1,length.out=breaks))
    
    ## Make a white to blue color ramp for NSV
    if (mnlv) col.mnlv <- hsv(h=240/360,v=1, s=seq(0,1,length.out=breaks))

    for (idVal in names(tumap)) {
    
        if (nsv) {
            for (ivgVal in grep("^nsv", names(tumap[[idVal]]@data), value=TRUE)) {
                ivgValNum <- substr(ivgVal, 5, nchar(ivgVal))
                nsv_classed <- cut(tumap[[idVal]][[ivgVal]], breaks=breaks)
                plot(tumap[[idVal]], col=col.nsv[as.numeric(nsv_classed)], axes=T, border="grey80", main=paste0("Number of separate visits\n", idVal, ", ivg=", ivgValNum), ...)
                if (legend!=FALSE) {
                    legend(x="topright", legend=levels(nsv_classed), fill=col.nsv, bg="white", cex=cex, title="nsv")
                }
                box()
            }
        }

        if (mnlv) {
            for (ivgVal in grep("^mnlv", names(tumap[[idVal]]@data), value=TRUE)) {
                ivgValNum <- substr(ivgVal, 5, nchar(ivgVal))
                mnlv_classed <- cut(tumap[[idVal]][[ivgVal]], breaks=breaks)
                plot(tumap[[idVal]], col=col.mnlv[as.numeric(mnlv_classed)], axes=T, border="grey80", main=paste0("Mean num locations per visits\n", idVal, ", ivg=", ivgValNum), ...)
                if (legend!=FALSE) {
                    legend(x="topright", legend=levels(mnlv_classed), fill=col.mnlv, bg="white", cex=cex, title="mnlv")
                }
                box()
            }
        }

    }
}