#' Plot time use maps
#' @export

plot.locoh.tumap <- function(tumap, breaks=10, nsv=TRUE, mnlv=TRUE, ivg=NULL) {

    if (!inherits(tumap, "locoh.tumap")) stop("tumap should be of class \"locoh.tumap\"")

    ## Make a white to red color ramp for NSV
    if (nsv) col.nsv <- hsv(h=360/360,v=1, s=seq(0,1,length.out=breaks))
    
    ## Make a white to blue color ramp for NSV
    if (mnlv) col.nsv <- hsv(h=240/360,v=1, s=seq(0,1,length.out=breaks))

    for (idVal in names(tumap)) {
    
        if (nsv) {
            for (ivgVal in grep("^nsv", names(tumap[[idVal]]@data), value=TRUE)) {
                ivgValNum <- substr(ivgVal, 5, nchar(ivgVal))
                nsv_classed <- cut(tumap[[idVal]][[ivgVal]], breaks=breaks)
                plot(tumap[[idVal]], col=col[as.numeric(nsv_classed)], axes=T, border="grey50", main=paste0("Number of separate visits\n", idVal, ", ivg=", ivgValNum))
            }
        }

        if (mnlv) {
            for (ivgVal in grep("^mnlv", names(tumap[[idVal]]@data), value=TRUE)) {
                ivgValNum <- substr(ivgVal, 5, nchar(ivgVal))
                mnlv_classed <- cut(tumap[[idVal]][[ivgVal]], breaks=breaks)
                plot(tumap[[idVal]], col=col[as.numeric(mnlv_classed)], axes=T, border="grey50", main=paste0("Mean num locations per visits\n", idVal, ", ivg=", ivgValNum))
            }
        }

    }
}