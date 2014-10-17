#' Add revisitation hull metrics to a LoCoH-hullset object
#'
#' Computes revisitation rate based on a minimum and maximum time away period
#'
#' @param lhs A LoCoH-hullset object
#' @param ta.min Minimum value(s) for time away in seconds (numeric vector)
#' @param ta.max Maximum value(s) for time away in seconds (numeric vector)
#' @param ta.cuts A numeric vector of time values in seconds that define the time-away intervals
#' @param status Show status messages. T/F
#'
#' @details
#' \code{ta.min} and \code{ta.max} define the minimum and maximum period of time (in seconds) which must pass for another occurence in the hull to 
#' be considered a 'revisit'. They should be the same lengths. If \code{ta.max} is NULL, no upper bound will be set
#'
#' \code{ta.cuts} is an alternative way of specifying the time-away intervals. The time values in \code{ta.cuts} will serve as the values between 
#' time away intervals. For example if \code{ta.cuts = c(1000,4000,7000,10000)}, three time-away intervals will be examined: 1000 to 4000 
#' seconds, 4000 to 7000 seconds, and 7000 to 10000 seconds. One way to get the values for \code{ta.cuts} is to plot the distribution
#' of revisit times using \code{\link{lhs.plot.revisit}}, and then use the \code{\link{get.vals}} function to select time values between clusters of revisit times.
#'
#' @return A LoCoH-hullset object
#'
#' @seealso \code{\link{lhs.plot.revisit}}, \code{\link{get.vals}}
#' @export


lhs.revisit.add <- function(lhs, ta.min=NULL, ta.max=NULL, ta.cuts=NULL, status=TRUE) {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (TRUE %in% sapply(lhs, function(hs) is.null(hs[["pts"]][["dt"]]))) stop("Date stamps not found, can't compute revisitation metrics")
    if (!is.null(ta.max)) {
        if (length(ta.max) != length(ta.min)) stop("ta.max should be the same length as ta.min")
        if (min(ta.max - ta.min) <= 0) stop("ta.max should be greater than ta.min") 
    }
    if (is.null(ta.cuts)) {
        if (is.null(ta.min)) stop("ta.min is a required parameter")        
    } else {
        if (!is.null(ta.min) || !is.null(ta.max)) stop("If ta.cuts is passed, do not pass ta.min and ta.max")
        if (length(ta.cuts) < 2) stop("ta.cuts must have at least two values")
        ta.cuts <- sort(ta.cuts)        
        ta.min <- ta.cuts[1:(length(ta.cuts)-1)]
        ta.max <- ta.cuts[2:length(ta.cuts)]
    }

    for (hs.idx in 1:length(lhs)) {
        if (status) cat(names(lhs)[hs.idx], "\n")
        
        dt.int <- as.numeric(lhs[[hs.idx]][["pts"]][["dt"]])
        tau <- lhs[[hs.idx]][["rw.params"]][["time.step.median"]]        
        
        for (ta.idx in 1:length(ta.min)) {
            taminVal <- ta.min[ta.idx]
            tamaxVal <- if (is.null(ta.max)) 0 else ta.max[ta.idx]
            
            if (status) cat(cw(paste("- ", ta.idx, " of ", length(ta.min), ". Computing revisits for ta.min=", taminVal, " (", secs.fmt(taminVal), ")",
                        if (tamaxVal==0) "" else paste(" to ta.max=", tamaxVal, " (", secs.fmt(tamaxVal), ")",  sep=""), ".", sep=""), indent=1, exdent=3, final.cr=F))
            
            if (tamaxVal==0) {
                nsr <- sapply(lhs[[hs.idx]][["enc.pts"]][["idx"]], function(x) sum(diff(dt.int[x]) >= taminVal))
            } else {
                nsr <- sapply(lhs[[hs.idx]][["enc.pts"]][["idx"]], function(x) {epdt <- diff(dt.int[x]); sum(epdt >= taminVal & epdt < tamaxVal)})
            }
            
            revisit.str <- paste("nsr.", taminVal, ".", tamaxVal, sep="")
            lhs[[hs.idx]][["hulls"]][[revisit.str]] <- nsr
            lhs[[hs.idx]][["hm"]][[revisit.str]] <- list(type="nsr", aux=list(ta.min=taminVal, ta.max=tamaxVal))
            
            if (status) cat(" Done. \n");flush.console()
                            
        }
        
        ## Add values of hm.params
        if (is.null(lhs[[hs.idx]][["hm.params"]])) lhs[[hs.idx]][["hm.params"]] <- list()
        lhs[[hs.idx]][["hm.params"]][["ta.min"]] <- unique(c(lhs[[hs.idx]][["hm.params"]][["ta.min"]], ta.min))
        lhs[[hs.idx]][["hm.params"]][["ta.max"]] <- unique(c(lhs[[hs.idx]][["hm.params"]][["ta.max"]], n2z(ta.max)))
        
    }    
    if (status) cat("Done.\n")
    return(lhs)

}
