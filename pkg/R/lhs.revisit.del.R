#' Delete all revisitation hull metrics in a LoCoH-hullset object
#'
#' @param lhs A LoCoH-hullset object
#' @param status Display summary, T/F
#'
#' @return A LoCoH-hullset object
#'
#' @seealso \code{\link{lhs.revisit.add}}
#' @export

lhs.revisit.del <- function(lhs, status=TRUE) {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")

    hm.deleted <- NULL
    for (hs.idx in 1:length(lhs)) {
        if (status) cat(names(lhs)[hs.idx], "\n")
        
        hm.nsr <- sapply(lhs[[hs.idx]][["hm"]], function(x) x$type == "nsr")
        hm.nsr.colnames <- names(lhs[[hs.idx]][["hm"]])[hm.nsr]
        
        ## Remove from lhs$hm
        lhs[[hs.idx]][["hm"]][hm.nsr] <- NULL
        
        ## Remove from lhs$hulls
        data.col.idx.del <- names(lhs[[hs.idx]][["hulls"]]@data) %in% hm.nsr.colnames
        
        lhs[[hs.idx]][["hulls"]]@data <- lhs[[hs.idx]][["hulls"]]@data[ , !data.col.idx.del]
        
        ## Delete from hm.params
        
        lhs[[hs.idx]][["hm.params"]][["ta.min"]] <- NULL
        lhs[[hs.idx]][["hm.params"]][["ta.max"]] <- NULL
        
        hm.deleted <- c(hm.deleted, paste(names(lhs)[hs.idx], ": ", hm.nsr.colnames, sep=""))
        
    }    
    
    
    if (status) {
        if (length(hm.deleted)>0) {
            cat("Revisit hull metrics deleted: \n") 
            print(hm.deleted)
        } else {
            cat("No reivist hull metrics founds \n") 
        }
    }
    
    return(lhs)

}
