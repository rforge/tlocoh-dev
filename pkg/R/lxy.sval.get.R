#' Get s values from a LoCoH-xy object
#'
#' Extracts the s-values from a LoCoH-xy object for a specified proportion of time-selected hulls
#'
#' @param lxy A \link[tlocoh]{LoCoH-xy} object
#' @param The desired proportion of time-select hulls (0..1)
#' @param id The name(s) of individuals to extract values from
#'
#' @return
#' A named list of \emph{s} values, with one element per id. 
#'
#' @details This function extracts the \emph{s} values from the individuals in a LoCoH-xy object.
#' This can be useful if you want to select \emph{s} values based on a consistent proportion
#' of time-selected hulls. 
#'
#' \emph{s} values are computed by \code{\link[tlocoh]{lxy.ptsh.add}} by iteratively
#' trying differet values of \emph{s} for a desired proportion of time selected hulls.
#' This function returns the corresponding value of \emph{s} that generates the desired 
#' ptsh within a given threshhold (see \code{\link[tlocoh]{lxy.ptsh.add}}). If no matching values
#' of \emph{s} were found, an empty vector will be returned for that list element. If
#' more than one set of ptsh tables are found (i.e, lxy.ptsh.add was run more than once),
#' only the first set of 's' values will be returned.
#'
#' @seealso \code{\link[tlocoh]{lxy.ptsh.add}}; Vignette on locoh.lxy data class
#'
#' @export


lxy.sval.get <- function(lxy, ptsh, id=NULL) {

    if (!inherits(lxy, "locoh.lxy")) stop("lxy should be of class \"locoh.lxy\"")
    if (length(ptsh)>1) stop("ptsh must be of length 1")
    if (is.null(lxy[["ptsh"]])) stop("ptsh table not found. Run lxy.ptsh.add and try again")
    
    id <- tlocoh::vectorize.parameter(id, type="character", sort=FALSE)
    if (is.null(id)) {
        id <- levels(lxy[["pts"]][["id"]])
    } else {
        if (FALSE %in% (id %in% levels(lxy[["pts"]][["id"]]))) stop("id value(s) not found")
    }

    if (max(sapply(id, function(x) length(lxy$ptsh[[x]]))) > 1) warning(cw("More than one ptsh table found. Only s-values from the first one will be returned"))
    
    svals <- lapply(id, function(x) lxy$ptsh[[x]][[1]][["target.s"]][  lxy$ptsh[[x]][[1]][["target.ptsh"]] == ptsh   ] )
    names(svals) <- id
    
    return(svals)


}