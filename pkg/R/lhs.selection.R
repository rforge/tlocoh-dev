#' Create a hullset selection object
#'
#' @param lhs A LoCoH-hullset object
#' @param k Optional initial value for \emph{k}
#' @param a Optional initial value for \emph{a}
#' @param r Optional initial value for \emph{r}
#'
#' @details This will create a new hullset selection object of class locoh.selection. You can then feed this
#' object into code{\link{lhs.shiny.select}} to visually select a hullset for each individual. This is typically
#' done when you have create a hullset containing multiple hullsets for each individual, and you want to pick one
#' per individual for the rest of the analysis.
#'
#' @seealso code{\link{lhs.shiny.select}}
#'
#' @export

lhs.selection <- function(lhs, k=NULL, a=NULL, r=NULL) {

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    lhs.mode <- unique(sapply(lhs, function(x) x$mode))
    if (length(lhs.mode) > 1) stop("More than one mode in this hullset")
    lhs.id <- unique(sapply(lhs, function(x) x$id))
    lhs.akr_vals <- sort(unique(sapply(lhs, function(x) x[[x$mode]])))

    if ((as.numeric(!is.null(k)) + as.numeric(!is.null(a)) + as.numeric(!is.null(r))) > 1) stop("Only pass one argument k, a, or r")
    if (!is.null(k)) {
        if (lhs.mode != "k") stop("This hullset does not use the k-method")
        if (!k %in% lhs.akr_vals) stop("k value not found")
        default_val <- k

    } else if (!is.null(a)) {
        if (lhs.mode != "a") stop("This hullset does not use the a-method")
        if (!a %in% lhs.akr_vals) stop("k value not found")
        default_val <- a

    } else if (!is.null(r)) {
        if (lhs.mode != "r") stop("This hullset does not use the r-method")
        if (!r %in% lhs.akr_vals) stop("k value not found")
        default_val <- r

    } else {
        default_val <- "NA"
    }

    hs.selection <- as.list(rep(as.character(default_val), length(lhs.id)))
    names(hs.selection) <- lhs.id
    attr(hs.selection, "mode") <- lhs.mode
    class(hs.selection) <- c("locoh.selection", "list")
    return(invisible(hs.selection))
}
