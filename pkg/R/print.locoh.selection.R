#' Prints a locoh.selection object
#'
#' @param x An object of class locoh.selection
#' @param ... Other arguments
#'
#' @export

print.locoh.selection <- function(x, ...) {

    for (i in 1:length(x)) {
        cat(names(x)[i], ": ", paste(x[[i]], collapse=",", sep=""), "\n", sep="")
    }

}
