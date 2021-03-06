% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lhs.shiny.select.R
\name{lhs.shiny.select}
\alias{lhs.shiny.select}
\title{Interactive selection of hullsets from a locoh-hullset}
\usage{
lhs.shiny.select(lhs, selection = NULL, gmap = "none")
}
\arguments{
\item{lhs}{A LoCoH-hullset object}

\item{selection}{An object of class \emph{locoh.selection} containing a named list (one per individual) of selected hullset parameter values}

\item{gmap}{The name of a background image that will be downloaded from Google: \code{"none"}, 
\code{"roadmap"}, \code{"satellite"}, \code{"hybrid"}, or \code{"terrain"}. May also be a object of type \code{locoh.gmap}, see Notes.}
}
\value{
An object of class \emph{locoh.selection} containing a list of selected hullset parameter value for each id (see Details)

with plots of isopleth area and isopleth edge:area ratio, and select the value of 'k', 'a', or 'r' that does the
best job balancing over- and under-estimation.

This function can be useful when you have a locoh-hullset object with multiple hullsets for each individual over a range of \emph{k} or \emph{a} values, and
you want to pick one hullset per individual for the rest of the analysis. Note the range of parameter values should be uniform. 

This function requires using RStudio and the Shiny package. To return the selection, be sure to click the 'Save and Return' button in the shiny app. The object returned can be passed to the function again as the 'selection' argument.

To display an image from Google in the background, set gmap to \code{"roadmap"}, \code{"satellite"}, \code{"hybrid"}, or \code{"terrain"}. 
This requires an internet connection. You may also set gmap to an object of type \code{locoh.gmap}, so the image(s) don't have to be 
downloaded each time. See \code{lhs.gmap} (in the tlocoh.dev package).
}
\description{
Visually select one hullset per individual from a \link[tlocoh]{LoCoH-hullset}
}
\seealso{
code{\link{lhs.selection}}, code{\link{lhs.gmap}}
}
