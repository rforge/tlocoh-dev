
lhs.exp.isodata <- function(fn, plot.iso=FALSE) {

    ## Given a series of *.RData file names containing LoCoH-hullset objects, this
    ## will extract the isopleth attributes and return a dataframe
    ## This function is used in combination with lxy.lhs.batch2disk for large datasets that
    ## can't fit into memory

    #if (!plotme %in% c("none","isoarea","isoear")) stop("Unknown value for plotme")

    isoc.info <- NULL

    for (lhs.fn in fn) {
        cat(" - loading ", lhs.fn, "...", sep=""); flush.console()
        obj <- load(lhs.fn)
        cat("Done.\n"); flush.console()
        lhs <- get(obj)
        for (hs.idx in 1:length(lhs)) {
            kar <- lhs[[hs.idx]][["mode"]]
            id <- lhs[[hs.idx]][["id"]]
            for (iso.idx in 1:length(lhs[[hs.idx]][["isos"]])) {
                iso.info <- data.frame(id=id, mode=kar, kar=lhs[[hs.idx]][[kar]], s=lhs[[hs.idx]][["s"]], iso.method=lhs[[hs.idx]][["isos"]][[iso.idx]][["iso.method"]], sort.metric=lhs[[hs.idx]][["isos"]][[iso.idx]][["sort.metric"]], lhs[[hs.idx]][["isos"]][[iso.idx]][["polys"]]@data)
                #print("Lets rbind");browser()
                isoc.info <- rbind(isoc.info, iso.info)
            }
            
        }

        if (plot.iso) plot(lhs, iso=T, png.dir="tmp", png.width=600, desc=0, shp.csv="gis_layers.csv", layers="t_border", axes.titles=F, axes.ticks=F, mar=c(0.5, 0.5, 2.8, 0.5))
        ## Change the first column name to the mode
        ## names(isoc.info)[1] <- kar

        rm(list=obj)
    }
    
    
        
    return(invisible(isoc.info))

}
