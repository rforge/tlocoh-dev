
lxy.lhs.batch <- function(lxy, a=NULL, k=NULL, r=NULL, s, dir=".", suf=NULL, ud=T, iso.levels=c(0.1,0.25,0.5,0.75,0.95), save.hulls=TRUE, save.enc.pts=TRUE) {

    ## Creates separate locoh hullsets for each parameter value and saves them to disk.
    ## For use with large datasets where memory limits make it difficult to create a hullset collection with multiple hullsets
    
    res <- NULL

    s <- vectorize.parameter(s, n2z=TRUE)

    for (sVal in s) {
    for (kVal in vectorize.parameter(k, n2z=TRUE)) {
    for (aVal in vectorize.parameter(a, n2z=TRUE)) {
    for (rVal in vectorize.parameter(r, n2z=TRUE)) {
    
        if (identical(kVal,0)) kVal <- NULL
        if (identical(aVal,0)) aVal <- NULL
        if (identical(rVal,0)) rVal <- NULL
    
        lhs <- lxy.lhs(lxy, k=kVal, a=aVal, r=rVal, s=sVal, ud=ud, iso.levels=iso.levels, save.hulls=save.hulls, save.enc.pts=save.enc.pts)
        
        #cat("Size of lhs: ", as.numeric(object.size(lhs)) / 2^20, " Mb \n", sep="")
        #if (ud) lhs <- lhs.iso.add(lhs, iso.levels=iso.levels)
        fn <- lhs.save(lhs, dir=dir, suf=NULL)
        rm(lhs)
        res <- c(res, fn)
        cat("\n");flush.console()
    }
    }
    }
    }

    cat("The following files were saved: \n")
    print(res)
    return(invisible(res))

}
