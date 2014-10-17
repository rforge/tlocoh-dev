#' Add hull metrics for association analysis
#'
#' @param lhs A \link{LoCoH-hullset} object
#' @param save.hso Whether to save the hull intersection list, T/F
#' @param skip.dups Skip duplicate hulls (faster)
#'
#' @param id A character vector of the hullset ids to compute metrics for. Can also be \code{'all'}.
#' @param hs2.id A character vector of the hullset ids to use as the comparison hullsets. Can also be \code{'all'}.
#' @param tbuff A temporal overlap threshhold (in seconds). See details.
#' @param ivg The intervisit gap period used to collapse intersecting hulls into discrete visits, see details
#' @param test A two-element numeric vector containing the number of hulls in hullset 1 and hullset 2 respectively to identify intersections
#' @param status Show status messages. T/F
#' @param piFUN The function to use to identify which pairs of hulls intersect: 'pIntersect' or 'pIntersectSat'
#'
#' @return A \link{LoCoH-hullset} object
#'
#' @details 
#' This function computes hull metrics for the spatially overlapping hulls from two ids. Typically this would be used
#' when you have hullset from two individuals (i.e., two animals) and you want to see the spatial and temporal patterns of shared space use.
#'
#' You can impose a temporal overlap requirement as well by passing a value for \code{tbuff}. Two hulls will be considered spatially overlapping only if their parent
#' points also were recorded within \code{tbuff} seconds of each other. This essentially produces metrics for spatially and temporally overlapping hulls.
#'
#' Hullset metrics are computed for each pair of ids. Thus if a hullset has hulls for three unique ids, 
#' each hull will have spatial overlap metrics computed for each of the other two hullsets. You can narrow which
#' id(s) to compute metrics for, and which hullset(s) to use as the comparison, with 
#' the \code{id} and \code{hs2.id} arguments.
#' 
#' Up to three spatial overlap metrics are computed. \code{so.count} is simply the number of hulls in hullset 2 that overlap. 
#' \code{so.dtmin} is the minimum amount of time (expressed in seconds) that passes between overlapping hulls. This reflects temporal partitioning 
#' of shared space - low values of \code{so.dtmin} suggest the two individuals don't mind being in the same area at the same time.
#' \code{so.nsv} (number of separate visits) is similiar to \code{so.count}, but collapses overlapping hulls 
#' into discrete visits based on an intervisit gap period \code{ivg}. \code{so.nsv}
#' is only computed if a value for \code{ivg} is passed.
#'
#' \code{pIntersect} and \code{pIntersectSat} are two functions that identify which pairs of hulls actually intersect. Neither are 
#' terribly fast, but \code{pIntersect} appears to work faster than \code{pIntersectSat}.
#'
#' @export
#' @import pbapply sp

lhs.so.add <- function(lhs, id="all", hs2.id="all", tbuff=0, ivg=NULL, test=0, skip.dups=TRUE, save.hso=TRUE,
                       status=TRUE, piFUN=c("pIntersect", "pIntersectSat")[1]) {


    ## In my tests, pIntersect is a lot faster than pIntersectSat

    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    
    start.time <- Sys.time()
    if (status) cat(" - start time: ", as.character(start.time), "\n", sep = "")
    
    if (!exists(piFUN)) stop(paste("Function ", piFUN, " not found", sep=""))
    pIntersectsFUN <- get(piFUN)

    ## Make a data frame of all of the hullsets and their properties
    lhs.info <- do.call(rbind, lapply(lhs, function(x) data.frame(id=x$id, s=x$s, k=n2z(x$k), a=n2z(x$a), r=n2z(x$r), stringsAsFactors=FALSE)))
    id.all <- unique(lhs.info[["id"]])
    if (length(id.all)==1) stop("This hull metric requires a hullset collection with more than one ID")

    if (identical(id, "all")) id <- id.all
    if (identical(hs2.id, "all")) hs2.id <- id.all
    if (length(test)==1) test <- rep(test,2)
    if (length(test) !=2 ) stop("Test must be a one or two-element numeric vector")

    ## Performance testing code
    #assign("pInter.log", NULL, envir=.GlobalEnv)
    ## fn.log="segment.intersections.txt", 
    #if (length(fn.log) > 0) assign("con.seg.test", file(fn.log, open = "a"), envir=.GlobalEnv)
    
    
    for (idVal in id) {
        ## Save the indices of the hullset(s) which have this id
        hsi.idVal1 <- which(lhs.info[["id"]] == idVal)
        
        for (h1 in hsi.idVal1) {
            ## Do I need to keep hulls1 as its own object?
            ## Is there a performance gain by coercing hulls to a SpatialPolygons object?
            ## Is there a performance gain if hullsx.bb is set by hullsx.coord.lst instead of @coords
            
            if (status) cat(names(lhs)[h1], "\n")
            
            if (save.hso) {if (is.null(lhs[[h1]][["hso"]])) lhs[[h1]][["hso"]] <- list()}

            hulls1.sp <- lhs[[h1]][["hulls"]]
            if (test[1] > 0) hulls1.sp <- hulls1.sp[1:min(length(hulls1.sp),test[1]), ]

            h1.pp.dt.int <- as.integer(lhs[[h1]][["pts"]][["dt"]][lhs[[h1]][["hulls"]][["pts.idx"]]])

            ## Extract the area and coordinates for all hulls. We will use this to identify duplicates and also passs to pIntersect
            ha1 <- sapply(hulls1.sp@polygons, function(p) p@area)
            hulls1.coord.lst <- lapply(hulls1.sp@polygons, function(x) x@Polygons[[1]]@coords)

            ## By default we will analyze all hulls
            h1.analyze.idx <- 1:length(hulls1.sp)

            ## To save time we can skip over duplicates and copy them back in later
            if (skip.dups) {

                ## Use hull area as a quick way to identify possible duplicates
                ## Make a list of the duplicate area hulls, the list element name will be the index of the first or 'master' polygon
                ha1.dups <- rep(NA, length(ha1))
                ha1.dups[duplicated(ha1)] <- match(ha1[duplicated(ha1)], ha1)
                ha1.dups.lst <- split(1:length(ha1), ha1.dups)
                dups1.master.idx <- as.numeric(names(ha1.dups.lst))

                ## Verify they are actually duplicates by comparing the coordinates
                ha1.dups.lst.verified <- lapply(1:length(ha1.dups.lst), function(i) ha1.dups.lst[[i]] [sapply(ha1.dups.lst[[i]], function(j) identical(hulls1.coord.lst[[ dups1.master.idx[i] ]], hulls1.coord.lst[[j]]))])
               
                h1.analyze.idx <- h1.analyze.idx[-unlist(ha1.dups.lst.verified, use.names=F)]
                cat(" - ", length(unlist(ha1.dups.lst.verified, use.names=F)), " duplicate hulls detected \n", sep="")
            }
            
            ## Extract the bounding boxes which we will use to create a short list
            hulls1.bb <- do.call(rbind, lapply(hulls1.sp@polygons[h1.analyze.idx], function(x) as.numeric(apply(x@Polygons[[1]]@coords,2,range))))

            for (idVal2 in hs2.id) {
                if (idVal2 != idVal) {
                    hsi.idVal2 <- which(lhs.info[["id"]] == idVal2)
                    
                    for (h2 in hsi.idVal2) {                    
                        #cat("h1=", h1, ". h2=", h2, "\n", sep="")
                        if (status) cat(" - finding intersections with ", names(lhs)[h2], "\n", sep="")

                        hulls2.sp <- lhs[[h2]][["hulls"]]
                        if (test[2]>0) hulls2.sp <- hulls2.sp[1:min(length(hulls2.sp),test[2]),]
                        
                        h2.tau <- lhs[[h2]][["rw.params"]][["time.step.median"]][1]
                        if (tbuff > 0 && tbuff < h2.tau) stop("tbuff should be larger than the median sampling interval")
                        
                        h2.pp.dt.int <- as.integer(lhs[[h2]][["pts"]][["dt"]][lhs[[h2]][["hulls"]][["pts.idx"]]])

                        #if (status) cat("   - preparing coordinate list, bounding boxes, and area..."); flush.console()
                        hulls2.coord.lst <- lapply(hulls2.sp@polygons, function(x) x@Polygons[[1]]@coords)
                        ha2 <- sapply(hulls2.sp@polygons, function(x) x@area)
                        #if (status) cat("Done.\n")

                        h2.analyze.idx <- 1:length(hulls2.sp)
                        
                        if (skip.dups) {

                            ## Extract hull area as a quick way to identify possible duplicates
                            ## ha <- sapply(hulls1.sp@polygons, function(p) p@Polygons[[1]]@area)

                            ## Make a list of the duplicate hulls, the list element name will be the index of the first or 'master' polygon
                            ha2.dups <- rep(NA, length(ha2))
                            ha2.dups[duplicated(ha2)] <- match(ha2[duplicated(ha2)], ha2)

                            ha2.dups.lst <- split(1:length(ha2), ha2.dups)
                            dups2.master.idx <- as.numeric(names(ha2.dups.lst))

                            #ha.dups.lst.verified <- lapply(1:length(ha.dups.lst), function(i) ha.dups.lst[[i]] [sapply(ha.dups.lst[[i]], function(j) identical(hulls1.sp@polygons[[dups.master.idx[i]]]@Polygons[[1]]@coords, hulls1.sp@polygons[[j]]@Polygons[[1]]@coords))])

                            ha2.dups.lst.verified <- lapply(1:length(ha2.dups.lst), function(i) ha2.dups.lst[[i]] [sapply(ha2.dups.lst[[i]], function(j) identical(hulls2.coord.lst[[ dups2.master.idx[i] ]], hulls2.coord.lst[[j]]))])

                            h2.analyze.idx <- h2.analyze.idx[-unlist(ha2.dups.lst.verified, use.names=F)]
                            cat("   - ", length(unlist(ha2.dups.lst.verified, use.names=F)), " duplicate hulls detected \n", sep="")
                        }

                        ## Compute the bounding boxes of h2 hulls
                        hulls2.bb <- do.call(rbind, lapply(hulls2.sp@polygons[h2.analyze.idx], function(x) as.numeric(apply(x@Polygons[[1]]@coords,2,range))))

                        ## Create a short list of h2 hulls that are *not* disjoint with the bounding boxes of hulls1
                        if (status) cat("   - creating short list based on bounding boxes..."); flush.console()
                        bb2sl <- lapply(1:length(h1.analyze.idx), function(i) which(!(hulls2.bb[,2] < hulls1.bb[i,1] | hulls2.bb[,1] > hulls1.bb[i,2] | hulls2.bb[,3] > hulls1.bb[i,4] | hulls2.bb[,4] < hulls1.bb[i,3])))
                        if (status) cat("Done.\n")

                        if (status) cat("   - identifying intersections \n")
                        #print("Just have to try");browser()
                        h1an.h2an.lst <- pblapply(1:length(h1.analyze.idx), function(i) bb2sl[[i]][unlist(sapply(bb2sl[[i]], function(j) pIntersectsFUN(p1=hulls1.coord.lst[[ h1.analyze.idx[i] ]], p2=hulls2.coord.lst[[ h2.analyze.idx[j] ]], p1.area=ha1[ h1.analyze.idx[i] ], p2.area=ha2[ h2.analyze.idx[j] ])))])
                        
                        if (skip.dups) {
                            ## Convert the hull2 indices from a reference to h2.analyze.idx to absolute
                            h1an.h2ab.lst <- lapply(h1an.h2an.lst, function(x) h2.analyze.idx[x])

                            ## Next, for each list element that contains a 'master' h2 index, append the indices of all duplicate h2 hulls
                            for (i in 1:length(dups2.master.idx)) {
                                have.master.idx <- which(sapply(h1an.h2ab.lst, function(x) dups2.master.idx[i] %in% x))
                                if (length(have.master.idx) > 0) {
                                    dup.idx.df <- data.frame(matrix(data=ha2.dups.lst.verified[[i]], ncol=length(have.master.idx), nrow=length(ha2.dups.lst.verified[[i]])))
                                    h1an.h2ab.lst[have.master.idx] <- mapply(c, h1an.h2ab.lst[have.master.idx], dup.idx.df, SIMPLIFY=FALSE)
                                }
                            }

                            
                            ## Save these to a list for all of hulls1
                            hulls1.hulls2.lst <- vector("list", length(hulls1.sp))
                            hulls1.hulls2.lst[h1.analyze.idx] <- h1an.h2ab.lst

                            ## Lastly, we copy the h1 elements that are master copies of duplicates to the duplicates
                            for (i in 1:length(dups1.master.idx)) {
                                for (j in ha1.dups.lst.verified[[i]]) {
                                    hulls1.hulls2.lst[[j]] <- hulls1.hulls2.lst[[dups1.master.idx[i]]]
                                }
                            }

                        } else {
                            ## We didn't skip duplicates
                            hulls1.hulls2.lst <- h1an.h2an.lst
                        }

                        if (length(hulls1.hulls2.lst) != length(hulls1.sp)) stop(cw("For some reason the length of the hull intersections list is different than the length of hulls analyzed", final.cr=F))

                        ## Sort the hulls by index number (which is equivalent to sorting by date-time)
                        hulls1.hulls2.lst <- lapply(hulls1.hulls2.lst, sort)

                        ## Construct a name for hullset 2 that will be used to save the spatial overlap hull metrics
                        hs2.name <- paste(idVal2, ".s", lhs[[h2]][["s"]], ".", lhs[[h2]][["mode"]], lhs[[h2]][[  lhs[[h2]][["mode"]]  ]], sep="")

                        ## Create a different list if there is a tbuff
                        if (tbuff==0) {
                            hin.lst.str <- "hulls1.hulls2.lst"
                        } else {
                            hulls1.hulls2.within.tbuff.lst <- lapply(1:length(hulls1.hulls2.lst), function(h1.idx) hulls1.hulls2.lst[[h1.idx]][sapply(hulls1.hulls2.lst[[i]], function(h2.idx) abs(h2.pp.dt.int[h2.idx] - h1.pp.dt.int[h1.idx]) <= tbuff)])
                            hin.lst.str <- "hulls1.hulls2.within.tbuff.lst"
                        }

                        if (max(test)==0) {

                            ## Count the number of hulls within tbuff that intersect
                            so.count <- sapply(get(hin.lst.str), length)
    
                            ## Add an item to hm containing the meta data for so.count
                            so.count.name <- paste("so.count.", hs2.name, sep="")
                            lhs[[h1]][["hm"]][[so.count.name]] <- list(type="so.count", aux=list(hs2.id=idVal2, hs2.s=lhs[[h2]][["s"]], hs2.k=lhs[[h2]][["k"]], hs2.a=lhs[[h2]][["a"]], hs2.r=lhs[[h2]][["r"]], tbuff=tbuff))

                            ## Add the hull metric to the SpatialPolygonsDataFrame
                            lhs[[h1]][["hulls"]]@data[[so.count.name]] <- so.count
    
                            ## Find the minimum time difference between the parent points of the intersecting hulls (dtmin)
                            so.dtmin <- sapply(1:length(hulls1.hulls2.lst), function(i) if (length(get(hin.lst.str)[[i]])==0) NA else min(abs(h2.pp.dt.int[get(hin.lst.str)[[i]]] - h1.pp.dt.int[i])))
    
                            ## Add an item to hm containing the meta data for so.dtmin
                            so.dtmin.name <- paste("so.dtmin.", hs2.name, sep="")
                            lhs[[h1]][["hm"]][[so.dtmin.name]] <- list(type="so.dtmin", aux=list(hs2.id=idVal2, hs2.s=lhs[[h2]][["s"]], hs2.k=lhs[[h2]][["k"]], hs2.a=lhs[[h2]][["a"]], hs2.r=lhs[[h2]][["r"]], tbuff=tbuff))
                            
                            ## Add the hull metric to the SpatialPolygonsDataFrame
                            lhs[[h1]][["hulls"]]@data[[so.dtmin.name]] <- so.dtmin
                            
                            ## Compute so.ivg if ivg is not null
                            if (!is.null(ivg)) {
                                for (ivgVal in ivg) {
                                    ivg.tab.lst <- lapply(get(hin.lst.str), function(x) as.numeric(table(cumsum(c(1, diff(h2.pp.dt.int[x]) >= ivgVal)))))
                                    so.nsv <- sapply(ivg.tab.lst, length)
                                    
                                    ## Add an item to hm containing the meta data for so.ivg
                                    so.ivg.name <- paste("so.ivg.", hs2.name, ".ivg", ivgVal, sep="")
                                    lhs[[h1]][["hm"]][[so.ivg.name]] <- list(type="so.ivg", aux=list(hs2.id=idVal2, hs2.s=lhs[[h2]][["s"]], hs2.k=lhs[[h2]][["k"]], hs2.a=lhs[[h2]][["a"]], hs2.r=lhs[[h2]][["r"]], tbuff=tbuff, ivg=ivgVal))

                                    ## Add the hull metric to the SpatialPolygonsDataFrame
                                    lhs[[h1]][["hulls"]]@data[[so.ivg.name]] <- so.nsv
                                }
                            }
                            

                        }

                        ## Save the hull intersection list
                        if (save.hso) {
                            lhs[[h1]][["hso"]][[names(lhs)[h2]]] <- list(id=idVal2, s=lhs[[h2]][["s"]], mode=lhs[[h2]][["mode"]], k=n2z(lhs[[h2]][["k"]]), a=n2z(lhs[[h2]][["a"]]), r=n2z(lhs[[h2]][["r"]]), hidx=hulls1.hulls2.lst)
                        }

                    }
                    
                    
                }
            }
        }
    }
    
    if (status) {
        time.taken = difftime(Sys.time(), start.time, units="auto")
        cat("Total time:", round(time.taken,1), units(time.taken), "\n", sep = " ")
    }

    return(lhs)

}
