#' Compute isopleth overlap
#'
#' Computes the area of isopleth overlap
#'
#' @param lhs A LoCoH-hullset object
#' @param id The id(s) of the hullsets to compare
#' @param k The k value of hullsets to compare
#' @param r The r value of hullsets to compare
#' @param a The a value of hullsets to compare
#' @param s The s value of hullsets to compare
#' @param hs.names The name(s) of saved hullsets to compare
#' @param iso.level A numeric vector of the isopleth level(s) of interest
#' @param hsnames_simplify If True, will simplify the hullset names to just the IDs for the row and column names of the overlap matrices (see Return)
#' @param status Show messages. T/F
#'
#' @details
#' This function computes the area intersection of isopleths among different hullsets. 
#' This might be done, for example, if the hullsets belong to different individuals, and you 
#' want to see which individuals share space. All pairs of hullsets in \code{lhs} will be compared,
#' and all isopleth levels will be compared. 
#' 
#' @return 
#' A list object with three named elements. The \emph{spdf} contains a SpatialPolygonsDataFrame,
#' with a data table that saves the areas of intersection expressed in map units and proportions of the isopleth area 
#' for each hullset. \emph{overlap_area} contains a square matrix whose values are the area of intersection in map units.
#' \emph{overlap_area} contains a square matrix whose values are the areas of intersection expressed as proportions
#' of the each the two isopleths.
#'
#' @seealso \code{\link[tlocoh]{isopleths}}
#'
#' @export

lhs.iso.overlap <- function(lhs, id=NULL, k=NULL, r=NULL, a=NULL, s=NULL, hs.names = NULL, 
                        iso.level=NULL, hsnames_simplify=TRUE, status=TRUE) {
     
    if (!inherits(lhs, "locoh.lhs")) stop("lhs should be of class \"locoh.lhs\"")
    if (!requireNamespace("rgeos")) stop("package rgeos required")
    
    if (is.null(id) && is.null(r) && is.null(k) && is.null(a) && is.null(s) && is.null(hs.names)) {
        hs <- lhs
    } else {    
        hs <- lhs.select(lhs, id=id, r=r, k=k, a=a, s=s, hs.names=hs.names)
    }
    
    if (length(hs)==0) stop("No hullsets found matching those criteria")
    rm(lhs)

    have.iso <- sapply(hs, function(x) !is.null(x$isos))
    if (sum(have.iso)==0) stop("No hullsets found with isopleths. See lhs.iso.add()")
    if (sum(have.iso)<length(hs)) {
        hs_with_iso <- names(hs)[have.iso]
        hs <- lhs.select(hs, hs.names=hs_with_iso)
    }
    
    ## Create a collection of the name(s) of run(s) to include
    if (is.null(hs.names)) hs.names <- names(hs)
    
    hsnames_matrix <- hs.names
    if (hsnames_simplify) {
        ids_all <- sapply(hs, function(x) x$id) 
        if (anyDuplicated(ids_all)) {
            stop("You can't simplify hullset names to their ids because ids are not unique. Set hsnames_simplify=FALSE")
        } else {
            hsnames_matrix <- as.character(ids_all)
        }
    }

    ## First, we make a list of all the hullsets pairs
    hs1_hs2 <- t(combn(hs.names, m=2, simplify=TRUE))
    
    if (status) {
        cat("Computing isopleth overlap\n")
        flush.console()
    }
    
    overlap_df <- NULL
    intersect_Polygons <- list()
    overlap_id <- 1
    iso_levels_with_overlap <- NULL
    compared_at_least_one_pair <- FALSE
    prj <- NULL
    slivers_cleaned <- NULL

    if (status) pb <- txtProgressBar(min=0, max=nrow(hs1_hs2), style=3)
    ## Loop through the pairs of hullsets
    for (i in 1:nrow(hs1_hs2)) {
        if (status) setTxtProgressBar(pb, i)
        
        hs1.name <- hs1_hs2[i,1]
        hs2.name <- hs1_hs2[i,2]
    
        ## Within a pairs of hullsets, loop through the isopleths objects in hs1 and h22
        for (hs1.iso.name in names(hs[[hs1.name]][["isos"]])) {
            
            ## Get the isopleth levels in hs1
            hs1.iso.levels <- hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]][["iso.level"]]
            if (is.null(iso.level)) {
                hs1.iso.levels.use <- hs1.iso.levels
            } else {
                hs1.iso.levels.use <- intersect(iso.level, hs1.iso.levels)
            }
            
            for (hs2.iso.name in names(hs[[hs2.name]][["isos"]])) {
                ## Get the isopleth levels in hs1
                hs2.iso.levels <- hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]][["iso.level"]]
                
                ## Intersect with passed iso.levels, if needed
                if (is.null(iso.level)) {
                    hs2.iso.levels.use <- hs2.iso.levels
                } else {
                    hs2.iso.levels.use <- intersect(iso.level, hs2.iso.levels)
                }
                
                ## Loop through the joint iso.levels
                for (this_iso_level in intersect(hs1.iso.levels.use, hs2.iso.levels.use)) {
                    compared_at_least_one_pair <- TRUE
                    
                    iso1 <- hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]][ hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]][["iso.level"]]==this_iso_level,]
                    iso2 <- hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]][ hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]][["iso.level"]]==this_iso_level,]
                    
                    ## Take the intersection of these two isopleths
                    #iso1_iso2 <- try(rgeos::gIntersection(iso1, iso2, id=sprintf("%04d", overlap_id)), silent=TRUE)
                    iso1_iso2 <- try(rgeos::gIntersection(iso1, iso2, id=as.character(overlap_id)), silent=TRUE)
                    
                    ## If one of the isopleth had a hole (or main piece) that had < 3 unit coordinates, it could trigger an error
                    if (class(iso1_iso2) == "try-error") {
                        #intersect_errors <- rbind(intersect_errors, data.frame(hs1=hs1.name, hs2=hs2.name, iso_level=this_iso_level))
                        
                        ## We will try to delete any invalid holes (or other polygons)
                        ## We clean out the slivers from the hs object, so we don't have to do it more than once
                        
                        iso1_all_cleaned_lst <- tlocoh::clean_slivers(hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]], status=FALSE)
                        if (!is.null(iso1_all_cleaned_lst$results)) {
                            hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]] <- iso1_all_cleaned_lst$sp
                            slivers_cleaned <- rbind(slivers_cleaned, data.frame(hs=hs1.name, iso=hs1.iso.name))
                        }
                        
                        iso2_all_cleaned_lst <- tlocoh:::clean_slivers(hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]], status=FALSE)
                        if (!is.null(iso2_all_cleaned_lst$results)) {
                            hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]] <- iso2_all_cleaned_lst$sp
                            slivers_cleaned <- rbind(slivers_cleaned, data.frame(hs=hs2.name, iso=hs2.iso.name))
                        }
                        
                        #hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]]  <- tlocoh:::clean_slivers(hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]], status=FALSE)
                        #hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]]  <- tlocoh:::clean_slivers(hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]], status=FALSE)
                        
                        ## Redefine iso1 and iso2
                        iso1 <- hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]][ hs[[hs1.name]][["isos"]][[hs1.iso.name]][["polys"]][["iso.level"]]==this_iso_level,]
                        iso2 <- hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]][ hs[[hs2.name]][["isos"]][[hs2.iso.name]][["polys"]][["iso.level"]]==this_iso_level,]
                        
                        ## Recompute the intersection
                        #iso1_iso2 <- rgeos::gIntersection(iso1, iso2, id=sprintf("%04d", overlap_id))
                        iso1_iso2 <- rgeos::gIntersection(iso1, iso2, id=as.character(overlap_id))
                    } 
                    
                    ## Check if rgeos returned something other than a SpatialPolygons object
                    ## If it returned a SpatialCollections object (e.g., because there was a shared node between the two ispleths)
                    ## Then just take the polyobj slot (lines and points have no area)
                    
                    if (is(iso1_iso2, "SpatialPoints")) iso1_iso2 <- NULL
                    if (is(iso1_iso2, "SpatialLines")) iso1_iso2 <- NULL
                    if (is(iso1_iso2, "SpatialCollections")) iso1_iso2 <- iso1_iso2@polyobj
                        
                    if (!is.null(iso1_iso2)) {
                        
                        ## If the intersection contains more than 1 polygon, union it into a single (multi-part) polygon object
                        if (length(iso1_iso2) > 1) iso1_iso2 <- unionSpatialPolygons(iso1_iso2, rep(as.character(overlap_id), length(iso1_iso2)))
                        
                        # Store the projection info (first pass only)
                        if (is.null(prj)) prj <- iso1_iso2@proj4string
                        
                        ## Compute the areas
                        iso1_area <- sum(sapply(iso1@polygons, function(x) x@area)) 
                        iso2_area <- sum(sapply(iso2@polygons, function(x) x@area))
                        iso1_iso2_area <- sum(sapply(iso1_iso2@polygons, function(x) x@area))
                        
                        ## Append the polygons[[1]] to a master list, and the dataframe
                        
                        if (hsnames_simplify) {
                            hs1name <- hs[[hs1.name]][["id"]]
                            hs2name <- hs[[hs2.name]][["id"]]
                        } else {
                            hs1name <- hs1.name
                            hs2name <- hs2.name
                        }
                        
                        overlap_df <- rbind(overlap_df, data.frame(hs1name=hs1name, hs2name=hs2name, iso_level=this_iso_level, area=iso1_iso2_area, area_prhs1=iso1_iso2_area/iso1_area, area_prhs2=iso1_iso2_area/iso2_area))
                        intersect_Polygons <- c(intersect_Polygons, iso1_iso2@polygons[[1]])

                        overlap_id <- overlap_id + 1
                        iso_levels_with_overlap <- c(iso_levels_with_overlap, this_iso_level) 
                        
                    }
                
                }
                
            }
            
            
        }
    
    }
    # Done with this
    if (status) close(pb)
    
    if (status && !is.null(slivers_cleaned)) {
        cat("Isopleth slivers deleted from: \n")
        print(slivers_cleaned)
    }
    
    if (!compared_at_least_one_pair) {
        if (status) cat("No pair of hullsets had matching isopleth levels \n")
        return(NULL)
    }
    
    if (length(intersect_Polygons) == 0) {
        cat("No intersections found \n")
        return(NULL)
    } else if (length(intersect_Polygons) != nrow(overlap_df)) {
        stop("Uh oh. The length of the SpatialPolygons object and the length of the data frame do not match.")
    }
    
    ## Create a SpatialPolygons object
    Sr <- sp::SpatialPolygons(intersect_Polygons, proj4string=prj)
    
    ## Create a SpatialPolygonsDataFrame object
    overlaps_spdf <- sp::SpatialPolygonsDataFrame(Sr, data=overlap_df, match.ID=FALSE)
    
    ## Create some blank lists that will be used to store Polygons objects
    proportion_overlap <- list()
    area_overlap <- list()
    
    # Make a vector of the values of hs1name and h2name in overlap_df
    if (hsnames_simplify) {
        hsnames_comp <- hsnames_matrix
    } else {
        hsnames_comp <- hs.names
    }

    # Loop through the iso.levels, and build the matricies for specific iso-levels
    for (isolev in unique(iso_levels_with_overlap)) {
        prp_overlap <- matrix(0, nrow=length(hsnames_matrix), ncol=length(hsnames_matrix))
        colnames(prp_overlap) <- hsnames_matrix
        rownames(prp_overlap) <- hsnames_matrix
        ar_overlap <- prp_overlap
        
        for (j in which(overlap_df[["iso_level"]] == isolev)) {
            
            hs1.idx <- which(hsnames_comp == overlap_df[j, "hs1name"])
            hs2.idx <- which(hsnames_comp == overlap_df[j, "hs2name"])
            
            prp_overlap[hs1.idx, hs2.idx] <- overlap_df[j, "area_prhs1"]
            prp_overlap[hs2.idx, hs1.idx] <- overlap_df[j, "area_prhs2"]
            ar_overlap[hs1.idx, hs2.idx] <- overlap_df[j, "area"]
            ar_overlap[hs2.idx, hs1.idx] <- overlap_df[j, "area"]            
        }
        
        ## Add to the list
        proportion_overlap[[paste0("iso", isolev * 100)]] <- prp_overlap
        area_overlap[[paste0("iso", isolev * 100)]] <- ar_overlap
    
    }
    
    return(list(spdf=overlaps_spdf, overlap_area=area_overlap, overlap_prop=proportion_overlap))  

}
