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
#' @param hsnames_simplify If True, will simplify the hullset names to just the IDs for the row and column names of the overlap matrix
#'
#' @details
#' This function computes the area intersection of isopleths among different hullsets. 
#' This might be done, for example, if the hullsets belong to different individuals, and you 
#' want to see which individuals share space. All pairs of hullsets in \code{lhs} will be compared,
#' and all isopleth levels will be compared. 
#' 
#' @return A list object with three named elements. The \emph{spdf} contains a SpatialPolygonsDataFrame
#' with the actual areas of intersection expressed in map units, and as the proportion of the isopleth area 
#' for each hullset.
#'
#' @seealso \code{\link{isopleths}}
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
        if (!anyDuplicated(ids_all)) {
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
                    
                    iso1_iso2 <- gIntersection(iso1, iso2, id=sprintf("%04d", overlap_id))    
                    
                    if (!is.null(iso1_iso2)) {
                    
                        # Save the projection info
                        if (is.null(prj)) prj <- iso1_iso2@proj4string
                        
                        # Check for an error
                        if (length(iso1_iso2@polygons) != 1) stop("Uh oh. length(iso1_iso2@polygons) != 1")
                        
                        iso1_area <- iso1@polygons[[1]]@area
                        iso2_area <- iso2@polygons[[1]]@area

                        ## Append the polygons[[1]]@Polygons slot to a master list, and the dataframe
                        
                        overlap_df <- rbind(overlap_df, data.frame(hs1name=hs1.name, hs2name=hs2.name, iso_level=this_iso_level, area=iso1_iso2@polygons[[1]]@area, area_prhs1=iso1_iso2@polygons[[1]]@area/iso1_area, area_prhs2=iso1_iso2@polygons[[1]]@area/iso2_area))
                        intersect_Polygons <- c(intersect_Polygons, iso1_iso2@polygons[[1]])

                        overlap_id <- overlap_id + 1
                        iso_levels_with_overlap <- c(iso_levels_with_overlap, this_iso_level) 
                        
                        ## Check for something more complex
                        if (length(iso1_iso2@polygons)>1) stop("Hold your horses - iso1_iso2@polygons has more than one element!")
                        
                    }
                
                }
                
            }
            
            
        }
    
    }
    if (status) close(pb)
    
    if (!compared_at_least_one_pair) {
        if (status) cat("No pair of hullsets had matching isopleth levels \n")
        return(NULL)
    }
    
    if (length(intersect_Polygons) != nrow(overlap_df)) stop("Uh oh. The length of the SpatialPolygons object and the length of the data frame do not match.")
    
    ## Create a SpatialPolygons object
    Sr <- SpatialPolygons(intersect_Polygons, proj4string=prj)
    
    ## Create a SpatialPolygonsDataFrame object
    overlaps_spdf <- SpatialPolygonsDataFrame(Sr, data=overlap_df, match.ID=FALSE)
    
    ## Create some blank lists that will be used to store Polygons objects
    proportion_overlap <- list()
    area_overlap <- list()
    
    # Loop through the iso.levels, and build the matricies for specific iso-levels
    for (isolev in unique(iso_levels_with_overlap)) {
        prp_overlap <- matrix(0, nrow=length(hsnames_matrix), ncol=length(hsnames_matrix))
        colnames(prp_overlap) <- hsnames_matrix
        rownames(prp_overlap) <- hsnames_matrix
        ar_overlap <- prp_overlap
        
        for (j in which(overlap_df[["iso_level"]] == isolev)) {
            hs1.idx <- which(hs.names == overlap_df[j, "hs1name"])
            hs2.idx <- which(hs.names == overlap_df[j, "hs2name"])
            prp_overlap[hs1.idx, hs2.idx] <-  overlap_df[j, "area_prhs1"]
            prp_overlap[hs2.idx, hs1.idx] <-  overlap_df[j, "area_prhs2"]
            ar_overlap[hs1.idx, hs2.idx] <-  overlap_df[j, "area"]
            ar_overlap[hs2.idx, hs1.idx] <-  overlap_df[j, "area"]
            
        }
        
        ## Add to the list
        proportion_overlap[[paste0("iso", isolev * 100)]] <- prp_overlap
        area_overlap[[paste0("iso", isolev * 100)]] <- ar_overlap
    
    }
    
    return(list(spdf=overlaps_spdf, overlap_area=area_overlap, overlap_prop=proportion_overlap))  

}
