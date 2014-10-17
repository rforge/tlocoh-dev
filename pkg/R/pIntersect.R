
pIntersect <- function(p1, p2, p1.area=0, p2.area=0, filter.bb=TRUE) {

    ## Returns TRUE if polygons intersect. Intersection is true if:
    ## 1) if any two ling segments cross, or
    ## 2) the smaller polygon is fully enclosed by the larger polygon
    
    ## p1 and p2 are two-column data frames or matrix containing the vertices of each polygon.
    ## This 'turbo' version of the function does not do any checking, it presumes that
    ## p1 and p2 are closed (first and last point the same). Also presumes package sp is loaded
    ## This turbo version also does not check if the polygon bounding boxes are disjoint
    
    ## Reformat the points into a 4 column line segment matrix (x1,y1,x2,y2)
    p1.fp <- 1:(nrow(p1)-1); p1.sp <- 2:nrow(p1)
    l1.mat <- cbind(x1=p1[p1.fp,1], y1=p1[p1.fp,2], x2=p1[p1.sp,1], y2=p1[p1.sp,2])
    p2.fp <- 1:(nrow(p2)-1); p2.sp <- 2:nrow(p2)
    l2.mat <- cbind(x1=p2[p2.fp,1], y1=p2[p2.fp,2], x2=p2[p2.sp,1], y2=p2[p2.sp,2])

    ## To save time, we only check intersection of those line segements whose bounding boxes are not disjoint
    
    ## Identify the bounding box for each line segment
    l1bb <- t(rbind(apply(l1.mat[,c(1,3)], 1, range), apply(l1.mat[,c(2,4)], 1, range)))
    l2bb <- t(rbind(apply(l2.mat[,c(1,3)], 1, range), apply(l2.mat[,c(2,4)], 1, range)))
    
    ## Get all combinations of line segments, and find out which are disjoint and could not possibly intersect
    l1l2.idx <- expand.grid(1:(nrow(p1)-1), 1:(nrow(p2)-1))
    #p1p2.bb <- cbind(p1p2.idx, p1bb[p1p2.idx[,1],], p2bb[p1p2.idx[,2],]) 
    l1l2.bb <- cbind(l1bb[l1l2.idx[,1],], l2bb[l1l2.idx[,2],]) 
    disjoint <- l1l2.bb[,1] > l1l2.bb[,6] | l1l2.bb[,2] < l1l2.bb[,5] | l1l2.bb[,3] > l1l2.bb[,8] | l1l2.bb[,4] < l1l2.bb[,7]
    
    row.idx.of.pairs.that.overlap <- (1:nrow(l1l2.idx))[!disjoint]
    for (i in row.idx.of.pairs.that.overlap) {
        #if (lines.intersect.turbo(p1[c(l1l2.idx[i,1],l1l2.idx[i,1]+1),], p2[c(l1l2.idx[i,2],l1l2.idx[i,2]+1),])) return(TRUE)    
        #num.sides.checked <<- c(num.sides.checked, i)
        if (lines.intersect.turbo.flat(l1.mat[l1l2.idx[i,1],], l2.mat[l1l2.idx[i,2],])) {
            
            #assign("pInter.log", rbind(get("pInter.log", envir=.GlobalEnv), data.frame(idx=match(i, row.idx.of.pairs.that.overlap), num=length(row.idx.of.pairs.that.overlap))), envir=.GlobalEnv)
            #writeLines(paste("\n======\n", "hulls2iso: ", Sys.time(), "\n========", sep=""), fdebug)
            
            return(TRUE)
        }
    }

    ## So no overlapping edges. Next we need to test for containment of the smaller poly in the larger poly
    ## Since there are no intersecting edges, either all of the points in the smaller poly
    ## lie within the larger poly (containment), or none of them do.
    ## So we only need to test one point. If we had tested for containment before testing
    ## for edge intersections, we would have had to test all points in the smaller polygon

    if (p1.area < p2.area) {
        return(point.in.polygon(point.x=p1[1,1], point.y=p1[1,2], pol.x=p2[,1], pol.y=p2[,2]) != 0)
    } else {
        return(point.in.polygon(point.x=p2[1,1], point.y=p2[1,2], pol.x=p1[,1], pol.y=p1[,2]) != 0)
    }
}
