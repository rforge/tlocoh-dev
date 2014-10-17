pIntersectSat <- function(p1, p2, p1.area=0, p2.area=0) {

    ## Returns TRUE if polygons intersect using the Separating Axis Theorem method

    ## p1 and p2 are two-column data frames or matrix containing the vertices of each polygon.
    ## This 'turbo' version of the function does not do any checking, it presumes that
    ## p1 and p2 are closed (first and last point the same).
    ## This turbo version also does not check if the polygon bounding boxes are disjoint
    ## which is a quick way to rule out that two polygons don't intersect

    #print("letspause");browser()

    ## Compute the edges (expressed as a data frame of delta.x and delta.y)
    edges.p1 <- data.frame(x=diff(p1[,1]), y=diff(p1[,2]))
    edges.p2 <- data.frame(x=diff(p2[,1]), y=diff(p2[,2]))
    edges.all <- rbind(edges.p1, edges.p2)

    ## Compute the length of each edge, then normalize and orthoganalize
    edges.mag <- sqrt(edges.all[,1]^2 + edges.all[,2]^2)
    edges.norm.orth <- data.frame(x = -1 * edges.all[,2] / edges.mag, y = edges.all[,1] / edges.mag)

    ## Loop through the normalized orthogonal edges, skipping those that are duplicates (parallel)
    ## for (i in (1:nrow(edges.norm.orth))[-which(duplicated(ifelse(is.infinite(edges.norm.orth$y / edges.norm.orth$x), Inf, edges.norm.orth$y / edges.norm.orth$x)))]) {

    ## this might be slower than its worth:
    #for (i in (1:nrow(edges.norm.orth))[!duplicated(ifelse(edges.norm.orth$x == 0, Inf, edges.norm.orth$y / edges.norm.orth$x))]) {
    
    for (i in (1:nrow(edges.norm.orth))) {

        ## Project points from both polygons onto the edge
        edge.vec <- as.numeric(edges.norm.orth[i,])
        p1.proj <- p1 %*% edge.vec
        p2.proj <- p2 %*% edge.vec

        p1.proj.range <- range(p1.proj)
        p2.proj.range <- range(p2.proj)
        if (p1.proj.range[2] <= p2.proj.range[1] || p1.proj.range[1] >= p2.proj.range[2]) return(FALSE)

    }
    
    return(TRUE)

    #print("No sep axis, going to look fro containment");browser()

    ## No separating axis, next we need to check for containment
    #if (p1.area < p2.area) {
    #    return(point.in.polygon(point.x=p1[1,1], point.y=p1[1,2], pol.x=p2[,1], pol.y=p2[,2]) != 0)
    #} else {
    #    return(point.in.polygon(point.x=p2[1,1], point.y=p2[1,2], pol.x=p1[,1], pol.y=p1[,2]) != 0)
    #}

}
