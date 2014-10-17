lines.intersect.turbo.flat <- function(L1, L2) {

    ## Returns a T/F

    ## L1 and L2 are each a four item numeric vector: x1,y1,x2,y2
    
    ## previously was 
    ## x1, y1,
    ## x2, y2
     
    
    ## Adapted from http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/Helpers.vb
    
    ## We presume that the bounding boxes of the two lines have already been checked for disjoint. If not, these lines can help:
    #L1.sort <- apply(L1,2,range)
    #L2.sort <- apply(L2,2,range)
    #if (L1.sort[1,1] > L2.sort[2,1] || L1.sort[2,1] < L2.sort[1,1] || L1.sort[1,2] > L2.sort[2,2] || L1.sort[2,2] < L2.sort[1,2]) {
    #    print("guess what, found a set of lines whose bounding box doesn't overlap");browser()
    #    return(FALSE)
    #}

    ## Dim d As Double = (L2.Y2 - L2.Y1) * (L1.X2 - L1.X1) - (L2.X2 - L2.X1) * (L1.Y2 - L1.Y1)
    d <- (L2[4] - L2[2]) * (L1[3] - L1[1]) - (L2[3] - L2[1]) * (L1[4] - L1[2])

    # Prevent a division by zero - this also indicates that the lines are parallel.
    if (d == 0) return(FALSE)

    ## n_a and n_b are calculated as seperate values for readability
    ## Dim n_a As Double = (L2.X2 - L2.X1) * (L1.Y1 - L2.Y1) - (L2.Y2 - L2.Y1) * (L1.X1 - L2.X1)
    n_a <- (L2[3] - L2[1]) * (L1[2] - L2[2]) - (L2[4] - L2[2]) * (L1[1] - L2[1])

    ## Dim n_b As Double = (L1.X2 - L1.X1) * (L1.Y1 - L2.Y1) - (L1.Y2 - L1.Y1) * (L1.X1 - L2.X1)
    n_b <- (L1[3] - L1[1]) * (L1[2] - L2[2]) - (L1[4] - L1[2]) * (L1[1] - L2[1])

    # If n_a and n_b were both equal to zero the lines would be on top of each
    # other (coincidental).  This check is not done because it is not
    # necessary for this implementation (the parallel check accounts for this).

    ## Calculate the intermediate fractional point that the lines potentially intersect.
    ## Dim ua As Double = n_a / d
    ua <- n_a / d

    ##Dim ub As Double = n_b / d
    ub <- n_b / d

    ## The fractional point will be between 0 and 1 inclusive if the lines
    ## intersect.  If the fractional calculation is larger than 1 or smaller
    ## than 0 the lines would need to be longer to intersect.

    return(ua >= 0 && ua <= 1 && ub > 0 && ub < 1)

}
