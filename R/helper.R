CheckClass <- function(object, classType = NULL) {
    # Check that object is of type classType
    #
    # Args:
    #   object: Object.
    #   classType: Class type of object.
    #
    # Returns:
    #    TRUE
    if (class(object) != classType) {
        objName <- deparse(substitute(object));
        stop(sprintf("%s is not an object of type %s.", object, classType));
    }
    return(TRUE);
}


CheckClassTxLocConsistency <- function(obj1, obj2) {
    # Check that the underlying properties of txloc
    # objects obj1 and obj2 are consistent.
    #
    # Args:
    #   obj1: txLoc object.
    #   obj2: txLoc object
    #
    # Returns:
    #   TRUE
    ref1 <- slot(obj1, "refGenome");
    ref2 <- slot(obj2, "refGenome");
    objName1 <- deparse(substitute(obj1));
    objName2 <- deparse(substitute(obj2));
    if (ref1 != ref2) {
        stop("%s and %s are not based on the same reference genome: %s != %s.",
             objName1, objName2, ref1, ref2);
    }
    if (!identical(names(slot(obj1, "loci")),
                   names(slot(obj2, "loci")))) {
        stop("Transcript regions in %s and %s do not match.",
             objName1, objName2);
    }
    return(TRUE);
}


EstimateCIFromBS <- function(x, breaks, nBS = 5000) {
    # Estimate 95% CI from empirical bootstrap
    #
    # Args:
    #   x: Data
    #   breaks: Vector of integer breaks for binning x.
    #   nBS: Number of bootstrap samples. Default is 5000.
    #
    # Returns:
    #    List of lower/upper 95% CI values for each bin.
    h0 <- hist(x, breaks = breaks, plot = FALSE);
    matBS <- matrix(0, ncol = length(breaks) - 1, nrow = nBS);
    for (j in 1:nBS) {
        xBS <- sample(x, size = length(x), replace = TRUE);
        h <- hist(xBS, breaks = breaks, plot = FALSE);
        matBS[j, ] <- h$counts;
    }
    sd <- apply(matBS, 2, sd);
    z1 <- apply(matBS, 2, function(x) {(x - mean(x)) / sd(x)});
    z1[is.na(z1)] <- 0;
    z1 <- apply(z1, 2, sort);
    y.low <- h0$counts - z1[round(0.025 * nBS), ] * sd;
    y.high <- h0$counts - z1[round(0.975 * nBS), ] * sd;
    return(list(x = h0$mids,
                y.low = y.low,
                y.high = y.high));
}
