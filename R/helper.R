#' Safe-loading an R package.
#'
#' Safe-loading an R package. See 'Details'.
#'
#' The function calls \code{require} to load a package.
#' If the package can be loaded it returns \code{TRUE},
#' else \code{FALSE}. 
#'
#' @param lib A character string; package to be loaded.
#'
#' @keywords internal
#' 
#' @return A logical scalar.
SafeLoad <- function(lib) {
    # Safe-loading an R package.
    #
    # Args:
    #   lib: Name of the package to be loaded
    #
    # Returns:
    #   A logical scalar; TRUE if lib was loaded successfully
    ret <- suppressMessages(require(lib,
                                    character.only = TRUE,
                                    quietly = TRUE));
    return(ret);
}


#' Check if object has class classType.
#'
#' Check if object has class classType.
#' 
#' @param object An R object.
#' @param classType String of class type.
#'
#' @return Logical. \code{TRUE} if object is of class classType.
#'
#' @keywords internal
#' 
#' @export
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


#' Check if two txLoc objects are consistent.
#'
#' Check if two txLoc objects are consistent, i.e. are based on the same
#' reference genome and contain the same transcript regions.
#' 
#' @param obj1 A txLoc object.
#' @param obj2 A txLoc object.
#'
#' @return Logical. \code{TRUE} if obj1 and obj2 are consistent.
#'
#' @keywords internal
#' 
#' @export
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


#' Calculate 95% confidence interval from data using empirical
#' bootstrap.
#'
#' Calculate 95% confidence interval from data using empirical
#' bootstrap.
#' 
#' @param x A data vector.
#' @param breaks Vector of integer breaks for binning x.
#' @param nBS Number of boostrap samples. Default is 5000.
#'
#' @return List of upper and lower 95% confidence interval
#' bounds for every bin value.
#'
#' @keywords internal
#' 
#' @export
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

#' Calculate mutual distances between entries for every transcript
#' region from two txLoc objects.
#'
#' Calculate mutual distances between entries for every transcript
#' region from two txLoc objects.
#' 
#' @param loc1 A txLoc object.
#' @param loc2 A txLoc object.
#' @param filter Only consider loci in transcript regions specified in filter. Default is NULL.
#' @param method Method to calculate relative distances. Possible arguments
#' are "ss" (start-start), "mm" (midpoint-midpoint), "se" (start-end), "es"
#' (end-start), "ee" (end-end).
#'
#' @return List of upper and lower 95% confidence interval
#' bounds for every bin value.
#' 
#' @export
GetRelativeDistance <- function(loc1,
                                loc2,
                                filter = NULL,
                                method = c("ss", "mm", "se", "es", "ee")) {
    CheckClassTxLocConsistency(loc1, loc2);
    method <- match.arg(method);
    id1 <- GetId(loc1);
    id2 <- GetId(loc2);
    refGenome <- GetRef(loc1);
    loc1 <- GetLoci(loc1);
    loc2 <- GetLoci(loc2);
    dist.list <- list();
    for (i in 1:length(loc1)) {
        txID <- intersect(loc1[[i]]$REFSEQ, loc2[[i]]$REFSEQ);
        if (length(txID) == 0) {
            dist.list[[length(dist.list)+1]] <- 0;
        } else {
            dist <- vector();
            for (j in 1:length(txID)) {
                loc1.sel <- loc1[[i]][which(loc1[[i]]$REFSEQ == txID[j]), ];
                loc2.sel <- loc2[[i]][which(loc2[[i]]$REFSEQ == txID[j]), ];
                if (method == "ss") {
                    dist <- c(dist,
                              as.vector(outer(loc1.sel$TXSTART,
                                              loc2.sel$TXSTART,
                                              "-")));
                } else if (method == "mm") {
                    dist <- c(dist,
                              as.vector(outer((loc1.sel$TXSTART + loc1.sel$TXEND) / 2,
                                              (loc2.sel$TXSTART + loc2.sel$TXEND) / 2,
                                              "-")));
                } else if (method == "se") {
                    dist <- c(dist,
                              as.vector(outer(loc1.sel$TXSTART,
                                              loc2.sel$TXEND,
                                              "-")));
                } else if (method == "es") {
                    dist <- c(dist,
                              as.vector(outer(loc1.sel$TXEND,
                                              loc2.sel$TXSTART,
                                              "-")));
                } else if (method == "ee") {
                    dist <- c(dist,
                              as.vector(outer(loc1.sel$TXEND,
                                              loc2.sel$TXEND,
                                              "-")));
                }
            }
            dist.list[[length(dist.list) + 1]] <- dist;
        }
    }
    names(dist.list) <- names(loc1);
}

