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


#' Load reference transcriptome.
#'
#' Load reference transcriptome. See 'Details'.
#'
#' The function loads transcriptome data stored in a \code{.RData}
#' file, and makes the objects assessible in the user's workspace.
#'
#' @param refGenome A character string; specifies a specific
#' reference genome assembly version based on which a transcriptome
#' is loaded; default is \code{"hg38"}.
#' @param env An \code{environment} object; default is the user's
#' workspace, i.e. \code{env = .GlobalEnv}.
#'
#' @keywords internal
#' 
#' @export
LoadRefTx <- function(refGenome = "hg38", env = .GlobalEnv) {
    refTx <- sprintf("tx_%s.RData", refGenome);
    if (!file.exists(refTx)) {
        ss <- sprintf("Reference transcriptome for %s not found.", refGenome);
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome);
        stop(ss);
    }
    load(refTx);
    requiredObj <- c("geneXID", "seqBySec", "txBySec");
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript objects not found.");
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "));
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome);
        stop(ss);
    }
    geneXID <- base::get("geneXID");
    seqBySec <- base::get("seqBySec");
    txBySec <- base::get("txBySec");
    assign("geneXID", geneXID, envir = env);
    assign("seqBySec", seqBySec, envir = env);
    assign("txBySec", txBySec, envir = env);
}


#' Filter sections of a \code{txLoc} object.
#'
#' Filter sections of a \code{txLoc} object.
#'
#' @param locus A \code{txLoc} object.
#' @param filter A character vector; only keep transcript sections
#' specified in \code{filter}; if \code{NULL} consider all sections.
#'
#' @return A \code{txLoc} object.
#' 
#' @export
FilterTxLoc <- function(locus, filter = NULL) {
    CheckClass(locus, "txLoc");
    id <- GetId(locus);
    refGenome <- GetRef(locus);
    version <- GetVersion(locus);
    locus <- GetLoci(locus);
    if (!is.null(filter)) {
        locus <- locus[which(names(locus) %in% filter)];
    }
    obj <- new("txLoc",
               loci = locus,
               id = id,
               refGenome = refGenome,
               version = version);
    return(obj);
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


#' Add transparency to a list of hex colors
#'
#' Add transparency to a list of colors
#'
#' @param hexList A vector or list of character strings.
#' @param alpha A real scalar; specifies the transparency; default
#' is \code{alpha = 0.5}.
#'
#' @keywords internal
#' 
#' @export
AddAlpha <- function(hexList, alpha = 0.5) {
    mat <- sapply(hexList, col2rgb, alpha = TRUE) / 255.0;
    mat[4, ] <- alpha;
    col <- vector();
    for (i in 1:ncol(mat)) {
        col <- c(col, rgb(mat[1, i],
                          mat[2, i],
                          mat[3, i],
                          mat[4, i]));
    }
    return(col);
}


#' Check if all entries in a character vector are empty.
#'
#' Check if all entries in a character vector are empty.
#' 
#' @param v A character vector.
#'
#' @return A logical scalar.
#'
#' @keywords internal
#'
#' @export
IsEmptyChar <- function(v) {
    return(all(nchar(v) == 0));
}


#' Return specific colour palette.
#'
#' Return specific colour palette.
#'
#' @param pal A character string.
#' @param n An integer scalar.
#' @param alpha A real scalar.
#'
#' @return A character vector.
#'
#' @keywords internal
#'
#' @export
GetColPal <- function(pal = c("apple", "google"), n = NULL, alpha = 1.0) {
    pal <- match.arg(pal);
    alpha <- alpha * 255;
    if (pal == "apple") {
        col <- c(
            rgb(95, 178, 51, alpha = alpha, maxColorValue = 255),
            rgb(106, 127, 147, alpha = alpha, maxColorValue = 255),
            rgb(245, 114, 6, alpha = alpha, maxColorValue = 255),
            rgb(235, 15, 19, alpha = alpha, maxColorValue = 255),
            rgb(143, 47, 139, alpha = alpha, maxColorValue = 255),
            rgb(19, 150, 219, alpha = alpha, maxColorValue = 255));
    } else if (pal == "google") {
        col <- c(
            rgb(61, 121, 243, alpha = alpha, maxColorValue = 255),
            rgb(230, 53, 47, alpha = alpha, maxColorValue = 255),
            rgb(249, 185, 10, alpha = alpha, maxColorValue = 255),
            rgb(52, 167, 75, alpha = alpha, maxColorValue = 255));
    }
    if (!is.null(n)) {
        col <- col[1:n];
    }
    return(col);
}
