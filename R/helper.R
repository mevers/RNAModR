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
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
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
                                    quietly = TRUE))
    return(ret)
}


#' Check if object has class classType.
#'
#' Check if object has class classType.
#'
#' @param object An R object.
#' @param classType String of class type for toplevel.
#' @param classType2 String of class type for level 2.
#'
#' @return Logical. \code{TRUE} if object is of class classType.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
CheckClass <- function(object, classType = NULL, classType2 = NULL) {
    # Check that object is of type classType
    #
    # Args:
    #   object: Object.
    #   classType: Class type of object.
    #
    # Returns:
    #    TRUE
    if (class(object) != classType) {
        objName <- deparse(substitute(object))
        stop(sprintf("%s is not an object of type %s.\n",
                     objName, classType),
             call. = FALSE)
    }
    if (classType == "list" &&
        !is.null(classType2) &&
        !all(lapply(object, class) == classType2)) {
        objName <- deparse(substitute(object))
        stop(sprintf("%s does not contain a list of objects of type %s.\n",
                     objName, classType2),
             call. = FALSE)
    }
    return(TRUE)
}


#' Check if entries of two \code{txLoc} objects are based on the
#' same reference genome.
#'
#' Check if entries of two \code{txLoc} objects are based on the
#' same reference genome.
#'
#' @param obj1 A \code{txLoc} object.
#' @param obj2 A \code{txLoc} object.
#'
#' @return A logical scalar.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
CheckClassTxLocRef <- function(obj1, obj2) {

    # Sanity checks
    CheckClass(obj1, "txLoc")
    CheckClass(obj2, "txLoc")

    # Check that reference genome versions match
    obj1Name <- deparse(substitute(obj1))
    obj2Name <- deparse(substitute(obj2))
    ref1 <- GetRef(obj1)
    ref2 <- GetRef(obj2)
    if (ref1 != ref2) {
        ss <- sprintf("%s and %s are not based on the same reference genome.",
                      obj1Name, obj2Name)
        ss <- sprintf("%s\n  %s: %s", ss, obj1Name, ref1)
        ss <- sprintf("%s\n  %s: %s", ss, obj2Name, ref2)
        stop(ss)
    } else TRUE

}


#' Check if entries of two \code{txLoc} objects are consistent.
#'
#' Check if entries of two \code{txLoc} objects are consistent.
#' See 'Details'.
#'
#' The function checks if the reference genome and transcript regions from
#' two \code{txLoc} objects match.
#'
#' @param obj1 A \code{txLoc} object.
#' @param obj2 A \code{txLoc} object.
#'
#' @return A logical scalar.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
CheckClassTxLocConsistency <- function(obj1, obj2) {

    # Check that `obj1` and `obj2` are based on the same reference genome
    CheckClassTxLocRef(obj1, obj2)

    # Check that transcript regions match
    obj1Name <- deparse(substitute(obj1))
    obj2Name <- deparse(substitute(obj2))
    if (!identical(GetRegions(obj1), GetRegions(obj2))) {
        ss <- sprintf(
            "Regions of %s and %s do not match:",
            obj1Name,
            obj2Name)
        ss <- sprintf(
            "%s\n  Regions in %s: %s",
            ss, obj1Name, paste(GetRegions(obj1), collapse = ", "))
        ss <- sprintf(
            "%s\n  Regions in %s: %s",
            ss, obj2Name, paste(GetRegions(obj2), collapse = ", "))
        stop(ss)
    } else TRUE

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
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
LoadRefTx <- function(refGenome = "hg38", env = .GlobalEnv) {
    refTx <- sprintf("tx_%s.RData", refGenome)
    if (!file.exists(refTx)) {
        ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    load(refTx)
    requiredObj <- c("geneXID", "seqBySec", "txBySec")
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript objects not found.")
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "))
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    geneXID <- base::get("geneXID")
    seqBySec <- base::get("seqBySec")
    txBySec <- base::get("txBySec")
    assign("geneXID", geneXID, envir = env)
    assign("seqBySec", seqBySec, envir = env)
    assign("txBySec", txBySec, envir = env)
}


#' Filter sections of a \code{txLoc} object.
#'
#' Filter sections of a \code{txLoc} object.
#'
#' @param txLoc A \code{txLoc} object.
#' @param filter A \code{character} vector; only keep transcript regions
#' that match entries from \code{filter}; if \code{NULL}, return the original
#' \code{txLoc} object.
#'
#' @return A \code{txLoc} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @export
FilterTxLoc <- function(txLoc, filter = NULL) {

    # Sanity check
    CheckClass(txLoc, "txLoc")

    # Get all slots from the `txLoc` object
    id <- GetId(txLoc)
    refGenome <- GetRef(txLoc)
    version <- GetVersion(txLoc)
    loci <- GetLoci(txLoc)

    # Filter regions
    if (!is.null(filter)) loci <- loci[which(names(loci) %in% filter)]

    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = id,
        refGenome = refGenome,
        version = version)

}


#' Downsample a \code{txLoc} object.
#'
#' Downsample a \code{txLoc1} object based on the number of sites per region
#' from a \code{txLoc2} object.
#'
#' @param txLoc1 A \code{txLoc} object; this is the \code{txLoc} object that
#' will be downsampled.
#' @param txLoc2 A \code{txLoc} object; this is the \code{txLoc} object that
#' will be used as a target for the downsampling.
#' @param seed A single value, interpreted as an \code{integer}, or \code{NULL};
#' this is to ensure reproducibility when subsampling \code{txLoc2} sites;
#' default is \code{NULL}.
#'
#' @return A \code{txLoc} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @export
DownsampleTxLoc <- function(txLoc1, txLoc2, seed = NULL) {

    # Sanity check
    CheckClassTxLocConsistency(txLoc1, txLoc2)
    objname1 <- deparse(substitute(txLoc1))
    objname2 <- deparse(substitute(txLoc2))
    if (any(GetNumberOfLoci(txLoc1) <= GetNumberOfLoci(txLoc2))) {
        ss <- sprintf(
            "Cannot downsample %s to %s",
            objname1,
            objname2)
        ss <- sprintf(
            "%s: %s",
            ss,
            sprintf(
                "There are more sites in %s than in %s!",
                objname2,
                objname1))
        stop(ss)
    }

    # Add "_downsampled" to id slot of `txLoc1`
    id <- GetId(txLoc1)
    id <- sprintf("%s_downsampled", id)

    # Get reference genome and version slots from `txLoc1`
    refGenome <- GetRef(txLoc1)
    version <- GetVersion(txLoc1)

    # If required, set fixed random seed
    if (!is.null(seed)) set.seed(seed)

    # Subsample
    loci <- Map(
        function(loci1, loci2) loci1[sample.int(nrow(loci1), nrow(loci2)), ],
        GetLoci(txLoc1),
        GetLoci(txLoc2))

    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = id,
        refGenome = refGenome,
        version = version)

}


#' Subsample a \code{txLoc} object.
#'
#' Subsample a \code{txLoc1} object based on a vector of fractions for every
#' transcript region.
#'
#' @param txLoc A \code{txLoc} object; this is the \code{txLoc} object from
#' which the subsampled \code{txLoc} will be created.
#' @param fractions A \code{numeric} vector, specifying the fraction of sites
#' that will be sampled from every transcript region. Note that the length
#' of \code{fractions} has to match the length of \code{GetRegions(txLoc)}.
#' @param seed A single value, interpreted as an \code{integer}, or \code{NULL};
#' this is to ensure reproducibility when subsampling \code{txLoc2} sites;
#' default is \code{NULL}.
#'
#' @return A \code{txLoc} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @export
SubsampleTxLoc <- function(txLoc, fractions, seed = NULL) {

    # Sanity check
    CheckClass(txLoc, "txLoc")
    objname <- deparse(substitute(txLoc))
    if (length(fractions) != length(GetRegions(txLoc))) {
        ss <- sprintf(
            "Need as many entries in `fractions` as there are regions in `%s`!",
            objname)
        stop(ss)
    }
    if (any(fractions > 1)) {
        ss <- "Entries in `fractions` need to be numbers from [0, 1]!"
        stop(ss)
    }

    # Add "_subsampled" to id slot of `txLoc`
    id <- GetId(txLoc)
    id <- sprintf("%s_subsampled", id)

    # Get reference genome and version slots from `txLoc`
    refGenome <- GetRef(txLoc)
    version <- GetVersion(txLoc)

    # If required, set fixed random seed
    if (!is.null(seed)) set.seed(seed)

    # Subsample
    loci <- Map(
        function(loci, frac) {
            idx <- sample.int(nrow(loci), round(frac * nrow(loci)))
            loci[idx, ]
        },
        GetLoci(txLoc),
        fractions)

    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = id,
        refGenome = refGenome,
        version = version)

}


#' Convert a \code{txLoc} object to a \code{GRangesList} object.
#'
#' Convert a \code{txLoc} object to a \code{GRangesList} object.
#' See 'Details'.
#'
#' The function converts a \code{txLoc} to a \code{GRangesList}
#' object. Coordinates can be either genomic coordinates
#' (\code{method = "genome"}) or transcript region coordinates
#' (\code{method = "tx_region"}).
#'
#' @param txLoc A \code{txLoc} object.
#' @param method A character string; specifies whether coordinates
#' are genome (\code{method = "genome"}) or transcriptome coordinates
#' (\code{method = "tx_region"}).
#'
#' @return A \code{list} of \code{GRanges} objects. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom S4Vectors DataFrame
TxLoc2GRangesList <- function(txLoc, method = c("tx_region", "genome")) {

    # Sanity checks
    CheckClass(txLoc, "txLoc")
    method <- match.arg(method)

    # Return a `list` of `GRanges`
    lapply(GetLoci(txLoc), function(locus) {
        gr <- switch(
            method,
            "genome" = locus$locus_in_genome,
            "tx_region" = locus$locus_in_tx_region)
        mcols(gr) <- DataFrame(
            source = GetId(txLoc),
            tx_refseq = locus$tx_refseq,
            tx_region = locus$tx_region,
            id = locus$id)
        gr
    })

}


#' Return list of nearest distances between entries from two \code{list}s of
#' \code{GRanges} objects.
#'
#' Return list of nearest distances between entries from two \code{list}s of
#' \code{GRanges} objects. See 'Details'.
#'
#' The function uses \code{GenomicRanges::distanceToNearest} to return the
#' nearest distances between the start positions of a \code{GRanges} object
#' from \code{lst1} and the corresponding \code{GRanges} object from
#' \code{lst2} with the same name (based on field \code{seqnames}).
#' Note that distances are given as signed \code{integer}s: Negative
#' distances correspond to pos(gr1) < pos(gr2), positive
#' distances correspond to pos(gr1) > pos(gr2).
#'
#'
#' @param lst1 A \code{list} of \code{GRanges} object.
#' @param lst2 A \code{list} of \code{GRanges} object.
#'
#' @return A list of \code{integer} vectors. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @import GenomicRanges
#' @importFrom S4Vectors queryHits subjectHits
GetRelDistNearest <- function(lst1, lst2) {

    # Sanity checks
    all(sapply(lst1, CheckClass, "GRanges"))
    all(sapply(lst2, CheckClass, "GRanges"))
    stopifnot(identical(names(lst1), names(lst2)))

    # Calculate nearest distances from start
    lst <- Map(
        function(gr1, gr2) {

            # Collapse range of gr1 and gr2 to the start coordinate
            end(gr1) <- start(gr1)
            end(gr2) <- start(gr2)

            # Calculate distance to nearest
            hits <- distanceToNearest(gr1, gr2, ignore.strand = TRUE)

            # Define and return distance d as
            #   d > 0 : if pos(gr1) > pos(gr2)
            #   d < 0 : if pos(gr1) < pos(gr2)
            # In words: Negative distances => gr1 is upstream of gr2
            #           Positive distances => gr1 is downstream of gr2
            ifelse(
                end(gr1)[queryHits(hits)] > start(gr2)[subjectHits(hits)],
                mcols(hits)$distance,
                -mcols(hits)$distance)
        },
        lst1, lst2)

    # Return lst
    lst

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
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom graphics hist
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
    h0 <- hist(x, breaks = breaks, plot = FALSE)
    matBS <- matrix(0, ncol = length(breaks) - 1, nrow = nBS)
    for (j in 1:nBS) {
        xBS <- sample(x, size = length(x), replace = TRUE)
        h <- hist(xBS, breaks = breaks, plot = FALSE)
        matBS[j, ] <- h$counts
    }
    sd <- apply(matBS, 2, sd)
    z1 <- apply(matBS, 2, function(x) {(x - mean(x)) / sd(x)})
    z1[is.na(z1)] <- 0
    z1 <- apply(z1, 2, sort)
    y.low <- h0$counts - z1[round(0.025 * nBS), ] * sd
    y.high <- h0$counts - z1[round(0.975 * nBS), ] * sd
    return(list(x = h0$mids,
                y.low = y.low,
                y.high = y.high))
}


#' Add transparency to a list of hex colors
#'
#' Add transparency to a list of colors
#'
#' @param hexList A vector or list of character strings.
#' @param alpha A float scalar; specifies the transparency; default
#' is \code{alpha = 0.5}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom grDevices col2rgb rgb
#'
#' @export
AddAlpha <- function(hexList, alpha = 0.5) {
    mat <- sapply(hexList, col2rgb, alpha = TRUE) / 255.0
    mat[4, ] <- alpha
    col <- vector()
    for (i in 1:ncol(mat)) {
        col <- c(col, rgb(mat[1, i],
                          mat[2, i],
                          mat[3, i],
                          mat[4, i]))
    }
    return(col)
}


#' Check if all entries in a character vector are empty.
#'
#' Check if all entries in a character vector are empty.
#'
#' @param v A character vector.
#'
#' @return A logical scalar.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
IsEmptyChar <- function(v) {
    return(all(nchar(v) == 0))
}


#' Return specific colour palette.
#'
#' Return specific colour palette.
#'
#' @param pal A character string.
#' @param n An integer scalar.
#' @param alpha A float scalar.
#'
#' @return A character vector.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
GetColPal <- function(pal = c("apple", "google"), n = NULL, alpha = 1.0) {
    pal <- match.arg(pal)
    alpha <- alpha * 255
    if (pal == "apple") {
        col <- c(
            rgb(95, 178, 51, alpha = alpha, maxColorValue = 255),
            rgb(106, 127, 147, alpha = alpha, maxColorValue = 255),
            rgb(245, 114, 6, alpha = alpha, maxColorValue = 255),
            rgb(235, 15, 19, alpha = alpha, maxColorValue = 255),
            rgb(143, 47, 139, alpha = alpha, maxColorValue = 255),
            rgb(19, 150, 219, alpha = alpha, maxColorValue = 255))
    } else if (pal == "google") {
        col <- c(
            rgb(61, 121, 243, alpha = alpha, maxColorValue = 255),
            rgb(230, 53, 47, alpha = alpha, maxColorValue = 255),
            rgb(249, 185, 10, alpha = alpha, maxColorValue = 255),
            rgb(52, 167, 75, alpha = alpha, maxColorValue = 255))
    }
    if (!is.null(n)) {
        col <- col[1:n]
    }
    return(col)
}


#' Unfactor entries in a \code{data.frame}.
#'
#' Unfactor entries in a \code{data.frame}.
#'
#' @param df A \code{data.frame} object.
#'
#' @return A \code{data.frame} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
Unfactor <- function(df) {
    idx <- sapply(df, is.factor)
    df[idx] <- lapply(df[idx], as.character)
    return(df)
}
