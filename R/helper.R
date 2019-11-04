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


#' Match `flank` function argument
#' 
#' @param flank An \code{integer} scalar or an \code{nteger} vector of length 2.
#'
#' @return An \code{integer} vector of length 2.
#' 
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
MatchFlank <- function(flank) {

    if (length(flank) == 1) {
        if (flank < 0) stop("`flank` cannot be negative!")
        return(c(-abs(flank), abs(flank)))
    } else if (length(flank) == 2) {
        if (any(flank < 0)) stop("`flank` cannot be negative!")
        return(c(-abs(flank[1]), abs(flank[2])))
    } else {
        stop("`flank` needs to be a scalar or a vector of length 2!")
    }

}


#' Load reference transcriptome.
#'
#' Load reference transcriptome. See 'Details'.
#'
#' The function loads transcriptome data generated from \code{BuildTx()} and
#' stored in a \code{.RData} file, and makes the objects assessible in the
#' parent enviroment. There should not be any need to call this function
#' directly.
#' Reference transcriptome data includes the following objects:
#' \enumerate{
#'     \item \code{txBySec}: A \code{list} of \code{GRangesList} objects
#'     \item \code{seqBySec}: A \code{list} of \code{DNAStringSet} objects
#'     \item \code{geneXID}: A \code{DataFrame} with transcript and gene IDs
#' }
#'
#' @param refGenome A \code{character} string; specifies a specific reference
#' genome assembly and gene annotation version.
#' @param verbose A \code{logical} scalar; determines the amount of output;
#' default is \code{FALSE}.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @importFrom utils packageVersion
#'
#' @export
LoadRefTx <- function(refGenome, verbose = FALSE) {

    # Check that RData file exists
    refTx <- sprintf("tx_%s.RData", refGenome)
    if (!file.exists(refTx)) {
        ss <- sprintf(
            "Reference transcriptome data for %s not found.",
            refGenome)
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }

    # Load and check objects
    load(refTx)
    requiredObj <- c("geneXID", "seqBySec", "txBySec")
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript data objects not found.")
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "))
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }

    # Get objets and check versions
    txBySec <- get("txBySec")
    seqBySec <- get("seqBySec")
    geneXID <- get("geneXID")
    if (is.null(attr(txBySec, "package_version"))) {
        ss <- "Transcriptome data are based on an unknown version of RNAModR:"
        ss <- sprintf(
            "%s\n    RNAModR version: %s, transcriptome data version: <NULL>", 
            ss, packageVersion("RNAModR"))
        ss <- sprintf(
            "%s\n  %s",
            ss,
            "There is not guarantee that RNAModR functions will work.")
        ss <- sprintf(
            "%s %s",
            ss,
            sprintf(
                "Consider re-running `BuildTx(\"%s\")`.",
                refGenome))
        warning(ss)
    } else {
        version <- attr(txBySec, "package_version")
        if (version$major != packageVersion("RNAModR")$major) {
            ss <- "Transcriptome data are based on an older version of RNAModR:"
            ss <- sprintf(
                "%s\n    RNAModR version: %s, transcriptome data version: %s", 
                ss, packageVersion("RNAModR"), version)
            ss <- sprintf(
                "%s\n  %s",
                ss,
                "Difference in major version. RNAModR may not work!")
            ss <- sprintf(
                "%s\n  %s",
                ss,
                sprintf(
                    "Consider re-running `BuildTx(\"%s\")`.",
                    refGenome))
            warning(ss)
            
        } else if (version$major == packageVersion("RNAModR")$major & 
                   version$minor != packageVersion("RNAModR")$minor) {
            ss <- "Transcriptome data are based on an older version of RNAModR:"
            ss <- sprintf(
                "%s\n    RNAModR version: %s, transcriptome data version: %s", 
                ss, packageVersion("RNAModR"), attr(txBySec, "package_version"))
            ss <- sprintf(
                "%s\n  %s",
                ss,
                "Difference in minor version. RNAModR might not work!")
            ss <- sprintf(
                "%s\n  %s",
                ss,
                sprintf(
                    "Consider re-running `BuildTx(\"%s\")`.",
                    refGenome))
            warning(ss)
        } else if (version$major == packageVersion("RNAModR")$major & 
                   version$minor == packageVersion("RNAModR")$minor & 
                   version$patchlevel != packageVersion("RNAModR")$patchlevel &
                   verbose == TRUE) {
            ss <- "Transcriptome data are based on an older version of RNAModR:"
            ss <- sprintf(
                "%s\n    RNAModR version: %s, transcriptome data version: %s", 
                ss, packageVersion("RNAModR"), attr(txBySec, "package_version"))
            ss <- sprintf(
                "%s\n  %s",
                ss,
                "Difference in patchlevel.")
            ss <- sprintf(
                "%s\n  %s",
                ss,
                sprintf(
                    "Consider re-running `BuildTx(\"%s\")`.",
                    refGenome))
            warning(ss)
        }
    }
    
    # Assign objects to parent environment
    assign("txBySec", txBySec, envir = parent.frame(1L))
    assign("seqBySec", seqBySec, envir = parent.frame(1L))
    assign("geneXID", geneXID, envir = parent.frame(1L))

    if (verbose) cat("Loaded reference transcriptome data.\n")

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


#' Calculate transcript region starting coordinates.
#'
#' Calculate transcript region starting coordinates along full transcript from
#' a reference transcriptome. See 'Details'.
#'
#' The function loads transcriptome data generated from \code{BuildTx()} and
#' stored in a \code{.RData} file, and calculates starting coordinates of the
#' transcript regions 5'UTR, CDS, 3'UTR along the full transcript. There should 
#' not be any need to call this function directly. The function returns a 
#' \code{list} with the following elements:
#' \enumerate{
#'     \item \code{start_region}: A \code{data.frame} of transcript region 
#'     start coordinates
#'     \item \code{gr.start}: A \code{GRanges} object of start codon transcript
#'     coordinates
#'     \item \code{gr.stop}: A \code{GRanges} object of stop codon transcript
#'     coordinates
#' } 
#'
#' @param refGenome A \code{character} string.
#'
#' @return A \code{list}. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
GetTxRegionCoordinates <- function(refGenome) {
    
    # Read transcriptome data
    seqBySec <- txBySec <- geneXID <- NULL
    LoadRefTx(refGenome)
    
    # Get transcript region starting coordinates
    regions <- c("5'UTR", "CDS", "3'UTR")
    df <- Reduce(
        function(x, y) merge(x, y, by = "seqnames"),
        Map(
            function(seq, reg) setNames(
                data.frame(names(seq), width(seq)),
                c("seqnames", reg)),
            seqBySec[regions],
            regions)
    )
    df[, regions] <- Reduce(`+`, df[, regions], accumulate = TRUE)
    
    # Return `list` consisting of the following elements
    #   - A `data.frame` of transcript region start coordinates
    #   - A `GRanges` object of start codon transcript coordinates
    #   - A `GRanges` object of stop codon transcript coordinates
    list(
        start_region = df,
        gr.start = GRanges(
            seqnames = df[, "seqnames"],
            IRanges(df[, "CDS"], df[, "CDS"] + 3),
            strand = "*"),
        gr.stop = GRanges(
            seqnames = df[, "seqnames"],
            IRanges(df[, "3'UTR"] - 4, df[, "3'UTR"] - 1),
            strand = "*"))

}
    

#' Get distance of sites to the nearest start/stop codon.
#'
#' Get distance of sites from a \code{txLoc} to the nearest start/stop codon.
#'
#' The function converts transcript region coordinates from a \code{txLoc} to 
#' transcript coordinates, and returns distances of sites to the nearest 
#' start/stop codon. A transcript is defined as the concatenation of
#' the following transcript regions: 5'UTR, CDS, 3'UTR. The return object is
#' a \code{list} of distances of sites from \code{txLoc} to the nearest start 
#' and stop codons.
#'
#' @param txLoc A \code{txLoc} object.
#'
#' @return A \code{list}. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#'
#' @keywords internal
GetDistNearestStartStop <- function(txLoc) {
    
    # Sanity checks
    CheckClass(txLoc, "txLoc")
    regions <- c("5'UTR", "CDS", "3'UTR")
    objName <- deparse(substitute(txLoc))
    if (!identical(GetRegions(txLoc), regions)) stop(sprintf(
        "%s must contain the regions %s!", 
        objName,
        paste(regions, collapse = ", ")))

    # Get transcript region start coordinates along transcript and
    # coordinates of start/stop codons
    lst <- GetTxRegionCoordinates(GetRef(txLoc))
    
    # Get transcript region coordinates of loci
    df.loci <- do.call(
        rbind, c(Map(
            function(locus, reg)
                setNames(cbind(
                    as.data.frame(locus$locus_in_tx_region)[, 1:3],
                    reg),
                    c("seqnames", "start", "end", "region")),
            GetLoci(txLoc),
            regions),
            make.row.names = FALSE)
    )

    # Merge transcript region coordinates and widths
    df.loci <- merge(
        df.loci, 
        setNames(
            cbind(
                lst[["start_region"]][, "seqnames"], 
                stack(lst[["start_region"]], select = -seqnames)),
            c("seqnames", "start_region", "region")),
        by = c("seqnames", "region"))
    
    # Convert transcript region to transcript coordinates
    df.loci[, "start"] <- df.loci[, "start"] + df.loci[, "start_region"]
    df.loci[, "end"] <- df.loci[, "end"] + df.loci[, "start_region"]

    # Store loci as `GRanges`
    gr.loci <- GRanges(
        seqnames = df.loci[, "seqnames"],
        IRanges(df.loci[, "start"], df.loci[, "end"]),
        strand = "*")
    
    # Calculate distances to nearest start/stop codons
    GetRelDistNearest(
        list(start = gr.loci, stop = gr.loci), 
        list(start = lst[["gr.start"]], stop = lst[["gr.stop"]]))
    
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
            
            # Make sure that seqlevels match
            lvls <- intersect(seqlevels(gr1), seqlevels(gr2))
            seqlevels(gr1, pruning.mode = "coarse") <- lvls
            seqlevels(gr2, pruning.mode = "coarse") <- lvls

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
