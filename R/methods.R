#' Map genome coordinates to transcript coordinates.
#'
#' Map genome coordinates to transcript coordinates. See
#' 'Details'.
#' This is a low-level function that is being called from
#' \code{SmartMap}.
#'
#' The function maps genomic coordinates from \code{locus} to
#' transcript section coordinates from \code{txBySec}. The
#' function returns a \code{list} of \code{data.frame} objects
#' with additional gene ID and sequence information.
#' 
#' @param locus A \code{GRanges} object; specifies the list of
#' genomic features to be mapped.
#' @param txBySec A \code{list} of \code{GRangesList} objects;
#' specifies the reference transcriptome; the object is usually
#' a result of running \code{BuildTx}.
#' @param seqBySec A \code{list} of \code{DNAStringSet} objects;
#' sequences of (some or all) segments from \code{txBySec};
#' the object is usually a result of running \code{BuildTx}.
#' @param geneXID A \code{data.frame}; specifies different gene
#' IDs; the object is usually a result of running \code{BuildTx}.
#' @param ignore.strand A logical scalar; if \code{TRUE} strand
#' information is ignored when mapping genome coordinates to
#' transcript coordinates; default is \code{FALSE}.
#' @param noFactors A logical scalar; if \code{TRUE}, the return
#' object contains no \code{factors}; default is \code{TRUE}.
#' @param showPb A logical scalar; if \code{TRUE} show a progress
#' bar; default is \code{FALSE}.
#'
#' @return A \code{list} of \code{data.frame} objects. See
#' 'Details'.
#'
#' @export
SmartMap.ToTx <- function(locus,
                          txBySec,
                          seqBySec,
                          geneXID,
                          ignore.strand = FALSE,
                          noFactors = TRUE,
                          showPb = FALSE) {
    # Map positions from locus to transcript regions.
    #
    # Args:
    #   locus: Loci of genomic features to be mapped as GRanges object.
    #   txBySec: List of GRangesList transcript features.
    #   seqBySec: List of DNAStringSet sequences.
    #   geneXID: Dataframe of cross-referencing gene IDs.
    #   ignore.strand: Ignore strand information during mapping.
    #                  Default is FALSE.
    #   noFactors: No (character) factors in output. Default is TRUE.
    #
    # Returns:
    #   List of dataframes with mapped transcript coordinates for all
    #   features in locus.
    locusInTx.list <- list();
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(txBySec), style = 3, width = 60);
    for (i in 1:length(txBySec)) {
        if (showPb) setTxtProgressBar(pb, i);
        gr <- mapToTranscripts(locus, txBySec[[i]], ignore.strand = ignore.strand);
        idxLoc <- mcols(gr)$xHits;
        idxTx  <- mcols(gr)$transcriptsHits;
        # THIS COULD DO WITH SOME TIDYING UP ...
        # Collate data from mapping, locus, txBySec, seqBySec and geneXID
#        idxSeq <- which(names(txBySec)[i] == names(seqBySec));
#        dataMap <- GenomicRanges::as.data.frame(gr);
#        dataLoc <- GenomicRanges::as.data.frame(locus[idxLoc]);
#        dataTx <-  GenomicRanges::as.data.frame(range(txBySec[[i]][idxTx]));
#        if (length(idxSeq) > 0) {
#            dataSeq <- seqBySec[[idxSeq]][match(
#                dataMap[, 1],
#                names(seqBySec[[idxSeq]]))];
#            dataSeq <- as.character(dataSeq);
#        } else {
#            dataSeq <- rep("", nrow(dataTx));
#        }
#        dataID <-  geneXID[match(dataMap[, 1], geneXID[, 1]), ];
        locusInTx <- GenomicRanges::as.data.frame(gr)[, 1:4];
        colnames(locusInTx) <- c("REFSEQ", "TXSTART", "TXEND", "TXWIDTH");
        dataFromQuery <- GenomicRanges::as.data.frame(locus[idxLoc]);
        colnames(dataFromQuery) <- c("CHR", "START", "STOP", "WIDTH",
                                     "STRAND", "SCORE", "ID");
        dataFromRef <- GenomicRanges::as.data.frame(
            range(txBySec[[i]][idxTx]))[, -1];
        idxSeq <- which(names(txBySec)[i] == names(seqBySec));
        if (length(idxSeq) > 0) {
            dataSeq <- as.character(seqBySec[[idxSeq]][match(
                dataFromRef[, 1],
                names(seqBySec[[idxSeq]]))]);
        } else {
            dataSeq <- rep("", nrow(dataFromRef));
        }
        dataFromRef <- cbind(rep(names(txBySec)[i], nrow(dataFromRef)),
                             dataFromRef[, 1],
                             geneXID[match(dataFromRef[, 1] ,geneXID[, 1]), ][, -1],
                             dataFromRef[, 2:ncol(dataFromRef)],
                             sum(width(txBySec[[i]][idxTx])),
                             dataSeq);
        colnames(dataFromRef) <- c("GENE_REGION", "GENE_REFSEQ", "GENE_ENTREZ",
                                   "GENE_SYMBOL", "GENE_ENSEMBL", "GENE_UNIGENE",
                                   "GENE_NAME", "GENE_CHR", "GENE_START",
                                   "GENE_STOP", "GENE_WIDTH", "GENE_STRAND",
                                   "REGION_TXWIDTH", "REGION_SEQ");
        locusInTx <- cbind(locusInTx,
                           dataFromQuery,
                           dataFromRef);
        if (noFactors) {
            idx <- sapply(locusInTx, is.factor);
            locusInTx[idx] <- lapply(locusInTx[idx], as.character);
        }
        locusInTx.list[[length(locusInTx.list) + 1]] <- locusInTx;
    }
    if (showPb) close(pb);
    names(locusInTx.list) <- names(txBySec);
    return(locusInTx.list);
}


#' Map genome coordinates to transcript coordinates.
#'
#' Map genome coordinates to transcript coordinates.
#'
#' The function maps genomic coordinates from \code{locus} to
#' transcript section coordinates. The function automatically
#' loads a reference transcriptome based on \code{refGenome}.
#' An error is thrown if a reference transcriptome could not
#' be found. This usually means that \code{BuildTx} was not yet
#' run successfully.
#' The function returns a \code{txLoc} object of mapped positions.
#'
#' @param locus A \code{GRanges} object; specifies the list of
#' of genomic features to be mapped.
#' @param id A character string; specifies a name for loci from
#' \code{locus}; if \code{NULL} then \code{id = ""}; default is
#' \code{NULL}.
#' @param refGenome A character string; specifies a specific
#' reference genome assembly version based on which a transcriptome
#' is loaded; default is \code{"hg38"}.
#' @param ignore.strand A logical scalar; if \code{TRUE} strand
#' information is ignored when mapping genome coordinates to
#' transcript coordinates; default is \code{FALSE}.
#'
#' @return A \code{txLoc} object. See 'Details'.
#' 
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' txSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' }

#' @export
SmartMap <- function(locus,
                     id = NULL,
                     refGenome = "hg38",
                     ignore.strand = FALSE) {
    # User function to map positions from locus to transcript regions.
    #
    # Args:
    #   locus: GRanges object of loci of genomic features.
    #   id: Identifier for loci. Default is NULL.
    #   refGenome: reference genome/transcriptome version
    #   ignore.strand: Ignore strand during mapping. Default is FALSE.
    #
    # Returns:
    #   txLoc object.
    #
    # Load reference transcriptome
    CheckClass(locus, "GRanges");
    refTx <- sprintf("tx_%s.RData", refGenome);
    if (!file.exists(refTx)) {
        ss <- sprintf("Reference transcriptome for %s not found.", refGenome);
        ss <- sprintf("%s\nRunning BuildTx(...) might fix that.", ss);
        stop(ss);
    }
    load(refTx);
    requiredObj <- c("geneXID", "seqBySec", "txBySec");
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript objects not found.");
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "));
        ss <- sprintf("%s\nRunning BuildTx(...) might fix that.", ss);
        stop(ss);
    }
    geneXID <- get("geneXID");
    seqBySec <- get("seqBySec");
    txBySec <- get("txBySec");
    # Map coordinates to transcript
    locusInTx.list <- SmartMap.ToTx(
        locus, txBySec, seqBySec, geneXID, 
        ignore.strand = ignore.strand,
        showPb = TRUE);
    if (is.null(id)) {
        id = "";
    }
    obj <- new("txLoc",
               loci = locusInTx.list,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    return(obj);
}


#' Generate a null distribution of transcript sites based on features from locus.
#'
#' Generate a null distribution of transcript sites based on features from locus.
#'
#' @param locus A txLoc object.
#' @param id Identifier for null sites. If NULL then id = \code{"null"}.
#' Default is NULL.
#' @param method Method used to generate null distribution.
#' If \code{method == "nucAbundance"} then the position of all nucleotides
#' specified by \code{nucleotide} will be used as null sites.
#' If \code{method == "permutation"} then null sites will be generated
#' by uniform-randomly shuffling candidate positions from locus within
#' the corresponding transcript region.
#' Default is \code{"nucAbundance"}.
#' @param nucleotide If \code{method == "nucAbundance"}, use \code{nucleotide}
#' to derive distribution of null sites.
#'
#' @return A \code{txLoc} object. Note that genome coordinates will not be
#' available for null sites.
#' 
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateSNMNull(posSites, method = "permutation");
#' }
#'
#' @export
GenerateSNMNull <- function(locus,
                            id = NULL,
                            method = c("ntAbund", "perm"),
                            nucleotide = "C")  {
    # Generate null distribuion of SNM's across different transcript regions.
    #
    # Args:
    #   locus: List of dataframes with mapped SNM's across different
    #          transcript regions.
    #   method: Method used to generate null distribution.
    #           If method == "nucAbundance" then the distribution of all
    #           nucleotides specified by nucleotide will be used as null sites.
    #           If method == "permutation" then null sites will be generated
    #           by uniformly randomly permuting candidate positions from locus
    #           within transcript region. Default is "nucAbundance".
    #   nucleotide: Nucleotide to be used for deriving a list of null sites,
    #               if method == "nucAbundance"
    #
    # Returns:
    #    A txLoc object. Note that genome coordinates are not available for
    #    null sites.
    CheckClass(locus, "txLoc");
    method <- match.arg(method);
    locusNull.list <- list();
    refGenome <- GetRef(locus);
    if (is.null(id)) {
        id <- sprintf("null_%s", GetId(locus));
    }
    locus <- GetLoci(locus);
    for (i in 1:length(locus)) {
        if (method == "ntAbund") {
            seqData <- locus[[i]][!duplicated(locus[[i]][, 1]),
                                  c("REFSEQ", "GENE_CHR", "GENE_START",
                                    "GENE_STOP", "STRAND", "REGION_SEQ")];
            txPos <- gregexpr(nucleotide, seqData$REGION_SEQ);
            names(txPos) <- seqData$REFSEQ;
            txPos <- data.frame(ID = rep(names(txPos), sapply(txPos, length)),
                              x = unlist(txPos), stringsAsFactors = FALSE);
# FIX THIS!!!
            gr <- GRanges(txPos$ID, IRanges(txPos$x, txPos$x));
#            load(sprintf("tx_%s.RData", refGenome));
#            genomePos <- as.data.frame(mapFromTranscripts(gr, txBySec[[i]]));
            rownames(txPos) <- seq(1, nrow(txPos));
            idxSeq <- match(txPos$ID, seqData$REFSEQ);
            idxLoc <- match(txPos$ID, locus[[i]]$GENE_REFSEQ);
            locusNull <- cbind.data.frame(
                txPos, txPos$x, rep(1, nrow(txPos)),
                seqData$GENE_CHR[idxSeq],
#                seqData$GENE_START[idxSeq] + txPos$x,
#                seqData$GENE_START[idxSeq] + txPos$x,
                rep("*", nrow(txPos)),
                rep("*", nrow(txPos)),
                rep(1, nrow(txPos)),
                seqData$STRAND[idxSeq],
                rep(0, nrow(txPos)),
                rep(sprintf("nucl_%s", nucleotide), nrow(txPos)),
                locus[[i]][idxLoc, 12:ncol(locus[[i]])]);      
            colnames(locusNull) <- colnames(locus[[i]]);
            idx <- sapply(locusNull, is.factor);
            locusNull[idx] <- lapply(locusNull[idx], as.character);
        } else if (method == "perm") {
            locusNull <- locus[[i]];
            locusNull$TXSTART <- apply(
                locusNull, 1, function(x) {
                    round(runif(1,
                                min = 1,
                                max = as.numeric(x["REGION_TXWIDTH"]) - 1))});
            locusNull$TXEND <- locusNull$TXSTART;
            locusNull[, c("CHR", "START", "STOP",
                          "WIDTH", "STRAND", "SCORE")] <- matrix(
                              "*",
                              nrow = nrow(locusNull),
                              ncol = 6);
            locusNull$ID <- "unif_perm";
        }
        locusNull.list[[length(locusNull.list) + 1]] <- locusNull;
    }
    names(locusNull.list) <- names(locus);
    obj <- new("txLoc",
               loci = locusNull.list,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    return(obj);
}


