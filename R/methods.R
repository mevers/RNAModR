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
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
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
        dataQuery <- GenomicRanges::as.data.frame(locus[idxLoc]);
        colnames(dataQuery) <- c("CHR", "START", "STOP", "WIDTH",
                                 "STRAND", "SCORE", "ID");
        dataRef <- GenomicRanges::as.data.frame(
            range(txBySec[[i]][idxTx]))[, -1];
        idxSeq <- which(names(txBySec)[i] == names(seqBySec));
        if (length(idxSeq) > 0) {
            dataSeq <- as.character(seqBySec[[idxSeq]][match(
                dataRef[, 1],
                names(seqBySec[[idxSeq]]))]);
        } else {
            dataSeq <- rep("", nrow(dataRef));
        }
        dataRef <- cbind(rep(names(txBySec)[i], nrow(dataRef)),
                             dataRef[, 1],
                             geneXID[match(dataRef[, 1], geneXID[, 1]), ][, -1],
                             dataRef[, 2:ncol(dataRef)],
                             sum(width(txBySec[[i]][idxTx])),
                             dataSeq);
        colnames(dataRef) <- c("GENE_REGION", "GENE_REFSEQ", "GENE_ENTREZ",
                               "GENE_SYMBOL", "GENE_ENSEMBL", "GENE_UNIGENE",
                               "GENE_NAME", "GENE_CHR", "GENE_START",
                               "GENE_STOP", "GENE_WIDTH", "GENE_STRAND",
                               "REGION_TXWIDTH", "REGION_SEQ");
        locusInTx <- cbind(locusInTx,
                           dataQuery,
                           dataRef);
        if (noFactors) {
            locusInTx <- Unfactor(locusInTx);
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
#' reference genome assembly version based on which the matching
#' transcriptome is loaded; default is \code{"hg38"}.
#' @param ignore.strand A logical scalar; if \code{TRUE} strand
#' information is ignored when mapping genome coordinates to
#' transcript coordinates; default is \code{FALSE}.
#'
#' @return A \code{txLoc} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicFeatures GenomicRanges methods
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

#' Construct \code{txLoc}-suitable \code{dataframe} object from
#' mapping transcript loci to genomic loci.
#' 
#' Construct \code{txLoc}-suitable \code{dataframe} object from
#' mapping transcript loci to genomic loci.
#'
#' @param gr A \code{GRanges} object; specifies transcript loci.
#' @param ref A \code{GRanges} object; the reference transcripts.
#' @param seq A \code{DNAStringSet}; the reference transcript
#' sequences.
#' @param section A character scalar; specifies the transcript
#' section.
#' @param geneXID A \code{data.frame}; output of function
#' \code{GetGeneIds} providing gene IDs from different gene
#' annotations.
#' 
#' @return A \code{dataframe} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import Biostrings GenomicRanges IRanges
#'
#' @keywords internal
#' 
#' @export
GetLocus.MapFromTranscripts <- function(gr, ref, seq, section, geneXID) {
    genomePos <- mapFromTranscripts(gr, ref);
    idxQuery <- mcols(genomePos)$xHits;
    idxRef <- mcols(genomePos)$transcriptsHits;
    dataRef <- ref[idxRef];
    idxSeq <- match(names(genomePos), names(seq));
    dataSeq <- seq[idxSeq];
    locus <- cbind.data.frame(
        names(genomePos),
        start(gr[idxQuery]), end(gr[idxQuery]), width(gr[idxQuery]),
        as.character(seqnames(genomePos)),
        start(genomePos), end(genomePos), width(genomePos),
        as.character(strand(genomePos)),
        rep(0, length(genomePos)),
        mcols(gr[idxQuery])$id,
        rep(section, length(genomePos)),
        geneXID[match(names(genomePos), geneXID$REFSEQ), ],
        as.character(seqnames(unlist(range(dataRef)))),
        start(unlist(range(dataRef))), end(unlist(range(dataRef))),
        width(unlist(range(dataRef))),
        as.character(strand(unlist(range(dataRef)))),
        sum(width(dataRef)),
        as.character(dataSeq));
    colnames(locus) <- c("REFSEQ", "TXSTART", "TXEND", "TXWIDTH",
                         "CHR", "START", "STOP", "WIDTH",
                         "STRAND", "SCORE", "ID",
                         "GENE_REGION", "GENE_REFSEQ", "GENE_ENTREZ",
                         "GENE_SYMBOL", "GENE_ENSEMBL", "GENE_UNIGENE",
                         "GENE_NAME", "GENE_CHR", "GENE_START",
                         "GENE_STOP", "GENE_WIDTH", "GENE_STRAND",
                         "REGION_TXWIDTH", "REGION_SEQ");
    return(locus);
}


#' Generate a list of null sites.
#'
#' Generate a list of null sites. See 'Details'.
#'
#' The function generates a null distribution of single-nucleotide
#' sites across different transcript sections, and returns a 
#' \code{txLoc} object.
#' Two different methods can be employed:
#' \enumerate{
#' \item \code{method = "ntAbund"}: Null sites are generated based
#' on the position of all non-modified nucleotides of type \code{nt} 
#' in those transcript sections that also contain a modified site of
#' the same type and as specified in \code{locus}. 
#' For example, if \code{locus} contains a list of m$^{6}$A sites, 
#' the list of null sites consists of all non-methylated adenosines 
#' in those transcripts that contain at least one m$^{6}$A site.
#' \item \code{method = "perm"}: Null sites are generated by 
#' permuting the position of sites from \code{locus} uniformly within
#' the corresponding transcript section. Note that this will generate
#' a list of null sites with the same abundance ratios across 
#' transcript sections as the list of sites from \code{locus}. It is
#' therefore not useful for assessing an enrichment of sites within
#' a particular transcript section.
#' }
#' It is import to emphasise that any downstream enrichment analysis
#' may depend critically on the choice of the null distribution. 
#' For example, a position permution-based null distribution may not 
#' be a valid null distribution, if the nucleotide position 
#' distribution is highly non-uniform across a transcript section. 
#' This is the case e.g. for the spatial distribution of cytosines 
#' within and across the 5'UTR, CDS and/or 3'UTR. In this case, a
#' better null distribution would be to consider all cytosines
#' in transcript sections that also contain a site in \code{locus}.
#' This can be achieved with \code{method = "ntAbund"}. A still 
#' better list of null sites would be based on all \emph{expressed} 
#' and non-modified cytosines based on the same sequencing data that
#' was used to identify modified cytosines. 
#'
#' @param locus A \code{txLoc} object.
#' @param id A character string; identifier for null sites; if
#' \code{NULL} then id = \code{"null"}; default is \code{NULL}.
#' @param method A character string; specifies the method used
#' to generate null distribution; if \code{method == "ntAbund"}
#' the position of all nucleotides specified by \code{nt} will
#' be used as null sites; if \code{method == "perm"} then null
#' sites will be generated by uniform-randomly shuffling candidate
#' positions from \code{locus} within the corresponding transcript
#' region; default is \code{"ntAbund"}.
#' @param nt A single character; if \code{method == "ntAbund"},
#' use \code{nt} to derive distribution of null sites.
#' @param showPb A logical scalar; if \code{TRUE} show a progress
#' bar; default is \code{TRUE}.
#'
#' @return A \code{txLoc} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#' 
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateNull(posSites, method = "ntAbund", nt = "A");
#' }
#'
#' @export
GenerateNull <- function(locus,
                         id = NULL,
                         method = c("ntAbund", "perm"),
                         nt = "C",
                         showPb = TRUE)  {
    # Generate null distribuion of SNM's across different transcript regions.
    #
    # Args:
    #   locus: List of dataframes with mapped SNM's across different
    #          transcript regions.
    #   method: Method used to generate null distribution.
    #           If method == "ntAbund" then the distribution of all
    #           nucleotides specified by nucleotide will be used as null sites.
    #           If method == "perm" then null sites will be generated
    #           by uniformly randomly permuting candidate positions from locus
    #           within transcript region. Default is "ntAbund".
    #   nt: Nucleotide to be used for deriving a list of null sites,
    #               if method == "ntAbund"
    #   showPb: If TRUE, show progress bar
    #
    # Returns:
    #    A txLoc object. Note that genome coordinates are not available for
    #    null sites.
    CheckClass(locus, "txLoc");
    method <- match.arg(method);
    refGenome <- GetRef(locus);
    if (is.null(id)) {
        id <- sprintf("null_%s", GetId(locus));
    }
    locus <- GetLoci(locus);
    locusNull.list <- list();
    if (method == "ntAbund") {
# FIX THIS!!!
#            LoadRefTx(refGenome, geneXID, seqBySec, txBySec);
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
        geneXID <- get("geneXID");
        seqBySec <- get("seqBySec");
        txBySec <- get("txBySec");
    }
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(locus), style = 3, width = 60);
    for (i in 1:length(locus)) {
        if (showPb) setTxtProgressBar(pb, i);
        if (method == "ntAbund") {
            seqData <- locus[[i]][!duplicated(locus[[i]][, 1]),
                                  c("REFSEQ", "GENE_CHR", "GENE_START",
                                    "GENE_STOP", "STRAND", "REGION_SEQ")];
            if (all(grepl("\\w+", seqData$REGION_SEQ))) {
                # Positive sites that should be excluded from null
                exclude <- locus[[i]][, c("REFSEQ", "TXSTART")];
                # All sites of a specific nucleotide
                txPos <- gregexpr(nt, seqData$REGION_SEQ);
                names(txPos) <- seqData$REFSEQ;
                txPos <- data.frame(
                    ID = rep(names(txPos), sapply(txPos, length)),
                    x = unlist(txPos), stringsAsFactors = FALSE);
                # Exclude positive sites
                txPos <- rbind(txPos, setNames(exclude, names(txPos)));
                txPos <- txPos[!duplicated(txPos, fromLast = FALSE) &
                               !duplicated(txPos, fromLast = TRUE), ];
                rownames(txPos) <- seq(1, nrow(txPos));
                # Map sites back to genome
                gr <- GRanges(txPos$ID, IRanges(txPos$x, txPos$x));
                idxSec <- which(names(locus)[i] == names(txBySec));
                genomePos <- mapFromTranscripts(gr, txBySec[[idxSec]]);
                genomePos <- as.data.frame(genomePos);           
            } else {
                sec <- names(locus)[i];
                ss <- sprintf("Cannot generate null for %s", sec);
                ss <- sprintf("%s: Missing sequence information.", ss);
                ss <- sprintf("%s\nEither exclude %s from downstream analysis",
                              ss, sec);
                ss <- sprintf("%s, or use different null.", ss);
                warning(ss);
                txPos <- cbind.data.frame(
                    ID = seqData$REFSEQ,
                    x = rep("*", nrow(seqData)));
                genomePos <- cbind.data.frame(
                    start = rep("*", nrow(seqData)),
                    end = rep("*", nrow(seqData)),
                    width = rep(1, nrow(seqData)));
            }
            idxSeq <- match(txPos$ID, seqData$REFSEQ);
            idxLoc <- match(txPos$ID, locus[[i]]$GENE_REFSEQ);
            locusNull <- cbind.data.frame(
                txPos, txPos$x, rep(1, nrow(txPos)),
                seqData$GENE_CHR[idxSeq],
                genomePos$start,
                genomePos$end,
                genomePos$width,
                seqData$STRAND[idxSeq],
                rep(0, nrow(txPos)),
                rep(sprintf("nucl_%s", nt), nrow(txPos)),
                locus[[i]][idxLoc, 12:ncol(locus[[i]])]);      
            colnames(locusNull) <- colnames(locus[[i]]);
            locusNull <- Unfactor(locusNull);
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
    if (showPb) close(pb);
    names(locusNull.list) <- names(locus);
    obj <- new("txLoc",
               loci = locusNull.list,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    return(obj);
}


GenerateNull.new <- function(locus,
                            id = NULL,
                            method = c("ntAbund", "perm"),
                            nt = "C",
                            showPb = TRUE)  {
    # Generate null distribuion of SNM's across different transcript regions.
    #
    # Args:
    #   locus: List of dataframes with mapped SNM's across different
    #          transcript regions.
    #   method: Method used to generate null distribution.
    #           If method == "ntAbund" then the distribution of all
    #           nucleotides specified by nucleotide will be used as null sites.
    #           If method == "perm" then null sites will be generated
    #           by uniformly randomly permuting candidate positions from locus
    #           within transcript region. Default is "ntAbund".
    #   nt: Nucleotide to be used for deriving a list of null sites,
    #               if method == "ntAbund"
    #   showPb: If TRUE, show progress bar
    #
    # Returns:
    #    A txLoc object. Note that genome coordinates are not available for
    #    null sites.
    CheckClass(locus, "txLoc");
    method <- match.arg(method);
    refGenome <- GetRef(locus);
    if (is.null(id)) {
        id <- sprintf("null_%s", GetId(locus));
    }
    locus <- GetLoci(locus);
    locusNull.list <- list();
    if (method == "ntAbund") {
        LoadRefTx(refGenome);
        geneXID <- get("geneXID");
        seqBySec <- get("seqBySec");
        txBySec <- get("txBySec");
    }
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(locus), style = 3, width = 60);
    for (i in 1:length(locus)) {
        if (showPb) setTxtProgressBar(pb, i);
        if (method == "ntAbund") {
            seqData <- locus[[i]][!duplicated(locus[[i]][, 1]),
                                  c("REFSEQ", "GENE_CHR", "GENE_START",
                                    "GENE_STOP", "STRAND", "REGION_SEQ")];
            if (all(grepl("\\w+", seqData$REGION_SEQ))) {
                # Positive sites that should be excluded from null
                exclude <- locus[[i]][, c("REFSEQ", "TXSTART")];
                # All sites of a specific nucleotide
                txPos <- gregexpr(nt, seqData$REGION_SEQ);
                names(txPos) <- seqData$REFSEQ;
                txPos <- data.frame(
                    ID = rep(names(txPos), sapply(txPos, length)),
                    x = unlist(txPos), stringsAsFactors = FALSE);
                # Exclude positive sites
                txPos <- rbind(txPos, setNames(exclude, names(txPos)));
                txPos <- txPos[!duplicated(txPos, fromLast = FALSE) &
                               !duplicated(txPos, fromLast = TRUE), ];
                rownames(txPos) <- seq(1, nrow(txPos));
                # Map sites back to genome
                gr <- GRanges(txPos$ID, IRanges(txPos$x, txPos$x));
                idxSec <- which(names(locus)[i] == names(txBySec));
                genomePos <- mapFromTranscripts(gr, txBySec[[idxSec]]);
                genomePos <- as.data.frame(genomePos);           
            } else {
                sec <- names(locus)[i];
                ss <- sprintf("Cannot generate null for %s", sec);
                ss <- sprintf("%s: Missing sequence information.", ss);
                ss <- sprintf("%s\nEither exclude %s from downstream analysis",
                              ss, sec);
                ss <- sprintf("%s, or use different null.", ss);
                warning(ss);
                txPos <- cbind.data.frame(
                    ID = seqData$REFSEQ,
                    x = rep("*", nrow(seqData)));
                genomePos <- cbind.data.frame(
                    start = rep("*", nrow(seqData)),
                    end = rep("*", nrow(seqData)),
                    width = rep(1, nrow(seqData)));
            }
            idxSeq <- match(txPos$ID, seqData$REFSEQ);
            idxLoc <- match(txPos$ID, locus[[i]]$GENE_REFSEQ);
            locusNull <- cbind.data.frame(
                txPos, txPos$x, rep(1, nrow(txPos)),
                seqData$GENE_CHR[idxSeq],
                genomePos$start,
                genomePos$end,
                genomePos$width,
                seqData$STRAND[idxSeq],
                rep(0, nrow(txPos)),
                rep(sprintf("nucl_%s", nt), nrow(txPos)),
                locus[[i]][idxLoc, 12:ncol(locus[[i]])]);      
            colnames(locusNull) <- colnames(locus[[i]]);
            locusNull <- Unfactor(locusNull);
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
    if (showPb) close(pb);
    names(locusNull.list) <- names(locus);
    obj <- new("txLoc",
               loci = locusNull.list,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    return(obj);
}


#' Calculate GC content within window around loci from a \code{txLoc}
#' object.
#'
#' Calculate GC content within window around loci from a \code{txLoc}
#' object. See 'Details'.
#'
#' The function calculates the GC content within a region around every
#' transcript locus from a \code{txLoc} object. The window is defined
#' by extending the position of every transcript locus upstream and
#' downstream by \code{flank} nucleotides (if possible).
#'
#' @param locus A \code{txLoc} object.
#' @param flank An integer scalar; see 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @import Biostrings
#'
#' @export
GetGC <- function(locus, flank = 10) {
    # Calculate GC contents of region and within window around site.
    #
    # Args:
    #   locus: List of dataframes with mapped features across different
    #          transcript regions
    #
    # Returns:
    #   locus with appended columns for GC content
    CheckClass(locus, "txLoc");
    id <- GetId(locus);
    refGenome <- GetRef(locus);
    locus <- GetLoci(locus);
    if (!SafeLoad("Biostrings")) {
        stop("Could not load library Biostrings.");
    }
    freq.list <- list();
    for (i in 1:length(locus)) {
        seq <- BStringSet(locus[[i]]$REGION_SEQ);
        nucFreq <- letterFrequency(
            seq,
            letters = c("A", "C", "G", "T"));
        sectionGC <- rowSums(nucFreq[, c(2, 3)]) / width(seq);
        if (is.numeric(locus[[i]]$TXSTART)) {
            x1 <- locus[[i]]$TXSTART - flank;
            x2 <- locus[[i]]$TXSTART + flank;
            subSeq <- BStringSet(substr(locus[[i]]$REGION_SEQ, x1, x2));
            nucFreq <- letterFrequency(subSeq,
                                       letters = c("A", "C", "G", "T"));
            siteGC <- rowSums(nucFreq[, c(2, 3)]) / nchar(subSeq);
        } else {
            siteGC <- rep(NA, nrow(locus[[i]]));
        }
        dfGC <- cbind(sectionGC, siteGC);
        rownames(dfGC) <- locus[[i]]$ID;
        freq.list[[length(freq.list) + 1]] <- dfGC;
    }
    names(freq.list) <- names(locus);
    return(freq.list);
}


#' Get exon-exon junctions.
#'
#' Get exon-exon junctions from transcriptome. See 'Details'.
#'
#' The function extracts exon-exon junction positions from a
#' reference transcriptome specified by \code{refGenome}, and 
#' returns a \code{txLoc} object of exon-exon junctions.
#' The position of an exon-exon junction is defined as the 
#' position of the 3'-last nucleotide of an exon followed by
#' an intron. 
#'
#' @param refGenome A character string; specifies a specific
#' reference genome assembly version based on which the matching
#' transcriptome is loaded; default is \code{"hg38"}.
#' @param filter A character vector; only consider transcript
#' sections specified in \code{filter}; default is 
#' \code{"CDS"}.
#'
#' @return A \code{txLoc} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#'
#' @export
GetEEJunct <- function(refGenome = "hg38", filter = "CDS") {
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
    geneXID <- get("geneXID");
    seqBySec <- get("seqBySec");
    txBySec <- get("txBySec");
    if (is.null(filter)) {
        filter <- c("5'UTR", "CDS", "3'UTR");
    }
    sel <- which(names(txBySec) %in% filter);
    locus <- GRanges();
    for (i in 1:length(sel)) {
        junct <- psetdiff(range(txBySec[[sel[i]]]), txBySec[[sel[i]]]);
        junct <- unlist(junct[elementLengths(junct) > 0]);
        eej <- GRanges(
            seqnames(junct),
            IRanges(ifelse(strand(junct) == "+",
                           start(junct) - 1,
                           end(junct) + 1),
                    ifelse(strand(junct) == "+",
                           start(junct) - 1,
                           end(junct) + 1)),
            strand(junct),
            score = 0,
            id = sprintf("eej|%s|%s",
                gene = names(junct),
                section = names(txBySec)[sel[i]]));
        locus <- append(locus, eej);
    }
    locusInTx.list <- SmartMap.ToTx(
        locus, txBySec, seqBySec, geneXID,
        ignore.strand = FALSE,
        showPb = TRUE);
    obj <- new("txLoc",
               loci = locusInTx.list,
               id = "eej",
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    obj <- FilterTxLoc(obj, filter);
    return(obj);
}


#' Get loci of motif(s).
#'
#' Get loci of motif(s) from transcriptome. See 'Details'.
#'
#' The function searches for one or multiple motifs within
#' sequences of a reference transcriptome specified by 
#' \code{refGenome}, and returns a \code{txLoc} object of the 
#' motif loci within different transcript sections. 
#' The maximum number of mismatches allowed in the motif search 
#' can be adjusted through \code{maxMM}. By default a text 
#' progressbar is shown \code{showPb = TRUE}.
#' Note that the motif search may take a few minutes, depending
#' on the size of the transcriptome and number of motifs.  
#'
#' @param motif A character vector; specifies the motif(s)
#' that will be matched against the transcriptome.
#' @param refGenome A character string; specifies the reference
#' genome version; default is \code{"hg38"}.
#' @param filter A character vector; only consider transcript
#' sections specified in \code{filter}; default is 
#' \code{c("5'UTR", "CDS", "3'UTR")}.
#' @param maxMM An integer scalar; specifies the maximum number
#' of mismatches that are allowed during the motif matching;
#' default is 0.
#' @param showPb A logical scalar; if \code{TRUE} show a progress
#' bar; default is \code{TRUE}.

#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#' @importFrom Biostrings vmatchPattern
#'
#' @export
GetMotifLoc <- function(motif, 
                        refGenome = "hg38", 
                        filter = c("5'UTR", "CDS", "3'UTR"), 
                        maxMM = 0, 
                        showPb = TRUE) {
    PASmotif <- c("AATAAA", "ATTAAA", "AGTAAA",
                  "TATAAA", "AAGAAA", "AATACA",
                  "AATATA", "CATAAA", "AATGAA",
                  "GATAAA", "ACTAAA", "AATAGA");
    if (is.null(motif)) {
      motif <- PASmotif;
    }
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
    geneXID <- get("geneXID");
    seqBySec <- get("seqBySec");
    txBySec <- get("txBySec");
    if (is.null(filter)) {
        filter <- c("5'UTR", "CDS", "3'UTR");
    }
    sel <- which(names(txBySec) %in% filter);
    locus.list <- list();
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(sel), style = 3, width = 60);
    for (i in 1:length(sel)) {
        if (showPb) setTxtProgressBar(pb, i);
        if (!IsEmptyChar(seqBySec[[sel[i]]])) {
# Disable warnings to get rid of messages of the form "Each of the
# 2 combined objects has sequence levels not in the other" when
# appending to GRanges object. Ugly! Is there a cleaner solution?
# Works for now, fix later.
            options(warn = -1);
            gr <- GRanges();
            for (j in 1:length(motif)) {
                m <- vmatchPattern(motif[j],
                                   seqBySec[[sel[i]]],
                                   max.mismatch = maxMM);
                m <- unlist(m);
                grMotif <- GRanges(names(m), m,
                                   id = rep(motif[j], length(m)));
                gr <- append(gr, grMotif);
            }
            options(warn = 0);
            locus <- GetLocus.MapFromTranscripts(
                gr,
                txBySec[[sel[i]]],
                seqBySec[[sel[i]]],
                names(txBySec)[sel[i]],
                geneXID);
            locus.list[[length(locus.list) + 1]] <- locus;
        }
    }
    if (showPb) close(pb);
    names(locus.list) <- names(seqBySec)[sel];
    obj <- new("txLoc",
               loci = locus.list,
               id = "motifs",
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    return(obj);
}

