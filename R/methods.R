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


#' Generate a null distribution of transcript sites.
#'
#' Generate a null distribution of transcript sites. See 'Details'.
#'
#' The function generates a null distribution of single-nucleotide
#' sites in transcript sections. Two different methods can be
#' employed:
#' \enumerate{
#' \item If \code{method == "ntAbund"} sites are generated based
#' on the position of all nucleotides \code{nt} in those transcript
#' sections that also contain a site in \code{locus}.
#' \item If \code{method == "perm"} sites are generated based on
#' sites in \code{locus} by permuting the position of those sites
#' uniformly across the associated transcript section.
#' }
#' It is import to emphasise that any downstream enrichment
#' analysis may depend critically on the choice of the null
#' distribution. For example, a position permution-based null
#' distribution may not be a valid null distribution, if the
#' position distribution is highly non-uniform across a transcript
#' section. This is the case e.g. for the spatial distribution of
#' cytosines across the 5'UTR, CDS and/or 3'UTR. In this case, a
#' better null distribution would be to consider all cytosines
#' in transcript regions that contain also a site in \code{locus}.
#' The return object is a \code{txLoc} of transcript loci in
#' different transcript sections.
#' 
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
#'
#' @return A \code{txLoc} object. See 'Details'.
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
#' negSites <- GenerateSNMNull(posSites, method = "permutation");
#' }
#'
#' @export
GenerateNull <- function(locus,
                         id = NULL,
                         method = c("ntAbund", "perm"),
                         nt = "C")  {
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
    for (i in 1:length(locus)) {
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
        freq.list[[length(freq.list) + 1]] <- cbind(sectionGC, siteGC);
    }
    names(freq.list) <- names(locus);
    return(freq.list);
}


#' Get exon-exon junctions from transcriptome.
#'
#' Get exon-exon junctions from transcriptome. See 'Details'.
#'
#' This function extracts exon-exon junction positions from a
#' reference transcriptome, and returns a \code{txLoc} object.
#'
#' @param refGenome A character string; specifies a specific
#' reference genome assembly version based on which the matching
#' transcriptome is loaded; default is \code{"hg38"}.
#' @param filter A character vector; only consider transcript
#' sections specified in \code{filter}; \code{filter = NULL}
#' corresponds to \code{filter = c("5'UTR", "CDS", "3'UTR")}.
#'
#' @return A \code{txLoc} object. See 'Details'.
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




