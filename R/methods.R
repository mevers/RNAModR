#' Map genome coordinates to transcript coordinates.
#'
#' Map genome coordinates to transcript coordinates. This function
#' should not be invoked by the end-user directly. It is called
#' from within \code{SmartMap}.
#'
#' @param locus GRanges object of genomic features.
#' @param txFeatures List of GRangesList transcript features.
#' @param txSequences List of DNAStringSet sequences.
#' @param geneXID Dataframe of cross-referencing gene IDs.
#' @param ignore.strand Ignore strand when mapping genome coordinates to
#' transcript coordinates. Default is FALSE.
#' @param noFactors No (character) factors in output. Default is TRUE.
#'
#' @return A list of dataframes with mapped features in transcript region.
#'
#' @keywords 
#'
#' @examples
#'
SmartMap.ToTx <- function(locus,
                          txFeatures,
                          txSequences,
                          geneXID,
                          ignore.strand = FALSE,
                          noFactors = TRUE) {
    # Map positions from locus to transcript regions.
    #
    # Args:
    #   locus: Loci of genomic features to be mapped as GRanges object.
    #   txFeatures: List of GRangesList transcript features.
    #   txSequences: List of DNAStringSet sequences.
    #   geneXID: Dataframe of cross-referencing gene IDs.
    #   ignore.strand: Ignore strand information during mapping.
    #                  Default is FALSE.
    #   noFactors: No (character) factors in output. Default is TRUE.
    #
    # Returns:
    #   List of dataframes with mapped transcript coordinates for all
    #   features in locus.
    if (!identical(names(txFeatures), names(txSequences))) {
        stop("Regions in transcript features and sequences do not match.");
    }
    locusInTx.list <- list();
    for (i in 1:length(txFeatures)) {
        gr <- mapToTranscripts(locus, txFeatures[[i]], ignore.strand = ignore.strand);
        locusInTx <- GenomicRanges::as.data.frame(gr)[, 1:4];
        colnames(locusInTx) <- c("REFSEQ", "TXSTART", "TXEND", "TXWIDTH");
        dataFromQuery <- GenomicRanges::as.data.frame(locus[mcols(gr)[, 1]]);
        colnames(dataFromQuery) <- c("CHR", "START", "STOP", "WIDTH",
                                     "STRAND", "SCORE", "ID");
        dataFromRef <- GenomicRanges::as.data.frame(
            range(txFeatures[[i]][mcols(gr)[, 2]]))[, -1];
        dataFromRef <- cbind(rep(names(txFeatures)[i], nrow(dataFromRef)),
                             dataFromRef[, 1],
                             geneXID[match(dataFromRef[, 1] ,geneXID[, 1]), ][, -1],
                             dataFromRef[, 2:ncol(dataFromRef)],
                             sum(width(txFeatures[[i]][mcols(gr)[, 2]])),
                             as.character(seqByFeat[[i]][match(dataFromRef[, 1],
                                                               names(seqByFeat[[i]]))]));
        colnames(dataFromRef) <- c("GENE_REGION", "GENE_REFSEQ", "GENE_ENTREZ",
                                   "GENE_SYMBOL", "GENE_ENSEMBL", "GENE_UNIGENE",
                                   "GENE_NAME", "GENE_GO", "GENE_CHR",
                                   "GENE_START", "GENE_STOP", "GENE_WIDTH",
                                   "GENE_STRAND", "REGION_TXWIDTH", "REGION_SEQ");
        locusInTx <- cbind(locusInTx,
                           dataFromQuery,
                           dataFromRef);
        if (noFactors) {
            idx <- sapply(locusInTx, is.factor);
            locusInTx[idx] <- lapply(locusInTx[idx], as.character);
        }
        locusInTx.list[[length(locusInTx.list) + 1]] <- locusInTx;
    }
    names(locusInTx.list) <- names(txFeatures);
    return(locusInTx.list);
}


#' Map genome coordinates to transcript coordinates.
#'
#' Map genome coordinates to transcript coordinates.
#'
#' @param gr GRanges object of genomic features
#' @param id Identifier for loci from \code{gr}. If NULL then id = "".
#' Default is NULL.
#' @param refGenome Reference transcriptome, defaults to hg38.
#' @param ignore.strand Ignore strand when mapping genome coordinates to
#' transcript coordinates. Default is FALSE.
#'
#' @return A \code{txLoc} object with coordinates from \code{gr} mapped to
#' the transcriptome based on assembly version \code{refGenome}.
#' 
#' @examples
#'
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' txSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");

#' @export
SmartMap <- function(gr,
                     id = NULL,
                     refGenome = "hg38",
                     ignore.strand = FALSE) {
    # User function to map positions from locus to transcript regions.
    #
    # Args:
    #   gr: GRanges object of loci of genomic features.
    #   id: Identifier for loci. Default is NULL.
    #   refGenome: reference genome/transcriptome version
    #   ignore.strand: Ignore strand during mapping. Default is FALSE.
    #
    # Returns:
    #   txLoc object.
    #
    # Load reference transcriptome
    refData <- sprintf("tx_%s", refGenome);
    try(data(list = refData));
    # Map coordinates to transcript
    locusInTx.list <- SmartMap.ToTx(gr, featByTx, seqByFeat, geneXID, 
                               ignore.strand = ignore.strand);
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
#' 
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateSNMNull(posSites, method = "permutation");
#'
#' @export
GenerateSNMNull <- function(locus,
                            id = NULL,
                            method = "nuclAbundance",
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
    locusInTx.list <- list();
    refGenome <- slot(locus, "refGenome");
    if (is.null(id)) {
        id <- slot(locus, "id");
        id <- sprintf("null_%s", id);
    }
    locus <- slot(locus, "loci");
    for (i in 1:length(locus)) {
        if (method == "nuclAbundance") {
            seqData <- locus[[i]][!duplicated(locus[[i]][, 1]),
                                  c(1, ncol(locus[[i]]))];
            pos <- reshape2::melt(lapply(strsplit(seqData$REGION_SEQ, ""),
                               function(x) grep(nucleotide, x)));
            locusInTx <- cbind.data.frame(seqData[pos[, 2], 1],
                                          pos[, 1], pos[, 1],
                                          rep(1, nrow(pos)),
                                          matrix("*",
                                                 nrow = nrow(pos),
                                                 ncol = 6),
                                          rep(sprintf("nucl_%s", nucleotide),
                                              nrow(pos)),
                                          locus[[i]][match(seqData[pos[, 2], 1],
                                                           locus[[i]][, 13]),
                                                     12:ncol(locus[[i]])]);
            colnames(locusInTx) <- colnames(locus[[i]])[1:ncol(locusInTx)];
        } else if (method == "permutation") {
            locusInTx <- locus[[i]];
            locusInTx$TXSTART <- apply(locusInTx, 1, function(x)
                                       round(runif(1, min = 1,
                                                   max = as.numeric(x["REGION_TXWIDTH"]) - 1)));
            locusInTx$TXEND <- locusInTx$TXSTART;
            locusInTx[c("CHR", "START", "STOP",
                        "WIDTH", "STRAND", "SCORE")] <- matrix("*",
                                                               nrow = nrow(locusInTx),
                                                               ncol = 6);
            locusInTx$ID <- "unif_perm";
        }
        locusInTx.list[[length(locusInTx.list) + 1]] <- locusInTx;
    }
    names(locusInTx.list) <- names(locus);
    obj <- new("txLoc",
               loci = locusInTx.list,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()));
    return(obj);
}

