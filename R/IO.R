#' Read BED-formatted file.
#'
#' Open and read a BED file, and store the annotation
#' features in a \code{GRanges} object.
#'
#' @param file A character string; specifies the input BED file.
#'
#' @return A \code{GRanges} object.
#' 
#' @examples
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' 
#' @export
ReadBED <- function(file) {
    # Read BED file and convert to GRanges object.
    #
    # Args:
    #   file: Input BED file
    #
    # Returns:
    #   GRanges object
    if (!file.exists(file)) {
        stop(sprintf("File %s doesn't exist.", file));
    }
    bed <- read.table(file, header = FALSE);
    if (ncol(bed) < 6) {
        stop("Need at least 6 columns in BED file.");
    }
    if (all(!grepl("chr", bed[, 1]))) {
        bed[, 1] <- sprintf("chr%s", bed[, 1]);
    }
    bed[, 1] <- gsub("chrMT", "chrM", bed[, 1]);
    bed <- bed[order(bed[, 1], bed[, 2]), ];
    if (length(grep("\\.", bed[, 6])) > 0) {
        bed[, 6] <- gsub("\\.", "*", bed[, 6]);
        warning(sprintf(
            "BED file %s contains entries without strand information.",
            file));
    }
    colnames(bed) <- c("chr", "start", "end", "id", "score", "strand");
    gr <- GRanges(bed$chr, IRanges(bed$start+1,bed$end), bed$strand,
                  score = bed$score, id = bed$id);
    return(gr);
}


#' Write a \code{list} of \code{GRangesList} to a BED file.
#'
#' Write a \code{list} of \code{GRangesList} to a BED file.
#' One file per transcript region is generated.
#'
#' @param txFeatures A \code{list} of \code{GRangesList} objects.
#'
#' @importFrom rtracklayer export
#' 
#' @export
WriteFeatToBED <- function(txFeatures) {
    # Write list of GRangesList transcript features to BED file.
    #
    # Args:
    #   txFeatures: List of GRangesList transcript features.
    #
    # Returns:
    #   NULL
    for (i in 1:length(txFeatures)) {
        rtracklayer::export(txFeatures[[i]], sprintf("%s.bed", names(txFeatures)[i]));
    }
}


#' Write txLoc object to BED file.
#'
#' Write a txLoc object to a/multiple BED file(s).
#' 
#' @param locus A \code{txLoc} object.
#' @param file Filename of output BED file. If NULL then file = "sites.bed".
#' @param formatChr Set output format of chromosome column. Default is
#' "noChrName".
#'
#' @export
WriteLocusToBED <- function(locus,
                            file = NULL,
                            formatChr = "noChrName") {
    # Save transcript features in BED file.
    #
    # Args:
    #   locus: List of dataframes with mapped features aross different
    #          transcript regions
    #   file: Filename. If file is NULL then output is written to site.bed.
    #         Default is NULL.
    #   formatChr: Format of chr column in BED file. Default is "noChrName".
    #
    # Returns:
    #    NULL
    CheckClass(locus, "txLoc");
    id <- GetId(locus);
    locus <- GetLoci(locus);
    BED <- vector();
    for (i in 1:length(locus)) {
        feat <- locus[[i]][, c("CHR", "START", "STOP", "ID", "SCORE", "STRAND")];
        if (is.numeric(feat[, 2]) && is.numeric(feat[, 3])) {
            if (formatChr == "noChrName") {
                feat[, 1] <- gsub("chr", "", feat[ ,1]);
            }
            feat[, 2] <- feat[ ,2] - 1;
            feat[, 4] <- sprintf("%s|%s|%s",
                                 feat[, 4],
                                 locus[[i]][, c("GENE_ENSEMBL")],
                                 names(locus)[i]);
            BED <- rbind(BED, feat);
        } else {
            ss <- sprintf("Skipping %s.\n", names(locus)[i]);
            warning(ss);
        }
        
    }
    BED <- BED[order(BED[, 1], BED[, 2]), ];
    if (is.null(file)) {
        if (nchar(id) > 0) {
            file <- sprintf("sites_%s.bed", id);
        } else {
            file <- "sites.bed";
            
        }
    }
    write.table(BED, file = file,
                sep = "\t", quote = FALSE,
                col.names = FALSE, row.names = FALSE);
}


#' Write txLoc object to CSV file.
#'
#' Write a txLoc object to a/multiple comma-separated-values file(s).
#' 
#' @param locus A \code{txLoc} object.
#' @param file Filename of output CSV file. If NULL then file = "sites.csv".
#' @param withSeq If TRUE then include full sequence. Default is FALSE.
#' @param withGC If TRUE then include GC content. Default is FALSE.
#'
#' @export
WriteLocusToCSV <- function(locus,
                            file = NULL, 
                            withSeq = FALSE, 
                            withGC = FALSE) {
    # Save list of dataframes with mapped features across different
    # transcript regions to csv file
    #
    # Args:
    #   locus: List of dataframes with mapped features across different
    #          transcript regions
    #   file: Filename. If file is NULL then output is written to site.csv.
    #         Default is NULL.
    #
    # Returns:
    #   NULL
    CheckClass(locus, "txLoc");
    locus <- GetLoci(locus);
    CSV <- vector();
    selCol <- c("GENE_REGION", "GENE_REFSEQ", "GENE_ENTREZ",
                "GENE_SYMBOL", "GENE_ENSEMBL", "GENE_UNIGENE",
                "ID", "CHR", "START", "STOP", "STRAND");
    if (withSeq) {
        selCol <- c(selCol, "REGION_SEQ");
    }
    for (i in 1:length(locus)) {
        data <- locus[[i]][, selCol];
        CSV <- rbind(CSV, data);
    }
    CSV <- CSV[order(CSV[, 8], CSV[, 9]), ];
    if (is.null(file)) {
        file <- "sites.csv";
    }
    write.csv(CSV, file = file,
              row.names = FALSE, quote = FALSE);
}
