#' Read BED-formatted file.
#'
#' Read BED-formatted file. See 'Details'.
#'
#' The function opens and reads in a BED6-formatted file, and 
#' stores the annotation features in a \code{GRanges} object.
#' Chromosome names must follow either UCSC (chr1, ..., chrM, chrX)
#' or Ensembl naming conventions (1, ..., MT, X).
#'
#' @param file A character string; specifies the input BED file.
#' @param collapseRange A logical scalar; should \code{TRUE} loci
#' spanning more than one nucleotide be collapsed to a single
#' nucleotide locus corresponding to the midpoint of the range;
#' default is \code{TRUE}.
#'
#' @return A \code{GRanges} object. See 'Details'.
#' 
#' @examples
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @import GenomicRanges IRanges
#' @importFrom utils read.table
#' 
#' @export
ReadBED <- function(file, collapseRange = TRUE) {

    # Sanity checks
    if (!file.exists(file)) {
        stop(sprintf("File %s doesn't exist.", file))
    }
    bed <- read.table(file, header = FALSE)
    if (ncol(bed) != 6) {
        stop("[ERROR] File must be a 6 column BED file.")
    }
    colnames(bed) <- c("chr", "start", "end", "id", "score", "strand")
    
    # Use UCSC chromosome names
    if (all(!grepl("chr", bed$chr, ignore.case = TRUE))) {
        bed$chr <- sprintf("chr%s", bed$chr)
    }
    bed$chr <- gsub("chrMT", "chrM", bed$chr)
    
    # Order BED file by chromosome then starting point
    bed <- bed[order(bed$chr, bed$start), ]
    
    # Replace strand = "." with strand = "*"
    if (length(grep("\\.", bed$strand)) > 0) {
        bed$strand <- gsub("\\.", "*", bed$strand)
        warning(sprintf(
            "BED file %s contains entries without strand information.",
            file))
    }
    
    # Collapse ranges
    if (collapseRange == TRUE) {
        bed$start <- floor(0.5 * (bed$start + bed$end))
        bed$end <- bed$start + 1
    }
    
    # Convert to GRanges
    GRanges(
        seqnames = bed$chr,
        ranges = IRanges(bed$start + 1, bed$end),
        strand = bed$strand,
        score = bed$score,
        id = bed$id)
    
}


#' Write a \code{list} of \code{GRangesList} to a BED file.
#'
#' Write a \code{list} of \code{GRangesList} to a BED file.
#' One file per transcript region is generated.
#'
#' @param txFeatures A \code{list} of \code{GRangesList} objects.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
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
        rtracklayer::export(txFeatures[[i]], sprintf("%s.bed", names(txFeatures)[i]))
    }
}


#' Write \code{txLoc} object to a BED file.
#'
#' Write \code{txLoc} object to a BED file. See 'Details'.
#'
#' The function writes entries from a \code{txLoc} object to a 6-column
#' BED file (BED6). Note that this process is not "splice-aware", i.e. 
#' if an entry spans an intron the BED entry gives the left and right-most 
#' genomic coordinate of the feature. If \code{file = NULL}, entries will 
#' be written to  \code{sites.bed}. If \code{noChrName = TRUE}, chromosome 
#' names in column 1 of the BED file will be written without "chr".
#' 
#' @param locus A \code{txLoc} object.
#' @param file A character string; specifies the filename of the output 
#' BED file. If \code{NULL}, then \code{file = "sites.bed"}; default is
#' \code{NULL}.
#' @param noChrName A logical scalar; if \code{TRUE}, chromosome names
#' will be written without "chr"; default is \code{FALSE}. 
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' txSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' WriteTxLocToBED(txSites);
#' }
#' 
#' @importFrom utils write.table
#' 
#' @export
WriteTxLocToBED <- function(locus,
                            file = NULL,
                            noChrName = FALSE) {
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
    CheckClass(locus, "txLoc")
    id <- GetId(locus)
    locus <- GetLoci(locus)
    BED <- vector()
    for (i in 1:length(locus)) {
        feat <- locus[[i]][, c("CHR", "START", "STOP", "ID", "SCORE", "STRAND")]
        if (is.numeric(feat[, 2]) && is.numeric(feat[, 3])) {
            if (noChrName == TRUE) {
                feat[, 1] <- gsub("chr", "", feat[ ,1])
            }
            feat[, 2] <- feat[ ,2] - 1
            feat[, 4] <- sprintf("%s|%s|%s",
                                 feat[, 4],
                                 locus[[i]][, c("GENE_ENSEMBL")],
                                 names(locus)[i])
            BED <- rbind(BED, feat)
        } else {
            ss <- sprintf("Skipping %s.\n", names(locus)[i])
            warning(ss)
        }
    }
    BED <- BED[order(BED[, 1], BED[, 2]), ]
    if (is.null(file)) {
        if (nchar(id) > 0) {
            file <- sprintf("sites_%s.bed", id)
        } else {
            file <- "sites.bed"
        }
    }
    write.table(BED, file = file,
                sep = "\t", quote = FALSE,
                col.names = FALSE, row.names = FALSE)
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
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom utils write.csv
#' 
#' @export
WriteTxLocToCSV <- function(locus,
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
    CheckClass(locus, "txLoc")
    locus <- GetLoci(locus)
    CSV <- vector()
    selCol <- c("GENE_REGION", "REGION_TXWIDTH", "GENE_REFSEQ",
                "GENE_ENTREZ", "GENE_SYMBOL", "GENE_ENSEMBL",
                "GENE_UNIGENE", "ID", "CHR", "START", "STOP", "STRAND")
    if (withSeq) {
        selCol <- c(selCol, "REGION_SEQ")
    }
    for (i in 1:length(locus)) {
        data <- locus[[i]][, selCol]
        CSV <- rbind(CSV, data)
    }
    CSV <- CSV[order(CSV[, 8], CSV[, 9]), ]
    if (is.null(file)) {
        file <- "sites.csv"
    }
    write.csv(CSV, file = file,
              row.names = FALSE, quote = FALSE)
}


#' Read a DBN file.
#'
#' Read a DBN file. See 'Details'.
#'
#' The function reads in a DBN (dot-bracket) structure file, and
#' returns a \code{dataframe} with the following data columns:
#' \enumerate{
#' \item Column 1: Sequence ID
#' \item Column 2: Length of the sequence (in nt)
#' \item Column 3: Mean free energy (MFE)
#' }
#' 
#' @param file A character string; specifies the input DBN file.
#'
#' @return A \code{dataframe} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @export
ReadDBN <- function(file) {
    if (!file.exists(file)) {
        ss <- sprintf("Could not open %s.", file)
        stop(ss)
    }
    d <- as.data.frame(matrix(readLines(file), ncol = 3, byrow = TRUE),
                       stringsAsFactors = FALSE)
    d[, 1] <- gsub("^>", "", d[, 1])
    d[, 2] <- nchar(d[, 2])
    d[, 3] <- gsub("^[\\(\\.\\)]+\\s", "", d[, 3])
    d[ ,3] <- as.numeric(gsub("[\\(\\)]", "", d[, 3]))
    colnames(d) <- c("id", "siteSeqLength", "siteMFE")
    return(d)
}
