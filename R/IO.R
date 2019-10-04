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
    bed <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
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
        score = as.numeric(bed$score),
        id = bed$id)
    
}


#' Write \code{txLoc} object to a BED file.
#'
#' Write \code{txLoc} object to a BED file. See 'Details'.
#'
#' The function writes entries from a \code{txLoc} object to a 6-column BED 
#' file (BED6). Note that this process is not "splice-aware", i.e. if an entry 
#' spans an intron the BED entry gives the left and right-most genomic 
#' coordinate of the feature. If \code{file = NULL}, entries will be written to  
#' \code{sites.bed}.
#' 
#' @param txLoc A \code{txLoc} object.
#' @param file A character string; specifies the filename of the output 
#' BED file. If \code{NULL}, then \code{file = "sites.bed"}; default is
#' \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' txSites <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' WriteBED(txSites)
#' }
#' 
#' @importFrom utils write.table
#' 
#' @export
WriteBED <- function(txLoc, file = NULL) {

    # Sanity checks
    CheckClass(txLoc, "txLoc")
    
    # Quieten R CMD check concerns regarding "no visible binding for global
    # variable ..."
    locus_in_genome <- id <- tx_region <- tx_refseq <- NULL

    # Get loci
    loci <- GetLoci(txLoc)
    
    # Extract BED6 entries from txLoc object and row-bind into `data.frame`
    lst <- lapply(loci, function(locus) 
        with(locus, data.frame(
            seqnames = seqnames(locus_in_genome),
            start = start(locus_in_genome) - 1L,
            end = end(locus_in_genome),
            id = paste(GetId(txLoc), id, tx_region, tx_refseq, sep = "|"),
            score = score,
            strand = strand(locus_in_genome))))
    df <- do.call(rbind, lst)
    
    # Sort entries
    df <- df[order(df[, 1], df[, 2]), ]
        
    # Output file name
    if (is.null(file)) {
        if (nchar(GetId(txLoc)) > 0) {
            file <- sprintf("sites_%s.bed", GetId(txLoc))
        } else {
            file <- "sites.bed"
        }
    }

    # Write to file
    write.table(
        df, 
        file = file,
        sep = "\t", 
        quote = FALSE,
        col.names = FALSE, 
        row.names = FALSE)
    cat(sprintf("Output in file %s.\n", file))

}


#' Write txLoc object to CSV file.
#'
#' Write a \code{txLoc} object to a comma-separated-values file.
#' 
#' @param txLoc A \code{txLoc} object.
#' @param file Filename of output CSV file. If NULL then file = "sites.csv".
#' @param withSeq If TRUE then include full sequence. Default is FALSE.
#' @param withGC If TRUE then include GC content. Default is FALSE. Unused.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @importFrom utils write.csv
#' 
#' @export
WriteCSV <- function(txLoc,
                     file = NULL, 
                     withSeq = FALSE, 
                     withGC = FALSE) {

    # Sanity checks
    CheckClass(txLoc, "txLoc")

    # Get loci
    loci <- GetLoci(txLoc)

    # Convert entries in columns to `data.frame`
    lst <- lapply(loci, function(locus) {
        if (!withSeq) 
            as.data.frame(locus[, -match("tx_region_sequence", names(locus))])
        else
            as.data.frame(locus)
    })
    df <- do.call(rbind, lst)    
   
    # Output file name
    if (is.null(file)) {
        if (nchar(GetId(txLoc)) > 0) {
            file <- sprintf("sites_%s.csv", GetId(txLoc))
        } else {
            file <- "sites.csv"
        }
    }
    
    # Write to file
    write.csv(
        df, 
        file = file,
        row.names = FALSE, 
        quote = FALSE)
    cat(sprintf("Output in file %s.\n", file))
    
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
