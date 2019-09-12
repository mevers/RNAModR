#' Check package dependencies.
#'
#' Check package dependencies. See 'Details'.
#'
#' This is a low-level function that is being called from
#' \code{BuildTx}. The function checks that the R packages
#' providing the relevant gene annotations and reference genome
#' for reference genome \code{genomeVersion} are available
#' in the current R environment. The function returns \code{TRUE}
#' if the necessary R packages are installed.
#'
#' @param genomeVersion A character string; refers to a specific
#' reference genome assembly version; default is \code{"hg38"}.
#'
#' @return A logical scalar. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom utils installed.packages
#'
#' @examples
#' \dontrun{
#' CheckPkgDependencies("hg38");
#' }
#'
#' @export
CheckPkgDependencies <- function(genomeVersion = "hg38") {
    pkgsInstalled <- installed.packages()
    pkgsReq <- list()
    
    #pkgsReq <- data.frame(
    #    genomeVersion = c(
    #        "hg38", "hg19", "hg18",
    #        "mm10", "mm9", "mm8",
    #        "dm6", "dm3", "dm2",
    #        "sacCer3", "sacCer2", "sacCer1"),
    #    annotation = c(
    #        rep("org.Hs.eg.db", 3),
    #        rep("org.Mm.eg.db", 3),
    #        rep("org.Dm.eg.db", 3),
    #        rep("org.Sc.sgd.db", 3)),
    #    genome = c(
    #        "BSgenome.Hsapiens.UCSC.hg38", 
    #        "BSgenome.Hsapiens.UCSC.hg19",
    #        "BSgenome.Hsapiens.UCSC.hg18",
    #        "BSgenome.Mmusculus.UCSC.mm10",
    #        "BSgenome.Mmusculus.UCSC.mm9",
    #        "BSgenome.Mmusculus.UCSC.mm8",
    #        "BSgenome.Dmelanogaster.UCSC.dm6",
    #        "BSgenome.Dmelanogaster.UCSC.dm3",
    #        "BSgenome.Dmelanogaster.UCSC.dm2",
    #        "BSgenome.Scerevisiae.UCSC.sacCer3",
    #        "BSgenome.Scerevisiae.UCSC.sacCer2",
    #        "BSgenome.Scerevisiae.UCSC.sacCer1"))
    #pkgsReq[match(genomeVersion, pkgsReq$genomeVersion), ]
        
    if (genomeVersion == "hg38") {
        pkgsReq[["Annot"]] <- "org.Hs.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Hsapiens.UCSC.hg38"
    } else if (genomeVersion == "hg19") {
        pkgsReq[["Annot"]] <- "org.Hs.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Hsapiens.UCSC.hg19"
    } else if (genomeVersion == "hg18") {
        pkgsReq[["Annot"]] <- "org.Hs.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Hsapiens.UCSC.hg18"
    } else if (genomeVersion == "mm10") {
        pkgsReq[["Annot"]] <- "org.Mm.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Mmusculus.UCSC.mm10"
    } else if (genomeVersion == "mm9") {
        pkgsReq[["Annot"]] <- "org.Mm.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Mmusculus.UCSC.mm9"
    } else if (genomeVersion == "mm8") {
        pkgsReq[["Annot"]] <- "org.Mm.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Mmusculus.UCSC.mm8"
    } else if (genomeVersion == "dm6") {
        pkgsReq[["Annot"]] <- "org.Dm.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Dmelanogaster.UCSC.dm6"
    } else if (genomeVersion == "dm3") {
        pkgsReq[["Annot"]] <- "org.Dm.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Dmelanogaster.UCSC.dm3"
    } else if (genomeVersion == "dm2") {
        pkgsReq[["Annot"]] <- "org.Dm.eg.db"
        pkgsReq[["Genome"]] <- "BSgenome.Dmelanogaster.UCSC.dm2"
    } else if (genomeVersion == "sacCer3") {
        pkgsReq[["Annot"]] <- "org.Sc.sgd.db"
        pkgsReq[["Genome"]] <- "BSgenome.Scerevisiae.UCSC.sacCer3"
    } else if (genomeVersion == "sacCer2") {
        pkgsReq[["Annot"]] <- "org.Sc.sgd.db"
        pkgsReq[["Genome"]] <- "BSgenome.Scerevisiae.UCSC.sacCer2"
    } else if (genomeVersion == "sacCer1") {
        pkgsReq[["Annot"]] <- "org.Sc.sgd.db"
        pkgsReq[["Genome"]] <- "BSgenome.Scerevisiae.UCSC.sacCer1"
    }
    ret <- sapply(pkgsReq, function(x) x %in% pkgsInstalled[, "Package"])
    df <- cbind(t(as.data.frame(pkgsReq)), ret)
    if (!all(ret)) {
        # Create meaningful error message
        ss <- "[ERROR] R/Bioconductor package dependency not met."
        pkgsMissing <- df[which(df[, 2] == FALSE), 1]
        for (i in 1:length(pkgsMissing)) {
            ss <- sprintf("%s\n  Not found: %s", ss, pkgsMissing[i])
        }
        ss <- sprintf("%s\n\n  Running", ss)
        biocSrc <- "source(\"http://www.bioconductor.org/biocLite.R\")"
        ss <- sprintf("%s\n   %s", ss, biocSrc)
        for (i in 1:length(pkgsMissing)) {
            ss <- sprintf("%s\n   biocLite(\"%s\")", ss, pkgsMissing[i])
        }
        ss <- sprintf("%s\n  might fix that.", ss)
        stop(ss)
    }
}


#' Download a new or load an existing reference transcriptome.
#'
#' Download a new or load an existing reference transcriptome.
#' See 'Details'.
#'
#' This is a low-level function that is being called from
#' \code{BuildTx}. If no local transcriptome db file exists or
#' \code{forceDownload} is \code{TRUE}, the function generates an
#' reference transcriptome by calling \code{makeTxDbFromUCSC} from
#' the \code{GenomicFeatures} package; the return object is a
#' \code{TxDb} object and is saved locally as an sqlite db.
#' Details on the \code{txdb} return object can be obtained from
#' \code{print(x)}, where \code{x} is the name of the \code{txdb}
#' return object.
#' The user may generate a custom transcriptome (e.g. based on a
#' different gene annotation) using \code{makeTxDb} or
#' \code{makeTxDb}-derived functions from the \code{GenomicFeatures}
#' package. Please see the \code{GenomicFeatures} reference manual.
#' In order to use a custom transcriptome, the \code{TxDb} object
#' needs to be saved locally as an sqlite db, e.g.
#' \code{saveDb(x), "txdb_y.sqlite")}, where \code{x} is the
#' \code{TxDb} object and \code{y} denotes the reference genome
#' version as in \code{genomeVersion}.
#'
#' @param genomeVersion A character string; refers to a specific
#' reference genome assembly version; default is \code{"hg38"}.
#' @param standardChrOnly A logical scalar; if \code{TRUE} and a
#' new db is created, keep only fully assembled chromosomes;
#' default is \code{TRUE}.
#' @param force A logical scalar; if \code{TRUE} force download
#' and overwrite existing db file; default is \code{FALSE}.
#' @param verbose A logical scalar; print additional information
#' about db; default is \code{FALSE}.
#'
#' @return A \code{TxDb} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' txdb <- GetTxDb();
#' print(txdb);
#' }
#'
#' @export
GetTxDb <- function(genomeVersion = "hg38",
                    standardChrOnly = TRUE,
                    force = FALSE,
                    verbose = FALSE) {
    # Download or load (if file exists) sqlite transcript database based on
    # RefSeq annotation from UCSC. Store database as sqlite file.
    #
    # Args:
    #   genomeVersion: Genome assembly version. Default is "hg38".
    #   standardChrOnly: As it says. Default is TRUE.
    #   force: Force re-downloading transcript database from UCSC.
    #          Default is FALSE.
    #   verbose: Print additional output. Default is FALSE.
    #
    # Returns:
    #   TxDb object.
    sqliteFile <- sprintf("txdb_%s.sqlite", genomeVersion)
    # RefSeq genes for human, mouse and fruitfly data
    # SGP genes for yeast
    geneRef <- ifelse(grepl("sacCer", genomeVersion),
                      "sgdGene",
                      "refGene")
    if ((!file.exists(sqliteFile)) || (force == TRUE)) {
        txdb <- makeTxDbFromUCSC(genome = genomeVersion,
                                 tablename = geneRef)
        saveDb(txdb, file = sqliteFile)
    } else {
        cat(sprintf("Found existing sqlite database %s.\n",
                    sqliteFile))
        txdb <- loadDb(file = sqliteFile)
    }
    if (standardChrOnly) {
        txdb <- keepStandardChromosomes(txdb)
    }
    if (verbose) {
        cat("List of chromosomes:\n")
        print(seqlevels(txdb))
        cat("List of columns:\n")
        print(columns(txdb))
    }
    return(txdb)
}


#' Get different gene IDs.
#'
#' Get different gene IDs for mapping between different gene
#' annotations (Ensembl, Entrez, etc.). See 'Details'.
#' This is a low-level function that is being called from
#' \code{BuildTx}.
#'
#' The function extracts the genome assembly version from
#' the \code{TxDb} object, and loads the suitable genome wide
#' annotation package. For example, if \code{TxDb} is based
#' on \code{"hg38"}, the function loads \code{org.Hs.eg.db}.
#' Various IDs (Ensembl, Unigene, gene symbols) are extracted
#' from the annotation package by mapping RefSeq IDs from
#' \code{TxDb} via \code{mapIds} from the \code{AnnotationDbi}
#' package. Please note that column entries are not necessarily
#' unique: For example, two different Ensembl IDs might have the
#' same Entrez ID. The return object is a \code{data.frame} with
#' the following columns:
#' \itemize{
#'   \item REFSEQ: RefSeq Accession number; these entries are unique
#'   and refer to the mRNA transcript
#'   \item ENTREZID: Entrez gene ID
#'   \item SYMBOL: Gene symbol
#'   \item ENSEMBL: Ensembl gene ID
#'   \item GENENAME: Gene name
#' }
#' Note that mapping between gene IDs is a many-to-many mapping process.
#' For example, different transcript RefSeq IDs can belong to the same
#' gene. Which transcripts are associated with which gene depends in
#' turn on the gene reference system (RefSeq, Ensembl, etc.). On top of
#' that, the same transcript may have >1 location in the reference
#' genome: NM_000451 (EntrezID 6473) is annotated on chromosomes X and Y,
#' NM_001001418 (EntrezID 414060) has two loci on chr17.
#'
#' @param txdb A \code{TxDb} object.
#'
#' @return A \code{data.frame}. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
GetGeneIds <- function(txdb) {
    # Get different cross-referencing gene IDs based on sqlite
    # transcript database.
    #
    # Args:
    #   txdb: A TxDb object.
    #
    # Returns:
    #   Dataframe of different gene IDs.
    genomeVersion <- unique(genome(txdb))
    if (length(genomeVersion)==0) {
        stop("Could not identify genome version from transcript database.")
    }
    # Load gene annotation package
    if (grepl("(hg38|hg19|hg18)", genomeVersion,
              ignore.case = TRUE)) {
        if (SafeLoad("org.Hs.eg.db")) {
            refDb <- get("org.Hs.eg.db")
        } else {
            stop(sprintf("[ERROR] Could not load gene annotation for %s.",
                         genomeVersion))
        }
    } else if (grepl("(mm10|mm9|mm8)", genomeVersion,
                     ignore.case = TRUE)) {
        if (SafeLoad("org.Mm.eg.db")) {
            refDb <- get("org.Mm.eg.db")
        } else {
            stop(sprintf("[ERROR] Could not load gene annotation for %s.",
                         genomeVersion))
        }
    } else if (grepl("(dm6|dm3|dm2)", genomeVersion,
                     ignore.case = TRUE)) {
        if (SafeLoad("org.Dm.eg.db")) {
            refDb <- get("org.Dm.eg.db")
        } else {
            stop(sprintf("[ERROR] Could not load gene annotation for %s.",
                         genomeVersion))
        }
    } else if (grepl("(sacCer3|sacCer2|sacCer1)", genomeVersion,
                     ignore.case = TRUE)) {
        if (SafeLoad("org.Sc.sgd.db")) {
            refDb <- get("org.Sc.sgd.db")
        } else {
            stop(sprintf("[ERROR] Could not load gene annotation for %s.",
                         genomeVersion))
        }
    } else {
        stop(sprintf("Unknown genome %s.", genomeVersion))
    }
    
    # Get transcripts from TxDb database; RefSeq transcript ID and Entrez gene ID
    # as metadata columns
    # Store REFSEQ and ENTREZID in DataFrame geneXID
    tx <- transcripts(txdb, c("TXNAME", "GENEID"))
    geneXID <- setNames(mcols(tx), c("REFSEQ", "ENTREZID"))
    geneXID$ENTREZID <- as.character(geneXID$ENTREZID)

    # Get various gene identifiers from OrgDb database and merge with
    # geneXID DataFrame
    # We suppress messages indicating matching of duplicate keys
    # That's ok, because duplicate keys will lead to duplicate
    # entries
    identifier <- c("ENTREZID", "SYMBOL", "ENSEMBL", "GENENAME")
    geneXID <- merge(
        geneXID, 
        suppressMessages(select(
            refDb, 
            keys = geneXID$ENTREZID, 
            columns = identifier, 
            keytype = "ENTREZID")),
        by = "ENTREZID",
        sort = FALSE)[c(2:1,3:5)]
    
    return(geneXID[!duplicated(geneXID), ])
}


#' Get transcript sections.
#'
#' Get transcript sections. See 'Details'.
#' This is a low-level function that is being called from
#' \code{BuildTx}.
#'
#' The function returns a \code{list} of \code{GRangesList}
#' objects. The number and names of list entries match the
#' number and names of \code{sections}. Entries in every
#' \code{GRangesList} are labelled according to their
#' transcript RefSeq ID. Note that multiple entries with
#' identical RefSeq ID exist.
#'
#' @param txdb A \code{TxDb} object.
#' @param sections A character vector; specifies which transcript
#' sections should be extracted from \code{txdb}; default is
#' \code{c("Promoter", "5'UTR", "CDS", "3'UTR", "Intron")}.
#' @param promUpstream An integer; specifies the number of
#' nucleotides to be included upstream of 5'UTR start as part of
#' the promoter.
#' @param promDownstream An integer; specifies the number of
#' nucleotides to be included downstream of 5'UTR start as part
#' of the promoter.
#' @param verbose A logical scalar; if \code{TRUE}, print additional
#' information; default is \code{FALSE}
#'
#' @return A \code{list} of \code{GRangesList} objects
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
GetTxBySec <- function(txdb,
                       sections = c(
                           "Promoter",
                           "5'UTR",
                           "CDS",
                           "3'UTR",
                           "Intron"),
                       promUpstream = 1000,
                       promDownstream = 0,
                       verbose = FALSE) {
    # Extract sections (UTR's, CDS etc.) from transcript database.
    #
    # Args:
    #   txdb: A TxdB object.
    #   features: Vector of transcript features to be extracted (in that order).
    #             Default is c("Promoter","5'UTR","CDS","3'UTR").
    #   promUpstream: Number of upstream nt's to be included in promoter.
    #   promDownstream: Number of downstream nt's to be included in promoter.
    #   verbose: Print additional information
    #
    # Returns:
    #   List of GRangesList transcript features.
    CheckClass(txdb, "TxDb")
    validFeat <- c("Promoter",
                   "CDS", "coding",
                   "5pUTR", "5'UTR", "UTR5",
                   "3pUTR", "3'UTR", "UTR3",
                   "Intron", "Intronic")
    isValidFeat <- grepl(sprintf("(%s)", paste0(validFeat, collapse = "|")),
                         sections, ignore.case = TRUE)
    if (!all(isValidFeat)) {
        stop(sprintf("%s is not a valid feature.", sections[!isValidFeat]))
    }
    txBySec <- list()
    for (i in 1:length(sections)) {
        # Suppress warning messages that arise from duplicate entries.
        # That's ok because we will deal with duplicates later
        if (grepl("(CDS|coding)", sections[i], ignore.case = TRUE)) {
            sec <- suppressWarnings(
                cdsBy(txdb, by="tx", use.names = TRUE))
        } else if (grepl("(5pUTR|5'UTR|UTR5)", sections[i], ignore.case = TRUE)) {
            sec <- suppressWarnings(
                fiveUTRsByTranscript(txdb, use.names = TRUE))
        } else if (grepl("(3pUTR|3'UTR|UTR3)", sections[i], ignore.case = TRUE)) {
            sec <- suppressWarnings(
                threeUTRsByTranscript(txdb, use.names = TRUE))
        } else if (grepl("(Intron|Intronic)", sections[i], ignore.case = TRUE)) {
            sec <- suppressWarnings(
                intronsByTranscript(txdb, use.names = TRUE))
        } else if (grepl("Promoter", sections[i], ignore.case = TRUE)) {
            promoter <- suppressWarnings(
                unname(promoters(txdb,
                          columns = c("tx_id", "tx_name"),
                          upstream = promUpstream,
                          downstream = promDownstream)))
            # Make sure we are still within valid transcript ranges
            promoter <- trim(promoter)
            sec <- GenomicRanges::split(promoter, promoter$tx_name)
        }
        txBySec[[length(txBySec) + 1]] <- sec
        if (verbose) {
            print(sprintf("%i entries for section %s.",
                          length(sec), sections[i]))
        }
    }
    names(txBySec) <- sections
    return(txBySec)
}


#' De-duplicate \code{GRangesList} entries based on distance.
#'
#' De-duplicate \code{GRangesList} entries based on distance.
#' See 'Details'.
#' This is a low-level function that is being called from
#' \code{CollapseTxBySec}.

#'
#' For two \code{GRangesList} objects select entries from
#' \code{query} that have an entry in \code{ref} and collapse
#' multiple entries with the same name from \code{query} based
#' on the minimum distance to the single entry from \code{ref}.
#'
#' @param query A \code{GRangesList} object; the query features
#' @param ref A \code{GRangesList} object; the reference features
#' @param showPb A logical scalar; if \code{TRUE} show a progress
#' bar.
#'
#' @return A \code{GRangesList} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#'
#' @export
DedupeBasedOnNearestRef <- function(query, ref, showPb = FALSE) {
    query <- query[which(elementNROWS(query) > 0)]
    query <- query[which(names(query) %in% names(ref))]
    query <- query[order(names(query), -sum(width(query)))]
    dupes <- which(duplicated(names(query)))
    if (length(dupes) > 0) {
        dupeID <- unique(names(query[dupes]))
        rem <- vector()
        if (showPb) {
            pb <- txtProgressBar(min = 0,
                                 max = length(dupeID),
                                 width = 80,
                                 style = 3)
        }
        for (i in 1:length(dupeID)) {
            idxQuery <- which(names(query) == dupeID[i])
            idxRef <- which(names(ref) == dupeID[i])
            if (length(idxQuery) > 1) {
                dist <- distance(
                    unlist(range(query[idxQuery])),
                    range(ref[[idxRef]]))
                idxMinDist <- which.min(dist)
                rem <- c(rem, idxQuery[-idxMinDist])
                if (showPb) {
                    setTxtProgressBar(pb, i)
                }
            }
        }
        if (showPb) {
            close(pb)
        }
        return(query[-rem])
    } else {
        return(query)
    }
}


#' Collapse entries from a \code{GetTxBySec} return object.
#'
#' Collapse entries from a \code{GetTxBySec} return object.
#' See 'Details'.
#' This is a low-level function that is being called from
#' \code{BuildTx}.
#'
#' Entries are collapsed using the following method:
#' \itemize{
#'   \item If multiple entries exist for one Entrez ID, the
#' transcript with the longest CDS is kept.
#'   \item For every CDS with a unique RefSeq ID, the
#' 5'/3'UTRs closest to the CDS (non-overlapping) is kept; in
#' the case of a distance tie, the longest UTR is kept.
#'   \item For every 5'UTR with a unique RefSeq ID, the
#' corresponding upstream promoter is chosen.
#' }
#'
#' @param txBySec A \code{list} of \code{GRangesList} objects;
#' output of function \code{GetTxBySec}.
#' @param geneXID A \code{data.frame}; output of function
#' \code{GetGeneIds} providing gene IDs from different gene
#' annotations.
#' @param verbose A logical scalar; if \code{TRUE} print
#' additional output; default it \code{FALSE}
#'
#' @return A \code{list} of \code{GRangesList} objects
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
CollapseTxBySec <- function(txBySec,
                            geneXID,
                            verbose = FALSE) {
    # Collapse list of transcripts split by sections to give one unique RefSeq
    # ID per transcript section per gene; employ the following method
    #  (1) Choose transcript with longest CDS for transcripts with multiple CDS
    #  (2) Choose 5'UTR/3'UTR closest to CDS (non-overlapping); choose the
    #      longest UTR if there is a distance tie
    #  (3) If promoters are given, choose promoter corresponding to selected
    #      5'UTR
    #
    # Args:
    #   txBySec: List of GRangesList transcript sections.
    #   geneXID: Gene annotation dataframe
    #   verbose: Print additional output. Default is FALSE.
    #
    # Returns:
    #   List of GRangesList transcript sections.
    # Error handling
    sections <- names(txBySec)
    whichCDS <- grep("(CDS|coding)", sections, ignore.case = TRUE)
    if (!length(whichCDS)) {
        stop("Transcript feature list does not contain CDS entries.")
    }
    whichUTR5 <- grep("(5pUTR|5'UTR|UTR5)", sections, ignore.case = TRUE)
    if (!length(whichUTR5)) {
        stop("Transcript feature list does not contain 5'UTR entries.")
    }
    whichUTR3 <- grep("(3pUTR|3'UTR|UTR3)", sections, ignore.case = TRUE)
    if (!length(whichUTR3)) {
        stop("Transcript feature list does not contain 3'UTR entries.")
    }
    whichPromoter <- grep("Promoter", sections, ignore.case = TRUE)
    whichIntron <- grep("(Intron|Intronic)", sections, ignore.case = TRUE)
    # (1) Collapse CDS:
    #      Step 1: Same RefSeq ID, multiple identical genomic copies
    #       Example: RefSeq ID = NM_000513
    #      Step 2: Same Entrez ID, multiple CDS isoforms
    #       Example: Entrez ID = 10001
    cds <- txBySec[[whichCDS]]
    cds <- cds[order(names(cds), -sum(GenomicRanges::width(cds)))]
    cds <- cds[!duplicated(names(cds))]
    geneData <- geneXID[which(geneXID$REFSEQ %in% names(cds)), ]
    geneData$lengthCDS <- sum(GenomicRanges::width(
        cds[match(geneData$REFSEQ, names(cds))]))
    geneData <- geneData[order(geneData$ENTREZID, -geneData$lengthCDS), ]
    geneData <- geneData[!duplicated(geneData$ENTREZID), ]
    cds <- cds[which(names(cds) %in% geneData$REFSEQ)]
    # (2) Collapse UTR's
    #      Step 1: Filter based on valid RefSeq ID from cds
    #       Note that UTR list will still contain duplicate RefSeq entries,
    #       due to multiple UTRs for one CDS, and multiple UTRs for multiple
    #       genomic copies of the same CDS
    #       Example: NM_000513 (multiple genomic copies)
    #                NM_000982 (multiple UTR lengths & multiple genomic loci)
    cat("Collapsing duplicate 5'UTR entries\n")
    utr5 <- txBySec[[whichUTR5]]
    utr5.match <- DedupeBasedOnNearestRef(utr5, cds, showPb = TRUE)
    cat("Collapsing duplicate 3'UTR entries\n")
    utr3 <- txBySec[[whichUTR3]]
    utr3.match <- DedupeBasedOnNearestRef(utr3, cds, showPb = TRUE)
    # (3) Collapse promoters (if included)
    #     Note that promoters have a different GRangesList structure:
    #     Duplicates are stored as multiple entries under the same ID
    #     Example: NM_000071
    if (length(whichPromoter) > 0) {
        cat("Collapsing duplicate promoter entries\n")
        promoter <- txBySec[[whichPromoter]]
        # Expand multiple entries under the same ID
        promoter <- as(unlist(promoter), "CompressedGRangesList")
        promoter.match <- DedupeBasedOnNearestRef(promoter, utr5.match, showPb = TRUE)
    }
    # (3) Collapse introns (if included)
    if (length(whichIntron) > 0) {
        cat("Collapsing duplicate introns entries\n")
        intron <- txBySec[[whichIntron]]
        intron.match <- DedupeBasedOnNearestRef(intron, cds, showPb = TRUE)
    }
    # Generate return object
    txBySec.collapsed <- list()
    for (i in 1:length(txBySec)) {
        if (grepl("(CDS|coding)", sections[i], ignore.case = TRUE)) {
            feat <- cds
        } else if (grepl("(5pUTR|5'UTR|UTR5)", sections[i], ignore.case = TRUE)) {
            feat <- utr5.match
        } else if (grepl("(3pUTR|3'UTR|UTR3)", sections[i], ignore.case = TRUE)) {
            feat <- utr3.match
        } else if (grepl("(Intron|Intronic)", sections[i], ignore.case = TRUE)) {
            feat <- intron.match
        } else if (grepl("Promoter", sections[i], ignore.case = TRUE)) {
            feat <- promoter.match
        }
        txBySec.collapsed[[length(txBySec.collapsed)+1]] <- feat
    }
    names(txBySec.collapsed) <- sections
    return(txBySec.collapsed)
}


#' Perform sanity check of collapsed transcriptome.
#'
#' Perform sanity check of collapsed transcriptome. See 'Details'.
#' This is a low-level function that is being called from
#' \code{BuildTx}.
#'
#' The function checks for entries with duplicate RefSeq IDs
#' across a list of transcripts sections. For example, the return
#' object of \code{CollapseTxBySec} should contain transcript
#' segments for every transcript section with a unique RefSeq ID.
#'
#' @param txBySec A \code{list} of \code{GRangesList} objects;
#' output of function \code{GetTxBySec} or \code{CollapseTxBySec}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
PerformSanityCheck <- function(txBySec) {
    # Perform sanity checks of list of transcript features
    #
    # Args:
    #   txBySec: List of GRangesList transcript features
    #
    # Returns:
    #   NULL
    dupesInSec <- vector()
    for (i in 1:length(txBySec)) {
        cat(sprintf("Sanity check for %s:\n", names(txBySec)[i]))
        cat(sprintf("(1) Number of entries = %i\n", length(txBySec[[i]])))
        cat("(2) Unique/duplicates?\n")
        t <- table(duplicated(names(txBySec[[i]])))
        print(t)
        cat("\n")
        if (length(grep("TRUE", names(t))) > 0) {
            dupesInSec <- c(dupesInSec, names(txBySec)[i])
        }
    }
    if (length(dupesInSec) > 0) {
        cat(sprintf("Found duplicate entries in sections %s.\n",
                    paste0(dupesInSec, collapse = ", ")))
    }
}

#' Get sequences for transcript segments.
#'
#' Get sequences for transcript segments.
#' This is a low-level function that is being called from
#' \code{BuildTx}.

#'
#' @param txBySec A \code{list} of \code{GRangesList} objects;
#' output of function \code{GetTxBySec} or \code{CollapseTxBySec}.
#' @param skipIntrons A logical scalar; if \code{TRUE} sequences
#' based on introns are not extracted; this is usually a good idea
#' as long intronic sequences can lead to large memory imprints
#' (and file sizes if objects are saved).
#'
#' @return A \code{list} of \code{DNAStringSet} objects.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @export
GetTxSeq <- function(txBySec,
                     skipIntrons = TRUE) {
    # Get a list of sequences for every feature from the list
    # of transcript sections
    #
    # Args:
    #   txBySec: List of GRangesList transcript sections.
    #   skipIntrons: If TRUE, skip intron sequences.
    #
    # Returns:
    #   List of DNAStringSet sequences
    genomeVersion <- unique(unlist(lapply(txBySec, function(x) genome(x))))
    if (grepl("^hg", genomeVersion)) {
        genomePkg <- sprintf("BSgenome.Hsapiens.UCSC.%s", genomeVersion)
        if (SafeLoad(genomePkg)) {
            genome <- get(genomePkg)
        } else {
            stop(sprintf("[ERROR] Could not load genome %s.", genomeVersion))
        }
    } else if (grepl("^mm",genomeVersion)) {
        genomePkg <- sprintf("BSgenome.Mmusculus.UCSC.%s", genomeVersion)
        if (SafeLoad(genomePkg)) {
            genome <- get(genomePkg)
        } else {
            stop(sprintf("[ERROR] Could not load genome %s.",genomeVersion))
        }
    } else if (grepl("^dm",genomeVersion)) {
        genomePkg <- sprintf("BSgenome.Dmelanogaster.UCSC.%s", genomeVersion)
        if (SafeLoad(genomePkg)) {
            genome <- get(genomePkg)
        } else {
            stop(sprintf("[ERROR] Could not load genome %s.",genomeVersion))
        }
    } else if (grepl("^sac",genomeVersion)) {
        genomePkg <- sprintf("BSgenome.Scerevisiae.UCSC.%s", genomeVersion)
        if (SafeLoad(genomePkg)) {
            genome <- get(genomePkg)
        } else {
            stop(sprintf("[ERROR] Could not load genome %s.",genomeVersion))
        }
    } else {
        stop(sprintf("[ERROR] Unknown genome %s.",genomeVersion))
    }
    sections <- names(txBySec)
    if (skipIntrons == TRUE) {
        whichIntrons <- grep("(Intron|Intronic)",
                             sections, ignore.case = TRUE)
        if (length(whichIntrons) > 0) {
            cat("Skipping introns.\n")
            sections <- sections[-whichIntrons]
        }
    }
    txSequences<-list()
    pb <- txtProgressBar(min = 0,
                         max = length(sections),
                         width = 80,
                         style = 3)
    for (i in 1:length(sections)) {
        seq <- extractTranscriptSeqs(genome, txBySec[[i]])
        txSequences[[length(txSequences)+1]] <- seq
        setTxtProgressBar(pb, i)
    }
    close(pb)
    names(txSequences) <- sections
    return(txSequences)
}


#' Build a custom transcriptome.
#'
#' Build a custom, organism-specific transcriptome. See 'Details'.
#'
#' The function builds an organism-specific transcriptome containing
#' one transcript per unique Entrez ID; the transcript is selected
#' from all UCSC RefSeq annotation-based isoforms as the transcript
#' with the longest CDS, and longest upstream/downstream adjoining
#' UTRs. Transcript segments are stored per transcript section, and
#' written into a \code{.RData} file.
#' For most operations, the user will run this function once, and
#' continue with further downstream analyses. Various RNAModR
#' routines will automatically load the transcriptome data to e.g.
#' map sites to and from the transcriptome.
#' Currently, RNAModR supports analyses of human, mouse, fruitfly
#' and yeast data, based on different reference genome versions:
#' \itemize{
#' \item Homo sapiens: hg38, hg19, hg18
#' \item Mus musculus: mm10, mm9, mm8
#' \item Drosophila melanogaster: dm6, dm3, dm2
#' \item Cerevisiae saccharomyces: sacCer3, sacCer2, sacCer1
#' }
#' Reconstruction of existing transcriptome data can be achieved
#' by running \code{BuildTx} with \code{force = TRUE}. Note that
#' this will overwrite the existing RData file.
#' Running \code{BuildTx} with \code{sanityCheck = TRUE} performs
#' additional checks of the various transcriptome components, and
#' is intended for debugging purposes only. It is usually safe to
#' run with the default \code{sanityCheck = FALSE}.
#'
#' @param genomeVersion A character string; refers to a specific
#' reference genome assembly version; default is \code{"hg38"}.
#' @param force A logical scalar; if \code{TRUE} force rebuild of
#' transcriptome; this will overwrite existing data.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @import AnnotationDbi GenomeInfoDb GenomicRanges GenomicFeatures
#' RSQLite
#'
#' @examples
#' \dontrun{
#' # Build the human hg38-based reference transcriptome
#' BuildTx("hg38");
#' }
#'
#' @export
BuildTx <- function(genomeVersion = c(
                        "hg38", "hg19", "hg18",
                        "mm10", "mm9", "mm8",
                        "dm6", "dm3", "dm2",
                        "sacCer3", "sacCer2", "sacCer1"),
                    force = FALSE) {
    # Build a custom transcriptome.
    #
    # Args:
    #   genomeVersion: Genome assembly version. Default is "hg38".
    #
    # Returns:
    #   NULL
    genomeVersion <- match.arg(genomeVersion)
    fn <- sprintf("tx_%s.RData", genomeVersion)
    if (file.exists(fn) && (force == FALSE)) {
        cat("Found existing transcriptome data. Nothing to do.\n")
        cat("To rebuild run with force = TRUE.\n")
    } else if (!file.exists(fn) || (force == TRUE)) {
        cat("Building the transcriptome. This will take a few minutes.\n")
        cat("Note that this step should only be done once.\n")
        # Stage 1 - Check package dependencies
        cat(sprintf("%s Stage 1/6: Checking package dependencies.\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
        CheckPkgDependencies(genomeVersion = genomeVersion)
        # Stage 2 - Create txdb object
        cat(sprintf("%s Stage 2/6: Getting transcripts and gene annotations.\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
        txdb <- GetTxDb(genomeVersion = genomeVersion,
                        standardChrOnly = TRUE,
                        force = force)
        geneXID <- GetGeneIds(txdb)
        # Stage 3 - Split txdb by transcript sections
        cat(sprintf("%s Stage 3/6: Splitting transcripts by section.\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
        txBySec <- GetTxBySec(txdb)
        # Stage 4 - Collapse isoforms
        cat(sprintf("%s Stage 4/6: Collapsing isoforms.\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
        txBySec <- CollapseTxBySec(txBySec, geneXID)
        # Stage 5 - Extract sequences
        cat(sprintf("%s Stage 5/6: Extracting sequences.\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
        seqBySec <- GetTxSeq(txBySec)
        # Stage 6 - Save objects
        cat(sprintf("%s Stage 6/6: Storing results in file.\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
        save(geneXID, txBySec, seqBySec,
             file = fn,
             compress = "gzip")
        cat(sprintf("%s [DONE]\n",
                    format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")))
    }
}


#' Build a test transcriptome.
#'
#' Build a test transcriptome. See 'Details'.
#'
#' The function builds a transcriptome for testing purposes.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @import GenomicRanges GenomicFeatures RSQLite
#'
#' @export
BuildTxTest <- function() {
# Generate txBySec
    sections <- c("5'UTR", "CDS", "3'UTR")
    tx <- data.frame(tx_id = seq(1,3),
                     tx_name=sprintf("tx%i",seq(1,3)),
                     tx_chrom="chr1",
                     tx_strand=c("-", "+", "+"),
                     tx_start=c(1, 2001, 3001),
                     tx_end=c(999, 2199, 5199))
    splice <-  data.frame(tx_id = c(1L, 2L, 2L, 2L, 3L, 3L),
                          exon_rank=c(1, 1, 2, 3, 1, 2),
                          exon_start=c(1, 2001, 2101, 2131, 3001, 4001),
                          exon_end=c(999, 2085, 2144, 2199, 3601, 5199),
                          cds_start=c(1, 2022, 2101, 2131, 3201, 4001),
                          cds_end=c(999, 2085, 2144, 2193, 3601, 4501))
    genes <- cbind.data.frame(tx_name = sprintf("tx%i", seq(1,3)),
                              gene_id = sprintf("gene%i", seq(1,3)))
    chrominfo <- cbind.data.frame(chrom = "chr1", length = 5199, is_circular = FALSE)
    txdb <- makeTxDb(tx, splice, genes = genes, chrominfo = chrominfo)
    txBySec <- list(fiveUTRsByTranscript(txdb, use.names = TRUE),
                    cdsBy(txdb, by = "tx", use.names = TRUE),
                    threeUTRsByTranscript(txdb, use.names = TRUE))
    names(txBySec) <- sections
# Generate genome
    bases <- c("A", "C", "G", "T")
    genome <- DNAString(paste(sample(bases, 5199, replace = TRUE), collapse = ""))
# Generate seqBySec
    seqBySec <- lapply(txBySec, function(x) extractTranscriptSeqs(genome, ranges(x)))
# Generate geneXID
    geneXID <- cbind.data.frame(REFSEQ = sprintf("tx%i", seq(1, 3)),
                                ENTREZID = seq(1000, 3000, length.out = 3),
                                SYMBOL = sprintf("gene%i", seq(1, 3)),
                                ENSEMBL = sprintf("ENSEMBL%i", seq(1, 3)),
                                UNIGENE = sprintf("UNIGENE%i", seq(1,3)),
                                GENENAME = sprintf("name%i", seq(1,3)))
    save(geneXID, txBySec, seqBySec, file = "tx_test.RData", compress = "gzip")
}
