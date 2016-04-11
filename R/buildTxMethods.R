#' Download a new or load an existing reference transcriptome.
#'
#' Download a new or load an existing reference transcriptome.
#' See 'Details'.
#' This is a low-level function that is being called from
#' \code{BuildTx}. 

#'
#' If no local transcriptome db file exists or \code{forceDownload}
#' is \code{TRUE}, the function generates an reference transcriptome
#' by calling \code{makeTxDbFromUCSC} from the \code{GenomicFeatures}
#' package; the return object is a \code{TxDb} object and is saved
#' locally as an sqlite db.
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
#' @param forceDownload A logical scalar; if \code{TRUE} force
#' download and overwrite existing db file; default is \code{FALSE}.
#' @param verbose A logical scalar; print additional information
#' about db; default is \code{FALSE}.
#'
#' @return A \code{TxDb} object. See 'Details'.
#'
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
                    forceDownload = FALSE,
                    verbose = FALSE) {
    # Download or load (if file exists) sqlite transcript database based on
    # RefSeq annotation from UCSC. Store database as sqlite file.
    #
    # Args:
    #   genomeVersion: Genome assembly version. Default is "hg38".
    #   standardChrOnly: As it says. Default is TRUE.
    #   forceDownload: Force re-downloading transcript database from UCSC.
    #                  Default is FALSE.
    #   verbose: Print additional output. Default is FALSE.
    #
    # Returns:
    #   TxDb object.
    sqliteFile <- sprintf("txdb_%s.sqlite", genomeVersion);
    if ((!file.exists(sqliteFile)) || (forceDownload)) {
        txdb <- makeTxDbFromUCSC(genome = genomeVersion,
                                 tablename = "refGene");
        saveDb(txdb, file = sqliteFile);
    } else {
        cat(sprintf("Found existing sqlite database. Using %s ...\n",
                    sqliteFile));
        txdb <- loadDb(file = sqliteFile);
    }
    if (standardChrOnly) {
        txdb <- keepStandardChromosomes(txdb);
    }
    if (verbose) {
        cat("List of chromosomes:\n");
        print(seqlevels(txdb));
        cat("List of columns:\n");
        print(columns(txdb));
    }
    return(txdb);
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
#'   \item UNIGENE: UniGene gene ID
#'   \item GENENAME: Gene name
#' }
#'
#' @param txdb A \code{TxDb} object.
#'
#' @return A \code{data.frame}. See 'Details'.
#' 
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
    genomeVersion <- unique(genome(txdb));
    if (length(genomeVersion)==0) {
        stop("Could not identify genome version from transcript database.");
    }
    if (grepl("(hg38|hg19)", genomeVersion, ignore.case = TRUE)) {
        if (SafeLoad("org.Hs.eg.db")) {
            refDb <- get("org.Hs.eg.db");
        } else {
            stop(sprintf("Could not load gene annotation for %s.",
                         genomeVersion));
        }
    } else if (grepl("(mm10|mm9)", genomeVersion, ignore.case = TRUE)) {
        if (SafeLoad("org.Mm.eg.db")) {
            refDb <- get("org.Mm.eg.db");
        } else {
            stop(sprintf("Could not load gene annotation for %s.",
                         genomeVersion));
        }
    } else {
        stop(sprintf("Unknown genome %s.", genomeVersion));
    }
    tx <- transcripts(txdb);
    geneXID <- select(txdb,
                      keys = tx$tx_name,
                      columns="GENEID",
                      keytype = "TXNAME");
    colnames(geneXID)[1:2] <- c("REFSEQ", "ENTREZID");
    identifier <- c("SYMBOL", "ENSEMBL", "UNIGENE", "GENENAME");
    for (i in 1:length(identifier)) {
        geneXID[, ncol(geneXID) + 1] <- mapIds(refDb,
                                               keys = geneXID[, 2],
                                               column = identifier[i],
                                               keytype = "ENTREZID",
                                               multiVals = "first");
    }
    colnames(geneXID)[3:6] <- identifier;
    geneXID<-geneXID[!duplicated(geneXID[, 1]), ];
    return(geneXID);
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
#' @param verbose A logical scalar; if \code{TRUE}, print more
#' output and additional warnings; default is \code{FALSE}.
#'
#' @return A \code{list} of \code{GRangesList} objects
#'
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
    #   verbose: Print additional output and warnings. Default is FALSE.
    #
    # Returns:
    #   List of GRangesList transcript features.
    CheckClass(txdb, "TxDb");
    oldWarningStat <- getOption("warn");
    if (!verbose) {
        options(warn = -1);
    }
    validFeat <- c("Promoter",
                   "CDS", "coding",
                   "5pUTR", "5'UTR", "UTR5",
                   "3pUTR", "3'UTR", "UTR3",
                   "Intron", "Intronic");
    isValidFeat <- grepl(sprintf("(%s)", paste0(validFeat, collapse = "|")),
                         sections, ignore.case = TRUE);
    if (!all(isValidFeat)) {
        stop(sprintf("%s is not a valid features.", sections[!isValidFeat]));
    }
    txBySec <- list();
    for (i in 1:length(sections)) {
        if (grepl("(CDS|coding)", sections[i], ignore.case = TRUE)) {
            sec <- cdsBy(txdb, by="tx", use.names = TRUE);
        } else if (grepl("(5pUTR|5'UTR|UTR5)", sections[i], ignore.case = TRUE)) {
            sec <- fiveUTRsByTranscript(txdb, use.names = TRUE);
        } else if (grepl("(3pUTR|3'UTR|UTR3)", sections[i], ignore.case = TRUE)) {
            sec <- threeUTRsByTranscript(txdb, use.names = TRUE);
        } else if (grepl("(Intron|Intronic)", sections[i], ignore.case = TRUE)) {
            sec <- intronsByTranscript(txdb, use.names = TRUE);
        } else if (grepl("Promoter", sections[i], ignore.case = TRUE)) {
            promoter <- promoters(txdb,
                                  columns = c("tx_id", "tx_name"),
                                  upstream = promUpstream,
                                  downstream = promDownstream);
            sec <- GenomicRanges::split(promoter, promoter$tx_name);
        }
        txBySec[[length(txBySec) + 1]] <- sec;
        if (verbose) {
            print(sprintf("%i entries for section %s.",
                          length(sec), sections[i]));
        }
    }
    names(txBySec) <- sections;
    options(warn = oldWarningStat);
    return(txBySec);
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
#' @keywords internal
#'
#' @export
DedupeBasedOnNearestRef <- function(query, ref, showPb = FALSE) {
    query <- query[which(elementLengths(query) > 0)];
    query <- query[which(names(query) %in% names(ref))];
    query <- query[order(names(query), -sum(width(query)))];
    dupes <- which(duplicated(names(query)));
    dupeID <- unique(names(query[dupes]));
    rem <- vector();
    if (showPb) {
        pb <- txtProgressBar(min = 0,
                             max = length(dupeID),
                             width = 80,
                             style = 3);
    }
    for (i in 1:length(dupeID)) {
        idxQuery <- which(names(query) == dupeID[i]);
        idxRef <- which(names(ref) == dupeID[i]);
        if (length(idxQuery) > 1) {
            dist <- distance(
                GenomicRanges::unlist(range(query[idxQuery])),
                range(ref[[idxRef]]));
            idxMinDist <- which.min(dist);
            rem <- c(rem, idxQuery[-idxMinDist]);
            if (showPb) {
                setTxtProgressBar(pb, i);
            }
        }
    }
    if (showPb) {
        close(pb);
    }
    return(query[-rem]);
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
    sections <- names(txBySec);
    whichCDS <- grep("(CDS|coding)", sections, ignore.case = TRUE);
    if (!length(whichCDS)) {
        stop("Transcript feature list does not contain CDS entries.");
    }
    whichUTR5 <- grep("(5pUTR|5'UTR|UTR5)", sections, ignore.case = TRUE);
    if (!length(whichUTR5)) {
        stop("Transcript feature list does not contain 5'UTR entries.");
    }
    whichUTR3 <- grep("(3pUTR|3'UTR|UTR3)", sections, ignore.case = TRUE);
    if (!length(whichUTR3)) {
        stop("Transcript feature list does not contain 3'UTR entries.");
    }
    whichPromoter <- grep("Promoter", sections, ignore.case = TRUE);
    whichIntron <- grep("(Intron|Intronic)", sections, ignore.case = TRUE);
    # (1) Collapse CDS:
    #      Step 1: Same RefSeq ID, multiple identical genomic copies
    #       Example: RefSeq ID = NM_000513
    #      Step 2: Same Entrez ID, multiple CDS isoforms
    #       Example: Entrez ID = 10001
    cds <- txBySec[[whichCDS]];
    cds <- cds[order(names(cds), -sum(GenomicRanges::width(cds)))];
    cds <- cds[!duplicated(names(cds))];
    geneData <- geneXID[which(geneXID[, 1] %in% names(cds)), ];
    geneData$lengthCDS <- sum(GenomicRanges::width(
        cds[match(geneData[, 1], names(cds))]));
    geneData <- geneData[order(geneData[, 2], -geneData[, ncol(geneData)]), ];
    geneData <- geneData[!duplicated(geneData[, 2]), ];
    cds <- cds[which(names(cds) %in% geneData[, 1])];
    # (2) Collapse UTR's
    #      Step 1: Filter based on valid RefSeq ID from cds
    #       Note that UTR list will still contain duplicate RefSeq entries,
    #       due to multiple UTRs for one CDS, and multiple UTRs for multiple
    #       genomic copies of the same CDS
    #       Example: NM_000513 (multiple genomic copies)
    #                NM_000982 (multiple UTR lengths & multiple genomic loci)
    cat("Collapsing duplicate 5'UTR entries\n");
    utr5 <- txBySec[[whichUTR5]];
    utr5.match <- DedupeBasedOnNearestRef(utr5, cds, showPb = TRUE);
    cat("Collapsing duplicate 3'UTR entries\n");
    utr3 <- txBySec[[whichUTR3]];
    utr3.match <- DedupeBasedOnNearestRef(utr3, cds, showPb = TRUE);
    # (3) Collapse promoters (if included)
    #     Note that promoters have a different GRangesList structure:
    #     Duplicates are stored as multiple entries under the same ID
    #     Example: NM_000071
    if (length(whichPromoter) > 0) {
        cat("Collapsing duplicate promoter entries\n");
        promoter <- txBySec[[whichPromoter]];
        # Expand multiple entries under the same ID
        promoter <- as(GenomicRanges::unlist(promoter), "GRangesList");
        promoter.match <- DedupeBasedOnNearestRef(promoter, utr5.match, showPb = TRUE);
    }
    # (3) Collapse introns (if included)
    if (length(whichIntron) > 0) {
        cat("Collapsing duplicate introns entries\n");
        intron <- txBySec[[whichIntron]];
        intron.match <- DedupeBasedOnNearestRef(intron, cds, showPb = TRUE);
    }
    # Generate return object
    txBySec.collapsed <- list();
    for (i in 1:length(txBySec)) {
        if (grepl("(CDS|coding)", sections[i], ignore.case = TRUE)) {
            feat <- cds;
        } else if (grepl("(5pUTR|5'UTR|UTR5)", sections[i], ignore.case = TRUE)) {
            feat <- utr5.match;
        } else if (grepl("(3pUTR|3'UTR|UTR3)", sections[i], ignore.case = TRUE)) {
            feat <- utr3.match;
        } else if (grepl("(Intron|Intronic)", sections[i], ignore.case = TRUE)) {
            feat <- intron.match;
        } else if (grepl("Promoter", sections[i], ignore.case = TRUE)) {
            feat <- promoter.match;
        }
        txBySec.collapsed[[length(txBySec.collapsed)+1]] <- feat;
    }
    names(txBySec.collapsed) <- sections;
    return(txBySec.collapsed);
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
    dupesInSec <- vector();
    for (i in 1:length(txBySec)) {
        cat(sprintf("Sanity check for %s:\n", names(txBySec)[i]));
        cat(sprintf("(1) Number of entries = %i\n", length(txBySec[[i]])));
        cat("(2) Unique/duplicates?\n");
        t <- table(duplicated(names(txBySec[[i]])));
        print(t);
        cat("\n");
        if (length(grep("TRUE", names(t))) > 0) {
            dupesInSec <- c(dupesInSec, names(txBySec)[i]);
        }
    }
    if (length(dupesInSec) > 0) {
        cat(sprintf("Found duplicate entries in sections %s.\n",
                    paste0(dupesInSec, collapse = ", ")));
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
#' @param verbose A logical scalar; if \code{TRUE} print
#' additional output; default is \code{FALSE}.
#'
#' @return A \code{list} of \code{DNAStringSet} objects.
#'
#' @keywords internal
#'
#' @export
GetTxSeq <- function(txBySec,
                     skipIntrons = TRUE,
                     verbose = FALSE) {
    # Get a list of sequences for every feature from the list
    # of transcript sections
    #
    # Args:
    #   txBySec: List of GRangesList transcript sections.
    #   verbose: Print additional output. Default is FALSE.
    #
    # Returns:
    #   List of DNAStringSet sequences
    genomeVersion <- unique(unlist(lapply(txBySec, function(x) genome(x))));
    validGenomes <- c("hg38", "hg19", "hg18", "mm10", "mm9", "mm8");
    if (length(genomeVersion) == 0) {
        stop("Could not identify genome version from list of transcript segments.");
    }
    if (!(genomeVersion %in% validGenomes)) {
        ss <- sprintf("%s is not a valid genome version.\n", genomeVersion);
        ss <- sprintf("%sValid genome versions are: %s",
                      ss,
                      paste0(validGenomes, collapse = ", "));
        stop(ss);
    }
    if (grepl("^hg", genomeVersion)) {
        genomePkg <- sprintf("BSgenome.Hsapiens.UCSC.%s", genomeVersion);
        if (SafeLoad(genomePkg)) {
            genome <- get(genomePkg);
        } else {
            stop(sprintf("Could not load genome %s.", genomeVersion));
        }
    } else if (grepl("^mm",genomeVersion)) {
        genomePkg <- sprintf("BSgenome.Mmusculus.UCSC.%s", genomeVersion);
        if (SafeLoad(genomePkg)) {
            genome <- get(genomePkg);
        } else {
            stop(sprintf("Could not load genome %s.",genomeVersion));
        }
    } else {
        stop(sprintf("Unknown genome %s.",genomeVersion));
    }
    if (verbose) {
        cat(sprintf("Loaded genome version %s.",genomeVersion));
    }
    sections <- names(txBySec);
    if (skipIntrons == TRUE) {
        whichIntrons <- grep("(Intron|Intronic)",
                             sections, ignore.case = TRUE);
        if (length(whichIntrons) > 0) {
            cat("Skipping introns.\n");
            sections <- sections[-whichIntrons];
        }
    }
    txSequences<-list();
    pb <- txtProgressBar(min = 0,
                         max = length(sections),
                         width = 80,
                         style = 3);
    for (i in 1:length(sections)) {
        seq <- extractTranscriptSeqs(genome, txBySec[[i]]);
        txSequences[[length(txSequences)+1]] <- seq;
        setTxtProgressBar(pb, i);
    }
    close(pb);
    names(txSequences) <- sections;
    return(txSequences);
}


#' Build a custom transcriptome.
#'
#' Build a custom transcriptome. See 'Details'.
#'
#' The function builds a custom transcriptome that consists of
#' one transcript per unique Entrez ID; the transcript is
#' determined from all RefSeq-based isoforms as the transcript
#' with the longest CDS, and adjoining upstream/downstream longest
#' UTRs. Transcript segments are stored per transcript section, and
#' written to a \code{.RData} file.
#' 
#' @param genomeVersion A character string; refers to a specific
#' reference genome assembly version; default is \code{"hg38"}.
#' @param sanityCheck A logical scalar; if \code{TRUE} perform
#' sanity checks.
#'
#' @import AnnotationDbi GenomeInfoDb GenomicRanges GenomicFeatures
#' IRanges RSQLite
#'
#' @export
BuildTx <- function(genomeVersion = "hg38", sanityCheck = FALSE) {
    # Build a custom transcriptome.
    #
    # Args:
    #   genomeVersion: Genome assembly version. Default is "hg38".
    #
    # Returns:
    #   NULL
    cat("Building the transcriptome. This will take a few minutes.\n");
    cat("This should only need to be done once.\n");
    cat(sprintf("%s Stage 1/5: Getting transcripts and gene annotations.\n",
                format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")));
    txdb <- GetTxDb(genomeVersion = genomeVersion);
    geneXID <- GetGeneIds(txdb);
    cat(sprintf("%s Stage 2/5: Splitting transcript by section.\n",
                format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")));
    txBySec <- GetTxBySec(txdb);
    if (sanityCheck) {
        PerformSanityCheck(txBySec);
    }
    cat(sprintf("%s Stage 3/5: Collapsing isoforms.\n",
                format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")));
    txBySec <- CollapseTxBySec(txBySec, geneXID);
    if (sanityCheck) {
        PerformSanityCheck(txBySec);
    }
    cat(sprintf("%s Stage 4/5: Obtaining sequences.\n",
                format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")));
    seqBySec <- GetTxSeq(txBySec);
    cat(sprintf("%s Stage 5/5: Storing results in file.\n",
                format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")));
    save(geneXID, txBySec, seqBySec,
         file = sprintf("tx_%s.RData", genomeVersion),
         compress = "gzip");
    cat(sprintf("%s [DONE]\n",
                format(Sys.time(), "[%a %b %d %Y %H:%M:%S]")));
}
