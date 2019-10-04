#' Map genome coordinates to transcript coordinates.
#'
#' Map genome coordinates to transcript coordinates. See 'Details'.
#' This is a low-level function that is being called from \code{SmartMap}. 
#' There is no guarantee that this function will get exported in future
#' releases of RNAModR. Use at your own risk.
#'
#' The function maps genomic coordinates from \code{gr} to transcript region
#' coordinates from \code{txBySec}. The function returns a \code{list} of
#' \code{DataFrame} objects, each with the following columns:
#' \enumerate{
#'     \item \code{locus_in_txx_region}: A \code{GRanges} object
#'     \item \code{locus_in_genome}: A \code{GRanges} object
#'     \item \code{score}: A \code{numeric} vector
#'     \item \code{id}: A \code{character} vector
#'     \item \code{tx_region}: A \code{character} vector
#'     \item \code{tx_region_width}: An \code{integer} vector
#'     \item \code{tx_region_sequence}: A \code{DNAStringSet} object
#'     \item \code{tx_refseq}: A \code{character} vector
#'     \item \code{gene_entrez}: A \code{character} vector
#'     \item \code{gene_symbol}: A \code{character} vector
#'     \item \code{gene_ensembl}: A \code{character} vector
#'     \item \code{gene_name}: A \code{character} vector
#' }
#' 
#' @param gr A \code{GRanges} object; specifies the list of genomic features to 
#' be mapped.
#' @param txBySec A \code{list} of \code{GRangesList} objects; specifies the 
#' reference transcriptome; the object is usually a result of running 
#' \code{BuildTx}.
#' @param seqBySec A \code{list} of \code{DNAStringSet} objects; sequences of 
#' (some or all) segments from \code{txBySec}; the object is usually a result 
#' of running \code{BuildTx}.
#' @param geneXID A \code{DataFrame}; specifies different gene IDs; the object 
#' is usually a result of running \code{BuildTx}.
#' @param ignore.strand A logical scalar; if \code{TRUE} strand information is 
#' ignored when mapping genome coordinates to transcript coordinates; default 
#' is \code{FALSE}.
#' @param showPb A logical scalar; if \code{TRUE} show a progress
#' bar; default is \code{FALSE}.
#'
#' @return A \code{list} of \code{data.frame} objects. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @importFrom S4Vectors DataFrame
#' 
#' @keywords internal
SmartMap.ToTx <- function(gr,
                          txBySec,
                          seqBySec,
                          geneXID,
                          ignore.strand = FALSE,
                          showPb = FALSE) {

    # Initialise progress bar if `showPb == TRUE`
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(txBySec), style = 3, width = 60)

    # Map to different `txBySec` regions and return `list` of `DataFrame`s
    lst <- lapply(
        setNames(1:length(txBySec), names(txBySec)), 
        function(i) {
    
            if (showPb) setTxtProgressBar(pb, i)
      
            region <- names(txBySec)[i]
      
            # [UPDATE September 2019] `mapToTranscripts` is now more stringent
            # in that it requires `seqlevels(gr)` to be a subset of 
            # `seqlevels(txBySec[[i]]@unlistData)`. In other words, if `gr` 
            # contains chromosomes that are _not_ included in `txBySec` this 
            # will  throw a  critical error. So we need to make sure that we 
            # only have  entries in `gr` on chromosomes that are included as 
            # seqlevels of transcripts in `txBySec`
            gr <- keepSeqlevels(
                gr, 
                intersect(
                  seqlevels(txBySec[[region]]@unlistData), 
                  seqlevels(gr)),
                pruning.mode = "coarse")
            
            # Map to transcripts; the output object has two metadata columns:
            #   xHits = index of the mapped `gr` object
            #   transcriptsHits = index of the `txBySec[[i]]` object
            hits <- mapToTranscripts(
              gr, 
              txBySec[[region]], 
              ignore.strand = ignore.strand)
            
            # Return DataFrame with columns
            #   locus_in_tx_region: <GRanges>
            #   locus_in_genome:    <GRanges>
            #   score:              <numeric>
            #   id:                 <character>
            #   tx_region:          <character>
            #   tx_region_width:    <integer>
            #   tx_region_sequence: <DNAStringSet>
            #   tx_refseq:          <character>
            #   gene_entrez:        <character>
            #   gene_symbol:        <character>
            #   gene_ensembl:       <character>
            #   gene_name:          <character>
            cbind(
              setNames(
                DataFrame(
                  granges(hits, use.names = F, use.mcols = F)), 
                "locus_in_tx_region"),
              setNames(
                DataFrame(
                  granges(gr[mcols(hits)$xHits], use.names = F, use.mcols = F)), 
                "locus_in_genome"),
              setNames(
                mcols(gr[mcols(hits)$xHits]),
                c("score", "id")),
              setNames(DataFrame(
                region, 
                sum(width(txBySec[[region]][mcols(hits)$transcriptsHits]))),
                c("tx_region", "tx_region_width")),
              setNames(DataFrame(
                seqBySec[[region]][match(
                  seqnames(hits), names(seqBySec[[region]]))]),
                "tx_region_sequence"),
              geneXID[match(seqnames(hits), geneXID$tx_refseq), ])
        
        }
    )
    
    if (showPb) close(pb)

    lst
    
}


#' Map genome coordinates to transcript coordinates.
#'
#' Map genome coordinates to transcript coordinates.
#'
#' The function maps genomic coordinates from \code{locus} to transcript 
#' section coordinates. The function automatically loads a reference 
#' transcriptome based on \code{refGenome}. An error is produced if a reference 
#' transcriptome could not be found. This usually means that \code{BuildTx} was 
#' not yet run successfully.
#' The function returns a \code{txLoc} object of mapped positions.
#'
#' @param gr A \code{GRanges} object; specifies the list of of genomic features 
#' to be mapped.
#' @param id A character string; specifies a name for loci from \code{gr}; if 
#' \code{NULL} then \code{id = ""}; default is \code{NULL}.
#' @param refGenome A character string; specifies a specific reference genome 
#' assembly version based on which the matching transcriptome is loaded; 
#' default is \code{"hg38"}.
#' @param ignore.strand A logical scalar; if \code{TRUE} strand information is 
#' ignored when mapping genome coordinates to transcript coordinates; default 
#' is \code{FALSE}.
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
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' txSites <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' }

#' @export
SmartMap <- function(gr,
                     id = NULL,
                     refGenome = "hg38",
                     ignore.strand = FALSE) {

    # Sanity check
    CheckClass(gr, "GRanges")
    
    # Load transcriptome data
    refTx <- sprintf("tx_%s.RData", refGenome)
    if (!file.exists(refTx)) {
        ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    load(refTx)
    requiredObj <- c("geneXID", "seqBySec", "txBySec")
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript objects not found.")
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "))
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    geneXID <- get("geneXID")
    seqBySec <- get("seqBySec")
    txBySec <- get("txBySec")
    
    # Map coordinates to transcript
    loci <- SmartMap.ToTx(
        gr, txBySec, seqBySec, geneXID, 
        ignore.strand = ignore.strand,
        showPb = TRUE)

    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = if (is.null(id)) "" else id,
        refGenome = refGenome,
        version = as.character(Sys.Date()))

}


#' Construct \code{txLoc}-suitable \code{dataframe} object from
#' mapping transcript loci to genomic loci. (CURRENTLY BROKEN)
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
#' @keywords internal
#' 
#' @import Biostrings GenomicRanges IRanges
GetLocus.MapFromTranscripts <- function(gr, ref, seq, section, geneXID) {
    genomePos <- mapFromTranscripts(gr, ref)
    idxQuery <- mcols(genomePos)$xHits
    idxRef <- mcols(genomePos)$transcriptsHits
    dataRef <- ref[idxRef]
    idxSeq <- match(names(genomePos), names(seq))
    dataSeq <- seq[idxSeq]
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
        as.character(dataSeq))
    colnames(locus) <- c("REFSEQ", "TXSTART", "TXEND", "TXWIDTH",
                         "CHR", "START", "STOP", "WIDTH",
                         "STRAND", "SCORE", "ID",
                         "GENE_REGION", "GENE_REFSEQ", "GENE_ENTREZ",
                         "GENE_SYMBOL", "GENE_ENSEMBL", "GENE_UNIGENE",
                         "GENE_NAME", "GENE_CHR", "GENE_START",
                         "GENE_STOP", "GENE_WIDTH", "GENE_STRAND",
                         "REGION_TXWIDTH", "REGION_SEQ")
    return(locus)
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
#' \item \code{method = "ntAbund"}: Null sites are generated based on the
#' position of all non-modified nucleotides of type \code{nt} in those
#' transcript sections that also contain a modified site of the same type and
#' as specified in \code{locus}. 
#' For example, if \code{locus} contains a list of m$^{6}$A sites, the list of
#' null sites consists of all non-methylated adenosines in those transcripts
#' that contain at least one m$^{6}$A site.
#' \item \code{method = "perm"}: Null sites are generated by randomly permuting
#' the position of sites from \code{locus} uniformly within the corresponding 
#' transcript section. Note that this will generate a list of null sites with 
#' the same abundance ratios across  transcript sections as the list of sites 
#' from \code{locus}. It is therefore not useful for assessing an enrichment of 
#' sites within a particular transcript section. In fact, this method should 
#' not be used and is included purely for paedagogical purposes (to demonstrate 
#' the importance of a sensible null distribution). It is likely that this 
#' method will be removed from future RNAModR versions.
#' 
#' It is import to emphasise that any downstream enrichment analysis may depend
#' critically on the choice of the null distribution. For example, a position 
#' permution-based null distribution may not be a valid null distribution, if 
#' the distribution of nucleotides is highly non-uniform across a transcript 
#' section. This is the case e.g. for the spatial distribution of cytosines 
#' within and across the 5'UTR, CDS and/or 3'UTR. In this case, a permutation-
#' based distribution of cytosines will not give a sensible null distribution.
#' Instead, a sensible null distribution can be derived from the position of 
#' all cytosines in the relevant transcript region containing the methylated 
#' cytosine site in \code{locus}. \code{method = "ntAbund"} generates a list of 
#' null sites using this approach.
#'
#' @param txLoc A \code{txLoc} object.
#' @param id A character string; identifier for null sites; if \code{NULL} then 
#' \code{id = "null"}; default is \code{NULL}.
#' @param method A character string; specifies the method used to the generate 
#' null distribution; if \code{method == "ntAbund"} the position of all 
#' nucleotides specified by \code{nt} will be used as null sites; if 
#' \code{method == "perm"} (DEPRECATED, see 'Details') then null sites will be 
#' generated by uniform-randomly shuffling candidatepositions from \code{locus} 
#' within the corresponding transcript region; default is \code{"ntAbund"}.
#' @param nt A single character; if \code{method == "ntAbund"}, use \code{nt} 
#' to derive distribution of null sites.
#' @param showPb A logical scalar; if \code{TRUE} show a progress
#' bar; default is \code{TRUE}.
#'
#' @return A \code{txLoc} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats runif setNames
#' @importFrom S4Vectors DataFrame
#' 
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' negSites <- GenerateNull(posSites, method = "ntAbund", nt = "A")
#' }
#'
#' @export
GenerateNull <- function(txLoc,
                         id = NULL,
                         method = c("ntAbund", "perm"),
                         nt = "C",
                         showPb = TRUE)  {

    # Sanity checks
    CheckClass(txLoc, "txLoc")
    method <- match.arg(method)
    refGenome <- GetRef(txLoc)
    if (is.null(id)) {
        id <- sprintf("null_%s", GetId(txLoc))
    }
    
    # Quieten R CMD check concerns regarding "no visible binding for global
    # variable ..."
    tx_region_sequence <- tx_refseq <- NULL
    
    # Get loci from `txLoc` object
    loci.pos <- GetLoci(txLoc)

    # Load transcriptome data for `method == "ntAbund"`
    if (method == "ntAbund") {
        
        # Load transcriptome data
        refTx <- sprintf("tx_%s.RData", refGenome)
        if (!file.exists(refTx)) {
            ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
            ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                          ss, refGenome)
            stop(ss)
        }
        load(refTx)
        requiredObj <- c("geneXID", "seqBySec", "txBySec")
        if (!all(requiredObj %in% ls())) {
            ss <- sprintf("Mandatory transcript objects not found.")
            ss <- sprintf("%s\nNeed all of the following: %s",
                          ss, paste0(requiredObj, collapse = ", "))
            ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                          ss, refGenome)
            stop(ss)
        }
        geneXID <- get("geneXID")
        seqBySec <- get("seqBySec")
        txBySec <- get("txBySec")
        
    }
    
    # Activate progressbar if `showPb == TRUE`
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(loci.pos), style = 3, width = 60)

    loci <- mapply(
        function(i, sec) {
        
          # Update progress bar
          if (showPb) setTxtProgressBar(pb, i)
        
          # Loci of positive sites for region i
          locus.pos <- loci.pos[[i]]
        
          # Get sequence data for unique RefSeq transcript regions
          # Returns a DNAStringSet object
          unique_tx_region_sequence <- with(
            locus.pos,
            tx_region_sequence[!duplicated(tx_refseq)])
        
          # Extract position of all `nt`s inside every sequence
          mindex <- vmatchPattern(nt, unique_tx_region_sequence)
        
          # Exclude positive sites from null sites
          gr <- subsetByOverlaps(
            GRanges(mindex), 
            locus.pos$locus_in_tx_region, 
            invert = TRUE)
        
          # Map sites from transcriptomic to genomic coordinates
          hits <- mapFromTranscripts(gr, sec)
          
          # Copy genomic strand information to tx data
          strand(gr) <- strand(hits)
        
          # Return DataFrame with columns
          #   locus_in_tx_region: <GRanges>
          #   locus_in_genome:    <GRanges>
          #   score:              <numeric>
          #   id:                 <character>
          #   tx_region:          <character>
          #   tx_region_width:    <integer>
          #   tx_region_sequence: <DNAStringSet>
          #   tx_refseq:          <character>
          #   gene_entrez:        <character>
          #   gene_symbol:        <character>
          #   gene_ensembl:       <character>
          #   gene_name:          <character>
          cbind(
            setNames(DataFrame(gr[mcols(hits)$xHits]), "locus_in_tx_region"),
            setNames(
              DataFrame(granges(hits, use.names = FALSE, use.mcols = FALSE)), 
              "locus_in_genome"),
            setNames(DataFrame(
              rep(0, length(gr)),
              paste0(seqnames(hits), ":", start(hits), "-", end(hits))), 
              c("score", "id")),
            locus.pos[match(
              seqnames(gr[mcols(hits)$xHits]), locus.pos$tx_refseq), -(1:4)])

      },
      setNames(seq_along(loci.pos), names(loci.pos)),
      txBySec[names(loci.pos)])
    
    if (showPb) close(pb)
    
    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = id,
        refGenome = refGenome,
        version = as.character(Sys.Date()))

}


#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats runif setNames
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
    CheckClass(locus, "txLoc")
    method <- match.arg(method)
    refGenome <- GetRef(locus)
    if (is.null(id)) {
        id <- sprintf("null_%s", GetId(locus))
    }
    locus <- GetLoci(locus)
    locusNull.list <- list()
    if (method == "ntAbund") {
        LoadRefTx(refGenome)
        geneXID <- get("geneXID")
        seqBySec <- get("seqBySec")
        txBySec <- get("txBySec")
    }
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(locus), style = 3, width = 60)
    for (i in 1:length(locus)) {
        if (showPb) setTxtProgressBar(pb, i)
        if (method == "ntAbund") {
            seqData <- locus[[i]][!duplicated(locus[[i]][, 1]),
                                  c("REFSEQ", "GENE_CHR", "GENE_START",
                                    "GENE_STOP", "STRAND", "REGION_SEQ")]
            if (all(grepl("\\w+", seqData$REGION_SEQ))) {
                # Positive sites that should be excluded from null
                exclude <- locus[[i]][, c("REFSEQ", "TXSTART")]
                # All sites of a specific nucleotide
                txPos <- gregexpr(nt, seqData$REGION_SEQ)
                names(txPos) <- seqData$REFSEQ
                txPos <- data.frame(
                    ID = rep(names(txPos), sapply(txPos, length)),
                    x = unlist(txPos), stringsAsFactors = FALSE)
                # Exclude positive sites
                txPos <- rbind(txPos, setNames(exclude, names(txPos)))
                txPos <- txPos[!duplicated(txPos, fromLast = FALSE) &
                               !duplicated(txPos, fromLast = TRUE), ]
                rownames(txPos) <- seq(1, nrow(txPos))
                # Map sites back to genome
                gr <- GRanges(txPos$ID, IRanges(txPos$x, txPos$x))
                idxSec <- which(names(locus)[i] == names(txBySec))
                genomePos <- mapFromTranscripts(gr, txBySec[[idxSec]])
                genomePos <- as.data.frame(genomePos);           
            } else {
                sec <- names(locus)[i]
                ss <- sprintf("Cannot generate null for %s", sec)
                ss <- sprintf("%s: Missing sequence information.", ss)
                ss <- sprintf("%s\nEither exclude %s from downstream analysis",
                              ss, sec)
                ss <- sprintf("%s, or use different null.", ss)
                warning(ss)
                txPos <- cbind.data.frame(
                    ID = seqData$REFSEQ,
                    x = rep("*", nrow(seqData)))
                genomePos <- cbind.data.frame(
                    start = rep("*", nrow(seqData)),
                    end = rep("*", nrow(seqData)),
                    width = rep(1, nrow(seqData)))
            }
            idxSeq <- match(txPos$ID, seqData$REFSEQ)
            idxLoc <- match(txPos$ID, locus[[i]]$GENE_REFSEQ)
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
            colnames(locusNull) <- colnames(locus[[i]])
            locusNull <- Unfactor(locusNull)
        } else if (method == "perm") {
            locusNull <- locus[[i]]
            locusNull$TXSTART <- apply(
                locusNull, 1, function(x) {
                    round(runif(1,
                                min = 1,
                                max = as.numeric(x["REGION_TXWIDTH"]) - 1))})
            locusNull$TXEND <- locusNull$TXSTART
            locusNull[, c("CHR", "START", "STOP",
                          "WIDTH", "STRAND", "SCORE")] <- matrix(
                              "*",
                              nrow = nrow(locusNull),
                              ncol = 6)
            locusNull$ID <- "unif_perm"
        }
        locusNull.list[[length(locusNull.list) + 1]] <- locusNull
    }
    if (showPb) close(pb)
    names(locusNull.list) <- names(locus)
    obj <- new("txLoc",
               loci = locusNull.list,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()))
    return(obj)
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
    CheckClass(locus, "txLoc")
    id <- GetId(locus)
    refGenome <- GetRef(locus)
    locus <- GetLoci(locus)
    if (!SafeLoad("Biostrings")) {
        stop("Could not load library Biostrings.")
    }
    freq.list <- list()
    for (i in 1:length(locus)) {
        seq <- BStringSet(locus[[i]]$REGION_SEQ)
        nucFreq <- letterFrequency(
            seq,
            letters = c("A", "C", "G", "T"))
        sectionGC <- rowSums(nucFreq[, c(2, 3)]) / width(seq)
        xLoc <- locus[[i]]$TXSTART
        if (is.numeric(xLoc)) {
            x1 <- xLoc - flank
            x2 <- xLoc + flank
            # Make sure we're still within transcript coordinates
            sel1 <- which(x1 <= 0)
            x1[sel1] <- 1
            x2[sel1] <- 2 * xLoc[sel1] - 1
            sel2 <- which(x2 > locus[[i]]$REGION_TXWIDTH)
            x2[sel2] <- locus[[i]]$REGION_TXWIDTH[sel2]
            x1[sel2] <- 2 * xLoc[sel2] - (locus[[i]]$REGION_TXWIDTH)[sel2]
            subSeq <- BStringSet(substr(locus[[i]]$REGION_SEQ, x1, x2))
            nucFreq <- letterFrequency(subSeq,
                                       letters = c("A", "C", "G", "T"))
            siteGC <- rowSums(nucFreq[, c(2, 3)]) / nchar(subSeq)
        } else {
            siteGC <- rep(NA, nrow(locus[[i]]))
        }
        dfGC <- cbind.data.frame("id" = locus[[i]]$ID,
                                 "sectionGC" = sectionGC,
                                 "siteGC" = siteGC,
                                 "siteSeq" = as.character(subSeq),
                                 stringsAsFactors = FALSE)
        freq.list[[length(freq.list) + 1]] <- dfGC
    }
    names(freq.list) <- names(locus)
    return(freq.list)
}


#' Get exon-exon boundary sites from reference transcriptome
#'
#' Get exon-exon boundary sites from transcriptome. See 'Details'.
#'
#' The function extracts exon-exon boundary sites (EEBS) from a 
#' reference trancriptome specified by \code{refGenome}, and returns 
#' a \code{GRanges} object. Boundary sites are defined as the 
#' location of those exonic single nucleotides that are closest 
#' upstream to the intronic donor site.
#'
#' @param refGenome A character string; specifies a specific
#' reference genome assembly version based on which the matching
#' transcriptome is loaded; default is \code{"hg38"}.
#' @param filter A character vector; only consider transcript
#' sections specified in \code{filter}; default is 
#' \code{c("CDS", "5'UTR")}.
#'
#' @return A \code{GRanges} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#'
#' @export
GetEEJunct <- function(refGenome = "hg38", filter = c("CDS", "5'UTR")) {
  
  # Read transcriptome data
  refTx <- sprintf("tx_%s.RData", refGenome)
  if (!file.exists(refTx)) {
    ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
    ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                  ss, refGenome)
    stop(ss)
  }
  load(refTx)
  requiredObj <- c("geneXID", "seqBySec", "txBySec")
  if (!all(requiredObj %in% ls())) {
    ss <- sprintf("Mandatory transcript objects not found.")
    ss <- sprintf("%s\nNeed all of the following: %s",
                  ss, paste0(requiredObj, collapse = ", "))
    ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                  ss, refGenome)
    stop(ss)
  }
  geneXID <- get("geneXID")
  seqBySec <- get("seqBySec")
  txBySec <- get("txBySec")
  
  # Filter regions
  if (is.null(filter)) {
    filter <- c("5'UTR", "CDS", "3'UTR")
  }
  
  # Create a list of GRanges for the acceptor/donor splice sites
  locus <- lapply(setNames(filter, filter), function(region) {
    
    # Get introns as the regions that are not annotated
    junct <- setdiff(range(txBySec[[region]]), txBySec[[region]])
    
    # Remove transcripts without introns and unlist
    junct <- unlist(junct[elementNROWS(junct) > 0])
    
    # Collapse ranges to start position (strand-specific)
    # corresponding to splicing donor sites
    junct <- resize(
      junct, width = 1, fix = "start", ignore.strand = F, use.names = T)
    
    # Shift position 1 nucleotide upstream (strand-specific) 
    # of splicing donor site and add score & id metadata columns
    # (This is to make sure that we are consistent with the BED6 format)
    junct <- promoters(junct, upstream = 1, downstream = 0)
    mcols(junct)$score = 0
    mcols(junct)$id = sprintf("upstream_donor|%s|%s", names(junct), region)

    # Remove names and return
    names(junct) <- NULL
    junct
  })      
  
  # Concatenate list of GRanges, unname GRanges and return
  gr <- unlist(as(locus, "GRangesList"))
  names(gr) <- NULL
  gr
  
}


#' Get splicing sites from reference transcriptome
#'
#' Get splicing sites from transcriptome. See 'Details'.
#'
#' The function extracts splicing sites from a reference
#' trancriptome specified by \code{refGenome}, and returns a
#' \code{GRanges} object. Splicing sites are identified as 
#' either donor (splicing site at the 5' end of the intron)
#' or acceptor (splicing site at the 3' end of the intron)
#' site.
#'
#' @param refGenome A character string; specifies a specific
#' reference genome assembly version based on which the matching
#' transcriptome is loaded; default is \code{"hg38"}.
#' @param filter A character vector; only consider transcript
#' sections specified in \code{filter}; default is 
#' \code{c("CDS", "5'UTR")}.
#'
#' @return A \code{GRanges} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @import GenomicRanges IRanges
#'
#' @export
GetSplicingSites <- function(refGenome = "hg38", filter = c("CDS", "5'UTR")) {
    
    # Read transcriptome data
    refTx <- sprintf("tx_%s.RData", refGenome)
    if (!file.exists(refTx)) {
        ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    load(refTx)
    requiredObj <- c("geneXID", "seqBySec", "txBySec")
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript objects not found.")
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "))
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    geneXID <- get("geneXID")
    seqBySec <- get("seqBySec")
    txBySec <- get("txBySec")
    
    # Filter regions
    if (is.null(filter)) {
        filter <- c("5'UTR", "CDS", "3'UTR")
    }

    # Create a list of GRanges for the acceptor/donor splice sites
    locus <- lapply(setNames(filter, filter), function(region) {
      
      # Get introns as the regions that are not annotated
      junct <- setdiff(range(txBySec[[region]]), txBySec[[region]])
      
      # Remove transcripts without introns and unlist
      junct <- unlist(junct[elementNROWS(junct) > 0])
      
      # Collapse ranges to start/end position (strand-specific)
      # corresponding to splicing donor and acceptor sites
      # Add score and id metadata columns to make sure that we are
      # consistent with the BED6 format
      junct <- unlist(as(mapply(
        function(fix, type) {
          eej <- resize(
            junct, width = 1, fix = fix, ignore.strand = F, use.names = T)
          mcols(eej)$id = sprintf("%s|%s|%s", type, names(eej), region)
          eej
        },
        c("start", "end"),
        c("donor", "acceptor") 
      ), "GRangesList"))
      
      # Remove names and return
      names(junct) <- NULL
      junct
    })      
    
    # Concatenate list of GRanges, unname GRanges and return
    gr <- unlist(as(locus, "GRangesList"))
    names(gr) <- NULL
    gr

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
GetMotifLoc <- function(motif, 
                        refGenome = "hg38", 
                        filter = c("5'UTR", "CDS", "3'UTR"), 
                        maxMM = 0, 
                        showPb = TRUE) {
    PASmotif <- c("AATAAA", "ATTAAA", "AGTAAA",
                  "TATAAA", "AAGAAA", "AATACA",
                  "AATATA", "CATAAA", "AATGAA",
                  "GATAAA", "ACTAAA", "AATAGA")
    if (is.null(motif)) {
      motif <- PASmotif
    }
    refTx <- sprintf("tx_%s.RData", refGenome)
    if (!file.exists(refTx)) {
        ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    load(refTx)
    requiredObj <- c("geneXID", "seqBySec", "txBySec")
    if (!all(requiredObj %in% ls())) {
        ss <- sprintf("Mandatory transcript objects not found.")
        ss <- sprintf("%s\nNeed all of the following: %s",
                      ss, paste0(requiredObj, collapse = ", "))
        ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
                      ss, refGenome)
        stop(ss)
    }
    geneXID <- get("geneXID")
    seqBySec <- get("seqBySec")
    txBySec <- get("txBySec")
    if (is.null(filter)) {
        filter <- c("5'UTR", "CDS", "3'UTR")
    }
    sel <- which(names(txBySec) %in% filter)
    locus.list <- list()
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(sel), style = 3, width = 60)
    for (i in 1:length(sel)) {
        if (showPb) setTxtProgressBar(pb, i)
        if (!IsEmptyChar(seqBySec[[sel[i]]])) {
# Disable warnings to get rid of messages of the form "Each of the
# 2 combined objects has sequence levels not in the other" when
# appending to GRanges object. Ugly! Is there a cleaner solution?
# Works for now, fix later.
            options(warn = -1)
            gr <- GRanges()
            for (j in 1:length(motif)) {
                m <- vmatchPattern(motif[j],
                                   seqBySec[[sel[i]]],
                                   max.mismatch = maxMM)
                m <- unlist(m)
                grMotif <- GRanges(names(m), m,
                                   id = rep(motif[j], length(m)))
                gr <- append(gr, grMotif)
            }
            options(warn = 0)
            locus <- GetLocus.MapFromTranscripts(
                gr,
                txBySec[[sel[i]]],
                seqBySec[[sel[i]]],
                names(txBySec)[sel[i]],
                geneXID)
            locus.list[[length(locus.list) + 1]] <- locus
        }
    }
    if (showPb) close(pb)
    names(locus.list) <- names(seqBySec)[sel]
    obj <- new("txLoc",
               loci = locus.list,
               id = "motifs",
               refGenome = refGenome,
               version = as.character(Sys.Date()))
    return(obj)
}

#' Fold sequences.
#'
#' Fold sequences. See 'Details'.
#'
#' The function takes a \code{dataframe}, extracts sequences from a
#' column specified by \code{colSeq}, and predicts secondary structures
#' using RNAfold \url{http://rna.tbi.univie.ac.at/}.
#' An optional column containing sequence IDs may be specified by
#' \code{colId}.
#' The function returns a \code{dataframe} with three columns:
#' \enumerate{
#' \item Column 1: Sequence ID
#' \item Column 2: Length of the sequence (in nt)
#' \item Column 3: Mean free energy (MFE)
#' }
#'
#' @param data A \code{dataframe} object. See 'Details'.
#' @param colSeq An integer scalar; specifies the column in
#' \code{data} containing the sequences.
#' @param colId An integer scalar; specifies the column in
#' \code{data} containing sequence IDs; if \code{NULL} IDs
#' are generated automatically; default is \code{NULL}.
#' 
#' @return A \code{dataframe} object. See 'Details'.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @import Biostrings GenomicRanges IRanges
#' 
GetMFE <- function(data, colSeq, colId = NULL) {
    sq <- DNAStringSet(data[, colSeq])
    if (is.null(colId) || length(data[, colId]) != length(sq)) {
        id <- sprintf("seq%i", seq(1, length(sq)))
    } else {
        id <- data[, colId]
    }
    names(sq) <- id
    writeXStringSet(sq, filepath = "tmp.fa", format = "fasta")
    cmd <- sprintf("RNAfold --noPS < tmp.fa > tmp.dbn")
    system(sprintf(cmd))
    str <- ReadDBN("tmp.dbn")
    file.remove(c("tmp.fa", "tmp.dbn"))
    return(str)
}
