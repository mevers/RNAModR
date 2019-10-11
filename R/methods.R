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
#' @keywords internal
#'
#' @importFrom S4Vectors DataFrame
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
#' @param id A \code{character} string; specifies a name for loci from
#' \code{gr}; if \code{NULL} then \code{id = ""}; default is \code{NULL}.
#' @param refGenome A \code{character} string; specifies a specific reference
#' genome assembly version based on which the matching transcriptome is loaded;
#' default is \code{"hg38"}.
#' @param ignore.strand A \code{logical} scalar; if \code{TRUE} strand
#' information is ignored when mapping genome coordinates to transcript
#' coordinates; default is \code{FALSE}.
#' @param showPb A \code{logical} scalar; if \code{TRUE} show a progress bar;
#' default is \code{TRUE}.
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
                     ignore.strand = FALSE,
                     showPb = TRUE) {

    # Sanity check
    CheckClass(gr, "GRanges")

    # Load transcriptome data
    seqBySec <- txBySec <- geneXID <- NULL
    LoadRefTx(refGenome)

    # Map coordinates to transcript
    loci <- SmartMap.ToTx(
        gr, txBySec, seqBySec, geneXID,
        ignore.strand = ignore.strand,
        showPb = showPb)

    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = if (is.null(id)) "" else id,
        refGenome = refGenome,
        version = as.character(Sys.Date()))

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
#' }
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
#' @param id A \code{character} string; identifier for null sites; if
#' \code{NULL} then \code{id = "null"}; default is \code{NULL}.
#' @param method A \code{character} string; specifies the method used to the
#' generate null distribution; if \code{method == "ntAbund"} the position of
#' all nucleotides specified by \code{nt} will be used as null sites; if
#' \code{method == "perm"} (DEPRECATED, see 'Details') then null sites will be
#' generated by uniform-randomly shuffling candidatepositions from \code{locus}
#' within the corresponding transcript region; default is \code{"ntAbund"}.
#' @param nt A single \code{character}; if \code{method == "ntAbund"}, use
#' \code{nt} to derive distribution of null sites.
#' @param showPb A \code{logical} scalar; if \code{TRUE} show a progress
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

    # Load transcriptome data
    seqBySec <- txBySec <- geneXID <- NULL
    if (method == "ntAbund") LoadRefTx(refGenome)

    # Activate progressbar if `showPb == TRUE`
    if (showPb == TRUE)
        pb <- txtProgressBar(
            max = length(GetRegions(txLoc)), style = 3, width = 60)

    loci <- Map(
        function(i, region) {

          # Update progress bar
          if (showPb) setTxtProgressBar(pb, i)

          # Loci of positive sites for region i
          locus.pos <- GetLoci(txLoc)[[i]]

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
          hits <- mapFromTranscripts(gr, region)

          # Make sure that we have one-to-one mapping between gr and hits
          gr <- gr[mcols(hits)$xHits]

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
            setNames(DataFrame(gr), "locus_in_tx_region"),
            setNames(
              DataFrame(granges(hits, use.names = FALSE, use.mcols = FALSE)),
              "locus_in_genome"),
            setNames(DataFrame(
              rep(0, length(gr)),
              paste0(seqnames(hits), ":", start(hits), "-", end(hits))),
              c("score", "id")),
            locus.pos[match(
              seqnames(gr), locus.pos$tx_refseq), -(1:4)])

      },
      setNames(seq_along(GetRegions(txLoc)), GetRegions(txLoc)),
      txBySec[GetRegions(txLoc)])

    if (showPb) close(pb)

    # Return `txLoc` object
    new("txLoc",
        loci = loci,
        id = id,
        refGenome = refGenome,
        version = as.character(Sys.Date()))

}


#' Calculate GC content within window around loci from a \code{txLoc}
#' object.
#'
#' Calculate GC content within window around loci from a \code{txLoc}
#' object. See 'Details'.
#'
#' The function calculates the GC content within a window around every locus
#' from a \code{txLoc} object. The window is defined by extending the position
#' of every \code{txLoc} site upstream and downstream by \code{flank}
#' nucleotides (if possible). This is the low-level function that gets called
#' from within \code{PlotGC}. There is no guarantee that this function will
#' be exported in future versions of RNAModR.
#'
#' @param txLoc A \code{txLoc} object.
#' @param flank An \code{integer} scalar; see 'Details'.
#'
#' @return A \code{list} of \code{DataFrame} objects, each with the following
#' columns:
#' \enumerate{
#'     \item \code{id}: A \code{character} vector
#'     \item \code{GC_tx_region}: A \code{numeric} vector
#'     \item \code{GC_window}: A \code{numeric} vector
#'     \item \code{sequence_window}: A \code{DNAStringSet} object
#' }
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom Biostrings letterFrequency
#'
#' @export
GetGC <- function(txLoc, flank = 10) {

    # Sanity check
    CheckClass(txLoc, "txLoc")

    # Create list of `Data.Frame`s
    lapply(GetLoci(txLoc), function(loci) {
        seqRegion <- loci$tx_region_sequence

        # GC content of transcript region
        nucFreq <- letterFrequency(seqRegion, letters = c("A", "C", "G", "T"))
        GCRegion <- rowSums(nucFreq[, 2:3]) / width(seqRegion)

        # Calculate windows for every site and make sure that we're still
        # within transcript region boundaries
        windows <- lapply(
            start(loci$locus_in_tx_region),
            function(x) x + c(-1, +1) * flank)
        windows <- Map(
            function(x, region_width) {
                x[x <= 0] <- 1
                x[x > region_width] <- region_width
                x
            },
            windows, loci$tx_region_width)

        # Extract subsequences
        # We use `subseq` instead of `substr` as it is more performant
        seqWindow <- as(Map(
            function(seq, x) subseq(seq, start = x[1], end = x[2]),
            seqRegion, windows), "DNAStringSet")

        # GC content of window
        nucFreq <- letterFrequency(seqWindow, letters = c("A", "C", "G", "T"))
        GCWindow <- rowSums(nucFreq[, 2:3]) / width(seqWindow)

        # Store in `DataFrame`
        cbind(
            setNames(
                DataFrame(loci$id, GCRegion, GCWindow, seqWindow),
                c("id", "GC_tx_region", "GC_window", "sequence_window")))

    })

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

    # Load transcriptome data
    seqBySec <- txBySec <- geneXID <- NULL
    LoadRefTx(refGenome)

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
    seqBySec <- txBySec <- geneXID <- NULL
    LoadRefTx(refGenome)

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
      junct <- unlist(as(Map(
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
#' Based on a set of query sequences \code{motif}, extract the loci of matching
#' sequences within different regions of a reference transcriptome
#' \code{refGenome}. The function returns a \code{txLoc} object including
#' transcriptomic and genomic loci of the query motif(s). The maximum number
#' of mismatches allowed in the motif search can be adjusted through
#' \code{maxMM}.
#' Note that the motif search may take a few minutes, depending on the number
#' of motif(s), hardware, and size of the transcriptome.
#'
#' @param motif A \code{character} vector; specifies the query sequences that
#' will be matched against the reference transcriptome.
#' @param refGenome A character string; specifies the reference genome version;
#' default is \code{"hg38"}.
#' @param id A character string; specifies an identifier for \code{motif};
#' default is \code{NULL}
#' @param maxMM An integer scalar; specifies the maximum number of mismatches
#' that are allowed during the motif matching; default is 0.
#' @param showPb A logical scalar; if \code{TRUE} show a progress bar; default
#' is \code{TRUE}.
#'
#' @return A \code{txLoc} object.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @import GenomicRanges IRanges
#' @importFrom Biostrings vmatchPattern
#'
#' @examples
#' \dontrun{
#' PAS <- GetMotifLoc(id = "PAS")
#' PAS <- FilterTxLoc(PAS, c("3'UTR", "CDS", "5'UTR"))
#' PlotSpatialDistribution(PAS, absolute = FALSE)
#' }
#'
#' @export
GetMotifLoc <- function(motif = NULL,
                        refGenome = "hg38",
                        id = NULL,
                        maxMM = 0,
                        showPb = TRUE) {


    # If motif == NULL, motif = PASmotif
    PASmotif <- c("AATAAA", "ATTAAA", "AGTAAA",
                  "TATAAA", "AAGAAA", "AATACA",
                  "AATATA", "CATAAA", "AATGAA",
                  "GATAAA", "ACTAAA", "AATAGA")
    if (is.null(motif)) {
      motif <- PASmotif
    }

    # If id == NULL, id = ""
    if (is.null(id)) id <- "motif"

    # Load transcriptome data
    seqBySec <- txBySec <- geneXID <- NULL
    LoadRefTx(refGenome)
    
    if (showPb == TRUE)
        pb <- txtProgressBar(max = length(seqBySec), style = 3, width = 60)

    lst <- list()
    for (i in 1:length(seqBySec)) {

        if (showPb) setTxtProgressBar(pb, i)

        # Extract position of all motif sites inside every sequence
        gr <- do.call(c, lapply(motif, function(x) {
            gr <- GRanges(vmatchPattern(x, seqBySec[[i]], max.mismatch = maxMM))
            mcols(gr)$motif <- x
            gr
        }))

        # Map sites from transcriptomic to genomic coordinates
        hits <- mapFromTranscripts(gr, txBySec[[i]])

        # Make sure that we have one-to-one mapping between gr and hits
        gr <- gr[mcols(hits)$xHits]

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
        lst[[length(lst) + 1]] <- cbind(
            setNames(DataFrame(
                granges(gr, use.names = F, use.mcols = F)),
                "locus_in_tx_region"),
            setNames(DataFrame(
                granges(hits, use.names = F, use.mcols = F)),
                "locus_in_genome"),
            setNames(DataFrame(
                rep(0, length(gr)),
                paste0(
                    id, "|", mcols(gr)$motif, "|",
                    seqnames(hits),":", start(hits), "-", end(hits))),
                c("score", "id")),
            setNames(DataFrame(
               rep(names(seqBySec)[i], length(gr)),
               width(seqBySec[[i]][match(
                   seqnames(gr), names(seqBySec[[i]]))]),
               seqBySec[[i]][match(
                   seqnames(gr), names(seqBySec[[i]]))]),
               c("tx_region", "tx_region_width", "tx_region_sequence")),
            geneXID[match(
                seqnames(gr), geneXID$tx_refseq), ])

    }

    if (showPb) close(pb)

    names(lst) <- names(seqBySec)

    new("txLoc",
               loci = lst,
               id = id,
               refGenome = refGenome,
               version = as.character(Sys.Date()))

}
