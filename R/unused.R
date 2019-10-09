##' Plot overlap of sites.
##'
##' Plot overlap of sites from two \code{txLoc} object.
##' See 'Details'.
##'
##' The function plots one or multiple Venn diagrams denoting the
##' spatial overlap between entries from two \code{txLoc} objects.
##' Two features are defined as overlapping, if they overlap by
##' at least one nucleotide. Overlaps are determined using the
##' function \code{GenomicRanges::countOverlaps}.
##'
##'
##' @param txLoc1 A \code{txLoc} object.
##' @param txLoc2 A \code{txLoc} object.
##'
##' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
##'
##' @import GenomicRanges IRanges
##' @importFrom gplots venn
#PlotOverlap <- function(loc1, loc2) {
#    CheckClassTxLocConsistency(loc1, loc2)
#    id1 <- GetId(loc1)
#    id2 <- GetId(loc2)
#    gr1 <- TxLoc2GRangesList(loc1)
#    gr2 <- TxLoc2GRangesList(loc2)
#    if (length(gr1) < 4) {
#        par(mfrow = c(1, length(gr1)))
#    } else {
#        par(mfrow = c(ceiling(length(gr1) / 2), 2))
#    }
#    for (i in 1:length(gr1)) {
#        # Supress warnings of sequences in gr1 not being in gr2
#        m <- suppressWarnings(countOverlaps(gr1[[i]], gr2[[i]]))
#        overlap <- length(m[m > 0])
#        grps <- list(
#            seq(1, length(gr1[[i]])),
#            seq(length(gr1[[i]]) - overlap + 1, length.out = length(gr2[[i]])))
#        names(grps) <- c(sprintf("%s (%3.2f%%)",
#                                 id1, overlap / length(gr1[[i]]) * 100),
#                         sprintf("%s (%3.2f%%)",
#                                 id2, overlap / length(gr2[[i]]) * 100))
#        venn(grps)
#        mtext(names(gr1)[i])
#    }
#}
#

##' Read a DBN file.
##'
##' Read a DBN file. See 'Details'.
##'
##' The function reads in a DBN (dot-bracket) structure file, and
##' returns a \code{dataframe} with the following data columns:
##' \enumerate{
##' \item Column 1: Sequence ID
##' \item Column 2: Length of the sequence (in nt)
##' \item Column 3: Mean free energy (MFE)
##' }
##' 
##' @param file A character string; specifies the input DBN file.
##'
##' @return A \code{dataframe} object. See 'Details'.
##'
##' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#ReadDBN <- function(file) {
#    if (!file.exists(file)) {
#        ss <- sprintf("Could not open %s.", file)
#        stop(ss)
#    }
#    d <- as.data.frame(matrix(readLines(file), ncol = 3, byrow = TRUE),
#                       stringsAsFactors = FALSE)
#    d[, 1] <- gsub("^>", "", d[, 1])
#    d[, 2] <- nchar(d[, 2])
#    d[, 3] <- gsub("^[\\(\\.\\)]+\\s", "", d[, 3])
#    d[ ,3] <- as.numeric(gsub("[\\(\\)]", "", d[, 3]))
#    colnames(d) <- c("id", "siteSeqLength", "siteMFE")
#    return(d)
#}



##' Fold sequences.
##'
##' Fold sequences. See 'Details'.
##'
##' The function takes a \code{dataframe}, extracts sequences from a
##' column specified by \code{colSeq}, and predicts secondary structures
##' using RNAfold \url{http://rna.tbi.univie.ac.at/}.
##' An optional column containing sequence IDs may be specified by
##' \code{colId}.
##' The function returns a \code{dataframe} with three columns:
##' \enumerate{
##' \item Column 1: Sequence ID
##' \item Column 2: Length of the sequence (in nt)
##' \item Column 3: Mean free energy (MFE)
##' }
##'
##' @param data A \code{dataframe} object. See 'Details'.
##' @param colSeq An integer scalar; specifies the column in
##' \code{data} containing the sequences.
##' @param colId An integer scalar; specifies the column in
##' \code{data} containing sequence IDs; if \code{NULL} IDs
##' are generated automatically; default is \code{NULL}.
##' 
##' @return A \code{dataframe} object. See 'Details'.
##'
##' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
##'
##' @import Biostrings GenomicRanges IRanges
#GetMFE <- function(data, colSeq, colId = NULL) {
#    sq <- DNAStringSet(data[, colSeq])
#    if (is.null(colId) || length(data[, colId]) != length(sq)) {
#        id <- sprintf("seq%i", seq(1, length(sq)))
#    } else {
#        id <- data[, colId]
#    }
#    names(sq) <- id
#    writeXStringSet(sq, filepath = "tmp.fa", format = "fasta")
#    cmd <- sprintf("RNAfold --noPS < tmp.fa > tmp.dbn")
#    system(sprintf(cmd))
#    str <- ReadDBN("tmp.dbn")
#    file.remove(c("tmp.fa", "tmp.dbn"))
#    return(str)
#}
#


##  #' Plot relative distance distribution of loci from \code{txLoc}
##  #' object to splice sites.
##  #'
##  #' Plot relative distance distribution of loci from \code{txLoc}
##  #' object to splice sites.
##  #'
##  #' @param locus A \code{txLoc} object.
##  #' @param ss A list of two \code{GRangesList} objects.
##  #' @param flank An integer scalar; specifies the absolute maximum
##  #' relative distance used as a cutoff; default is 1000.
##  #' @param binWidth An integer scalar; specifies the spatial width
##  #' by which distances will be binned; default is 20.
##  #'
##  #' @export
##  PlotRelDistSSDistribution <- function(locus, ss, flank = 1000, binWidth = 20) {
##      dist <- GetRelDistSS(locus, ss, flank)
##      breaks <- seq(-flank, flank, by = binWidth)
##      bwString <- sprintf("bw = %i nt", binWidth)
##      h0 <- hist(dist, breaks = breaks, plot = FALSE)
##      plot(h0$mids, h0$counts,
##           type = "s",
##           lwd = 2,
##           xlab = "Relative distance from 1st CDS splice site [nt]",
##           ylab = "Abundance",
##           xlim = c(-1000, 1000),
##           font.main = 1)
##      CIFromBS <- EstimateCIFromBS(dist,
##                                   breaks = breaks,
##                                   nBS = 5000)
##      x1 <- c(CIFromBS$x[1],
##              rep(CIFromBS$x[-1], each = 2),
##              CIFromBS$x[length(CIFromBS$x)])
##      CI <- cbind(c(x1,
##                    rev(x1)),
##                  c(rep(CIFromBS$y.low, each = 2),
##                    rev(rep(CIFromBS$y.high, each = 2))))
##      polygon(CI[, 1], CI[, 2], col = rgb(1, 0, 0, 0.2),
##              lwd = 1, border = NA, lty = 1)
##      # Loess smoothing of boostrap CI
##      lines(lowess(CIFromBS$x, CIFromBS$y.low, f = 1/5),
##            col = "red", lty = 2, lwd = 1)
##      lines(lowess(CIFromBS$x, CIFromBS$y.high, f = 1/5),
##            col = "red", lty = 2, lwd = 1)
##      legend("topleft",
##             c(sprintf("Abundance (%s)", bwString),
##               "95%CI (empirical bootstrap)",
##               "Lowess-smoothed 95%CI"),
##             lwd = c(2, 5, 1),
##             col = c("black", rgb(1, 0, 0, 0.2), "red"),
##             lty = c(1, 1, 2),
##             bty = "n")
##  }
##  
##  #' Plot relative distance enrichment of sites from \code{locPos}
##  #' and \code{locPos} to splice sites.
##  #' 
##  #' Plot relative distance enrichment of sites from \code{locPos}
##  #' and \code{locPos} to splice sites.
##  #'
##  #' @param locPos A \code{txLoc} object; specifies the positive
##  #' sites.
##  #' @param locNeg A \code{txLoc} object; specifies the negative
##  #' sites.
##  #' @param ss A list of two \code{GRangesList} objects.
##  #' @param flank An integer scalar; specifies the absolute maximum
##  #' relative distance used as a cutoff; default is 1000.
##  #' @param binWidth An integer scalar; specifies the spatial width
##  #' by which distances will be binned; default is 20.
##  #'
##  #' @export
##  PlotRelDistSSEnrichment <- function(locPos, locNeg, ss, flank = 1000, binWidth = 20) {
##      idPos <- GetId(locPos)
##      idNeg <- GetId(locNeg)
##      distPos <- GetRelDistSS(locPos, ss, flank)
##      distNeg <- GetRelDistSS(locNeg, ss, flank)
##      breaks <- seq(-flank, flank, by = binWidth)
##      ctsPos <- table(cut(distPos, breaks = breaks))
##      ctsNeg <- table(cut(distNeg, breaks = breaks))
##      ctsMat <- as.matrix(rbind(ctsPos, ctsNeg))
##      rownames(ctsMat) <- c("pos", "neg")
##      title <- sprintf("N(%s) = %i, N(%s) = %i\n(bw = %i nt)",
##                       idPos, sum(ctsPos),
##                       idNeg, sum(ctsNeg),
##                       binWidth)
##      tmp <- PlotEnrichment.Generic(ctsMat,
##                                    title = title,
##                                    x.las = 2, x.cex = 0.8, x.padj = 0.8)
##  }


##  #' Plot distribution of relative distances to exon-exon junctions.
##  #' 
##  #' Plot distribution of relative distances to exon-exon junctions.
##  #'
##  #' @param dist A list of integer vectors.
##  #' @param flank An integer scalar; specifies the absolute maximum
##  #' relative distance used as a cutoff; default is 1000.
##  #' @param binWidth An integer scalar; specifies the spatial width
##  #' by which distances will be binned; default is 20.
##  #' @param doBootstrap A logical scalar; if \code{YES} calculate
##  #' 95% CI based on empirical bootstrap of sites within transcript
##  #' region; default is \code{TRUE}.
##  #' 
##  #' @export
##  PlotRelDistEEJDistribution <- function(dist,
##                                         flank = 1000,
##                                         binWidth = 20,
##                                         doBootstrap = TRUE) {
##      if (length(dist) < 4) {
##          par(mfrow = c(1, length(dist)))
##      } else {
##          par(mfrow = c(ceiling(length(dist) / 2), 2))
##      }
##      breaks <- seq(-flank, flank, by = binWidth)
##      bwString <- sprintf("bw = %3.2f", binWidth)
##      for (i in 1:length(dist)) {
##          if (flank > 0) {
##              dist[[i]] <- dist[[i]][abs(dist[[i]]) <= flank]
##          }
##          title <- sprintf("%s (N=%i)",
##                           names(dist)[i],
##                           length(dist[[i]]))
##          xlab <- "Relative distance to exon-exon junction [nt]"
##          PlotAbundance.generic(dist[[i]],
##                                xmin = -flank, xmax = flank,
##                                binWidth = binWidth,
##                                title = title,
##                                xlab = xlab,
##                                doBootstrap = doBootstrap)
##      }
##  }


##  #' Perform and plot enrichment analysis of relative distances of
##  #' "positive" and "negative" control sites to exon-exon junctions.
##  #' 
##  #' Perform and plot enrichment analysis of relative distances of
##  #' positive and negative control sites to exon-exon junctions.
##  #'
##  #' @param distPos A list of integer vectors; specifies distances
##  #' of "positive" sites to exon-exon junctions.
##  #' @param distNeg A list of integer vectors; specifies distances
##  #' of "negative" sites to exon-exon junctions.
##  #' @param flank An integer scalar; specifies the absolute maximum
##  #' relative distance used as a cutoff; default is 1000.
##  #' @param binWidth An integer scalar; specifies the spatial width
##  #' by which distances will be binned; default is 20.
##  #' 
##  #' @export
##  PlotRelDistEEJEnrichment <- function(distPos,
##                                       distNeg,
##                                       flank = 1000,
##                                       binWidth = 20) {
##      if (length(distPos) < 4) {
##          par(mfrow = c(1, length(distPos)))
##      } else {
##          par(mfrow = c(ceiling(length(distPos) / 2), 2))
##      }
##      breaks <- seq(-flank, flank, by = binWidth)
##      bwString <- sprintf("bw = %3.2f", binWidth)
##      for (i in 1:length(distPos)) {
##          if (flank > 0) {
##              distPos[[i]] <- distPos[[i]][abs(distPos[[i]]) <= flank]
##              distNeg[[i]] <- distNeg[[i]][abs(distNeg[[i]]) <= flank]
##          }
##          ctsPos <- table(cut(distPos[[i]], breaks = breaks))
##          ctsNeg <- table(cut(distNeg[[i]], breaks = breaks))
##          ctsMat <- as.matrix(rbind(ctsPos, ctsNeg))
##          rownames(ctsMat) <- c("pos", "neg")
##          title <- sprintf("N(%s) = %i, N(%s) = %i\n(bw = %i nt)",
##                           "pos", sum(ctsPos),
##                           "neg", sum(ctsNeg),
##                           binWidth)
##          tmp <- PlotEnrichment.Generic(ctsMat,
##                                        title = title,
##                                        x.las = 2, x.cex = 0.8, x.padj = 0.8)
##      }
##  }


## #' Calculate relative distances between loci from \code{txLoc}
## #' object to splice sites.
## #'
## #' Calculate relative distances between loci from \code{txLoc}
## #' and splice sites. See 'Details'.
## #' This is a low-level function that is being called from
## #' \code{PlotRelDistSSDistribution} and
## #' \code{PlotRelDistSSEnrichment}.
## #'
## #' The function calculates relative distances between 5' splice
## #' sites from \code{ss} and loci from \code{locus}. Distances
## #' are taken to be negative (positive) if a locus from
## #' \code{txLoc} is upstream (downstream) of a splice site.
## #' Distances flank <= d <= flank are considered. The fun
## #'
## #' @param locus A \code{txLoc} object.
## #' @param ss A list of two \code{GRangesList} objects.
## #' @param flank An integer scalar; specifies the absolute maximum
## #' relative distance used as a cutoff; default is 1000.
## #'
## #' @return An integer vector.
## #'
## #' @keywords internal
## #'
## #' @export
## GetRelDistSS <- function(locus,
##                          ss,
##                          flank = 1000) {
##     CheckClass(locus, "txLoc")
##     CheckClass(ss, "list", "GRangesList")
##     grLoc <- unlist(
##         TxLoc2GRangesList(locus,
##                           filter = c("5'UTR", "CDS", "3'UTR")))
##     # Select 5' splice sites in CDS
##     grSS <- unlist(ss[[grep("5p", names(ss))]])
##     grSS <- grSS[grep("(CDS|coding)",
##                       grSS$section,
##                       ignore.case = TRUE)]
##     # Filter genes with ss _and_ modification site
##     genes <- intersect(grSS$gene, grLoc$gene)
##     grSS <- grSS[which(grSS$gene %in% genes)]
##     grLoc <- grLoc[which(grLoc$gene %in% genes)]
##     # Select first 5' splice site upstream of start codon in CDS
##     grSSPos <- grSS[which(strand(grSS) == "+")]
##     grSSNeg <- grSS[which(strand(grSS) == "-")]
##     dupIDPos <- duplicated(grSSPos$gene,
##                           fromLast = FALSE)
##     dupIDNeg <- duplicated(grSSNeg$gene,
##                           fromLast = TRUE)
##     grSSPos <- grSSPos[!dupIDPos]
##     grSSNeg <- grSSNeg[!dupIDNeg]
##     grSS <- sort(append(grSSPos, grSSNeg))
## #    rtracklayer::export(grSS, "ss5p_firstCDS.bed")
##     # Get distances
##     d <- as.data.frame(distanceToNearest(grLoc, grSS))
##     idx <- d$subjectHits
##     dist <- ifelse(end(grLoc) > start(grSS[idx]), d$distance, -d$distance)
##     dist <- ifelse(strand(grLoc) == "+", dist, -dist)
##     dist <- dist[abs(dist) <= flank]
##     return(dist)
## }


##  #' Calculate mutual distances between entries for every transcript
##  #' region from two txLoc objects.
##  #'
##  #' Calculate mutual distances between entries for every transcript
##  #' region from two txLoc objects.
##  #'
##  #' The function calculates relative distances between loci from
##  #' \code{loc1} and \code{loc2} within the same transcript section.
##  #' Positive distances refer to loc1 > loc2, negative distances to
##  #' loc1 < loc2. By default only the smallest distance for every
##  #' locus from \code{loc1} is kept (\code{method = "nearest"}).
##  #' By default relative distances correspond to distances between
##  #' start coordinates (\code{refPoint = "ss"}). 
##  #' 
##  #' @param loc1 A \code{txLoc} object.
##  #' @param loc2 A \code{txLoc} object.
##  #' @param filter A character vector; only consider transcript sections
##  #' specified in \code{filter}; if \code{NULL} consider all sections.
##  #' @param method A character strings; specifies the distance
##  #' calculation method; possible arguments are "all", "nearest";
##  #' default is "nearest".
##  #' @param refPoint A character string; specifies the reference point
##  #' for calculating relative distances; possible arguments are "ss"
##  #' (start-start), "mm" (midpoint-midpoint), "se" (start-end), "es"
##  #' (end-start), "ee" (end-end); default is "ss".
##  #'
##  #' @return List of upper and lower 95% confidence interval
##  #' bounds for every bin value.
##  #' 
##  #' @export
##  GetRelDist <- function(loc1,
##                         loc2,
##                         filter = NULL,
##                         method = c("nearest", "all"),
##                         refPoint = c("ss", "mm", "se", "es", "ee")) {
##      CheckClass(loc1, "txLoc")
##      CheckClass(loc2, "txLoc")
##      CheckClassTxLocRef(loc1, loc2)
##      refPoint <- match.arg(refPoint)
##      method <- match.arg(method)
##      id1 <- GetId(loc1)
##      id2 <- GetId(loc2)
##      refGenome <- GetRef(loc1)
##      sec <- intersect(names(GetLoci(loc1)), names(GetLoci(loc2)))
##      loc1 <- FilterTxLoc(loc1, sec)
##      loc2 <- FilterTxLoc(loc2, sec)
##      loc1 <- GetLoci(loc1)
##      loc2 <- GetLoci(loc2)
##      dist.list <- list()
##      for (i in 1:length(loc1)) {
##          txID <- intersect(loc1[[i]]$REFSEQ, loc2[[i]]$REFSEQ)
##          if (length(txID) == 0) {
##              dist.list[[length(dist.list)+1]] <- 0
##          } else {
##              dist <- vector()
##              for (j in 1:length(txID)) {
##                  loc1.sel <- loc1[[i]][which(loc1[[i]]$REFSEQ == txID[j]), ]
##                  loc2.sel <- loc2[[i]][which(loc2[[i]]$REFSEQ == txID[j]), ]
##                  distMat <- switch(refPoint,
##                                    "ss" = outer(loc1.sel$TXSTART,
##                                        loc2.sel$TXSTART,
##                                        "-"),
##                                    "mm" = outer((loc1.sel$TXSTART + loc1.sel$TXEND) / 2,
##                                        (loc2.sel$TXSTART + loc2.sel$TXEND) / 2,
##                                        "-"),
##                                    "se" = outer(loc1.sel$TXSTART,
##                                        loc2.sel$TXEND,
##                                        "-"),
##                                    "es" = outer(loc1.sel$TXEND,
##                                        loc2.sel$TXSTART,
##                                        "-"),
##                                    "ee" = outer(loc1.sel$TXEND,
##                                        loc2.sel$TXEND,
##                                        "-"))
##                  if (method == "nearest") {
##                      dist <- c(dist, apply(distMat, 1, function(x) x[which.min(abs(x))]))
##                  } else {
##                      dist <- c(dist, as.vector(distMat))
##                  }
##              }
##              dist.list[[length(dist.list) + 1]] <- dist
##          }
##      }
##      names(dist.list) <- names(loc1)
##      return(dist.list)
##  }


## #' Get 5'/3' splice sites from transcriptome.
## #'
## #' Get 5'/3' splice sites from transcriptome. See 'Details'.
## #'
## #' This function extracts splice site positions from a reference
## #' transcriptome, and returns a list of two \code{GRangesList}
## #' objects of all 5' and 3' splice sites per gene.
## #'
## #' @param refGenome A character string; specifies a specific
## #' reference genome assembly version based on which the matching
## #' transcriptome is loaded; default is \code{"hg38"}.
## #' @param writeBED A logical scalar; if \code{TRUE} write 5'/3'
## #' splice sites to two BED files; default is \code{FALSE}.
## #'
## #' @import GenomicRanges IRanges
## #' @importFrom rtracklayer export
## #' 
## #' @return A list of two \code{GRangesList} objects. See 'Details'.
## #'
## #' @export
## GetSpliceSites <- function(refGenome = "hg38", writeBED = FALSE) {
##     refTx <- sprintf("tx_%s.RData", refGenome)
##     if (!file.exists(refTx)) {
##         ss <- sprintf("Reference transcriptome for %s not found.", refGenome)
##         ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
##                       ss, refGenome)
##         stop(ss)
##     }
##     load(refTx)
##     requiredObj <- c("geneXID", "seqBySec", "txBySec")
##     if (!all(requiredObj %in% ls())) {
##         ss <- sprintf("Mandatory transcript objects not found.")
##         ss <- sprintf("%s\nNeed all of the following: %s",
##                       ss, paste0(requiredObj, collapse = ", "))
##         ss <- sprintf("%s\nRunning BuildTx(\"%s\") might fix that.",
##                       ss, refGenome)
##         stop(ss)
##     }
##     geneXID <- get("geneXID")
##     seqBySec <- get("seqBySec")
##     txBySec <- get("txBySec")
##     secWithSpliceSites <- c("5'UTR", "CDS", "3'UTR")
##     sel <- which(names(txBySec) %in% secWithSpliceSites)
##     ss5p.all <- GRanges()
##     ss3p.all <- GRanges()
##     for (i in 1:length(sel)) {
##         junct <- psetdiff(range(txBySec[[sel[i]]]), txBySec[[sel[i]]])
##         junct <- unlist(junct[elementLengths(junct) > 0])
##         ss5p <- GRanges(
##             seqnames(junct),
##             IRanges(ifelse(strand(junct) == "+",
##                            start(junct),
##                            end(junct)),
##                     ifelse(strand(junct) == "+",
##                            start(junct),
##                            end(junct))),
##             strand(junct),
##             type = "ss5p",
##             gene = names(junct),
##             section = names(txBySec)[sel[i]])
##         ss3p <- GRanges(
##             seqnames(junct),
##             IRanges(ifelse(strand(junct) == "+",
##                            end(junct),
##                            start(junct)),
##                     ifelse(strand(junct) == "+",
##                            end(junct),
##                            start(junct))),
##             strand(junct),
##             type = "ss3p",
##             gene = names(junct),
##             section = names(txBySec)[sel[i]])
##         ss5p.all <- append(ss5p.all, ss5p)
##         ss3p.all <- append(ss3p.all, ss3p)
##     }
##     ss5p.all <- sort(ss5p.all)
##     ss3p.all <- sort(ss3p.all)
##     ss5p.all$name <- paste(ss5p.all$type, ss5p.all$gene, ss5p.all$section, sep = "|")
##     ss3p.all$name <- paste(ss3p.all$type, ss3p.all$gene, ss3p.all$section, sep = "|")
##     if (writeBED) {
##         rtracklayer::export(ss5p.all, "ss5p.bed")
##         rtracklayer::export(ss3p.all, "ss3p.bed")
##     }
##     ret <- list(GenomicRanges::split(ss5p.all, ss5p.all$gene),
##                 GenomicRanges::split(ss3p.all, ss3p.all$gene))
##     names(ret) <- c("ss5p", "ss3p")
##     return(ret)
## }

