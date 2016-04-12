#' Plot length distribution of transcript sections.
#'
#' Plot length distribution of transcript sections from list of
#' GRangesList transcript features.
#'
#' @param txBySec A \code{list} of \code{GRangesList} objects;
#' specifies the list of transcript sections.
#' @param asLog A logical scalar; if \code{TRUE} plot length
#' distribution on a log scale; Default is \code{TRUE}.
#' @param printMetrics A character string; specifies which metrics
#' should be printed as part of the box labels; default is \code{"median"}.
#' @param filter A character vector; only plot transcript sections
#' specified in \code{filter}; if \code{NULL} plot all sections; default
#' is \code{NULL}.
#' @param outliers A logical scalar; if \code{TRUE} show outliers in
#' boxplot; default is \code{FALSE}.
#' @param ... Any additional parameters passed to \code{axis}.
#'
#' @export
PlotTxSecLength <- function(txBySec,
                            asLog = TRUE,
                            printMetrics = c("median", "mean"),
                            filter = NULL,
                            outliers = FALSE,
                            ...) {
    # Plot length distribution of transcript features.
    #
    # Args:
    #   txBySec: List of GRangesList transcript features.
    #   asLog: Plot length distribution on log scale. Default is TRUE.
    #   printMetrics: Output median and mad of length (TRUE),
    #          or mean and sd (FALSE). Default is TRUE.
    #   filter: Only plot transcript regions specified in filter.
    #   outliers: Show outliers in boxplot. Default is FALSE.
    #   ...: Additional parameters passed to x-axis.
    #
    # Returns:
    #   NULL
    par(mfrow = c(1, 1));
    len <- lapply(txBySec, function(x) sum(width(x)));
    cdsIdx <- grep("CDS", names(txBySec), ignore.case = TRUE);
    if (length(cdsIdx) > 0) {
        metaData <- slot(txBySec[[cdsIdx]], "metadata")$genomeInfo;
        refOrganism <- metaData[[grep("Organism", names(metaData))]];
        refGenome <- metaData[[grep("Genome", names(metaData))]];
        refSource <- metaData[[grep("Data source", names(metaData))]];
        len[["CDS_exon"]] <- IRanges::unlist(width(txBySec[[cdsIdx]]));
    } else {
        cat("No CDS entry in %s. Could not infer organism meta information",
            deparse(substitute(txBySec)));
    }
    intronIdx <- grep("(Intron|Intronic)",
                      names(txBySec),
                      ignore.case = TRUE);
    if (length(intronIdx) > 0) {
        len[[intronIdx]] <- IRanges::unlist(width(txBySec[[intronIdx]]));
    }
    if (!is.null(filter)) {
        len <- len[which(names(len) %in% filter)];
    }
    labels <- sprintf("%s\nmedian = %i nt\nmad = %i nt",
                      names(len),
                      round(sapply(len, stats::median)),
                      round(sapply(len, stats::mad)));
    if (printMetrics == "mean") {
        labels <- sprintf("%s\nmean = %i nt\nsd = %i nt",
                          names(len),
                          round(sapply(len, mean)),
                          round(sapply(len, stats::sd)));
    } else if (printMetrics == "none") {
        labels <- sprintf("%s", names(len));
    }
    ylab <- "Length [nt]";
    if (asLog) {
        len <- lapply(len,log10);
        ylab <- "log10 Length";
    }
    par(font.main = 1);
    boxplot(len, names = rep("", length(len)),
            xaxt = "n", ylab = ylab,
            main = sprintf("%s (version = %s, source = %s)",
                refOrganism,
                refGenome,
                refSource),
            font.main = 1,
            outline = outliers);
    axis(1, at = 1:length(len),
         labels = labels,
         padj = 1, cex.axis = 0.7, ...);
}


#' Plot length distribution of transcript sections.
#'
#' Plot length distribution of transcript sections.
#'
#' @param txBySec A \code{list} of \code{GRangesList} objects;
#' specifies the list of transcript sections.
#' @param asLog A logical scalar; if \code{TRUE} plot length
#' distribution on a log scale; Default is \code{TRUE}.
#' @param printMetrics A character string; specifies which metrics
#' should be printed as part of the box labels; default is \code{"median"}.
#' @param filter A character vector; only plot transcript sections
#' specified in \code{filter}; if \code{NULL} plot all sections; default
#' is \code{NULL}.
#' @param ... Any additional parameters passed to \code{axis}.
#'
#' @keywords internal
#'
#' @export
PlotTxSecLength.bean <- function(txBySec,
                                 asLog = TRUE,
                                 printMetrics = c("median", "mean"),
                                 filter = NULL,
                                 ...) {
    # Plot length distribution of transcript features.
    #
    # Args:
    #   txBySec: List of GRangesList transcript features.
    #   asLog: Plot length distribution on log scale. Default is TRUE.
    #   printMetrics: Output median and mad of length (TRUE),
    #          or mean and sd (FALSE). Default is TRUE.
    #   filter: Only plot transcript regions specified in filter.
    #   outliers: Show outliers in boxplot. Default is FALSE.
    #   ...: Additional parameters passed to x-axis.
    #
    # Returns:
    #   NULL
    par(mfrow = c(1, 1));
    printMetrics <- match.arg(printMetrics);
    len <- lapply(txBySec, function(x) sum(width(x)));
    cdsIdx <- grep("CDS", names(txBySec), ignore.case = TRUE);
    if (length(cdsIdx) > 0) {
        metaData <- slot(txBySec[[cdsIdx]], "metadata")$genomeInfo;
        refOrganism <- metaData[[grep("Organism", names(metaData))]];
        refGenome <- metaData[[grep("Genome", names(metaData))]];
        refSource <- metaData[[grep("Data source", names(metaData))]];
        len[["CDS_exon"]] <- IRanges::unlist(width(txBySec[[cdsIdx]]));
    } else {
        cat("No CDS entry in %s. Could not infer organism meta information",
            deparse(substitute(txBySec)));
    }
    intronIdx <- grep("(Intron|Intronic)",
                      names(txBySec),
                      ignore.case = TRUE);
    if (length(intronIdx) > 0) {
        len[[intronIdx]] <- IRanges::unlist(width(txBySec[[intronIdx]]));
    }
    if (!is.null(filter)) {
        len <- len[which(names(len) %in% filter)];
    }
    labels <- sprintf("%s\nmedian = %i nt\nmad = %i nt",
                      names(len),
                      round(sapply(len, stats::median)),
                      round(sapply(len, stats::mad)));
    if (printMetrics == "mean") {
        labels <- sprintf("%s\nmean = %i nt\nsd = %i nt",
                          names(len),
                          round(sapply(len, mean)),
                          round(sapply(len, stats::sd)));
    } else if (printMetrics == "none") {
        labels <- sprintf("%s", names(len));
    }
    ylab <- "Length [nt]";
    if (asLog) {
        len <- lapply(len,log10);
        ylab <- "log10 Length";
    }
    par(font.main = 1);
    col <- addAlpha(rainbow(length(len)));
    col <- split(cbind(col, rgb(0, 0, 0, 0.2)), col);
    beanplot(len,
             ll = 0.02,
             bw = "nrd0",
             border = NA,
             col = as.list(col),
             show.names = FALSE,
             ylab = ylab,
             main = sprintf("%s (version = %s, source = %s)",
                 refOrganism,
                 refGenome,
                 refSource),
             font.main = 1,
             method = "jitter");
    axis(1, at = 1:length(len),
         labels = labels,
         padj = 1, cex.axis = 0.7, ...);
}


#' Plot piechart of the number of loci in every transcript section.
#'
#' Plot piechart of the number of loci in every transcript section.
#'
#' @param locus A \code{txLoc} object.
#' @param filter A character vector; only consider transcript sections
#' specified in \code{filter}; if \code{NULL} consider all sections. 
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' PlotPieNumberOfLoci(posSites);
#' }
#'
#' @export
PlotPieNumberOfLoci <- function(locus, filter = NULL) {
    # Plot piechart of the mumber of loci in every transcript section.
    #
    # Args:
    #   locus: A txLoc object.
    #
    # Returns:
    #   NULL
    id <- GetId(locus);
    locus <- FilterTxLoc(locus, filter);
    N <- GetNumberOfLoci(locus);
    labels <- names(GetLoci(locus));
    percentage <- N / sum(N) * 100.0;
    labels <- sprintf("%s %2.1f%% (%i)", labels, percentage, N);
    pie(N, labels = labels, col = rainbow(length(labels)),
        main = sprintf("Distribution of %i %s sites across transcript sections",
            sum(N), id),
        font.main = 1);
}


#' Plot spatial distribution of loci from \code{txLoc} object.
#'
#' Plot spatial distribution of loci from \code{txLoc} object within every
#' transcript section.
#'
#' @param locus A \code{txLoc} object.
#' @param filter Only plot loci in transcript regions specified in filter.
#' @param nbreaks Number of spatial bins. Default is 100.
#' @param absolute Plot spatial distribution in absolute coordinates.
#' Default is \code{FALSE}.
#' @param binWidth Spatial bin width. Overrides \code{nbreaks} if not
#' \code{NULL}.
#' @param posMax If \code{absolute == TRUE}, show spatial distribution
#' within a window given by \code{posMax} from the 5'/3' position of
#' the transcript feature. Default is 1000 nt.
#' @param doBootstrap Calculate 95% CI based on empirical bootstrap of
#' sites within transcript region. Default is \code{TRUE}.
#' @param ... Additional parameters passed to plot.
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' PlotSpatialDistribution(posSites);
#' PlotSpatialDistribution(posSites,
#'                         absolute = TRUE,
#'                         filter = c("5'UTR", "CDS", "3'UTR"),
#'                         ylim = c(0, 200));
#' }
#' 
#' @export
PlotSpatialDistribution <- function(locus,
                                    filter = NULL,
                                    nbreaks = 100,
                                    absolute = FALSE,
                                    binWidth = NULL,
                                    posMax = 1000,
                                    doBootstrap = TRUE,
                                    ...) {
    # Plot the spatial distribution of features across different
    # transcript regions.
    #
    # Args:
    #   locus: List of dataframes with mapped features across different
    #          transcript regions.
    #   filter: Only plot loci in transcript regions specified in filter.
    #   nbreaks: Number of spatial bins.
    #   absolute: Plot as a function of absolute coordinates.
    #             Default is FALSE.
    #   binWidth: Spatial bin width. Overrides nbreaks if not NULL.
    #   posMax: If absolute == TRUE, plot up distribution up to
    #                   this distance. Default is 1000 nt.
    #   doBootstrap: Calculate bootstrap CI (only if absolute == FALSE).
    #                Default is TRUE.
    #   ...: Additional parameters passed to plot.
    #
    # Returns:
    #   NULL
    CheckClass(locus, "txLoc");
    id <- GetId(locus);
    refGenome <- GetRef(locus);
    locus <- GetLoci(locus);
    if (!is.null(filter)) {
        locus <- locus[which(names(locus) %in% filter)];
    }
    if (!absolute) {
        if (length(locus) < 4) {
            par(mfrow = c(1, length(locus)));
        } else {
            par(mfrow = c(ceiling(length(locus) / 2), 2));
        }
        breaks <- seq(0.0, 1.0, length.out = nbreaks);
        bw <- 1.0 / nbreaks;
#        if (!is.null(binWidth)) {
#            breaks <- seq(0.0, 1.0, by = binWidth);
#            bw <- binWidth;
#        }
        bwString <- sprintf("bw = %3.2f", bw);
    } else {
        par(mfrow = c(length(locus), 2));
        breaks <- seq(1, posMax, length.out = nbreaks);
        bw <- round((posMax - 1) / nbreaks);
#        if (!is.null(binWidth)) {
#            breaks <- seq(1, posMax, by = binWidth);
#            bw <- binWidth;
#        }
        bwString <- sprintf("bw = %i nt", bw);
    }
    for (i in 1:length(locus)) {
        # Store site positions
        #  (1) relative to 5'start, and
        #  (2) relative to 3'end.
        pos <- list("5p" = locus[[i]]$TXSTART,
                    "3p" = locus[[i]]$REGION_TXWIDTH - locus[[i]]$TXSTART + 1);
        if (!absolute) {
            pos <- lapply(pos, function(x) x / locus[[i]]$REGION_TXWIDTH);
            if (grepl("(5'UTR|UTR5|5pUTR|Promoter)", names(locus)[i],
                      ignore.case = TRUE)) {
                pos <- pos["3p"];
                xrange <- list(c(1.0, 0.0));
                xlab <- list("Relative position (relative to 3' end)");
            } else {
                pos <- pos["5p"];
                xrange <- list(c(0.0, 1.0));
                xlab <- list("Relative position (relative to 5' start)");
            }
        } else {
            pos <- lapply(pos, function(x) x[x <= posMax]);
            xrange <- list(c(1, posMax), c(posMax, 1));
            xlab <- list("Absolute position (relative to 5' start) [nt]",
                         "Absolute position (relative to 3' end) [nt]");
        }
        for (j in 1:length(pos)) {
            h0 <- hist(pos[[j]], breaks = breaks, plot = FALSE);
            plot(h0$mids, h0$counts,
                 type = "s",
                 lwd = 2,
                 xlab = xlab[[j]],
                 ylab = "Abundance",
                 xlim = xrange[[j]],
                 main = sprintf("%s in %s\nN = %i",
                     id,
                     names(locus)[i],
                     length(pos[[j]])),
                 font.main = 1,
                 ...);
            if (doBootstrap) {
                CIFromBS <- EstimateCIFromBS(pos[[j]],
                                             breaks = breaks,
                                             nBS = 5000);
                x1 <- c(CIFromBS$x[1],
                        rep(CIFromBS$x[-1], each = 2),
                        CIFromBS$x[length(CIFromBS$x)]);
                CI <- cbind(c(x1,
                              rev(x1)),
                            c(rep(CIFromBS$y.low, each = 2),
                              rev(rep(CIFromBS$y.high, each = 2))));
                polygon(CI[, 1], CI[, 2], col = rgb(1, 0, 0, 0.2),
                        lwd = 1, border = NA, lty = 1);
                # Loess smoothing of boostrap CI
                lines(lowess(CIFromBS$x, CIFromBS$y.low, f = 1/5),
                      col = "red", lty = 2, lwd = 1);
                lines(lowess(CIFromBS$x, CIFromBS$y.high, f = 1/5),
                      col = "red", lty = 2, lwd = 1);
                legend("topleft",
                       c(sprintf("Abundance (%s)", bwString),
                         "95%CI (empirical bootstrap)",
                         "Lowess-smoothed 95%CI"),
                       lwd = c(2, 5, 1),
                       col = c("black", rgb(1, 0, 0, 0.2), "red"),
                       lty = c(1, 1, 2),
                       bty = "n");
            }
        }
    }
}


#' Generic function to perform enrichment analysis and plot results.
#'
#' Generic function to perform enrichment analysis and plot results.
#' Enrichment/depletion is evaluated using (multiple) Fisher's exact test(s).
#' Multiple hypothesis testing correction is applied following the method of
#' Bejamini and Hochberg.
#' Note: This function should not be invoked directly by the user.
#'
#' @param mat The data matrix.
#' @param title An identifier for \code{mat}.
#' @param x.las Orientation of x-axis labels. Default is 1.
#' @param x.cex Scaling factor for x-axis labels. Default is 1.
#' @param x.padj Vertical adjustment of x-axis labels. Default is 1.
#' @param plotType Plot style. Default is "l".
#' @param reverseXaxis Reverse the order of values from \code{mat}.
#' @param withExtendedAxisLabel Print extended axis label.
#'
#' @keywords internal
#' 
#' @return A list of \code{fisher.test} return objects and \code{mat}
PlotEnrichment.Generic <- function(mat, title = "",
                                   x.las = 1, x.cex = 1, x.padj = 1,
                                   plotType = "l",
                                   reverseXaxis = FALSE,
                                   withExtendedAxisLabel = 0) {
    # Generic function to perform enrichment analysis and plot results.
    #
    # Args:
    #   mat: The data matrix.
    #   title: Plot title.
    #   x.las: Orientation of x-axis labels. Default is 1.
    #   x.cex: Scaling factor for x-axis labels. Default is 1.
    #   x.padj: Vertical adjustment of x-axis labels. Default is 1.
    #   plotType: Plot style. Default is "l".
    #   reverseXaxis: Reverse the order of values from mat.
    #   withExtendedAxisLabel: Print extended axis label.
    #
    # Returns:
    #   A list of fisher.test return objects.
    log10pval.limit <- -12;
    ft <- list();
    for (i in 1:ncol(mat)) {
        redMat <- cbind(mat[, i], rowSums(mat[, -i]));
        ft[[i]]<-fisher.test(redMat);
    }
    names(ft) <- colnames(mat);
    OR <- sapply(ft, function(x) x$estimate);
    OR[is.infinite(OR)] <- 1;
    OR <- log10(OR);
    names(OR) <- colnames(mat);
    pval <- p.adjust(sapply(ft, function(x) x$p.value), method="BH");
    names(pval) <- colnames(mat);
    pval[is.na(pval)] <- 1;
    pval <- log10(pval);
    pval.uncapped <- pval;
    pval[pval < log10pval.limit] <- log10pval.limit;
    CI <- sapply(ft, function(x) x$conf.int);
    CI[CI == 0] <- 1.e-12;
    CI[is.infinite(CI)] <- 1e12;
    CI <- log10(CI);
    if (reverseXaxis) {
        OR <- rev(OR);
        pval <- rev(pval);
        CI[1, ] <- rev(CI[1, ]);
        CI[2, ] <- rev(CI[2, ]);
    }
    #ymin <- floor(min(-abs(pval), -2.0));
    #ymax <- round(max(max(OR), 2.0));
    ymin <- log10pval.limit;
    ymax <- 1;
    # Plot log10(p-value)'s
    par(mar = c(7, 4, 4, 4) + 0.1);
    xrange <- 1.2 * length(pval);
    xmin <- -0.04 * xrange;
    xmin <- 0.2;
    xmax <- xrange + 0.04 * xrange;
    xmax <- 0.96 * xmax;
    mp <- barplot(pval,
                  xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax),
                  col = rgb(0.2, 0.2, 1, 0.1),
                  axes = FALSE, names = "",
                  border = NA,
                  main = title, font.main = 1);
    abline(h = -2, col = "blue", lty = 3, lwd = 2);
    if (withExtendedAxisLabel > 0) {
        if (withExtendedAxisLabel == 1) {
            labels = sprintf("%s\nOR=%4.3f,p=%4.3e",
                names(OR), 10^OR, 10^pval.uncapped);
        } else {
            labels = sprintf("%s\nN(%s)=%i\nN(%s)=%i\nOR=%4.3f\np=%4.3e",
                names(OR),
                rownames(mat)[1], mat[1, ],
                rownames(mat)[2], mat[2, ],
                10^OR, 10^pval.uncapped);
        }
    } else {
        labels = names(OR);
    }
    x <- mp;
    axis(1, at = x,
         labels = labels,
         las = x.las, padj = x.padj, cex.axis = x.cex);
    axis(4, at = seq(ymin, ymax));
    mtext("log10(p-value)", side = 4, line = 2.5);
    # Draw significance labels
    sig <- pval;
    sig[pval <= -4] <- "****";
    sig[pval <= -3 & pval > -4] <- "***";
    sig[pval <= -2 & pval > -3] <- "**";
    sig[pval <= -1.30103 & pval > -2] <- "*";
    sig[pval > -1.30103] <- "ns";
    text(x, y = 0.5, labels = sig, srt = 90, col = "blue");
    legend("bottomleft",
           c("Odds-ratio (OR)",
             "95% CI OR",
             "p-value"),
           lwd = c(2, 5, 5),
           col = c("red", rgb(1, 0, 0, 0.2), rgb(0.2, 0.2, 1, 0.1)),
           lty = c(1, 1, 1),
           bty = "n");
    # Plot odds-ratios
    par(new = TRUE);
    plot(x, OR, col = rgb(1, 0.2, 0.2, 1),
         lwd = 2, type = plotType,
         xlim = c(xmin, xmax),
         ylim = c(-1.0, 1.0),
         axes = FALSE, xlab = "", ylab = "");
    CI <- cbind(c(x, rev(x)),
                c(CI[1 ,], rev(CI[2, ])));
    polygon(CI[,1], CI[,2], col = rgb(1, 0, 0, 0.2),
            lwd = 1, border = NA, lty = 1);
    axis(2, at = seq(-1.0, 1.0));
    mtext("log10(OR)", side = 2, line = 2.5);
    abline(h = 0.0, col = "red", lty = 3, lwd = 2);
    return(list(ft = ft,
                mat = mat));
}


#' Perform transcript section enrichment analysis and plot results.
#'
#' Perform transcript section enrichment analysis and plot results.
#' Enrichment/depletion is evaluated using (multiple) Fisher's exact test(s).
#' Multiple hypothesis testing correction is applied following the method of
#' Bejamini and Hochberg.
#'
#' @param locPos A \code{txLoc} object. These should be the positive control sites.
#' @param locNeg A \code{txLoc} object. These should be the negative control sites.
#' @param filter Only plot loci in transcript regions specified in filter.  Default is NULL.
#' @param withExtendedAxisLabel Plot extended axis labels. Default is 2. See ??? for details.
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateSNMNull(posSites, method = "nuclAbundance");
#' PlotSectionEnrichment(posSites, negSites,
#'                       filter = c("5'UTR", "CDS", "3'UTR"));
#' }
#'
#' @export
PlotSectionEnrichment <- function(locPos,
                                  locNeg,
                                  filter = NULL,
                                  withExtendedAxisLabel = 2) {
    # Perform transcript section enrichment analysis and plot results.
    #
    # Args:
    #   locPos: A txLoc object of the positive control sites.
    #   locNeg: A txLoc object of the negative control sites.
    #   filter: Only plot loci in transcript regions specified in filter.
    #
    # Returns:
    #   NULL
    CheckClassTxLocConsistency(locPos, locNeg);
    idPos <- GetId(locPos);
    idNeg <- GetId(locNeg);
    refGenome <- GetRef(locPos);
    locPos <- GetLoci(locPos);
    locNeg <- GetLoci(locNeg);
    if (!is.null(filter)) {
        locPos <- locPos[which(names(locPos) %in% filter)];
        locNeg <- locNeg[which(names(locNeg) %in% filter)];
    }
    ctsPos <- sapply(locPos, nrow);
    ctsNeg <- sapply(locNeg, nrow);
    ctsMat <- rbind(ctsPos, ctsNeg);
    rownames(ctsMat) <- c(idPos, idNeg);
    title <- sprintf("N(%s) = %i, N(%s) = %i",
                     idPos, sum(ctsPos),
                     idNeg, sum(ctsNeg));
    par(mfrow = c(1,1));
    ret <- PlotEnrichment.Generic(ctsMat,
                                  title = title,
                                  x.las = 1, x.cex = 0.8, x.padj = 0.8,
                                  withExtendedAxisLabel = withExtendedAxisLabel);
}

#' Perform spatial enrichment analysis and plot results.
#'
#' Perform spatial enrichment analysis and plot results.
#' Enrichment/depletion is evaluated using (multiple) Fisher's exact test(s).
#' Multiple hypothesis testing correction is applied following the method of
#' Bejamini and Hochberg.
#'
#' @param locPos A \code{txLoc} object. These should be the positive control sites.
#' @param locNeg A \code{txLoc} object. These should be the negative control sites.
#' @param filter Only plot loci in transcript regions specified in filter.
#' @param binWidth Spatial bin width. Default is 20 nt.
#' @param posMax Evaluate enrichment within a window given by \code{posMax}.
#' Default is 1000 nt.
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateSNMNull(posSites, method = "permutation");
#' PlotSpatialEnrichment(posSites, negSites,
#'                       filter = c("5'UTR", "CDS", "3'UTR"));
#' }
#'
#' @export
PlotSpatialEnrichment <- function(locPos,
                                  locNeg,
                                  filter = NULL,
                                  binWidth = 20,
                                  posMax = 1000) {
    # Perform spatial enrichment analysis and plot results.
    #
    # Args:
    #   locPos: A txLoc object of the positive control sites.
    #   locNeg: A txLoc object of the negative control sites.
    #   filter: Only plot loci in transcript regions specified in filter.
    #   binWidth: Spatial bin width. Default is 20 nt.
    #   posMax: Evaluate enrichment within a window given by posMax.
    #           Default is 1kb.
    #
    # Returns:
    #   NULL
    CheckClassTxLocConsistency(locPos, locNeg);
    idPos <- GetId(locPos);
    idNeg <- GetId(locNeg);
    refGenome <- GetRef(locPos);
    locPos <- GetLoci(locPos);
    locNeg <- GetLoci(locNeg);
    if (!is.null(filter)) {
        locPos <- locPos[which(names(locPos) %in% filter)];
        locNeg <- locNeg[which(names(locNeg) %in% filter)];
    }
    par(mfrow = c(length(locPos), 2));
    breaks <- seq(0, posMax, by = binWidth);
    for (i in 1:length(locPos)) {
        posPos <- list("5p" = locPos[[i]]$TXSTART,
                       "3p" = locPos[[i]]$REGION_TXWIDTH - locPos[[i]]$TXSTART + 1);
        posNeg <- list("5p" = locNeg[[i]]$TXSTART,
                       "3p" = locNeg[[i]]$REGION_TXWIDTH - locNeg[[i]]$TXSTART + 1);
        posPos <- lapply(posPos, function(x) x[x <= posMax]);
        posNeg <- lapply(posNeg, function(x) x[x <= posMax]);
        revAxis <- list(FALSE, TRUE);
        xlab <- list("Absolute position (relative to 5' start) [nt]",
                     "Absolute position (relative to 3' end) [nt]");
        for (j in 1:length(posPos)) {
            ctsPos <- table(cut(posPos[[j]], breaks = breaks));
            ctsNeg <- table(cut(posNeg[[j]], breaks = breaks));
            ctsMat <- as.matrix(rbind(ctsPos, ctsNeg));
            rownames(ctsMat) <- c(idPos, idNeg);
            title <- sprintf("%s\nN(%s) = %i, N(%s) = %i\n%s (bw = %i nt)",
                             names(locPos)[i],
                             idPos,
                             sum(ctsPos),
                             idNeg,
                             sum(ctsNeg),
                             xlab[[j]],
                             binWidth);
            tmp <- PlotEnrichment.Generic(ctsMat,
                                          title = title,
                                          x.las = 2, x.cex = 0.8, x.padj = 0.8,
                                          reverseXaxis = revAxis[[j]]);
        }
    }
}


#' Plot ratio of two spatial distributions.
#'
#' Plot ratio of two spatial distributions.
#'
#' @param locPos A \code{txLoc} object. These should be the positive control sites.
#' @param locNeg A \code{txLoc} object. These should be the negative control sites.
#' @param filter Only plot loci in transcript regions specified in filter.
#' @param binWidth Spatial bin width. Overrides \code{nbreaks} if not
#' \code{NULL}.
#' @param posMax If \code{absolute == TRUE}, show spatial distribution
#' within a window given by \code{posMax} from the 5'/3' position of
#' the transcript feature. Default is 1000 nt.
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateSNMNull(posSites, method = "permutation");
#' PlotSpatialRatio(posSites, negSites, c("3'UTR", "CDS", "5'UTR"));
#' }
#' 
#' @export
PlotSpatialRatio <- function(locPos, locNeg,
                             filter = NULL,
                             binWidth = 20,
                             posMax = 1000) {
    # Plot the per-bin ratio of spatial distributions from locPos
    # and locNeg.
    #
    # Args:
    #   locPos: A txLoc object of the positive control sites.
    #   locNeg: A txLoc object of the negative control sites.
    #   filter: Only plot loci in transcript regions specified in filter.
    #   binWidth: Spatial bin width. Overrides nbreaks if not NULL.
    #   posMax: If absolute == TRUE, plot up distribution up to
    #                   this distance. Default is 1000 nt.
    #
    # Returns:
    #   NULL
    CheckClassTxLocConsistency(locPos, locNeg);
    idPos <- GetId(locPos);
    idNeg <- GetId(locNeg);
    refGenome <- GetRef(locPos);
    if (abs(log10(sum(GetNumberOfLoci(locPos)) /
                  sum(GetNumberOfLoci(locNeg)))) > 0.1) {
        stop("Number of positive and negative sites are too different.");
    }
    locPos <- GetLoci(locPos);
    locNeg <- GetLoci(locNeg);
    if (!is.null(filter)) {
        locPos <- locPos[which(names(locPos) %in% filter)];
        locNeg <- locNeg[which(names(locNeg) %in% filter)];
    }
    par(mfrow = c(length(locPos), 2));
    breaks <- seq(0, posMax, by = binWidth);
    mid <- breaks[-1] - binWidth / 2;
    bwString <- sprintf("bw = %i nt", binWidth);
    for (i in 1:length(locPos)) {
        posPos <- list("5p" = locPos[[i]]$TXSTART,
                       "3p" = locPos[[i]]$REGION_TXWIDTH - locPos[[i]]$TXSTART + 1);
        posNeg <- list("5p" = locNeg[[i]]$TXSTART,
                       "3p" = locNeg[[i]]$REGION_TXWIDTH - locNeg[[i]]$TXSTART + 1);
        posPos <- lapply(posPos, function(x) x[x <= posMax]);
        posNeg <- lapply(posNeg, function(x) x[x <= posMax]);
        xlim <- list(c(mid[1], mid[length(mid)]), c(mid[length(mid)], mid[1]));
        revAxis <- list(FALSE, TRUE);
        xlab <- list("Absolute position (relative to 5' start) [nt]",
                     "Absolute position (relative to 3' end) [nt]");
        for (j in 1:length(posPos)) {
            ctsPos <- table(cut(posPos[[j]], breaks = breaks));
            ctsNeg <- table(cut(posNeg[[j]], breaks = breaks));
            CIFromBS.Pos <- EstimateCIFromBS(posPos[[j]],
                                         breaks = breaks,
                                         nBS = 5000);
            CIFromBS.Neg <- EstimateCIFromBS(posNeg[[j]],
                                            breaks = breaks,
                                            nBS = 5000);
            r <- ctsPos / ctsNeg;
            r <- log10(r);
            r[is.na(r)] <- 1;
            r[is.infinite(r)] <- 1;
            plot(mid, as.numeric(r), type = "s",
                 xlab = xlab[[j]],
                 ylab = "log10(Ratio)",
                 xlim = xlim[[j]],
                 ylim = c(-2, 2),
                 main = sprintf("Ratio of occurances in %s\nN(%s)=%i w.r.t. N(%s)=%i",
                     names(locPos)[i],
                     idPos,sum(ctsPos),
                     idNeg, sum(ctsNeg)),
                 font.main = 1);
            abline(h = 0, col = "black", lty = 3);
            # Plot boostrap-based standard errors of ratio
            sePos <- abs(CIFromBS.Pos$y.high - CIFromBS.Pos$y.low) / (2 * 1.96);
            seNeg <- abs(CIFromBS.Neg$y.high - CIFromBS.Neg$y.low) / (2 * 1.96);
            dr <- 1.0 / log(10) * sqrt((sePos / ctsPos) ^ 2 + (seNeg / ctsNeg) ^ 2);
            dr[is.na(dr)] <- 1;
            dr[is.infinite(dr)] <- 1;
            x1 <- c(mid[1],
                    rep(mid[-1], each = 2),
                    mid[length(mid)]);
            CI <- cbind(c(x1,
                          rev(x1)),
                        c(rep(r - dr, each = 2),
                          rev(rep(r + dr, each = 2))));
            polygon(CI[, 1], CI[, 2], col = rgb(1, 0, 0, 0.2),
                    lwd = 1, border = NA, lty = 1);
            # Loess smoothing of standard errors
            lines(lowess(mid, r + dr, f = 1/5),
                  col = "red", lty = 2, lwd = 1);
            lines(lowess(mid, r - dr, f = 1/5),
                  col = "red", lty = 2, lwd = 1);
            legend("topleft",
                   c(sprintf("Ratio of number of sites (%s)", bwString),
                     "Boostrap-based standard error",
                     "Lowess-smoothed standard error"),
                   lwd = c(2, 5, 1),
                   col = c("black", rgb(1, 0, 0, 0.2), "red"),
                   lty = c(1, 1, 2),
                   bty = "n");
        }
    }
}


#' Plot GC content.
#'
#' Plot GC content.
#'
#' @param locPos A \code{txLoc} object.
#' @param locNeg A \code{txLoc} object.
#' @param filter A logical scalar; only consider loci in transcript
#' regions specified in filter; Default is \code{NULL}.
#' @param geneNorm A logical scalar; if \code{TRUE} normalise GC
#' content in window to GC content of transcript section; default
#' is \code{FALSE}.
#'
#' @importFrom beanplot beanplot
#' 
#' @export
PlotGC <- function(locPos, locNeg,
                   filter = NULL,
                   geneNorm = FALSE) {
    # Plot and compare GC content.
    #
    # Args:
    #   locPos: A txLoc object of the positive control sites
    #   locNeg: A txLoc object of the negative control sites.
    #
    # Returns:
    #   NULL
    CheckClassTxLocConsistency(locPos, locNeg);
    idPos <- GetId(locPos);
    idNeg <- GetId(locNeg);
    refGenome <- GetRef(locPos);
    locPos <- GetLoci(locPos);
    locNeg <- GetLoci(locNeg);
    if (!is.null(filter)) {
        locPos <- locPos[which(names(locPos) %in% filter)];
        locNeg <- locNeg[which(names(locNeg) %in% filter)];
    }
    if (length(grep("GC", colnames(locPos[[1]]))) < 1) {
        ss <- sprintf("Could not find GC content data in %s.",
                      idPos);
        ss <- sprintf("%s\nRunning CalculateGC(...) might fix that.",
                      ss);
        stop(ss);
        
    }
    if (length(grep("GC", colnames(locNeg[[1]]))) < 1) {
        ss <- sprintf("Could not find GC content data in %s.",
                      idNeg);
        ss <- sprintf("%s\nRunning CalculateGC(...) might fix that.",
                      ss);
        stop(ss);
    }
    df <- data.frame();
    namesBean <- vector();
    for (i in 1:length(locPos)) {
        GC1 <- locPos[[i]]$SITE_GC;
        GC2 <- locNeg[[i]]$SITE_GC;
        if (geneNorm) {
            GC1 <- GC1 / locPos[[i]]$REGION_GC;
            GC2 <- GC2 / locNeg[[i]]$REGION_GC;
        }
        ttest <- t.test(GC1, GC2);
        wtest <- wilcox.test(GC1, GC2);
        lbl <- sprintf("%s (%i,%i)\ndiff=%4.3f\n95%%CI=(%4.3f,%4.3f)\np=%4.3e",
                       names(locPos)[i],
                       length(GC1), length(GC2),
                       ttest$estimate[1] - ttest$estimate[2],
                       min(ttest$conf.int),
                       max(ttest$conf.int),
                       ttest$p.value,
                       wtest$p.value);
        namesBean <- c(namesBean, lbl);
        df <- rbind(df, cbind.data.frame(
            c(GC1, GC2),
            c(rep(sprintf("%s %s", names(locPos)[i], idPos), length(GC1)),
              rep(sprintf("%s %s", names(locNeg)[i], idNeg), length(GC2))),
            stringsAsFactors = FALSE));
    }
    levels <- unique(df[, 2]);
    levels <- levels[c(grep("Promoter", levels),
                       grep("5'UTR", levels),
                       grep("CDS", levels),
                       grep("3'UTR", levels),
                       grep("Introns", levels))];
    df[, 2] <- factor(df[, 2],
                      levels = levels); 
    col <- list(c(rgb(1,0,0,0.5), rgb(0.1,0.1,0.1,0.2), rgb(0.1,0.1,0.1,0.2)),
                c(rgb(0,0,1,0.5), rgb(0.1,0.1,0.1,0.2), rgb(0.1,0.1,0.1,0.2)));
    ylab <- "GC content";
    if (geneNorm) {
        ylab <- "GC content / transcript section GC content";
    }
    par(mar = c(7, 4, 4, 4) + 0.1, font.main = 1);
    beanplot(df[,1] ~ df[,2],
             bw = "nrd0",
             side = "both",
             border = NA,
             col = col,
             ylab = ylab,
             show.names = FALSE,
             main = ylab,
             font.main = 1,
             method = "jitter");
    axis(1,
         at = seq(1, length(locPos)),
         labels = namesBean,
         padj = 1,
         las = 1);
    legend("bottomleft",
           fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)),
           legend = c(idPos, idNeg),
           bty = "n");
}

