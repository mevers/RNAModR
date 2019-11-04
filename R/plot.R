#' Plot piechart of the number of loci in every transcript region
#'
#' Plot piechart of the number of loci in every transcript region.
#'
#' @param txLoc A \code{txLoc} object.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @importFrom graphics pie
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' PlotSectionDistribution(posSites)
#' }
#'
#' @export
PlotSectionDistribution <- function(txLoc) {

    # Sanity checks
    CheckClass(txLoc, "txLoc")

    # Get numbers and plot
    N <- GetNumberOfLoci(txLoc)
    labels <- names(N)
    percentage <- N / sum(N) * 100.0
    labels <- sprintf("%s: %2.1f%% (%i)", labels, percentage, N)
    pie(
        N,
        labels = labels,
        col = GetColPal("apple", length(N)),
        main = sprintf(
          "Distribution of %i %s sites across transcript sections",
          sum(N), GetId(txLoc)),
        font.main = 1)

}


#' Plot spatial distribution of loci from a \code{txLoc} object.
#'
#' Plot spatial distribution of loci from a \code{txLoc} object within every
#' transcript region.
#'
#' @param txLoc A \code{txLoc} object.
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
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' PlotSpatialDistribution(posSites)
#' PlotSpatialDistribution(posSites,
#'                         absolute = TRUE,
#'                         filter = c("5'UTR", "CDS", "3'UTR"),
#'                         ylim = c(0, 200))
#' }
#'
#' @export
PlotSpatialDistribution <- function(txLoc,
                                    nbreaks = 100,
                                    absolute = FALSE,
                                    binWidth = NULL,
                                    posMax = 1000,
                                    doBootstrap = TRUE,
                                    ...) {

    # Sanity check
    CheckClass(txLoc, "txLoc")

    # Quieten R CMD check concerns regarding "no visible binding for global
    # variable ..."
    locus_in_tx_region <- tx_region_width <- NULL

    # Get id and loci
    id <- GetId(txLoc)
    loci <- GetLoci(txLoc)

    # Determine figure panel layout and binwidth
    if (!absolute) {
        if (length(loci) < 4) {
            par(mfrow = c(1, length(loci)))
        } else {
            par(mfrow = c(ceiling(length(loci) / 2), 2))
        }
        breaks <- seq(0.0, 1.0, length.out = nbreaks)
        bw <- 1.0 / nbreaks
#        if (!is.null(binWidth)) {
#            breaks <- seq(0.0, 1.0, by = binWidth)
#            bw <- binWidth
#        }
        bwString <- sprintf("bw = %3.2f", bw)
    } else {
        par(mfrow = c(length(loci), 2))
        breaks <- seq(1, posMax, length.out = nbreaks)
        bw <- round((posMax - 1) / nbreaks)
#        if (!is.null(binWidth)) {
#            breaks <- seq(1, posMax, by = binWidth)
#            bw <- binWidth
#        }
        bwString <- sprintf("bw = %i nt", bw)
    }

    for (i in 1:length(loci)) {

        # Store site positions
        #  (1) relative to 5'start, and
        #  (2) relative to 3'end.
        pos <- with(loci[[i]], list(
          "5p" = start(locus_in_tx_region),
          "3p" = tx_region_width - end(locus_in_tx_region) + 1))

        if (!absolute) {

            pos <- lapply(pos, function(x) x / loci[[i]]$tx_region_width)
            if (grepl("(5'UTR|UTR5|5pUTR|Promoter)", names(loci)[i],
                      ignore.case = TRUE)) {
                pos <- pos["3p"]
                xrange <- list(c(1.0, 0.0))
                xlab <- list("Relative position (relative to 3' end)")
            } else {
                pos <- pos["5p"]
                xrange <- list(c(0.0, 1.0))
                xlab <- list("Relative position (relative to 5' start)")
            }
        } else {
            pos <- lapply(pos, function(x) x[x <= posMax])
            xrange <- list(c(1, posMax), c(posMax, 1))
            xlab <- list("Absolute position (relative to 5' start) [nt]",
                         "Absolute position (relative to 3' end) [nt]")
        }
        for (j in 1:length(pos)) {
            h0 <- hist(pos[[j]], breaks = breaks, plot = FALSE)
            plot(h0$mids, h0$counts,
                 type = "s",
                 lwd = 2,
                 xlab = xlab[[j]],
                 ylab = "Abundance",
                 xlim = xrange[[j]],
                 main = sprintf("%s in %s\nN = %i",
                     id,
                     names(loci)[i],
                     length(pos[[j]])),
                 font.main = 1,
                 ...)
            lblLegend <- c(sprintf("Abundance (%s)", bwString))
            lwdLegend <- c(2)
            colLegend <- c("black")
            ltyLegend <- c(1)
            if (doBootstrap) {
                CIFromBS <- EstimateCIFromBS(pos[[j]],
                                             breaks = breaks,
                                             nBS = 5000)
                x1 <- c(CIFromBS$x[1],
                        rep(CIFromBS$x[-1], each = 2),
                        CIFromBS$x[length(CIFromBS$x)])
                CI <- cbind(c(x1,
                              rev(x1)),
                            c(rep(CIFromBS$y.low, each = 2),
                              rev(rep(CIFromBS$y.high, each = 2))))
                polygon(CI[, 1], CI[, 2], col = rgb(1, 0, 0, 0.2),
                        lwd = 1, border = NA, lty = 1)
                # Loess smoothing of boostrap CI
                lines(lowess(CIFromBS$x, CIFromBS$y.low, f = 1/5),
                      col = "red", lty = 2, lwd = 1)
                lines(lowess(CIFromBS$x, CIFromBS$y.high, f = 1/5),
                      col = "red", lty = 2, lwd = 1)
                lblLegend <- c(lblLegend,
                               "95%CI (empirical bootstrap)",
                               "Lowess-smoothed 95%CI")
                lwdLegend <- c(lwdLegend, 5, 1)
                colLegend <- c(colLegend, rgb(1, 0, 0, 0.2), "red")
                ltyLegend <- c(ltyLegend, 1, 2)
            }
            legend("topleft",
                   lblLegend,
                   lwd = lwdLegend,
                   col = colLegend,
                   lty = ltyLegend,
                   bty = "n")
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
#' @param mat A \code{dataframe}; specifies the contingency table.
#' @param title A character string; specifies plot title.
#' @param xlab A character string; specifies the x-axis label.
#' @param x.las An integer scalar; specifies the orientation of
#' x-axis labels; Default is 1.
#' @param x.cex A float scalar; specifies a scaling factor for
#' x-axis labels; default is 1.
#' @param x.padj A float scalar; specifies the vertical adjustment
#' of x-axis labels; default is 1.
#' @param plotType A character string; specifies the plot style;
#' default is "l" (for line).
#' @param revXaxis A logical scalar; if \code{TRUE}, the order
#' of columns from \code{mat} is reversed; default is \code{FALSE}.
#' @param xAxisLblFmt An integer scalar; specifies the
#' format of the x-axis label.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom graphics abline axis barplot legend lines mtext
#' par plot polygon text
#' @importFrom stats fisher.test lowess p.adjust
#'
#' @return A list of \code{fisher.test} return objects and \code{mat}.
PlotEnrichment.Generic <- function(mat,
                                   title = "",
                                   xlab = "",
                                   x.las = 1, x.cex = 1, x.padj = 1,
                                   plotType = "l",
                                   revXaxis = FALSE,
                                   xAxisLblFmt = 0) {
    # Generic function to perform enrichment analysis and plot results.
    #
    # Args:
    #   mat: The data matrix.
    #   title: Plot title.
    #   x.las: Orientation of x-axis labels. Default is 1.
    #   x.cex: Scaling factor for x-axis labels. Default is 1.
    #   x.padj: Vertical adjustment of x-axis labels. Default is 1.
    #   plotType: Plot style. Default is "l".
    #   revXaxis: Reverse the order of values from mat.
    #   xAxisLblFmt: Print extended axis label.
    #
    # Returns:
    #   A list of fisher.test return objects.
    log10pval.limit <- -12
    ft <- list()
    if (ncol(mat) > 2) {
        for (i in 1:ncol(mat)) {
            redMat <- cbind(mat[, i], rowSums(mat[, -i]))
            ft[[i]] <- fisher.test(redMat)
        }
        names(ft) <- colnames(mat)
    } else {
        redMat <- mat
        ft[[1]] <- fisher.test(redMat)
        names(ft) <- colnames(mat)[1]
    }
    OR <- sapply(ft, function(x) x$estimate)
    OR[is.infinite(OR)] <- 1
    OR <- log10(OR)
    pval <- p.adjust(sapply(ft, function(x) x$p.value), method="BH")
    pval[is.na(pval)] <- 1
    pval <- log10(pval)
    pval.uncapped <- pval
    pval[pval < log10pval.limit] <- log10pval.limit
    if (ncol(mat) > 2) {
        names(OR) <- colnames(mat)
        names(pval) <- colnames(mat)
    } else {
        names(OR) <- colnames(mat)[1]
        names(pval) <- colnames(mat)[1]
    }
    CI <- sapply(ft, function(x) x$conf.int)
    CI[CI == 0] <- 1.e-12
    CI[is.infinite(CI)] <- 1e12
    CI <- log10(CI)
    if (revXaxis) {
        OR <- rev(OR)
        pval <- rev(pval)
        CI[1, ] <- rev(CI[1, ])
        CI[2, ] <- rev(CI[2, ])
    }
    #ymin <- floor(min(-abs(pval), -2.0))
    #ymax <- round(max(max(OR), 2.0))
    ymin <- log10pval.limit
    ymax <- 1
    # Plot log10(p-value)'s
    par(mar = c(5, 4, 4, 4) + 0.1)
    xrange <- 1.2 * length(pval)
    xmin <- -0.04 * xrange
    xmin <- 0.2
    xmax <- xrange + 0.04 * xrange
    xmax <- 0.96 * xmax
    mp <- barplot(pval,
                  xlim = c(xmin, xmax),
                  ylim = c(ymin, ymax),
                  col = rgb(0.2, 0.2, 1, 0.1),
                  axes = FALSE, names = "",
                  border = NA,
                  main = title, font.main = 1)
    abline(h = -2, col = "blue", lty = 3, lwd = 2)
    # Draw x axis
    idxLabels <- seq(1, length(OR))
    labels <- names(OR)
    if (xAxisLblFmt == 1) {
        labels = sprintf("%s\nOR=%4.3f,p=%4.3e",
            names(OR),
            10^OR,
            10^pval.uncapped)
    } else if (xAxisLblFmt == 2) {
        labels = sprintf("%s\nN(%s)=%i\nN(%s)=%i\nOR=%4.3f\np=%4.3e",
            names(OR),
            rownames(mat)[1], mat[1, ],
            rownames(mat)[2], mat[2, ],
            10^OR,
            10^pval.uncapped)
    } else if (xAxisLblFmt >= 3) {
        x1 <- as.numeric(
            gsub("(\\(|,[+-]*\\d+\\.*\\d*e*\\+*\\d*])", "", labels))
        x2 <- as.numeric(
            gsub("(\\([+-]*\\d+\\.*\\d*e*\\+*\\d*,|])", "", labels))
        labels <- (x2 + x1) / 2
        # Zero-crossing?
        idxZero <- which(x2 == 0 | (x1 < 0 & x2 > 0))
        if (length(idxZero) > 0) {
            abline(v = mp[idxZero],
                   col = rgb(0, 0, 0, 0.2),
                   lwd = 2,
                   lty = 3)
        }
        if (xAxisLblFmt == 3) {
            deltaIdx <- floor(length(x1) / 10 + 0.5)
            if (length(idxZero) == 0) {
                idxLabels <- seq(1, length(x1), deltaIdx)
                if (revXaxis == TRUE) {
                    idxLabels <- sort(seq(length(x1), 1, -deltaIdx))
                }
            } else {
                deltaIdx <- round(length(labels) / 10 + 0.5)
                idxLabels <- sort(
                    c(seq(idxZero, by = -deltaIdx),
                      seq(idxZero + deltaIdx, length(labels), by = deltaIdx)))
            }
        }
    }
    axis(1,
         at = mp[idxLabels],
         labels = labels[idxLabels],
         las = x.las,
         padj = x.padj,
         cex.axis = x.cex)
    mtext(xlab, side = 1, line = 2.7)
    # Draw right y axis
    axis(4, at = seq(ymin, ymax))
    mtext("log10(p-value)", side = 4, line = 2.5)
    # Draw significance labels
    sig <- pval
    sig[pval <= -4] <- "****"
    sig[pval <= -3 & pval > -4] <- "***"
    sig[pval <= -2 & pval > -3] <- "**"
    sig[pval <= -1.30103 & pval > -2] <- "*"
    sig[pval > -1.30103] <- "ns"
    text(mp, y = 0.5, labels = sig, srt = 90, col = "blue")
    # Draw legend
    legend("bottomleft",
           c("Odds-ratio (OR)",
             "95% CI OR",
             "p-value"),
           lwd = c(2, 5, 5),
           col = c("red", rgb(1, 0, 0, 0.2), rgb(0.2, 0.2, 1, 0.1)),
           lty = c(1, 1, 1),
           bty = "n")
    # Plot odds-ratios
    par(new = TRUE)
    plot(mp, OR, col = rgb(1, 0.2, 0.2, 1),
         lwd = 2, type = plotType,
         xlim = c(xmin, xmax),
         ylim = c(-1.0, 1.0),
         axes = FALSE, xlab = "", ylab = "")
    # Draw 95% confidence intervals
    CI <- cbind(c(mp, rev(mp)),
                c(CI[1 ,], rev(CI[2, ])))
    polygon(CI[,1], CI[,2], col = rgb(1, 0, 0, 0.2),
            lwd = 1, border = NA, lty = 1)
    axis(2, at = seq(-1.0, 1.0))
    mtext("log10(OR)", side = 2, line = 2.5)
    abline(h = 0.0, col = "red", lty = 3, lwd = 2)
    return(list(ft = ft,
                mat = mat))
}


#' Perform enrichment analysis of sites per transcript region and plot results.
#'
#' Perform enrichment analysis of the number of positive sites in
#' \code{txLoc.pos} relative to the number of null sites in \code{txLoc.neg}
#' per transript region and plot results.
#' Enrichment/depletion is evaluated using (multiple) Fisher's exact test(s).
#' Multiple hypothesis testing correction is applied following the method of
#' Bejamini and Hochberg.
#'
#' @param txLoc.pos A \code{txLoc} object. These correspond to the positive
#' sites.
#' @param txLoc.neg A \code{txLoc} object. These correspond to the negative
#' (null) sites.
#' @param xAxisLblFmt Plot extended axis labels. Default is 2. See ??? for details.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' pos <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' neg <- GenerateNull(pos, method = "nuclAbundance")
#' PlotSectionEnrichment(
#'     FilterTxLoc(pos, c("5'UTR", "CDS", "3'UTR")),
#'     FilterTxLoc(neg, c("5'UTR", "CDS", "3'UTR")))
#' }
#'
#' @export
PlotSectionEnrichment <- function(txLoc.pos,
                                  txLoc.neg,
                                  xAxisLblFmt = 2) {

    # Sanity check
    CheckClassTxLocConsistency(txLoc.pos, txLoc.neg)

    # Get counts and plot
    ctsPos <- sapply(GetLoci(txLoc.pos), nrow)
    ctsNeg <- sapply(GetLoci(txLoc.neg), nrow)
    ctsMat <- rbind(ctsPos, ctsNeg)
    rownames(ctsMat) <- c(GetId(txLoc.pos), GetId(txLoc.neg))
    title <- sprintf("N(%s) = %i, N(%s) = %i",
                     GetId(txLoc.pos), sum(ctsPos),
                     GetId(txLoc.neg), sum(ctsNeg))
    par(mfrow = c(1,1))
    ret <- PlotEnrichment.Generic(
        ctsMat,
        title = title,
        xlab = "",
        x.las = 1, x.cex = 0.8, x.padj = 0.8,
        xAxisLblFmt = xAxisLblFmt)
}


#' Perform spatial enrichment analysis and plot results.
#'
#' Perform enrichment analysis of the spatial distribution of positive sites in
#' \code{txLoc.pos} relative to the distribution of null sites in
#' \code{txLoc.neg} per transcript region and plot results.
#' Enrichment/depletion is evaluated using (multiple) Fisher's exact test(s).
#' Multiple hypothesis testing correction is applied following the method of
#' Bejamini and Hochberg.
#'
#' @param txLoc.pos A \code{txLoc} object. These correspond to the positive
#' sites.
#' @param txLoc.neg A \code{txLoc} object. These correspond to the negative
#' (null) sites.
#' @param binWidth Spatial bin width. Default is 20 nt.
#' @param posMax Evaluate enrichment within a window given by \code{posMax}.
#' Default is 1000 nt.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR")
#' sites <- ReadBED(bedFile)
#' pos <- SmartMap(sites, id = "m6A", refGenome = "hg38")
#' null <- GenerateNull(posSites, method = "permutation")
#' PlotSpatialEnrichment(
#'     FilterTxLoc(pos, c("5'UTR", "CDS", "3'UTR")),
#'     FilterTxLoc(neg, c("5'UTR", "CDS", "3'UTR")))
#' }
#'
#' @export
PlotSpatialEnrichment <- function(txLoc.pos,
                                  txLoc.neg,
                                  binWidth = 20,
                                  posMax = 1000) {

    # Sanity check
    CheckClassTxLocConsistency(txLoc.pos, txLoc.neg)

    # Quieten R CMD check concerns regarding "no visible binding for global
    # variable ..."
    locus_in_tx_region <- tx_region_width <- NULL

    # Determine figure panel layout and number of breaks
    par(mfrow = c(length(GetLoci(txLoc.pos)), 2))
    breaks <- seq(0, posMax, by = binWidth)

    # Plot
    invisible(Map(
        function(loci.pos, loci.neg, region) {

            # Store site positions
            #  (1) relative to 5'start, and
            #  (2) relative to 3'end.
            pos.pos <- with(loci.pos, list(
              "5p" = start(locus_in_tx_region),
              "3p" = tx_region_width - end(locus_in_tx_region) + 1))
            pos.neg <- with(loci.neg, list(
              "5p" = start(locus_in_tx_region),
              "3p" = tx_region_width - end(locus_in_tx_region) + 1))

            # Limit distance to window [0, posMax]
            pos.pos <- lapply(pos.pos, function(x) x[x <= posMax])
            pos.neg <- lapply(pos.neg, function(x) x[x <= posMax])

            # Set axis properties
            revAxis <- list(FALSE, TRUE)
            xlab <- list("Absolute position (relative to 5' start) [nt]",
                         "Absolute position (relative to 3' end) [nt]")

            # Plot
            for (j in 1:length(pos.pos)) {
                ctsPos <- table(cut(pos.pos[[j]], breaks = breaks))
                ctsNeg <- table(cut(pos.neg[[j]], breaks = breaks))
                ctsMat <- as.matrix(rbind(ctsPos, ctsNeg))
                rownames(ctsMat) <- c(GetId(txLoc.pos), GetId(txLoc.neg))
                title <- sprintf(
                    "%s\nN(%s) = %i, N(%s) = %i\nbw = %i nt, window = %i nt",
                    region,
                    GetId(txLoc.pos),
                    sum(ctsPos),
                    GetId(txLoc.neg),
                    sum(ctsNeg),
                    binWidth,
                    posMax)
                tmp <- PlotEnrichment.Generic(
                    ctsMat,
                    title = title,
                    xlab = xlab[[j]],
                    x.las = 1, x.cex = 0.8, x.padj = 0,
                    revXaxis = revAxis[[j]],
                    xAxisLblFmt = 3)
            }
        },
        GetLoci(txLoc.pos),
        GetLoci(txLoc.neg),
        GetRegions(txLoc.pos)))

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
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' negSites <- GenerateNull(posSites, method = "permutation");
#' PlotSpatialRatio(posSites, negSites, c("3'UTR", "CDS", "5'UTR"));
#' }
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
    CheckClassTxLocConsistency(locPos, locNeg)
    idPos <- GetId(locPos)
    idNeg <- GetId(locNeg)
    refGenome <- GetRef(locPos)
    if (abs(log10(sum(GetNumberOfLoci(locPos)) /
                  sum(GetNumberOfLoci(locNeg)))) > 0.1) {
        stop("Number of positive and negative sites are too different.")
    }
    locPos <- GetLoci(locPos)
    locNeg <- GetLoci(locNeg)
    if (!is.null(filter)) {
        locPos <- locPos[which(names(locPos) %in% filter)]
        locNeg <- locNeg[which(names(locNeg) %in% filter)]
    }
    par(mfrow = c(length(locPos), 2))
    breaks <- seq(0, posMax, by = binWidth)
    mid <- breaks[-1] - binWidth / 2
    bwString <- sprintf("bw = %i nt", binWidth)
    for (i in 1:length(locPos)) {
        posPos <- list("5p" = locPos[[i]]$TXSTART,
                       "3p" = locPos[[i]]$REGION_TXWIDTH - locPos[[i]]$TXSTART + 1)
        posNeg <- list("5p" = locNeg[[i]]$TXSTART,
                       "3p" = locNeg[[i]]$REGION_TXWIDTH - locNeg[[i]]$TXSTART + 1)
        posPos <- lapply(posPos, function(x) x[x <= posMax])
        posNeg <- lapply(posNeg, function(x) x[x <= posMax])
        xlim <- list(c(mid[1], mid[length(mid)]), c(mid[length(mid)], mid[1]))
        revAxis <- list(FALSE, TRUE)
        xlab <- list("Absolute position (relative to 5' start) [nt]",
                     "Absolute position (relative to 3' end) [nt]")
        for (j in 1:length(posPos)) {
            ctsPos <- table(cut(posPos[[j]], breaks = breaks))
            ctsNeg <- table(cut(posNeg[[j]], breaks = breaks))
            CIFromBS.Pos <- EstimateCIFromBS(posPos[[j]],
                                         breaks = breaks,
                                         nBS = 5000)
            CIFromBS.Neg <- EstimateCIFromBS(posNeg[[j]],
                                            breaks = breaks,
                                            nBS = 5000)
            r <- ctsPos / ctsNeg
            r <- log10(r)
            r[is.na(r)] <- 1
            r[is.infinite(r)] <- 1
            plot(mid, as.numeric(r), type = "s",
                 xlab = xlab[[j]],
                 ylab = "log10(Ratio)",
                 xlim = xlim[[j]],
                 ylim = c(-2, 2),
                 main = sprintf("Ratio of occurances in %s\nN(%s)=%i w.r.t. N(%s)=%i",
                     names(locPos)[i],
                     idPos,sum(ctsPos),
                     idNeg, sum(ctsNeg)),
                 font.main = 1)
            abline(h = 0, col = "black", lty = 3)
            # Plot boostrap-based standard errors of ratio
            sePos <- abs(CIFromBS.Pos$y.high - CIFromBS.Pos$y.low) / (2 * 1.96)
            seNeg <- abs(CIFromBS.Neg$y.high - CIFromBS.Neg$y.low) / (2 * 1.96)
            dr <- 1.0 / log(10) * sqrt((sePos / ctsPos) ^ 2 + (seNeg / ctsNeg) ^ 2)
            dr[is.na(dr)] <- 1
            dr[is.infinite(dr)] <- 1
            x1 <- c(mid[1],
                    rep(mid[-1], each = 2),
                    mid[length(mid)])
            CI <- cbind(c(x1,
                          rev(x1)),
                        c(rep(r - dr, each = 2),
                          rev(rep(r + dr, each = 2))))
            polygon(CI[, 1], CI[, 2], col = rgb(1, 0, 0, 0.2),
                    lwd = 1, border = NA, lty = 1)
            # Loess smoothing of standard errors
            lines(lowess(mid, r + dr, f = 1/5),
                  col = "red", lty = 2, lwd = 1)
            lines(lowess(mid, r - dr, f = 1/5),
                  col = "red", lty = 2, lwd = 1)
            legend("topleft",
                   c(sprintf("Ratio of number of sites (%s)", bwString),
                     "Boostrap-based standard error",
                     "Lowess-smoothed standard error"),
                   lwd = c(2, 5, 1),
                   col = c("black", rgb(1, 0, 0, 0.2), "red"),
                   lty = c(1, 1, 2),
                   bty = "n")
        }
    }
}


#' Plot GC content.
#'
#' Plot and assess the difference in the distributions of GC content within
#' a window around sites from two \code{txLoc} objects. See 'Details'.
#'
#' The function calculates the GC content within a window around every site
#' from two \code{txLoc} objects. The window is defined by extending the
#' position of every \code{txLoc} site upstream and downstream by \code{flank}
#' nucleotides (if possible).
#' The means of the resulting GC content distributions are assessed using a
#' two-sample two-tailed t-test. If \code{norm_region = TRUE}, the GC content
#' in the window is normalised to the GC content of the entire transcript
#' region. If \code{downsample = TRUE}, the number of windows/sites from the
#' \emph{second} \code{txLoc2} object is downsampled to match the number of
#' windows/sites from the \code{txLoc1} object. This is useful (and therefore
#' the default setting) when comparing the GC distribution around positive
#' and null sites, as the list of null sites is often significantly larger than
#' that of the positive sites.
#' Note that this function calls \code{GetGC}, which performs the sequence
#' extraction and GC calculation. See \code{?GetGC} for details.
#'
#' @param txLoc1 A \code{txLoc} object.
#' @param txLoc2 A \code{txLoc} object.
#' @param flank An \code{integer} scalar; see 'Details'.
#' @param norm_region A \code{logical} scalar; if \code{TRUE} normalise the GC
#' content in the window to the GC content of the corresponding transcript
#' region; default is \code{FALSE}.
#' @param downsample A \code{logical} scalar; if \code{TRUE}, subsample sites
#' from \code{txLoc2} to match the number of sites per transcript region from
#' \code{txLoc1}; default is \code{TRUE}.
#' @param seed A single value, interpreted as an \code{integer}, or \code{NULL};
#' this is to ensure reproducibility when subsampling \code{txLoc2} sites;
#' ignored when \code{downsample == FALSE}; default is \code{NULL}.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @importFrom beanplot beanplot
#' @importFrom stats t.test wilcox.test
#'
#' @export
PlotGC <- function(txLoc1, txLoc2,
                   flank = 10,
                   norm_region = FALSE,
                   downsample = TRUE,
                   seed = NULL) {

    # Sanity check
    CheckClassTxLocConsistency(txLoc1, txLoc2)

    id1 <- GetId(txLoc1)
    id2 <- GetId(txLoc2)
    refGenome <- GetRef(txLoc1)

    # Downsample txLoc2 to txLoc1
    if (downsample == TRUE)
        txLoc2 <- DownsampleTxLoc(txLoc2, txLoc1, seed = seed)

    # Get GC content of window around sites and transcript region
    dataGC1 <- GetGC(txLoc1, flank = flank)
    dataGC2 <- GetGC(txLoc2, flank = flank)

    df <- data.frame()

    namesBean <- vector()
    for (i in 1:length(dataGC1)) {
        GC1 <- dataGC1[[i]][, "GC_window"]
        GC2 <- dataGC2[[i]][, "GC_window"]
        if (norm_region) {
            GC1 <- GC1 / dataGC1[[i]][, "GC_tx_region"]
            GC2 <- GC2 / dataGC2[[i]][, "GC_tx_region"]
        }
        ttest <- t.test(GC1, GC2)
        wtest <- wilcox.test(GC1, GC2)
        lbl <- sprintf("%s (%i,%i)\ndiff=%4.3f\n95%%CI=(%4.3f,%4.3f)\np=%4.3e",
                       names(dataGC1)[i],
                       length(GC1), length(GC2),
                       ttest$estimate[1] - ttest$estimate[2],
                       min(ttest$conf.int),
                       max(ttest$conf.int),
                       ttest$p.value,
                       wtest$p.value)
        namesBean <- c(namesBean, lbl)
        df <- rbind(df, cbind.data.frame(
            c(GC1, GC2),
            c(rep(sprintf("%s %s", names(dataGC1)[i], id1), length(GC1)),
              rep(sprintf("%s %s", names(dataGC2)[i], id2), length(GC2))),
            stringsAsFactors = FALSE))
    }

    levels <- unique(df[, 2])
    levels <- levels[c(grep("Promoter", levels),
                       grep("5'UTR", levels),
                       grep("CDS", levels),
                       grep("3'UTR", levels),
                       grep("Introns", levels))]
    df[, 2] <- factor(df[, 2],
                      levels = levels);
    col <- list(c(rgb(1,0,0,0.5), rgb(0.1,0.1,0.1,0.2), rgb(0.1,0.1,0.1,0.2)),
                c(rgb(0,0,1,0.5), rgb(0.1,0.1,0.1,0.2), rgb(0.1,0.1,0.1,0.2)))
    ylab <- "GC content"
    if (norm_region) {
        ylab <- "GC content / transcript section GC content"
    }
    par(mar = c(7, 4, 4, 4) + 0.1, font.main = 1)
    beanplot(df[,1] ~ df[,2],
             ll = 0.04,
             bw = 0.05,
             side = "both",
             border = NA,
             col = col,
             ylab = ylab,
             show.names = FALSE,
             main = ylab,
             font.main = 1,
             method = "jitter")
    axis(1,
         at = seq(1, length(dataGC1)),
         labels = namesBean,
         padj = 1,
         las = 1)
    legend("bottomleft",
           fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)),
           legend = c(id1, id2),
           bty = "n")
}


#' Generic function for plotting abundances (histograms).
#'
#' Generic function for plotting abundances (histograms).
#'
#' @param data A vector of integers.
#' @param xmin An integer scalar; default is \code{NULL}.
#' @param xmax An integer scalar; default is \code{NULL}.
#' @param binWidth An integer scalar; default is \code{NULL}.
#' @param plotType A single character; default is \code{"s"}.
#' @param lwd An integer scalar; default is 2.
#' @param title A character string; default is \code{""}.
#' @param xlab A character string; default is \code{""}.
#' @param ylab A character string; default is \code{""}.
#' @param doBootstrap A logical scalar; default is \code{TRUE}.
#' @param nBS An integer scalar; default is 5000.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#'
#' @importFrom graphics legend lines plot polygon
#' @importFrom stats lowess
#'
#' @export
PlotAbundance.generic <- function(data,
                                  xmin = NULL,
                                  xmax = NULL,
                                  binWidth = NULL,
                                  plotType = "s",
                                  lwd = 2,
                                  title = "",
                                  xlab = "",
                                  ylab = "Abundance",
                                  doBootstrap = TRUE,
                                  nBS = 5000,
                                  doLowess = TRUE) {
    if (!is.null(xmin) && !is.null(xmax) && !is.null(binWidth)) {
        breaks <- seq(xmin, xmax, by = binWidth)
    } else {
        breaks <- seq(min(data), max(data), length.out = 50)
    }
    bwString <- sprintf("bw = %3.2f", binWidth)
    h0 <- hist(data,
               breaks = breaks,
               plot = FALSE)
    plot(h0$mids, h0$counts,
         type = plotType,
         lwd = lwd,
         xlab = xlab,
         ylab = ylab,
         main = title,
         font.main = 1)
    if (doBootstrap) {
        CIFromBS <- EstimateCIFromBS(data,
                                     breaks = breaks,
                                     nBS = nBS)
        x1 <- c(CIFromBS$x[1],
                rep(CIFromBS$x[-1], each = 2),
                CIFromBS$x[length(CIFromBS$x)])
        CI <- cbind(c(x1,
                      rev(x1)),
                    c(rep(CIFromBS$y.low, each = 2),
                      rev(rep(CIFromBS$y.high, each = 2))))
        polygon(CI[, 1], CI[, 2], col = rgb(1, 0, 0, 0.2),
                lwd = 1, border = NA, lty = 1)
        lblLegend <- c(sprintf("Abundance (%s)", bwString),
                      "95%CI (empirical bootstrap)")
        lwdLegend <- c(2, 5)
        colLegend <- c("black", rgb(1, 0, 0, 0.2))
        ltyLegend <- c(1, 1)
        # Loess smoothing of boostrap CI
        if (doLowess) {
            lines(lowess(CIFromBS$x, CIFromBS$y.low, f = 1/5),
                  col = "red", lty = 2, lwd = 1)
            lines(lowess(CIFromBS$x, CIFromBS$y.high, f = 1/5),
                  col = "red", lty = 2, lwd = 1)
            lblLegend <- c(lblLegend,
                           "Lowess-smoothed 95%CI")
            lwdLegend <- c(lwdLegend, 1)
            colLegend <- c(colLegend, "red")
            ltyLegend <- c(ltyLegend, 2)
        }
        legend("topleft",
               lblLegend,
               lwd = lwdLegend,
               col = colLegend,
               lty = ltyLegend,
               bty = "n")
    }
}


#' Plot distribution of relative distances.
#'
#' Plot distribution of relative distances between sites
#' from two \code{txLoc} objects. See 'Details'.
#'
#' The function calculates minimum distances per transcript region, between
#' entries from \code{txLoc} relative to \code{txLocRef}. Relative distances
#' are shown within a window (-\code{flank}, \code{flank}), where negative
#' distances correspond to a feature from \code{txLoc} that is upstream of a
#' site from \code{txLocRef}, and positive distances indicate a feature from
#' \code{txLoc} that is downstream of a site from \code{txLocRef}. Relative
#' distances are binned in bins of \code{binWidth} nt, and shown as an
#' abundance histogram. If \code{doBootstrap = TRUE}, 95% confidence intervals
#' are calculated and shown, based on an empirical bootstrap of relative
#' distances.
#'
#' @param txLoc A \code{txLoc} object.
#' @param txLocRef A \code{txLoc} object.
#' @param flank An integer scalar or an integer vector of length 2; specifies 
#' the absolute maximum relative distance(s) used as a cutoff; by specifying 
#' \code{flank} as a vector, a non-symmetric window can be  defined; e.g. 
#' `flank = c(1000, 2000)` corresponds to a window defined by 1000 nt 
#' downstream and 2000 nt upstream of the reference site. Default is 1000.
#' @param binWidth An \code{integer} scalar; specifies the spatial width in nt
#' by which distances will be binned; default is 20.
#' @param doBootstrap A \code{logical} scalar; if \code{YES} calculate
#' 95% CI based on empirical bootstrap of relative distances within
#' transcript region; default is \code{TRUE}.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @export
PlotRelDistDistribution <- function(txLoc,
                                    txLocRef,
                                    flank = 1000,
                                    binWidth = 20,
                                    doBootstrap = TRUE) {

    # Sanity checks
    CheckClass(txLoc, "txLoc")
    CheckClass(txLocRef, "txLoc")
    CheckClassTxLocConsistency(txLoc, txLocRef)
    
    # Allow for variable window sizes
    if (length(flank) == 1) {
        flank <- c(-abs(flank), abs(flank))
    } else if (length(flank) == 2) {
        flank <- c(-abs(flank[1]), abs(flank[2]))
    } else {
        stop("`flank` needs to be a scalar or vector of length 2!")
    }

    # Get ids
    id <- GetId(txLoc)
    idRef <- GetId(txLocRef)

    # Convert sites to `list` of `GRanges`
    lst <- TxLoc2GRangesList(txLoc, method = "tx_region")
    lstRef <- TxLoc2GRangesList(txLocRef, method = "tx_region")

    # Calculate distances
    lstDist <- GetRelDistNearest(lst, lstRef)

    # Determine figure panel layout
    if (length(lstDist) < 4) {
        par(mfrow = c(1, length(lstDist)))
    } else {
        par(mfrow = c(ceiling(length(lstDist) / 2), 2))
    }

    # Set breaks & binwidth
    breaks <- seq(flank[1], flank[2], by = binWidth)
    bwString <- sprintf("bw = %3.2f", binWidth)


    invisible(Map(
        function(dist, region) {

            # Filter distances that are within window [flank[1], flank[2]]
            dist <- dist[dist >= flank[1] & dist <= flank[2]]

            # Plot
            title <- sprintf(
                "d(%s,%s) in %s (N=%i)\nbw = %i nt",
                id,
                idRef,
                region,
                length(dist),
                binWidth)
            xlab <- sprintf("Relative distance to %s [nt]", idRef)
            PlotAbundance.generic(
                dist,
                xmin = flank[1], xmax = flank[2],
                binWidth = binWidth,
                title = title,
                xlab = xlab,
                doBootstrap = doBootstrap)

        },
        lstDist, names(lstDist)))

}


#' Perform enrichment analysis of relative distances.
#'
#' Perform enrichment analysis and plot results of two relative
#' distance distributions. See 'Details'.
#'
#' The function calculates minimum distances per transcript region, between
#' entries from \code{txLoc1} relative to \code{txLocRef}, and \code{txLoc2}
#' relative to \code{txLocRef}.
#' Enrichment/depletion is assessed using multiple Fisher's exact tests on the
#' counts per distance bin relative to the counts in all other bins within the
#' window defined by (-\code{flank}, \code{flank}). Resulting enrichment plots
#' show odds-ratios (including 95\% confidence intervals) and associated
#' p-values as a function of relative distance bins. Negative distances
#' indicate sites from \code{txLoc1} and \code{txLoc2} that are \emph{upstream}
#' of sites from \code{txLocRef}; positive distances correspond to sites from
#' \code{txLoc1} and \code{txLoc2} that are \emph{downstream} of
#' \code{txLocRef}. The bin width and window size can be adjusted with
#' \code{flank} and \code{binWidth}.
#'
#' @param txLoc1 A \code{txLoc} object.
#' @param txLoc2 A \code{txLoc} object.
#' @param txLocRef A \code{txLoc} object.
#' @param flank An integer scalar or an integer vector of length 2; specifies 
#' the absolute maximum relative distance(s) used as a cutoff; by specifying 
#' \code{flank} as a vector, a non-symmetric window can be  defined; e.g. 
#' `flank = c(1000, 2000)` corresponds to a window defined by 1000 nt 
#' downstream and 2000 nt upstream of the reference site. Default is 1000.
#' @param binWidth An \code{integer} scalar; specifies the spatial width in nt 
#' by which distances will be binned; default is \code{20}.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @export
PlotRelDistEnrichment <- function(txLoc1,
                                  txLoc2,
                                  txLocRef,
                                  flank = 1000,
                                  binWidth = 20) {
    # Sanity checks
    CheckClass(txLoc1, "txLoc")
    CheckClass(txLoc2, "txLoc")
    CheckClass(txLocRef, "txLoc")
    CheckClassTxLocConsistency(txLoc1, txLocRef)
    CheckClassTxLocConsistency(txLoc2, txLocRef)

    # Allow for variable window sizes
    if (length(flank) == 1) {
        flank <- c(-abs(flank), abs(flank))
    } else if (length(flank) == 2) {
        flank <- c(-abs(flank[1]), abs(flank[2]))
    } else {
        stop("`flank` needs to be a scalar or vector of length 2!")
    }

    # Get ids
    id1 <- GetId(txLoc1)
    id2 <- GetId(txLoc2)
    idRef <- GetId(txLocRef)

    # Convert sites to `list` of `GRanges`
    lst1 <- TxLoc2GRangesList(txLoc1, method = "tx_region")
    lst2 <- TxLoc2GRangesList(txLoc2, method = "tx_region")
    lstRef <- TxLoc2GRangesList(txLocRef, method = "tx_region")

    # Calculate distances
    lstDist1 <- GetRelDistNearest(lst1, lstRef)
    lstDist2 <- GetRelDistNearest(lst2, lstRef)

    # Determine figure panel layout
    if (length(lstDist1) < 4) {
        par(mfrow = c(1, length(lstDist1)))
    } else {
        par(mfrow = c(ceiling(length(lstDist1) / 2), 2))
    }

    # Set breaks & binwidth
    breaks <- seq(flank[1], flank[2], by = binWidth)
    bwString <- sprintf("bw = %3.2f", binWidth)

    invisible(Map(
        function(dist1, dist2, region) {

            # Filter distances that are within window [flank[1], flank[2]]
            dist1 <- dist1[dist1 >= flank[1] & dist1 <= flank[2]]
            dist2 <- dist2[dist2 >= flank[1] & dist2 <= flank[2]]

            # Bin distances and count matrix
            cts1 <- table(cut(dist1, breaks = breaks))
            cts2 <- table(cut(dist2, breaks = breaks))
            ctsMat <- as.matrix(rbind(cts1, cts2))
            rownames(ctsMat) <- c("pos", "neg")

            # Plot
            title <- sprintf(
                "%s\nN(d(%s,%s)) = %i, N(d(%s,%s)) = %i\nbw = %i nt",
                region,
                id1, idRef, sum(cts1),
                id2, idRef, sum(cts2),
                binWidth)
            PlotEnrichment.Generic(
                ctsMat,
                title = title,
                xlab = sprintf("Distance relative to %s [nt]", idRef),
                x.las = 1, x.cex = 0.8, x.padj = 0,
                xAxisLblFmt = 3)

        },
        lstDist1, lstDist2, names(lstDist1)))

}


#' Plot sequence logo.
#'
#' Plot sequence logo.
#'
#' The function determines the sequence logo within a window defined by
#' extending sites from \code{txLoc} upstream and downstream by \code{flank}
#' nucleotides.
#'
#' @param txLoc A \code{txLoc} object.
#' @param flank An \code{integer} scalar; see 'Details'.
#' @param ylim An \code{integer} vector; specifies limits for the y-axis;
#' automatically determined if \code{ymin = NULL}; default is \code{c(0, 2)}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @import Biostrings
#'
#' @export
PlotSeqLogo <- function(txLoc, flank = 5, ylim = c(0, 2)) {

    # Sanity check
    CheckClass(txLoc, "txLoc")

    # Determine figure panel layout
    if (length(GetRegions(txLoc)) < 4) {
        par(mfrow = c(1, length(GetRegions(txLoc))))
    } else {
        par(mfrow = c(ceiling(length(GetRegions(txLoc)) / 2), 2))
    }

    invisible(Map(
        function(locus, region) {

            # Define window coordinates
            x1 <- start(locus$locus_in_tx_region) - flank
            x2 <- start(locus$locus_in_tx_region) + flank

            # Extract subsequences
            # We use substr here because it's vectorised in all arguments
            # and we don't have to worry about start < 0 arguments
            seq_window <- substr(
                locus$tx_region_sequence, start = x1, stop = x2)

            # Only keep sequences that have the width 2 * flank + 1
            seq_window <- seq_window[nchar(seq_window) == 2 * flank + 1L]

            # Concensus matrix
            # We only keep the first 4 rows to avoid any non-ACGT bases which
            # would show up in rows > 4
            mat <- consensusMatrix(seq_window, as.prob = TRUE)[1:4, ]

            # Convert to `data.frame` and calculate Shannon entropy
            df <- as.data.frame(t(mat))
            df$pos <- seq(-flank, flank, by = 1)
            df$height <- apply(df[, 1:4], 1, function(x) {
                x[x == 0] <- 1.e-9
                2 + sum(x * log2(x))
                })
            df <- data.frame(
                A = df$A * df$height,
                C = df$C * df$height,
                G = df$G * df$height,
                T = df$T * df$height,
                pos = df$pos)

            # Plot
            title <- sprintf("%s, N(%s)=%i\nSequence logo in %i nt window",
                             region,
                             GetId(txLoc),
                             nrow(locus),
                             2 * flank + 1)
            mp <- barplot(t(df[ , 1:4]),
                          col = GetColPal("google", 4),
                          ylim = ylim,
                          ylab = "Information content [bits]",
                          main = title,
                          font.main = 1)
            axis(1, at = mp, labels = df[, "pos"])
            mtext("Relative position [nt]", 1, padj = 4)
            legend("topright",
                   fill = GetColPal("google", 4),
                   legend = c("A", "C", "G", "T"),
                   bty = "n")

        },
        GetLoci(txLoc), GetRegions(txLoc)))

}


#' Enrichment analysis of sites relative to start/stop codon.
#'
#' Perform and visualise the enrichment analysis of sites from a \code{txLoc} 
#' object relative to the start/stop codons from a reference transcriptome. 
#' See 'Details'.
#'
#' The function calculates relative distances of sites from a \code{txLoc}
#' object to the corresponding transcript's start and stop codons.
#' Enrichment/depletion is assessed using multiple Fisher's exact tests on the
#' counts per distance bin relative to the counts in all other bins within the
#' window defined by (-\code{flank}, \code{flank}). Resulting enrichment plots
#' show odds-ratios (including 95\% confidence intervals) and associated
#' p-values as a function of relative distance bins. Negative distances
#' indicate sites from \code{txLoc} that are \emph{upstream} of the start/stop
#' codon; positive distances correspond to sites from \code{txLoc} that are 
#' \emph{downstream} of \code{txLocRef}. The bin width and window size can be 
#' adjusted with \code{flank} and \code{binWidth}.
#' Note that this function can be quite slow if \code{txLoc2} is based on all 
#' null sites. If this is the case, consider downsampling \code{txLoc2} using
#' \code{DownsampleTxLoc}.
#'
#' @param txLoc1 A \code{txLoc} object.
#' @param txLoc2 A \code{txLoc} object.
#' @param flank An integer scalar or an integer vector of length 2; specifies 
#' the absolute maximum relative distance(s) used as a cutoff; by specifying 
#' \code{flank} as a vector, a non-symmetric window can be  defined; e.g. 
#' `flank = c(1000, 2000)` corresponds to a window defined by 1000 nt 
#' downstream and 2000 nt upstream of the reference site. Default is 1000.
#' @param binWidth An \code{integer} scalar; specifies the spatial width in nt
#' by which distances will be binned; default is 20.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @keywords internal
#'
#' @export
PlotRelStartStopEnrichment <- function(txLoc1, 
                                       txLoc2,
                                       flank = 550, 
                                       binWidth = 20) {
    
    # Get distance of sites to start/stop codon
    lstDist1 <- GetDistNearestStartStop(txLoc1)
    lstDist2 <- GetDistNearestStartStop(txLoc2)
    
    # Allow for variable window sizes
    if (length(flank) == 1) {
        flank <- c(-abs(flank), abs(flank))
    } else if (length(flank) == 2) {
        flank <- c(-abs(flank[1]), abs(flank[2]))
    } else {
        stop("`flank` needs to be a scalar or vector of length 2!")
    }
    
    # Determine figure panel layout
    par(mfrow = c(1, 2))

    # Set breaks & binwidth
    breaks <- seq(flank[1], flank[2], by = binWidth)
    bwString <- sprintf("bw = %3.2f", binWidth)
    
    invisible(Map(
        function(dist1, dist2, region) {
            
            # Filter distances that are within window [flank[1], flank[2]]
            dist1 <- dist1[dist1 >= flank[1] & dist1 <= flank[2]]
            dist2 <- dist2[dist2 >= flank[1] & dist2 <= flank[2]]

            # Bin distances and count matrix
            cts1 <- table(cut(dist1, breaks = breaks))
            cts2 <- table(cut(dist2, breaks = breaks))
            ctsMat <- as.matrix(rbind(cts1, cts2))
            rownames(ctsMat) <- c("pos", "neg")
            
            # Plot
            title <- sprintf(
                "Relative to %s\nN(%s) = %i, N(%s) = %i\nbw = %i nt ",
                region, 
                GetId(txLoc1), sum(cts1),
                GetId(txLoc2), sum(cts2),
                binWidth)
            PlotEnrichment.Generic(
                ctsMat,
                title = title,
                xlab = sprintf("Distance relative to %s [nt]", region),
                x.las = 1, x.cex = 0.8, x.padj = 0,
                xAxisLblFmt = 3)
            
        },
        lstDist1, lstDist2, names(lstDist1)))

}


#' Plot overlap of sites.
#'
#' Plot overlap of sites from two \code{txLoc} objects. See 'Details'.
#'
#' The function plots one or multiple Venn diagrams showing the spatial 
#' overlap between entries from two \code{txLoc} objects. Two features are 
#' defined as overlapping, if they overlap by at least one nucleotide. 
#' Overlaps are determined using the function 
#' \code{GenomicRanges::countOverlaps}.
#'
#' @param txLoc1 A \code{txLoc} object.
#' @param txLoc2 A \code{txLoc} object.
#'
#' @return \code{NULL}.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#'
#' @import GenomicRanges IRanges
#' @importFrom gplots venn
#' @importFrom graphics title
PlotOverlap <- function(txLoc1, txLoc2) {
    
    # Sanity check
    CheckClassTxLocConsistency(txLoc1, txLoc2)
    
    # Determine figure panel layout
    if (length(GetRegions(txLoc1)) < 4) {
        par(mfrow = c(1, length(GetRegions(txLoc1))))
    } else {
        par(mfrow = c(ceiling(length(GetRegions(txLoc1)) / 2), 2))
    }
    
    invisible(Map(
        function(loci1, loci2, region) {
            
            # Get transcriptome coordinates
            gr1 <- loci1$locus_in_tx_region
            gr2 <- loci2$locus_in_tx_region
            
            # Make sure that seqlevels match
            lvls <- intersect(seqlevels(gr1), seqlevels(gr2))
            seqlevels(gr1, pruning.mode = "coarse") <- lvls
            seqlevels(gr2, pruning.mode = "coarse") <- lvls
            
            # Count overlap
            m <- countOverlaps(gr1, gr2)
            overlap <- sum(m > 0)
            
            # Plot
            grps <- list(
                seq(1, length(gr1)),
                seq(length(gr1) - overlap + 1, length.out = length(gr2)))
            names(grps) <- c(
                sprintf(
                    "%s (%3.2f%%)",
                    GetId(txLoc1), overlap / length(gr1) * 100),
                sprintf(
                    "%s (%3.2f%%)",
                    GetId(txLoc2), overlap / length(gr2) * 100))
            venn(grps)
            title(region)
            mtext(names(gr1))
            
        },
        GetLoci(txLoc1),
        GetLoci(txLoc2),
        GetRegions(txLoc1)))
    
}


