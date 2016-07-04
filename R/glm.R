#' @importFrom stats binomial confint.default glm poisson
#' @importFrom graphics points
testGLM <- function(locPos,
                    locNeg,
                    filter = NULL,
                    binWidth = 20,
                    posMax = 1000) {
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
        for (j in 1:length(posPos)) {
            # Generate bootstrap samples
            nBS <- 1000;
            mat1 <- matrix(sample(posPos[[j]], 
                                  size = nBS * length(posPos[[j]]), 
                                  replace = TRUE),
                           ncol = nBS);
            mat1 <- apply(mat1, 2, function(x) {table(cut(x, breaks = breaks))});
            mat2 <- matrix(sample(posNeg[[j]], 
                                  size = nBS * length(posPos[[j]]), 
                                  replace = TRUE),
                            ncol = nBS);
            mat2 <- apply(mat2, 2, function(x) {table(cut(x, breaks = breaks))});
            # Poisson regression
            df <- data.frame(matrix(0, ncol = 5, nrow = nrow(mat1)));
            df[ ,1] <- breaks[-length(breaks)] + binWidth / 2;
            for (k in 1:nrow(mat1)) {
                x <- factor(c(rep(idPos, nBS), rep(idNeg, nBS)),
                            levels = c(idNeg, idPos));
                y <- c(mat1[k, ], mat2[k, ]);
                fit <- glm(y ~ x, family = poisson(link = "log"));
                df[k, 2] <- summary(fit)$coefficients[2, 1];
                df[k, c(3, 4)] <- confint.default(fit, level = 0.95)[2, ];
                df[k, 5] <- log10(summary(fit)$coefficients[2, 4]);
            }
            if (j == 1) {
                xmin <- 0;
                xmax <- posMax;
            } else {
                xmin <- posMax;
                xmax <- 0;
            }
            plot(df[, 1], df[, 2],
                 xlim = c(xmin, xmax),
                 ylim = c(-round(max(abs(df[, 2]))), round(max(abs(df[, 2])))),
                 col = rgb(1, 0.2, 0.2, 1), type = "l",
                 lwd = 2,
                 xlab = "Absolute position [nt]",
                 ylab = "log(Fold change)");
            CI <- cbind(c(df[, 1], rev(df[, 1])),
                        c(df[, 3], rev(df[, 4])));
            polygon(CI[,1], CI[,2], col = rgb(1, 0, 0, 0.2),
                    lwd = 1, border = NA, lty = 1);
            abline(h = 0.0, col = "red", lty = 3, lwd = 2);
            points(df[, 1], df[, 5],
                   type = "h",
                   lwd = 10,
                   lend = "butt",
                   col = rgb(0.2, 0.2, 1, 0.1));
        }
    }
}
