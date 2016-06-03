#' @importFrom stats binomial glm
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
    breaks <- seq(0, posMax, by = binWidth);
    for (i in 1:length(locPos)) {
        posPos <- list("5p" = locPos[[i]]$TXSTART,
                       "3p" = locPos[[i]]$REGION_TXWIDTH - locPos[[i]]$TXSTART + 1);
        posNeg <- list("5p" = locNeg[[i]]$TXSTART,
                       "3p" = locNeg[[i]]$REGION_TXWIDTH - locNeg[[i]]$TXSTART + 1);
        posPos <- lapply(posPos, function(x) x[x <= posMax]);
        posNeg <- lapply(posNeg, function(x) x[x <= posMax]);
        for (j in 1:length(posPos)) {
            ctsPos <- table(cut(posPos[[j]], breaks = breaks));
            ctsNeg <- table(cut(posNeg[[j]], breaks = breaks));
            ctsMat <- as.matrix(rbind(ctsPos, ctsNeg));
            rownames(ctsMat) <- c(idPos, idNeg);
            # Generate bootstrap samples
            nBS <- 1000;
            mat1 <- matrix(sample(posPos[[j]], 
                                  size = nBS * length(posPos[[j]]), 
                                  replace = TRUE),
                            ncol = nBS);
            mat1 <- apply(mat1, 2, function(x) {table(cut(x, breaks = breaks))});
            mat2 <- matrix(sample(posNeg[[j]], 
                                  size = 1000 * length(posPos[[j]]), 
                                  replace = TRUE),
                            ncol = 1000);
            mat2 <- apply(mat2, 2, function(x) {table(cut(x, breaks = breaks))});
            #plot(rowMeans(mat2), col = "blue", ylim = c(0, 200));
            #points(rowMeans(mat1), col = "red");
            # Logistic regression
            for (k in 1:nrow(mat1)) {
              y <- cbind(mat1[k, ], mat2[k, ]);
              x <- factor(c(rep(idPos, nBS), rep(idNeg, nBS)), 
                          levels = c(idPos, idNeg));
              y <- factor(c(rep(idPos, nBS), rep(idNeg, nBS)), 
                          levels = c(idPos, idNeg));
              y <- c(rep(0, nBS), rep(1, nBS))
              x <- c(mat1[k, ], mat2[k, ]);
              fit <- glm(y ~ x, family = binomial("logit"));
            }
        }
    }
}
