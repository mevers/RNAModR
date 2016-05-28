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
            mat1 <- matrix(sample(posPos[[j]], 
                                  size = 1000 * length(posPos[[j]]), 
                                  replace = TRUE),
                            ncol = 1000);
            mat1 <- apply(mat1, 2, function(x) {table(cut(x, breaks = breaks))});
            colnames(mat1) <- sprintf("%s_%i", idPos, seq(1, 1000));
            mat2 <- matrix(sample(posNeg[[j]], 
                                  size = 1000 * length(posNeg[[j]]), 
                                  replace = TRUE),
                            ncol = 1000);
            mat2 <- apply(mat2, 2, function(x) {table(cut(x, breaks = breaks))});
            colnames(mat2) <- sprintf("%s_%i", idNeg, seq(1, 1000));
            # Count data and design matrix
            countData <- cbind(mat1, mat2);
            nSamples <- ncol(countData);
            colData <- data.frame(condition = c(rep(idPos, 1000), 
                                                rep(idNeg, 1000)));
            dds <- DESeqDataSetFromMatrix(countData = countData,
                                          colData = colData,
                                          design = ~condition);
            dds$condition <- relevel(dds$condition, ref = idNeg);
            # This is very very slow.
            # Need to find a better and faster method.
            # Logistic regression ...
            dds <- DESeq(dds);
            res <- results(dds);
        }
    }
}
