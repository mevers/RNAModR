#' txLoc object definition
#' 
#' @rdname txLoc
#'
#' @slot loci A list of loci within a transcript section.
#' @slot id Identifier for sites in \code{txLoc} object.
#' @slot refGenome Reference genome upon which the transcriptome
#' and site positions are based.
#' @slot version Timestamp for version tracking.
#' 
#' @export
setClass("txLoc",
         representation(loci = "list",
                        id = "character",
                        refGenome = "character",
                        version = "character"),
         prototype = prototype(
             loci = list(),
             id = "",
             refGenome = "",
             version = ""),
         );


#' Generic method "info" for S4 object txLoc.
#' 
#' @rdname txLoc
#'
#' @param x A \code{txLoc} object.
#' 
#' @export
setGeneric(name = "info",
           def = function(x) {
               standardGeneric("info");
           });


#' Method "print" for S4 object txLoc.
#'
#' Print general information of a \code{txLoc} object.
#'
#' @rdname txLoc
#'
#' @param x A \code{txLoc} object.
#'
#' @examples
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' info(posSites);
#' 
#' @export
setMethod("info",
          signature = "txLoc",
          definition = function(x) {
              cat("Object of class \"txLoc\".\n");
              cat("\n");
              cat(sprintf("ID               = %s\n", slot(x, "id")));
              cat(sprintf("Reference genome = %s\n", slot(x, "refGenome")));
              cat(sprintf("Version          = %s\n", slot(x, "version")));
              cat(sprintf("Package          = RNAModR"));
              cat("\n");
              loc <- slot(x, "loci");
              nSec <- length(loc);
              cat(sprintf("%i transcript sections: %s\n",
                          nSec,
                          paste0(names(loc), collapse = ", ")));
              for (i in 1:nSec) {
                  cat(sprintf(" %10s: Number of loci = %i\n", names(loc)[i],
                              nrow(loc[[i]])));
              }
          });


#' Method "head" for S4 object txLoc.
#'
#' @rdname txLoc
#'
#' @param x A \code{txLoc} object.
#' @param n Number of rows to be printed.
#' 
#' @export
setMethod("head",
          signature = "txLoc",
          definition = function(x, n = 6L) {
              loc <- slot(x, "loci");
              nSec <- length(loc);
              for (i in 1:nSec) {
                  cat(sprintf("$%s\n", names(loc)[i]));
                  locPrint <- loc[[i]];
                  locPrint[, ncol(locPrint)] <- sprintf(
                      "%s... (truncated, total length = %i)",
                      substr(locPrint[, ncol(locPrint)], 1, 5),
                      nchar(locPrint[, ncol(locPrint)]));
                  print(head(locPrint, n = n));
                  cat(sprintf("\n", names(loc)[i]));
              }
          });


#' Generic method "GetNumberOfLoci" for S4 object txLoc.
#' 
#' @rdname txLoc
#'
#' @param x A \code{txLoc} object.
#' 
#' @export
setGeneric(name = "GetNumberOfLoci",
           def = function(x) {
               standardGeneric("GetNumberOfLoci");
           });


#' Method "GetNumberOfLoci" for S4 object txLoc.
#'
#' Get the number of loci of a \code{txLoc} object in every
#' transcript section.
#'
#' @rdname txLoc
#'
#' @param x A \code{txLoc} object.
#'
#' @return A named integer vector with the number of sites per
#' transcript section.
#'
#' @export
setMethod("GetNumberOfLoci",
          signature = "txLoc",
          definition = function(x) {
              loc <- slot(x, "loci");
              return(sapply(loc, nrow));
          });
