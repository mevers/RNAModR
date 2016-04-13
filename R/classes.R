###########################################################################
########################### CLASS DEFINITION ##############################
###########################################################################

#' txLoc object.

#' The S4 object stores information about mapped loci per transcript
#' section. Meta-data such as an identifier and the reference genome
#' are stored in separate slots.
#' 
#' @slot loci A \code{list} of \code{GRangesList} objects; specifies
#' the list of loci per transcript section.
#' @slot id A character string; identifier for loci in \code{txLoc}
#' object.
#' @slot refGenome A character string; gives the reference genome
#' upon which transcriptome-derived positions are based.
#' @slot version A character string; can be used to flag specific
#' version, e.g. using the current system time & date.
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

###########################################################################
#################### GENERAL CLASS SLOT ACCESSORS #########################
###########################################################################

#' Generic method "GetLoci" for S4 object txLoc.
#' 
#' @keywords internal
#' 
#' @export
setGeneric(name = "GetLoci",
           def = function(x) {
               standardGeneric("GetLoci");
           });


#' Method "GetLoci" for S4 object txLoc.
#'
#' Get the loci of a \code{txLoc} object
#' as a list of dataframes.
#'
#' @return A list of dataframes with loci for
#' every transcript section.
#'
#' @keywords internal
#' 
#' @export
setMethod("GetLoci",
          signature = "txLoc",
          definition = function(x) {
              loc <- slot(x, "loci");
              return(loc);
          });


#' Generic method "GetId" for S4 object txLoc.
#' 
#' @keywords internal
#' 
#' @export
setGeneric(name = "GetId",
           def = function(x) {
               standardGeneric("GetId");
           });


#' Method "GetId" for S4 object txLoc.
#'
#' Get the id string of a \code{txLoc} object.
#'
#' @return A string.
#'
#' @keywords internal
#' 
#' @export
setMethod("GetId",
          signature = "txLoc",
          definition = function(x) {
              id <- slot(x, "id");
              return(id);
          });


#' Generic method "GetRef" for S4 object txLoc.
#' 
#' @keywords internal
#' 
#' @export
setGeneric(name = "GetRef",
           def = function(x) {
               standardGeneric("GetRef");
           });


#' Method "GetRef" for S4 object txLoc.
#'
#' Get the reference genome string of a \code{txLoc} object.
#'
#' @return A string.
#'
#' @keywords internal
#' 
#' @export
setMethod("GetRef",
          signature = "txLoc",
          definition = function(x) {
              ref <- slot(x, "refGenome");
              return(ref);
          });


#' Generic method "GetVersion" for S4 object txLoc.
#' 
#' @keywords internal
#' 
#' @export
setGeneric(name = "GetVersion",
           def = function(x) {
               standardGeneric("GetVersion");
           });


#' Method "GetVersion" for S4 object txLoc.
#'
#' Get the version string of a \code{txLoc} object.
#'
#' @return A string.
#'
#' @keywords internal
#' 
#' @export
setMethod("GetVersion",
          signature = "txLoc",
          definition = function(x) {
              version <- slot(x, "version");
              return(version);
          });


###########################################################################
###################### SPECIFIC CLASS ACCESSORS ###########################
###########################################################################

#' Generic method "info" for S4 object txLoc.
#' 
#' @rdname txLoc-class
#'
#' @export
setGeneric(name = "info",
           def = function(x) {
               standardGeneric("info");
           });


#' Method "info" for S4 object txLoc.
#'
#' Print general information of a \code{txLoc} object.
#'
#' @rdname txLoc-class
#'
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' info(posSites);
#' }
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
              cat(sprintf("Total # of sites = %i\n",
                          sum(sapply(slot(x, "loci"), nrow))));
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


#' Generic method "head" for S4 object \code{txLoc}.
setGeneric("head");


#' Method "head" for S4 object txLoc.
#'
#' @rdname txLoc-class
#'
#' @param x A \code{txLoc} object.
#' @param n Number of rows to be printed.
#' @param ... Additional parameters passed to \code{head}.
#' 
#' @export
setMethod("head",
          signature = "txLoc",
          definition = function(x, n = 6L, ...) {
              loc <- slot(x, "loci");
              nSec <- length(loc);
              for (i in 1:nSec) {
                  cat(sprintf("$%s\n", names(loc)[i]));
                  locPrint <- loc[[i]];
                  locPrint[, "REGION_SEQ"] <- sprintf(
                      "%s... (truncated, total length = %i)",
                      substr(locPrint[, ncol(locPrint)], 1, 5),
                      nchar(locPrint[, ncol(locPrint)]));
                  print(head(locPrint, n = n));
                  cat(sprintf("\n", names(loc)[i]));
              }
          });


#' Generic method "GetNumberOfLoci" for S4 object txLoc.
#'
#' @rdname txLoc-class
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
#' @rdname txLoc-class
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
