###########################################################################
########################### CLASS DEFINITION ##############################
###########################################################################

#' txLoc object.

#' The S4 object stores information about mapped loci per transcript
#' section. Meta-data such as an identifier and the reference genome
#' are stored in separate slots.
#' 
#' @slot loci A \code{list} of \code{dataframe} objects; specifies
#' the list of loci per transcript section.
#' @slot id A character string; identifier for loci in \code{txLoc}
#' object.
#' @slot refGenome A character string; gives the reference genome
#' upon which transcriptome-derived positions are based.
#' @slot version A character string; can be used to flag specific
#' version, e.g. using the current system time & date.
#'
#' @name txLoc-class
#' @rdname txLoc-class
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportClass txLoc
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


#' Method "GetLoci" for S4 object txLoc.
#'
#' Get the loci of a \code{txLoc} object
#' as a list of dataframes.
#'
#' @return A list of dataframes with loci for
#' every transcript section.
#'
#' @param x A \code{txLoc} object.
#' 
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod GetLoci
setGeneric("GetLoci", function(x) standardGeneric("GetLoci"));

#' @rdname GetLoci
setMethod("GetLoci",
          signature = "txLoc",
          definition = function(x) {
              loc <- slot(x, "loci");
              return(loc);
          });


#' Method "GetId" for S4 object txLoc.
#'
#' Get the id string of a \code{txLoc} object.
#'
#' @param x A \code{txLoc} object.
#' 
#' @return A string.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod GetId
setGeneric("GetId", function(x) standardGeneric("GetId"));

#' @rdname GetId
setMethod("GetId",
          signature = "txLoc",
          definition = function(x) {
              id <- slot(x, "id");
              return(id);
          });


#' Method "GetRef" for S4 object txLoc.
#'
#' Get the reference genome string of a \code{txLoc} object.
#'
#' @param x A \code{txLoc} object.
#' 
#' @return A string.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod GetRef
setGeneric("GetRef", function(x) standardGeneric("GetRef"));

#' @rdname GetRef
setMethod("GetRef",
          signature = "txLoc",
          definition = function(x) {
              ref <- slot(x, "refGenome");
              return(ref);
          });


#' Method "GetVersion" for S4 object txLoc.
#'
#' Get the version string of a \code{txLoc} object.
#'
#' @param x A \code{txLoc} object.
#' 
#' @return A string.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod GetVersion
setGeneric("GetVersion", function(x) standardGeneric("GetVersion"));

#' @rdname GetVersion
setMethod("GetVersion",
          signature = "txLoc",
          definition = function(x) {
              version <- slot(x, "version");
              return(version);
          });


###########################################################################
###################### SPECIFIC CLASS ACCESSORS ###########################
###########################################################################

#' Method "Info" for S4 object txLoc.
#'
#' Print general information of a \code{txLoc} object.
#'
#' @param x A \code{txLoc} object.
#' 
#' @examples
#' \dontrun{
#' bedFile <- system.file("extdata",
#'                        "miCLIP_m6A_Linder2015_hg38.bed",
#'                        package = "RNAModR");
#' sites <- ReadBED(bedFile);
#' posSites <- SmartMap(sites, id = "m6A", refGenome = "hg38");
#' Info(posSites);
#' }
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod Info
setGeneric("Info", function(x) standardGeneric("Info"));

#' @rdname Info
setMethod("Info",
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


#' Method "GetNumberOfLoci" for S4 object txLoc.
#'
#' Get the number of loci of a \code{txLoc} object in every
#' transcript section.
#'
#' @param x A \code{txLoc} object.
#' 
#' @return A named integer vector with the number of sites per
#' transcript section.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod GetNumberOfLoci
setGeneric("GetNumberOfLoci", function(x) standardGeneric("GetNumberOfLoci"));

#' @rdname GetNumberOfLoci
setMethod("GetNumberOfLoci",
          signature = "txLoc",
          definition = function(x) {
              loc <- slot(x, "loci");
              return(sapply(loc, nrow));
          });


###########################################################################
#################### GENERAL CLASS SLOT MUTATORS ##########################
###########################################################################

##   #' Method "SetId" for S4 object txLoc.
##   #'
##   #' Set the id string of a \code{txLoc} object.
##   #'
##   #' @param x A \code{txLoc} object.
##   #' @param id A character string.
##   #'
##   #' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
##   #' 
##   #' @exportMethod SetId
##   setGeneric("SetId", function(x, id) standardGeneric("SetId"));
##   
##   #' @rdname SetId
##   setMethod("SetId",
##             signature = "txLoc",
##             definition = function(x, id) {
##                 slot(x, "id") <- id;
##             });
##   
##   
##   #' Method "SetVersion" for S4 object txLoc.
##   #'
##   #' Set the version string of a \code{txLoc} object.
##   #'
##   #' @param x A \code{txLoc} object.
##   #' @param version A character string.
##   #'
##   #' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
##   #' 
##   #' @exportMethod SetVersion
##   setGeneric("SetVersion", function(x, version) standardGeneric("SetVersion"));
##   
##   #' @rdname SetVersion
##   setMethod("SetVersion",
##             signature = "txLoc",
##             definition = function(x, version) {
##                 slot(x, "version") <- version;
##             });
