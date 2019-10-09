###########################################################################
########################### CLASS DEFINITION ##############################
###########################################################################

#' \code{txLoc} object.
#' 
#' The S4 object stores information about mapped loci per transcript region.
#' 
#' The \code{txLoc} S4 object contains the following four slots:
#' \enumerate{
#'     \item \code{loci}: A \code{list} of \code{DataFrame} objects, each with 
#'     the following columns:
#'     \itemize{
#'         \item \code{locus_in_txx_region}: A \code{GRanges} object
#'         \item \code{locus_in_genome}: A \code{GRanges} object
#'         \item \code{score}: A \code{numeric} vector
#'         \item \code{id}: A \code{character} vector
#'         \item \code{tx_region}: A \code{character} vector
#'         \item \code{tx_region_width}: An \code{integer} vector
#'         \item \code{tx_region_sequence}: A \code{DNAStringSet} object
#'         \item \code{tx_refseq}: A \code{character} vector
#'         \item \code{gene_entrez}: A \code{character} vector
#'         \item \code{gene_symbol}: A \code{character} vector
#'         \item \code{gene_ensembl}: A \code{character} vector
#'         \item \code{gene_name}: A \code{character} vector
#'     }
#'     The \code{loci} slot can be accessed using \code{GetLoci(txLoc)}. 
#'     \item \code{id}: An identifier specified by the user. The \code{id} slot
#'     can be accessed using \code{GetId(txLoc)}.
#'     \item \code{refGenome}: The reference genome version (e.g. 
#'     \code{"hg38"}), which determines the mapping between genomic and 
#'     transcriptomic coordinates. The \code{refGenome} slot can be accessed 
#'     using \code{GetRef(txLoc)}.
#'     \item `version`: A version identifier; currently this slot is used to
#'     to store the current system time & date. The \code{version} slot can be
#'     accessed using \code{GetVersion(txLoc)}.
#' }
#' 
#' @slot loci A \code{list} of \code{DataFrame} objects; see 'Details'.
#' @slot id A \code{character} string; see 'Details.
#' @slot refGenome A \code{character} string; see 'Details'.
#' @slot version A \code{character} string; see 'Details'.
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
             id = NA_character_,
             refGenome = NA_character_,
             version = NA_character_)
)


###########################################################################
########################### CLASS SLOT ACCESSORS ##########################
###########################################################################

#' \code{txLoc} accessors.
#' 
#' Various \code{txLoc} accessors.
#' 
#' \itemize{
#'     \item \code{show(object)}: Print summary information about the 
#'     \code{txLoc} object.
#'     \item \code{GetLoci(object)}: Get a \code{list} of \code{DataFrame} 
#'     objects of the loci of the \code{txLoc} object.
#'     \item \code{GetId(object)}: Get the identifier string of the \code{txLoc} 
#'     object.
#'     \item \code{GetRef(object)}: Get the reference genome string of the
#'     \code{txLoc} object.
#'     \item \code{GetVersion(object)}: Get the version string of a \code{txLoc} 
#'     object.
#'     \item \code{GetNumberOfLoci(object)}: Get the number of loci of a 
#'     \code{txLoc} object in every transcript section.
#'     \item \code{GetRegions(objects)}: Character vector of the transcript 
#'     regions of a \code{txLoc} object.
#' }
#' 
#' @param object A \code{txLoc} object.
#' @name txLoc-accessors
NULL


#' Method \code{show} for S4 object \code{txLoc}.
#' 
#' @author Maurits Evers, \email{maurits.evers@anu.edu.au}
#' 
#' @importFrom utils packageVersion
#' 
#' @exportMethod show
#' 
#' @rdname txLoc-accessors
setMethod(
    "show", 
    signature = "txLoc", 
    function(object) {
        cat("Object of class \"txLoc\".\n")
        cat("\n")
        cat(sprintf("ID               = %s\n", slot(object, "id")))
        cat(sprintf("Reference genome = %s\n", slot(object, "refGenome")))
        cat(sprintf("Version          = %s\n", slot(object, "version")))
        cat(sprintf("Total # of sites = %i\n",
                    sum(sapply(slot(object, "loci"), nrow))))
        cat(sprintf("Package          = RNAModR (%s)", utils::packageVersion("RNAModR")))
        cat("\n")
        loc <- slot(object, "loci")
        nSec <- length(loc)
        cat(sprintf("%i transcript sections: %s\n",
                    nSec,
                    paste0(names(loc), collapse = ", ")))
        for (i in 1:nSec) {
            cat(sprintf(" %10s: Number of loci = %i\n", names(loc)[i],
                        nrow(loc[[i]])))
        }
    }
)


#' Method \code{GetLoci} for S4 object \code{txLoc}.
#'
#' Get a \code{list} of \code{DataFrame} objects of the loci of a \code{txLoc} 
#' object.
#'
#' @return A \code{list} of \code{DataFrame} objects with loci for every
#' transcript region.
#'
#' @param object A \code{txLoc} object.
#' @keywords internal
#' 
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' 
#' @exportMethod GetLoci
setGeneric("GetLoci", function(object) standardGeneric("GetLoci"))

#' @rdname txLoc-accessors
setMethod("GetLoci",
          signature = "txLoc",
          definition = function(object) slot(object, "loci"))


#' Method \code{GetId} for S4 object \code{txLoc}.
#'
#' Get the identifier string of a \code{txLoc} object.
#'
#' @param object A \code{txLoc} object.
#' 
#' @return A string.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @exportMethod GetId
setGeneric("GetId", function(object) standardGeneric("GetId"))

#' @rdname txLoc-accessors
setMethod("GetId",
          signature = "txLoc",
          definition = function(object) slot(object, "id"))


#' Method \code{GetRef} for S4 object \code{txLoc}.
#'
#' Get the reference genome string of a \code{txLoc} object.
#'
#' @param object A \code{txLoc} object.
#' 
#' @return A string.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @exportMethod GetRef
setGeneric("GetRef", function(object) standardGeneric("GetRef"))

#' @rdname txLoc-accessors
setMethod("GetRef",
          signature = "txLoc",
          definition = function(object) slot(object, "refGenome"))


#' Method \code{GetVersion} for S4 object \code{txLoc}.
#'
#' Get the version string of a \code{txLoc} object.
#'
#' @param object A \code{txLoc} object.
#' 
#' @return A string.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @exportMethod GetVersion
setGeneric("GetVersion", function(object) standardGeneric("GetVersion"))

#' @rdname txLoc-accessors
setMethod("GetVersion",
          signature = "txLoc",
          definition = function(object) slot(object, "version"))


#' Method \code{GetNumberOfLoci} for S4 object \code{txLoc}.
#'
#' Get the number of loci of a \code{txLoc} object in every transcript section.
#'
#' @param object A \code{txLoc} object.
#' 
#' @return A named integer vector with the number of sites per
#' transcript section.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @exportMethod GetNumberOfLoci
setGeneric("GetNumberOfLoci", function(object) 
    standardGeneric("GetNumberOfLoci"))

#' @rdname txLoc-accessors
setMethod("GetNumberOfLoci",
          signature = "txLoc",
          definition = function(object) sapply(slot(object, "loci"), nrow))


#' Method \code{GetRegions} for S4 object \code{txLoc}.
#'
#' Character vector of the transcript regions of a \code{txLoc} object.
#'
#' @param object A \code{txLoc} object.
#' 
#' @return A \code{character} vector.
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @keywords internal
#' 
#' @exportMethod GetRegions
setGeneric("GetRegions", function(object) standardGeneric("GetRegions"))

#' @rdname txLoc-accessors
setMethod("GetRegions",
          signature = "txLoc",
          definition = function(object) names(GetLoci(object)))

