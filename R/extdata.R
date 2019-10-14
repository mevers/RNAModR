#' Internal raw datasets.
#'
#' Internal raw datasets in BED6 format. All genome positions where converted 
#' to the most recent genome assembly version; i.e. for human data, all
#' positions are given in GRCh38/hg38 coordinates. For mouse data, positions
#' are given in mm10 coordinates.
#' 
#' \describe{
#'     \item{\code{bsRNAseq_m5C_Squires2012_hg38.bed}}{
#'     bsRNA-seq 5-methylcytosine (m5C) sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/22344696}}
#'     \item{\code{HITSCLIP_eIF4A3_Sauliere2012_hg38.bed}}{
#'     HITSCLIP eIF4AIII binding sites. 
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/23085716}}
#'     \item{\code{MeRIPseq_m1A_Dominissini2016_hg38.bed}}{
#'     MeRIP-seq N(1)-methyladenosine (m1A) sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/26863196}}
#'     \item{\code{miCLIP_m6A_Linder2015_hg38.bed}}{
#'     miCLIP N(6)-methyladenosine (m6A) sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/26121403}}
#'     \item{\code{PARCLIP_AGO1234_Hafner2010_hg38.bed}}{
#'     PAR-CLIP AGO1-4 binding sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/20371350}}
#'     \item{\code{PARCLIP_eIF3_Meyer2015_hg38.bed}}{
#'     PAR-CLIP eIF3 binding sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/26593424}}
#'     \item{\code{PARCLIP_eIF3_Lee2015_hg38.bed}}{
#'     PAR-CLIP eIF3 binding sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/25849773}}
#'     \item{\code{PARCLIP_YTHDF2_Wang2014_hg38.bed}}{
#'     PAR-CLIP YTHDF2 binding sites.
#'     Reference: \url{https://www.ncbi.nlm.nih.gov/pubmed/24284625}}
#'     \item{\code{PseudoU_Schwartz2014_Li2015_hg38.bed}}{
#'     Merged Psi-seq and CeU-seq pseudouridine (Psi) sites.
#'     References: \url{https://www.ncbi.nlm.nih.gov/pubmed/25219674},
#'     \url{https://www.ncbi.nlm.nih.gov/pubmed/26075521}}
#'     \item{\code{TargetScan_miRNAtargets_hg38.bed}}{
#'     TargetScan-based miRNA target sites.
#'     Reference: \url{http://www.targetscan.org/vert_72/}}
#' }
#'
#' @name Datasets
#' @rdname Datasets
#'
#' @author Maurits Evers, \email{maurits.evers@@anu.edu.au}
#' @examples
#' \dontrun{
#' bedFile <- ReadBED(system.file(
#'     "extdata",
#'     "miCLIP_m6A_Linder2015_hg38.bed",
#'     package = "RNAModR")
#' }
NULL
