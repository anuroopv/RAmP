#' File containing sample information
#'
#' .xlsx file containing all the required information about the all samples that
#' will be used for the analysis. The file contains details of both
#' input proteome and modified proteome data
#'
#' @docType data
#'
#' @usage data(sampleTable)
#'
#' @format An object of class \code{"data.frame"}
#' \describe{
#'  \item{label}{Mandatory column. Exact column names of LFQ (or iBAQ) and Intensity columns from MaxQuant output files}
#'  \item{sample_name}{Contains sample name}
#'  \item{condition}{Contains complete details of different conditions}
#'  \item{replicate}{Mandatory column. Replicate information of samples from different conditions}
#'  \item{batch}{Mandatory column. Batch information of samples from different conditions. Especially important for batch correction (if samples were prepared or analyzed in different batches)}
#'  \item{fraction}{Mandatory column. Contains information if the label corresponds to "Proteome" or "Enriched"}
#'  \item{genotype}{Optional column. Genotype information if required for PCA (Additional columns can be added, for example, time)}
#' }
#' @references This data set was manually created RAmP package
#' @keywords sample
#' @examples
#'
#' data(sampleTable)
#' head(sampleTable)
#'
"sampleTable"
