#' Input proteome data
#'
#' This file is obtained from MaxQuant as proteinGroups.txt. Column information is not explained as these can be
#' obtained in the summary.pdf file from MaxQuant output. This corresponding input proteome data is obtained from a
#' recently published study comparing the whole cell extracts of proteomes between wild-type and chameau RNAi flies
#' (Further details can be found in the link provided below)
#'
#' @docType data
#'
#' @usage data(prot.Data)
#'
#' @format An object of class \code{"data.frame"}
#' @references Venkatasubramani AV, et al., (2023), EMBO Rep, 24(10):e57023
#' (\href{https://www.embopress.org/doi/full/10.15252/embr.202357023}{EMBO Reports})
#'
#' @source (\href{http://www.ebi.ac.uk/pride/archive/projects/PXD042471})
#'
#' @keywords input
#' @examples
#'
#' data(prot.Data)
#' head(prot.Data)
#'
"prot.Data"
