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
#' @source \href{https://www.embopress.org/doi/full/10.15252/embr.202357023}{EMBO Reports}
#'
#'
#' @keywords input
#' @examples
#'
#' data(prot.Data)
#' head(prot.Data)
#'
"prot.Data"

#' Modified proteome data
#'
#' This file is obtained from MaxQuant as Acetyl(K)Sites.txt. Column information is not explained as these can be
#' obtained in the summary.pdf file from MaxQuant output This corresponding acetylome data is obtained from a
#' recently published study comparing the whole cell extracts of acetylomes between wild-type and chameau RNAi flies
#' (Further details can be found in the link provided below)
#'
#' @docType data
#'
#' @usage data(enrich.Data)
#'
#' @format An object of class \code{"data.frame"}
#' @references Venkatasubramani AV, et al., (2023), EMBO Rep Oct 9;24(10):e57023
#' @source \href{https://www.embopress.org/doi/full/10.15252/embr.202357023}{EMBO Reports}
#'
#'
#' @keywords enriched
#' @examples
#'
#' data(enrich.Data)
#' head(enrich.Data)
#'
"enrich.Data"


#' File from UNIPROT containing protein sequences from Drosophila melanogaster. Exact fasta file used also in MaxQuant search
#' Contains sequence information of all proteins from Drosophila melanogaster in fasta format. Obtained from UNIPROT
#' @docType data
#'
#' @usage data(fasta)
#'
#' @format An object of class \code{"Large character"}
#'
#' @source \href{https://www.uniprot.org/proteomes/UP000000803}{UNIPROT}
#'
#'
#' @keywords fasta
#' @examples
#'
#' data(fasta)
#' head(fasta)
#'
"fasta"
