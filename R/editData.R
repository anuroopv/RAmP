#' @title editData
#'
#' @description Cleans the raw data from MaxQuant for downstream analysis
#'
#' @param data Input or modified proteome file from MaxQuant (.txt)
#' @param probability Numeric value between 0-1. Filters out modified peptides with probabilities less than the given value (Only used if Fraction = "Enriched")
#' @param Fraction Can be either "Proteome" or "Enriched". Indicates the type of input data used

#' @return A data frame that will be cleaner than the original raw data to be used for differential expression analysis
#'
#' @examples
#' lfqdata <- editData(data = proteome, probability = 0.9)
#'
#' @export
#'
#' @importFrom AnnotationDbi mapIds
#' @import org.Dm.eg.db
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @import org.Sc.sgd.db
################# Data preparation for differential expression analysis ###################

editData <- function(data, Fraction, probability = NULL){

  if(Fraction == "Proteome"){
    dir.create(paste(path1,"/",Fraction,sep = ""), showWarnings = FALSE)
    # Remove contaminants
    #data <- data[!(data$Only.identified.by.site=="+" | data$Reverse=="+" | data$Potential.contaminant=="+"),]
      data.sub <- subset(data, data$Only.identified.by.site=="+" | data$Reverse=="+" | data$Potential.contaminant=="+")
      data <- data[!data$Protein.IDs %in% data.sub$Protein.IDs,]

    # Edit data to obtain symbols from Uniprot
    data$Uniprot <- gsub("((sp\\|)|(tr\\|))", "", data$Fasta.headers)
    data$Uniprot  <- gsub("[|].*", "", data$Uniprot)

    # Obtain LFQ information
    if(quantification == "LFQ"){
      lfq.data <- data[,grep("LFQ.", colnames(data))]
      lfq.data$Uniprot <- data$Uniprot
      if(org != "sce"){
        lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
      }else {
        lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
      }

    } else if(quantification ==  "iBAQ"){
      lfq.data <- data[,grep("iBAQ.", colnames(data))]
      lfq.data[,1] <- NULL
      lfq.data$Uniprot <- data$Uniprot
      if(org != "sce"){
        lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
      }else {
        lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
      }

    }else {
      stop("Values can be either LFQ or iBAQ")
    }

    lfq.data[lfq.data==0] <- NA # Convert zeros to NA
    # if(org != "sce"){
    #   lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
    # }else {
    #   lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
    # }

  } else if(Fraction == "Enriched"){
    dir.create(paste(path1,"/",Fraction,sep = ""), showWarnings = FALSE)
    # Remove contaminants
    #data <- data[!(data$Only.identified.by.site=="+" | data$Reverse=="+" | data$Potential.contaminant=="+"),]
    data.sub <- subset(data, data$Reverse=="+" | data$Potential.contaminant=="+")
    data <- data[!data$Protein %in% data.sub$Protein,]

    # Edit data to obtain symbols from Uniprot
    data$Uniprot <- gsub("((sp\\|)|(tr\\|))", "", data$Fasta.headers)
    data$Uniprot  <- gsub("([|].*)|([;].*)", "", data$Uniprot)

    # Filter based on localization probability
    data.filter <- subset(data, data$Localization.prob >= probability)

    # Obtain intensities
    lfq.data <- data.filter[,grep("Intensity.", colnames(data.filter))]
    lfq.data <- lfq.data[,-grep("__", colnames(lfq.data))]
    lfq.data$Uniprot <- data.filter$Uniprot
    lfq.data$Probability <- data.filter$Localization.prob
    lfq.data$siteScore <- data.filter$Score.diff
    lfq.data$Sequence <- data.filter[grep("..Probabilities", colnames(data.filter))]
    lfq.data[lfq.data==0] <- NA
    if(org != "sce"){
      lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
    }else {
      lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
    }

  }else{
    stop("Accepted values are Proteome or Enriched")
  }
  dir.create(paste(path1,"/",Fraction,"/Processed_data",sep = ""), showWarnings = FALSE)
  writexl::write_xlsx(x = lfq.data, path = paste(path1,"/",Fraction,"/Processed_data/",Fraction,"_processedData.xlsx", sep = ""), col_names = TRUE, format_headers = TRUE)
  return(lfq.data)
}
