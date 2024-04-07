#' @title editData
#'
#' @description Cleans the raw data from MaxQuant for downstream analysis
#'
#' @param data Input or modified proteome file from MaxQuant (.txt)
#' @param fraction Can be either "Proteome" or "Enriched". Indicates the type of input data used
#' @param probability Numeric value between 0-1. Filters out modified peptides with probabilities less than the given value (Only used if fraction = "Enriched")
#' @param org Database of the organism. Drosophila melanogaster = "dme", Mus muscuslus ' "mmu", Homo sapiens = "hsa", Saccharomyces cerevisae = "sce". Default is "dme"
#' @param quantification Default is LFQ. Can be either "LFQ" or "iBAQ"
#'
#' @return A data frame that will be cleaner than the original raw data to be used for differential expression analysis
#'
#' @examples
#' lfqdata <- editData(data = proteome, fraction = "Enriched", probability = 0.9)
#' @export
#' @importFrom dplyr "%>%"

################# Data preparation for differential expression analysis ###################

editData <- function(data, fraction = c("Proteome", "Enriched"), probability = NULL, org = "dme", quantification = "LFQ"){

  # Decide the organism database

  if(org == "dme"){
    orgDB = org.Dm.eg.db
  }else if(org == "hsa"){
    orgDB = org.Hs.eg.db
  }else if(org == "mmu"){
    orgDB = org.Mm.eg.db
  }else if(org == "sce"){
    orgDB = org.Sc.sgd.db
  }else{
    stop("Only drosophila, human, mouse and yeast databases are supported")
  }

  if(fraction == "Proteome"){

    # Remove contaminants
    #data <- data[!(data$Only.identified.by.site=="+" | data$Reverse=="+" | data$Potential.contaminant=="+"),]
    if(quantification != "DIA-NN"){
      data.sub <- subset(data, data$Only.identified.by.site=="+" | data$Reverse=="+" | data$Potential.contaminant=="+")
      data <- data[!data$Protein.IDs %in% data.sub$Protein.IDs,]

    # Edit data to obtain symbols from Uniprot
    data$Uniprot <- gsub("((sp\\|)|(tr\\|))", "", data$Fasta.headers)
    data$Uniprot  <- gsub("[|].*", "", data$Uniprot)
    }

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

    } else if(quantification == "DIA-NN"){
      lfq.data <- data[,c(6:ncol(data))]
      lfq.data$Uniprot <- data$Protein.Ids
      lfq.data$Uniprot <- gsub(";.*", "", lfq.data$Uniprot)
      if(org != "sce"){
        lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
      }else {
        lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
      }

      }else {
      stop("Values can be either LFQ or iBAQ or DIA-NN")
    }

    lfq.data[lfq.data==0] <- NA # Convert zeros to NA
    # if(org != "sce"){
    #   lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
    # }else {
    #   lfq.data$symbol <- mapIds(x = orgDB, keys =  as.character(lfq.data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
    # }

  } else if(fraction == "Enriched" & quantification != "DIA-NN"){

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

  }else if(fraction == "Enriched" & quantification == "DIA-NN"){
    stop("Analysis of modified proteomes with DIA-NN output is currently not available")
  }
    else{
    stop("Accepted values are Proteome or Enriched")
  }
  return(lfq.data)
}
