#' @title enrich_normalization
#'
#' @description Normalizes modified proteomes to the corresponding input proteomes
#'
#' @param protein.data Input proteome file from MaxQuant (.txt)
#' @param enrich.data Modified proteome file from MaxQuant (.txt)
#' @param enrich.batch If the normalization should be done based on paired replicates or avergae of replicates
#' @param probability Numeric value between 0-1. Filters out modified peptides with probabilities less than the given value
#' @param org Database of the organism. Drosophila melanogaster = "dme", Mus muscuslus ' "mmu", Homo sapiens = "hsa", Saccharomyces cerevisae = "sce". Default is "dme"
#' @param sampleTable .xlsx file containing information about the samples. Three columns are mandatory (label, condition and replicate)
#'
#' @return Filtered and vsn normalized data
#'
#' @export
############### Obtain normalized data for differential expression analysis of enriched fraction ###############

enrich_normalization <- function(protein.data, enrich.data, enrich.batch = c(TRUE, FALSE), probability, quantification,
                                 org, sampleTable){

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

  protein.data <- editData(data = protein.data, Fraction = "Proteome", org = org, quantification = quantification)
  enrich.data <- editData(data = enrich.data, Fraction = "Enriched", probability = probability, org = org, quantification = quantification)

  if(enrich.batch == TRUE){
    # Combine proteome and enrichment data using Uniprot IDs
    merge.data <- merge(enrich.data, protein.data, by = "Uniprot")

    # Obtain column names for each
    enrich.colnames <- as.character(sampleTable[grep("Enriched", sampleTable$fraction),]$label)
    prot.colnames <- as.character(sampleTable[grep("Proteome", sampleTable$fraction),]$label)

    # Subset the whole proteome and enriched proteome
    data.ProtSub <- merge.data[,c(prot.colnames)]
    data.EnrichSub <- merge.data[,c(enrich.colnames)]

    # Get the sequence of column
    col_num <- seq(from = length(unique(sampleTable$replicate)), to = (round(nrow(sampleTable)/2, digits = 0)+1), by = length(unique(sampleTable$replicate)))

    # Initialize a data frame for normalizing enriched intensity from each condition with corresponding proteome via the for loop
    k=1
    norm.acetData <- as.data.frame(matrix(nrow = nrow(merge.data), ncol = length(enrich.colnames)))

    while (k < length(col_num)){
      for(l in seq(from=1, to=nrow(sampleTable)/2, by=length(unique(sampleTable$replicate)))){
        norm.acetData[,(l:col_num[k])] <- data.EnrichSub[,(l:col_num[k])]/data.ProtSub[,(l:col_num[k])]
        k=k+1
      }
    }

    colnames(norm.acetData) <- enrich.colnames
    norm.acetData$Uniprot <- merge.data$Uniprot
    norm.acetData$Probability <- merge.data$Probability
    norm.acetData$siteScore <- merge.data$siteScore
    norm.acetData$Sequence <- merge.data$Sequence

    # Convert UNIPROT to symbols
    if(org != "sce"){
      norm.acetData$symbol <- mapIds(x = orgDB, keys =  as.character(norm.acetData$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
    }else {
      norm.acetData$symbol <- mapIds(x = orgDB, keys =  as.character(norm.acetData$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
    }
  }else if(enrich.batch == FALSE){
    # Combine proteome and enrichment data using Uniprot IDs
    merge.data <- merge(enrich.data, protein.data, by = "Uniprot")

    # Obtain column names for each
    enrich.colnames <- as.character(sampleTable[grep("Enriched", sampleTable$fraction),]$label)
    prot.colnames <- as.character(sampleTable[grep("Proteome", sampleTable$fraction),]$label)

    # Subset the whole proteome and enriched proteome
    data.ProtSub <- merge.data[,c(prot.colnames)]
    data.EnrichSub <- merge.data[,c(enrich.colnames)]

    # Get the sequence of column
    col_num_p <- seq(from = length(unique(sampleTable[sampleTable$fraction %in% "Proteome",]$replicate)), to = (round(nrow(sampleTable[sampleTable$fraction %in% "Proteome",]), digits = 0)+1), by = length(unique(sampleTable[sampleTable$fraction %in% "Proteome",]$replicate)))
    col_num_e <- seq(from = length(unique(sampleTable[sampleTable$fraction %in% "Enriched",]$replicate)), to = (round(nrow(sampleTable[sampleTable$fraction %in% "Enriched",]), digits = 0)+1), by = length(unique(sampleTable[sampleTable$fraction %in% "Enriched",]$replicate)))

    # Initialize a data frame for calculating average intensity from each condition via the for loop
    avg.intensity <- as.data.frame(matrix(nrow = nrow(merge.data), ncol = length(unique(sampleTable$condition))))
    i=1
    j=1
    while (i <= length(col_num_p)){
      avg.intensity[,i] <- as.data.frame(rowMeans(data.ProtSub[,j:col_num_p[i]]))
      i=i+1
      j=j+length(unique(sampleTable[sampleTable$fraction %in% "Proteome",]$replicate))
    }
    for (i in 1:length(col_num_p)){
      colnames(avg.intensity)[i] <- colnames(data.ProtSub[col_num_p[i]])
    }

    # Initialize a data frame for normalizing enriched intensity from each condition with corresponding proteome via the for loop
    k=1
    norm.acetData <- as.data.frame(matrix(nrow = nrow(merge.data), ncol = length(enrich.colnames)))

    while (k < length(col_num_e)){
      for(l in seq(from=1, to=round(nrow(sampleTable[sampleTable$fraction %in% "Enriched",])), by=length(unique(sampleTable[sampleTable$fraction %in% "Enriched",]$replicate)))){
        norm.acetData[,(l:col_num_e[k])] <- data.EnrichSub[,(l:col_num_e[k])]/avg.intensity[,colnames(data.ProtSub[col_num_p[k]])]
        k=k+1
      }
    }
    colnames(norm.acetData) <- enrich.colnames
    norm.acetData$Uniprot <- merge.data$Uniprot
    norm.acetData$Probability <- merge.data$Probability
    norm.acetData$siteScore <- merge.data$siteScore
    norm.acetData$Sequence <- merge.data$Sequence

    # Remove infinite and NaN values (not zero values)
    norm.acetData <- norm.acetData %>% filter_if(~is.numeric(.), all_vars(!is.infinite(.)))
    norm.acetData <- norm.acetData %>% filter_if(~is.numeric(.), all_vars(!is.nan(.)))

    # Convert UNIPROT to symbols
    if(org != "sce"){
      norm.acetData$symbol <- mapIds(x = orgDB, keys =  as.character(norm.acetData$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
    }else {
      norm.acetData$symbol <- mapIds(x = orgDB, keys =  as.character(norm.acetData$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
    }

  }else{
    stop("Accepted values for enrich.batch can be either TRUE or FALSE")
  }

  return(norm.acetData)
}
