#' @title obtain_exclusive
#'
#' @description Get exclusive proteins/sites (only used when filter.protein.type = "fraction")
#'
#' @param data Input or modified proteome file from MaxQuant (.txt)
#' @param Fraction Can be either "Proteome" or "Enriched". Indicates the type of input data used
#' @param sampleTable .xlsx file containing information about the samples. Three columns are mandatory (label, condition and replicate)
#' @param contrasts Mentions the conditions to be compared. Ex. MUTANT_vs_WILDTYPE or (MUTANT-A_vs_WILDTYPE-A)_vs_(MUTANT-B_vs_WILDTYPE-B). This how the contrasts should be provided and the values should be the same as the one given in condition column of sampleTable. (MUTANT-A_vs_WILDTYPE-A)_vs_(MUTANT-B_vs_WILDTYPE-B): This type can be used for complex data when comparing the interaction between two conditions such genotype and time
#'
#' @return List of exclsive proteins/sites from all comparisons (except difference of difference comparisons)
#'
#' @export
#'
##################### Obtain exclusive proteins and sites ########################

obtain_exclusive <- function(data, Fraction, sampleTable, contrasts){

  data.na <- data

  if(Fraction=="Proteome"){
    data.na <- make_unique(proteins = data.na, names = "Uniprot", ids = "symbol", delim = ";")
    bin_data <- data.na[,1:(ncol(data.na)-4)]
    idx <- is.na(bin_data)
    bin_data[!idx] <- 1
    bin_data[idx] <- 0
    rownames(bin_data) <- data.na$name
  }else if(Fraction=="Enriched"){
    data.na$symbol <- NULL
    data.na <- data.na %>%
      mutate_all(~replace(., . == 0, NA))
    data.na$symbol <- data$symbol
    data.na <- make_unique(proteins = data.na, names = "Uniprot", ids = "symbol", delim = ";")
    bin_data <- data.na[,1:(ncol(data.na)-7)]
    idx <- is.na(bin_data)
    bin_data[!idx] <- 1
    bin_data[idx] <- 0
    rownames(bin_data) <- data.na$name
  }else{
    stop("Accepted values are Proteome or Enriched")
  }

  keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(label, value, -rowname) %>%
    left_join(., sampleTable, by = "label") %>%
    group_by(rowname, condition) %>%
    dplyr::summarize(miss_val = n() - sum(value)) %>%
    spread(condition, miss_val)
  keep$`<NA>` <- NULL

  exclusive.list <- list()
  contrasts.sep <- list()

  new.contrasts <- contrasts[str_count(contrasts, "vs")==1]

  for(i in 1:length(new.contrasts)){
    contrasts.sep[i] = strsplit(gsub("_vs_|[()]", ",", new.contrasts[i]), ",")
    contrasts.sep[i] = lapply(contrasts.sep[i], function(z){ z[!is.na(z) & z != ""]})
  }

  for (i in 1:length(new.contrasts)){ # Exclusive list of difference of difference is removed as they might not exactly give the exclusive set of proteins
    # if(length(unlist(contrasts.sep[[i]])) == 2){
    keep.sub <- keep[,colnames(keep) %in% unlist(contrasts.sep[[i]])]
    keep.sub <- as.data.frame(keep.sub)
    rownames(keep.sub) <- keep$rowname
    exclusive.list[[i]] <- keep.sub[keep.sub[,1] - keep.sub[,2] <= -2 | keep.sub[,1] - keep.sub[,2] >= 2,]
    # } else if (length(unlist(contrasts.sep[[i]])) == 4){
    #   keep.sub <- keep[,colnames(keep) %in% unlist(contrasts.sep[[i]])]
    #   keep.sub <- as.data.frame(keep.sub)
    #   rownames(keep.sub) <- keep$rowname
    #   keep.sub <- keep.sub[(keep.sub[,1] - keep.sub[,2] <=-2 | keep.sub[,1] - keep.sub[,2] >= 2) & (keep.sub[,3] - keep.sub[,4] <=-2 | keep.sub[,3] - keep.sub[,4] >= 2),]
    #
    #   keep.sub$diff1 <- sapply(keep.sub[,2] - keep.sub[,1], function(a){
    #     if(a < 0){-1
    #     } else{1}})
    #   keep.sub$diff2 <- sapply(keep.sub[,4] - keep.sub[,3], function(a){
    #     if(a < 0){-1
    #     } else{1}})
    #
    #   keep.sub <- keep.sub[!(keep.sub$diff2 - keep.sub$diff1)==0,]
    #   keep.sub$diff1 <- NULL
    #   keep.sub$diff2 <- NULL
    #   exclusive.list[[i]] <- keep.sub
    # } else{
    #   print("Error in contrasts")
    # }
  }

  if(Fraction=="Proteome"){
    exclusive.list <- lapply(exclusive.list, dplyr::add_rownames, "Uniprot")
    for(i in 1:length(new.contrasts)){
      exclusive.list[[i]] <- merge(x = exclusive.list[[i]], y = data.na[ , c("name", "Uniprot", "symbol")], by = "Uniprot")
      colnames(exclusive.list[[i]]) <- gsub("Symbol", "symbol", colnames(exclusive.list[[i]]))
      exclusive.list[[i]] <- select(exclusive.list[[i]], unlist(c("Uniprot", contrasts.sep[i], "name", "symbol")))
      stopifnot(identical(unlist(c("Uniprot", contrasts.sep[i], "name", "symbol")), colnames(exclusive.list[[i]])))
    }
  }else if(Fraction=="Enriched"){
    exclusive.list <- lapply(exclusive.list, dplyr::add_rownames, "name")
    for(i in 1:length(new.contrasts)){
      exclusive.list[[i]] <- merge(x = exclusive.list[[i]], y = data.na[ , c("name", "Uniprot", "symbol", "Sequence")], by = "name")
      tmp <- exclusive.list[[i]]$Sequence
      exclusive.list[[i]]$Sequence <- tmp[,1]
      exclusive.list[[i]] <- select(exclusive.list[[i]], unlist(c("name", contrasts.sep[i], "Uniprot", "symbol", "Sequence")))
      stopifnot(identical(unlist(c("name", contrasts.sep[i], "Uniprot", "symbol", "Sequence")), colnames(exclusive.list[[i]])))
    }
  } else{
    stop("Accepted values are Proteome or Enriched")
  }
  names(exclusive.list) <- new.contrasts
  dir.create(paste(getwd(),"/Results/",Fraction,"/Exclusive_files",sep = ""), showWarnings = TRUE)
  writexl::write_xlsx(x = exclusive.list, path = paste(getwd(),"/Results/",Fraction,"/Exclusive_files/",Fraction,"_exclusive.xlsx",sep = ""), col_names = TRUE, format_headers = TRUE)
  return(exclusive.list)
}
