#' @title filter.identical
#'
#' @description Remove identical proteins/sites from exclusive and non-exclusive data (only used when filter.protein.type = "fraction")
#'
#' @param nonexclusive.data Data frame(s) (as list) obtained after limma differential expression analysis
#' @param exlusive.data Output from obtain_exclusive
#' @param fraction Can be either "Proteome" or "Enriched". Indicates the type of input data used
#' @param contrasts Mentions the conditions to be compared. Ex. MUTANT_vs_WILDTYPE or (MUTANT-A_vs_WILDTYPE-A)_vs_(MUTANT-B_vs_WILDTYPE-B). This how the contrasts should be provided and the values should be the same as the one given in condition column of sampleTable. (MUTANT-A_vs_WILDTYPE-A)_vs_(MUTANT-B_vs_WILDTYPE-B): This type can be used for complex data when comparing the interaction between two conditions such genotype and time
#'
#' @return Data frame(s) (as list) for volcano plots
#'
#' @export
################################################################################
# Remove identical values from non-exclusive data that are also present in exclusive data

filter.identical <- function(nonexclusive.data, exlusive.data, fraction, contrasts){
  if(fraction == "Enriched"){
    for (i in 1:length(contrasts)){
      if(str_count(contrasts[i], "vs")==1){
        test.ex <- exlusive.data[[i]]
        test.NonEx <- nonexclusive.data[[i]]
        colnames(test.NonEx)[5] <- "Sequence"

        duplicate <- semi_join(test.NonEx, test.ex, by = c("name", "Sequence"))
        duplicate <- subset(duplicate, duplicate$imputed == TRUE)

        nonexclusive.data[[i]] <- nonexclusive.data[[i]][!nonexclusive.data[[i]]$name %in% duplicate$name,]
      }
    }
  }else{
    for (i in 1:length(contrasts)){
      if(str_count(contrasts[i], "vs")==1){
        test.ex <- exlusive.data[[i]]
        test.NonEx <- nonexclusive.data[[i]]

        duplicate <- semi_join(test.NonEx, test.ex, by = "name")
        duplicate <- subset(duplicate, duplicate$imputed == TRUE)

        nonexclusive.data[[i]] <- nonexclusive.data[[i]][!nonexclusive.data[[i]]$name %in% duplicate$name,]
      }
    }
  }
  return(nonexclusive.data)
}
