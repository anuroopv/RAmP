#' @title filter.identical
#'
#' @description Remove identical proteins/sites from exclusive and non-exclusive data (only used when filter.protein.type = "fraction")
#'
#' @param nonexclusive.data Data frame(s) (as list) obtained after limma differential expression analysis
#' @param exlusive.data Output from obtain_exclusive
#'
#' @return Data frame(s) (as list) for volcano plots
#'
#' @export
################################################################################
# Remove identical values from non-exclusive data that are also present in exclusive data

filter.identical <- function(nonexclusive.data, exlusive.data){
  if(Fraction == "Enriched"){
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
