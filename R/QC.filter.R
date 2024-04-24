#' @title QC.filter
#'
#' @description Performs quality control and filters the missing values based on parameters provided
#'
#' @param data output file (or output) from editData function
#' @param filter.thr Only if filter.protein.type = condition. Numerical value less than the number of relicates in either condition (Ex. 0 indicates the protein should have no NAs in all replicates of atleast one condition while 1 indicates they can have one NAs)
#' @param filter.protein.min Only if filter.protein.type = fraction. Any value between 0-1 Any value between 0-1 (Ex. 0.75 indicates the protein should not have NAs in 75\% of all samples)
#'
#' @return Filtered and vsn normalized data
#'
#' @examples
#' norm_data <- QC.filter(data = lfqdata, filter.protein.type = "fraction", filter.protein.min = 0.75)
#' @export
###################################################### Filtering and quality control of the data ##############################################

QC.filter <- function(data, filter.thr = NA, filter.protein.min = NULL){

  if(org != "sce"){
    data$symbol <- mapIds(x = orgDB, keys =  as.character(data$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")
  }else {
    data$symbol <- mapIds(x = orgDB, keys =  as.character(data$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")
  }

  # Are there any duplicated gene names?
  data$symbol %>% duplicated() %>% any()

  # Make a table of duplicated gene names
  data %>% dplyr::group_by(symbol) %>% dplyr::summarize(frequency = n()) %>% dplyr::arrange(desc(frequency)) %>% filter(frequency > 1)

  # Make unique names using the annotation in the "Symbol" column as primary names and the annotation in "Uniprot" as name for those that do not have an gene name.
  data_unique <- make_unique(proteins = data, names = "Uniprot", ids = "symbol", delim = ";")

  # Are there any duplicated names?
  data_unique$name %>% duplicated() %>% any()

  if(Fraction == "Proteome"){
    LFQ_columns <- grep(quantification, colnames(data_unique))
    experimental_design <- sampleTable[grep(Fraction, sampleTable$fraction),]
    } else if(Fraction == "Enriched"){
    colnames(data_unique) <- gsub("Normalized.", "", colnames(data_unique))
    LFQ_columns <- grep("Intensity", colnames(data_unique))
    experimental_design <- sampleTable[grep(Fraction, sampleTable$fraction),]
  } else{
    stop("Mention one of these: Proteome or Enriched")
  }

  # Converting factor to character
  experimental_design$label <- as.character(experimental_design$label)
  experimental_design$condition <- as.character(experimental_design$condition)
  data_se <- make_se_edit(proteins_unique = data_unique, columns = LFQ_columns, expdesign = experimental_design)

  if(filter.protein.type == "complete"){
    data_filt <- filter_proteins(data_se, type = filter.protein.type)
    data_norm <- normalize_vsn(data_filt)
  }else if(filter.protein.type == "condition"){
    data_filt <- filter_proteins(data_se, type = filter.protein.type, thr = filter.thr)
    data_norm <- normalize_vsn(data_filt)
  }else if(filter.protein.type == "fraction"){
    data_filt <- filter_proteins(data_se, type = filter.protein.type, min = filter.protein.min)
    data_norm <- normalize_vsn(data_filt)
  } else{
    stop("Mention one of these: complete, condition or fraction")
  }

  dir.create(paste(path1,"/",Fraction,"/Initial_filtering",sep = ""), showWarnings = FALSE)
  pdf(file = paste(path1,"/",Fraction,"/Initial_filtering/",Fraction,"_InitialFiltering-plots.pdf",sep = ""))

  print(plot_frequency(data_se))
  print(plot_numbers(data_filt))
  print(plot_coverage(data_filt))
  print(plot_normalization(data_filt, data_norm))
  meanSdPlot(data_norm)
  if(filter.protein.type!="complete"){
    plot_missval(data_norm) # Plot a heatmap of proteins with missing values
    plot_detect(data_norm)
  }
  dev.off()
  return(data_norm)
}
