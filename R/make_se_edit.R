#' @title make_se_edit
#'
#' @description Generate summarized experiment object (Obtained and edited from DEP package)
#'
#' @param proteins_unique data containing the list of proteins with unique protein names (output from make_unique function)
#' @param columns LFQ or iBAQ
#' @param expdesign Sample table containing at least label, condition and replicate column
#'
#' @return An object of SummarizedExperiment class
#'
#' @export
#'
#' @importFrom tibble rownames_to_column is_tibble
#' @import dplyr
# @importFrom dplyr semi_join select add_rownames mutate %>% left_join summarize group_by filter_if
#' @importFrom tidyr gather spread unite
#' @importFrom stringr str_count str_extract_all str_match_all
##########################################################
# Make summarized experiment (modified from DEP package)

make_se_edit <- function (proteins_unique, columns, expdesign)
{
  assertthat::assert_that(is.data.frame(proteins_unique),
                          is.integer(columns), is.data.frame(expdesign))
  if (any(!c("name", "ID") %in% colnames(proteins_unique))) {
    stop("'name' and/or 'ID' columns are not present in '",
         deparse(substitute(proteins_unique)), "'.\nRun make_unique() to obtain the required columns",
         call. = FALSE)
  }
  if (any(!c("label", "condition", "replicate") %in% colnames(expdesign))) {
    stop("'label', 'condition' and/or 'replicate' columns",
         "are not present in the experimental design", call. = FALSE)
  }
  if (any(!apply(proteins_unique[, columns], 2, is.numeric))) {
    stop("specified 'columns' should be numeric", "\nRun make_se_parse() with the appropriate columns as argument",
         call. = FALSE)
  }
  if (tibble::is_tibble(proteins_unique))
    proteins_unique <- as.data.frame(proteins_unique)
  if (tibble::is_tibble(expdesign))
    expdesign <- as.data.frame(expdesign)
  rownames(proteins_unique) <- proteins_unique$name
  raw <- proteins_unique[, columns]
  raw[raw == 0] <- NA
  raw <- log2(raw)
  expdesign <- mutate(expdesign, condition = make.names(condition)) %>%
    unite(ID, condition, replicate, remove = FALSE)
  rownames(expdesign) <- expdesign$ID
  matched <- match(make.names(expdesign$label),
                   make.names(colnames(raw)))
  if (any(is.na(matched))) {
    stop("None of the labels in the experimental design match ",
         "with column names in 'proteins_unique'", "\nRun make_se() with the correct labels in the experimental design",
         "and/or correct columns specification")
  }
  colnames(raw)[matched] <- expdesign$ID
  raw <- raw[, !is.na(colnames(raw))][rownames(expdesign)]
  row_data <- proteins_unique[, -columns]
  rownames(row_data) <- row_data$name
  se <- SummarizedExperiment(assays = as.matrix(raw), colData = expdesign,
                             rowData = row_data)
  return(se)
}
