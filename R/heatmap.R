#' @title heatmap
#'
#' @description Generates heatmap from pheatmap
#'
#' @param data_impute Output from plot_pca function
#' @param lfq.data Title of the plot
#' @param distance.matrix One of "spearman", "pearson", "uncentered correlation", "absolute pearson", "sqrt", "weird"
#' @param exclusive.data Data frame(s) (as list) containing exclusive proteins or sites
#' @param clustering.method One of "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
#' @param title Title of the plot
#'
#' @return Generates heatmap and a file containing the order of protein(s) or site(s) in the heatmap
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom seriation seriate get_order
#' @importFrom ClassDiscovery distanceMatrix
#' @importFrom grid convertX convertY gpar
################################################################################

# Function for generating heatmaps for exclusive and significant proteins/sites

heatmap <- function(data_impute, lfq.data = lfq.data,
                    distance.matrix = c("spearman", "pearson", "uncentered correlation", "absolute pearson", "sqrt", "weird"), exclusive.data = exclusive.data,
                    clustering.method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
                    title = NA){

  sign.data <- list()
  protlist <- list()

  for(i in 1:length(contrasts)){
    conditions <- unlist(strsplit(gsub("_vs_|[()]", ",", contrasts[i]), ","))
    conditions <- conditions[nzchar(conditions)]
    if(length(conditions) == 2){conditions = conditions}else{conditions = conditions[c(1,3)]}

    if(filter.protein.type == "fraction" & str_count(contrasts[i], "vs") == 1){
      title <- "Heat map of exclusive and significant features"
    }else if(filter.protein.type == "condition" & str_count(contrasts[i], "vs") == 1){
      title <- "Heat map of significant features"
    }else if(filter.protein.type == "complete" & str_count(contrasts[i], "vs") == 1){
      title <- "Heat map of significant features"
    }else if(str_count(contrasts[i], "vs") != 1){
      title <- "Heat map of significant features \n(Note: Exclusive data is not available for difference of difference))"
      }

    data.heatmap <- as.data.frame(rowData(data_impute))
    colnames(data.heatmap) <- gsub("X.", "", colnames(data.heatmap))
    colnames(data.heatmap) <- gsub("\\._vs_\\.", "_vs_", colnames(data.heatmap))
    colnames(data.heatmap) <- gsub("\\.\\.", ".", colnames(data.heatmap))
    sign.data[[i]] <- data.heatmap[,grep(contrasts[i], colnames(data.heatmap))]
    sign.data[[i]]$Uniprot <- data.heatmap$Uniprot
    sign.data[[i]]$symbol <- data.heatmap$symbol
    if(Fraction == "Enriched"){sign.data[[i]]$Sequence <- data.heatmap[,5]}
    sign.data[[i]]$UniqueNames <- data.heatmap$name
    sign.data[[i]] <- subset(sign.data[[i]], abs(sign.data[[i]][,1]) > lfcCutOff & (sign.data[[i]][,9] < pvalCutOff | sign.data[[i]][,11] < sigmaCutOff))

    if(Fraction == "Enriched" & filter.protein.type=="fraction" & str_count(contrasts[i], "vs") == 1){rownames(exclusive.data[[i]]) <- exclusive.data[[i]]$name}

    protlist[[i]] <- if(Fraction == "Enriched"){
      if(filter.protein.type=="fraction" & str_count(contrasts[i], "vs") == 1){
        distinct(rbind(sign.data[[i]][,c(12,14)], exclusive.data[[i]][,c(4,6)]))
      }else{sign.data[[i]][,c(12,14)]
      }
    }else{
      if(filter.protein.type=="fraction" & str_count(contrasts[i], "vs") == 1){
        unique(c(sign.data[[i]]$Uniprot, exclusive.data[[i]]$Uniprot))
      }else{
        unique(sign.data[[i]]$Uniprot)
      }
    }

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

    sampleTable.sub <- sampleTable[sampleTable$condition %in% conditions,]

    if(Fraction == "Proteome"){
      protIntensityData <- lfq.data[,colnames(lfq.data) %in% sampleTable.sub$label]
      protIntensityData$symbol <- lfq.data$symbol
      protIntensityData$Uniprot <- lfq.data$Uniprot
      protIntensityData <- protIntensityData[protIntensityData$Uniprot %in% protlist[[i]],]

      protIntensityData2 <- protIntensityData[,c(1:(ncol(protIntensityData)-2))]
      rownames(protIntensityData2) <- protIntensityData$Uniprot
      protIntensityData2 <- log2(protIntensityData2)

    }else if(Fraction == "Enriched"){
      protIntensityData <- lfq.data[,colnames(lfq.data) %in% sampleTable.sub$label]
      protIntensityData$symbol <- lfq.data$symbol
      protIntensityData$Uniprot <- lfq.data$Uniprot

      protIntensityData$Sequence <- lfq.data$Sequence[[1]]
      protIntensityData <- protIntensityData[protIntensityData$Sequence %in% protlist[[i]]$Sequence,]

      protIntensityData2 <- protIntensityData[,c(1:(ncol(protIntensityData)-3))]
      rownames(protIntensityData2) <- make.unique(paste(protIntensityData$Uniprot, protIntensityData$Sequence, sep="_"))
      protIntensityData2 <- log2(protIntensityData2)

    }else{
      stop("Check the value provided for Fraction")
    }

    hc.features <- as.matrix(protIntensityData2) %>% t() %>%
      distanceMatrix(metric=distance.matrix) %>%
      hclust(method=clustering.method)

    hc.samples <- as.matrix(protIntensityData2) %>%
      distanceMatrix(metric=distance.matrix) %>%
      hclust(method=clustering.method)

    names(protIntensityData2) <- sampleTable[grep(Fraction, sampleTable$fraction),]$condition[match(names(protIntensityData2), sampleTable[grep(Fraction, sampleTable$fraction),]$label)]

    calc_ht_size = function(ht, unit = "inch") {
      pdf(NULL)
      ht = ComplexHeatmap::draw(ht)
      w = ComplexHeatmap:::width(ht)
      w = convertX(w, unit, valueOnly = TRUE)
      h = ComplexHeatmap:::height(ht)
      h = convertY(h, unit, valueOnly = TRUE)
      dev.off()

      c(w, h)
    }

    lgd = ComplexHeatmap::Legend(at = 1, labels = "NA", legend_gp = gpar(fill = "grey"))
    x <- ComplexHeatmap::Heatmap(matrix = t(scale(t(as.matrix(protIntensityData2)))), border = "grey", na_col = "darkgrey", column_title = title, cluster_rows = hc.features,
                                 cluster_columns = hc.samples, width = ncol(protIntensityData2)*unit(5, "mm"), height = nrow(protIntensityData2)*unit(0.2, "mm"),
                                 show_row_names = FALSE, show_column_names = TRUE, heatmap_legend_param = list(title = "Scaled values"))

    size = calc_ht_size(ht = x, unit = "inch")

    x <- ComplexHeatmap::draw(x)

    dir.create(paste(path1,"/",Fraction,"/Heatmap",sep = ""), showWarnings = FALSE)
    pdf(file = paste(path1,"/",Fraction,"/Heatmap/",Fraction,"_",contrasts[i],"_heatMaps.pdf",sep = ""), height = size[2]*1.25, width = size[1]*1.5)
    ComplexHeatmap::draw(x, heatmap_legend_list = lgd)
    dev.off()

    geneOrder <- data.frame("Uniprot" = rownames(protIntensityData2)[ComplexHeatmap::row_order(x)])
    if(Fraction=="Enriched"){geneOrder$Sequence = gsub(".*_", "", geneOrder$Uniprot)
    geneOrder$Uniprot = gsub("_.*", "", geneOrder$Uniprot)}

    geneOrder$symbol <- if(org != "sce"){mapIds(x = orgDB, keys =  as.character(geneOrder$Uniprot), column = "SYMBOL", keytype="UNIPROT", multiVals="first")}else{mapIds(x = orgDB, keys =  as.character(geneOrder$Uniprot), column = "GENENAME", keytype="UNIPROT", multiVals="first")}
    write.csv(geneOrder, paste(path1,"/",Fraction,"/Heatmap/",Fraction,"_",contrasts[i],"_proteinOrder-heatMaps.csv",sep = ""), row.names = FALSE)
  }

}
