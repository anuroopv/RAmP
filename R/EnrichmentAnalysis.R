#' @title EnrichmentAnalysis
#'
#' @description Performs Gene Set Enrichment Analysis and Over Representation Analaysis using clusterProfiler package
#'
#' @param enrich Can be either "gsea" or "enrich"
#' @param nonExclusive.data Default is NULL. If not, the data should be given as list containing data frame(s)
#' @param rankBy Can be either stat or logFC. Genes are ranked accordingly for GSEA
#' @param KEGG Logical, default is FALSE
#' @param ont Can be "BP", "CC", "MF" or "ALL"
#' @param padjustMethod Method for adjusting p-value (Refer clusterProfiler for available options)
#' @param exclusive.data Default is NULL. If not, the data should be given as list containing data frame(s)
#' @param background Logical, default is TRUE. If FALSE, background proteins will not be considered for over representation analysis
#' @param minGS minimal size of annotated genes in geneset
#' @param maxGS maximal size of annotated genes in geneset
#' @param simplify Logical, default is FALSE. If TRUE, redundant GO terms will be simplified based on the simplify_cutoff
#' @param simplify_cutoff Cutoff for simplify method. Numerical value from 0-1. Default is 0.7
#'
#' @return Generates list of significant GO terms
#'
#' @export
#'
#' @importFrom DOSE setReadable
#' @importFrom clusterProfiler gseGO gseKEGG simplify enrichGO
################################################################################
# GO term analysis

##### Gene Set Enrichement analysis ####

EnrichmentAnalysis <- function(enrich = c('gsea', 'ora'), nonExclusive.data = NULL,
                        rankBy = "stat", KEGG = FALSE,
                        ont = "BP",
                        padjustMethod = "fdr",
                        exclusive.data = NULL,
                        background = TRUE,
                        minGS = 50, maxGS = 500,
                        simplify = FALSE, simplify_cutoff = 0.7){

  if(enrich == "gsea"){
  if(Fraction == "Enriched"){
    stop(print("GSEA cannot be perfromed for enriched data as they cannot be ranked due to repetitiveness of symbol. Proceed with over-representation analysis"))
  }else if(is.null(exclusive.data) == FALSE){
    print("GSEA cannot be perfromed for exclusive data but the analysis will be proceeded for non-exclusive data. To use both, perform over-representation analysis")
  }else{
    print("Performing Gene Set Enrichment Analysis")
  }
  gse.list <- list()
  pathways.list <- list()

  for (i in 1:length(contrasts)){
    data.sub <- nonExclusive.data[[i]]
    rownames(data.sub) <- data.sub$name

    if(rankBy == "stat"){
      allGenes = data.sub[,grep("[.]t$", colnames(data.sub))]
    } else if(rankBy == "logFC"){
      allGenes = data.sub[,grep("logFC", colnames(data.sub))]
    } else{
      stop("Accepted values are stat or logFC")
    }

    names(allGenes) = data.sub$Uniprot
    allGenes <- na.omit(allGenes)
    allGenes = sort(allGenes, decreasing = TRUE)
    allGenes = allGenes[!duplicated(names(allGenes))]

    if(KEGG == FALSE){
      gse <- gseGO(geneList=allGenes,
                   ont = ont,
                   keyType = "UNIPROT",
                   minGSSize = 50,
                   maxGSSize = 600,
                   pvalueCutoff = pvalCutOff,
                   verbose = TRUE,
                   OrgDb = orgDB,
                   pAdjustMethod = padjustMethod, by = 'fgsea', eps = 0)

      if(org!="sce"){gse <- setReadable(x = gse, OrgDb = orgDB, keyType = "UNIPROT")
      }else{
        gse <- gse
        print("Uniprot IDs for yeast cannot be translated to corresponding gene names")}

    } else {
      gse <- gseKEGG(geneList=allGenes,
                     organism = org,
                     keyType = 'uniprot',
                     minGSSize = 50,
                     maxGSSize = 600,
                     pvalueCutoff = pvalCutOff, pAdjustMethod = padjustMethod,
                     verbose = TRUE, eps = 0)

      if(org!="sce"){gse <- setReadable(x = gse, OrgDb = orgDB, keyType = "UNIPROT")
      }else{
        gse <- gse
        print("Uniprot IDs for yeast cannot be translated to corresponding gene names")}

    }

    if(simplify == TRUE & KEGG == FALSE){
      gse <- clusterProfiler::simplify(gse, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min, measure = "Wang", semData = NULL)
    } else if(simplify == TRUE & KEGG == TRUE){
      stop("Simplification cannot be performed for GO terms from KEGG pathways")
    } else {
      gse <- gse
    }

    pathways <- as.data.frame(gse)
    if(nrow(pathways) == 0){
      print("No significant go-terms for the provided cutoff parameters")
      # next
    }
    gse.list[[i]] <- gse
    pathways.list[[i]] <- pathways

    names(gse.list)[[i]] <- contrasts[i]
    names(pathways.list)[[i]] <- contrasts[i]
  }
  dir.create(paste(path1,"/",Fraction,"/Enrichment_analysis",sep = ""), showWarnings = FALSE)
  writexl::write_xlsx(x = pathways.list, path = paste(path1,"/",Fraction,"/Enrichment_analysis/",Fraction,"_GSEA.xlsx",sep = ""), col_names = TRUE, format_headers = TRUE)
  return(gse.list)
}else if(enrich == "ora"){

  ##### Over representation analysis ####
  go_enrich.list <- list()
  pathways.list <- list()
  for(i in 1:length(contrasts)){

    nonExdata.sub <- nonExclusive.data[[i]]
    if(is.null(exclusive.data)==FALSE & str_count(contrasts[i], "vs") == 1){Exdata.sub <- exclusive.data[[i]]
    background.genes <- unique(c(Exdata.sub$Uniprot, nonExdata.sub$Uniprot))
    }else{Exdata.sub <-NULL
    background.genes <- unique(nonExdata.sub$Uniprot)
    }

    conditions <- unlist(strsplit(gsub("_vs_|[()]", ",", contrasts[i]), ","))
    conditions <- conditions[nzchar(conditions)]
    if(length(conditions) == 2){conditions <- conditions}else{conditions <- conditions[c(1,3)]}
    j=1

    genelist.1 <- unique(c(subset(nonExdata.sub, (nonExdata.sub[,grep("q.val", colnames(nonExdata.sub))] < pvalCutOff | nonExdata.sub[,grep("SIGMA", colnames(nonExdata.sub))] < sigmaCutOff) & nonExdata.sub[,grep("logFC", colnames(nonExdata.sub))] > lfcCutOff)[,2],
                           subset(Exdata.sub, Exdata.sub[,grep(conditions[j], colnames(Exdata.sub))] <= (length(unique(sampleTable$replicate))-2) & Exdata.sub[,grep(conditions[j+1], colnames(Exdata.sub))] >= (length(unique(sampleTable$replicate))-1))[,1]))
    genelist.1 <- gsub("[.].*","", genelist.1)

    genelist.2 <- unique(c(subset(nonExdata.sub, (nonExdata.sub[,grep("q.val", colnames(nonExdata.sub))] < pvalCutOff | nonExdata.sub[,grep("SIGMA", colnames(nonExdata.sub))] < sigmaCutOff) & nonExdata.sub[,grep("logFC", colnames(nonExdata.sub))] < lfcCutOff)[,2],
                           subset(Exdata.sub, Exdata.sub[,grep(conditions[j+1], colnames(Exdata.sub))] <= (length(unique(sampleTable$replicate))-2) & Exdata.sub[,grep(conditions[j], colnames(Exdata.sub))] >= (length(unique(sampleTable$replicate))-1))[,1]))
    genelist.2 <- gsub("[.].*","", genelist.2)

    k = i-1
    if(background == TRUE){
      go_enrich <- enrichGO(gene = genelist.1,
                            universe = background.genes,
                            OrgDb = orgDB,
                            keyType = "UNIPROT",
                            readable = if(org!="sce"){T}else{F},
                            ont = ont,
                            pvalueCutoff = pvalCutOff,
                            pAdjustMethod = padjustMethod, minGSSize = minGS, maxGSSize = maxGS)

      if(org!="sce"){go_enrich <-setReadable(x = go_enrich, OrgDb = orgDB, keyType = "UNIPROT")
      }else{
        go_enrich <- go_enrich
        print("Uniprot IDs for yeast cannot be translated to corresponding symbols")}

      go_enrich.list[[i+k]] <- go_enrich
      pathways.list[[i+k]] <- as.data.frame(go_enrich)
      k = k+1
      go_enrich <- enrichGO(gene = genelist.2,
                            universe = background.genes,
                            OrgDb = orgDB,
                            keyType = "UNIPROT",
                            readable = if(org!="sce"){T}else{F},
                            ont = ont,
                            pvalueCutoff = pvalCutOff,
                            pAdjustMethod = padjustMethod, minGSSize = minGS, maxGSSize = maxGS)

      if(org!="sce"){go_enrich <-setReadable(x = go_enrich, OrgDb = orgDB, keyType = "UNIPROT")
      }else{
        go_enrich <- go_enrich
        print("Uniprot IDs for yeast cannot be translated to corresponding symbols")}

      go_enrich.list[[i+k]] <- go_enrich
      pathways.list[[i+k]] <- as.data.frame(go_enrich)

    }else{
      go_enrich <- enrichGO(gene = genelist.1,
                            OrgDb = orgDB,
                            keyType = "UNIPROT",
                            readable = if(org!="sce"){T}else{F},
                            ont = ont,
                            pvalueCutoff = pvalCutOff,
                            pAdjustMethod = padjustMethod, minGSSize = minGS, maxGSSize = maxGS)

      if(org!="sce"){go_enrich <-setReadable(x = go_enrich, OrgDb = orgDB, keyType = "UNIPROT")
      }else{
        go_enrich <- go_enrich
        print("Uniprot IDs for yeast cannot be translated to corresponding symbols")}

      go_enrich.list[[i+k]] <- go_enrich
      pathways.list[[i+k]] <- as.data.frame(go_enrich)
      k=k+1
      go_enrich <- enrichGO(gene = genelist.2,
                            OrgDb = orgDB,
                            keyType = "UNIPROT",
                            readable = if(org!="sce"){T}else{F},
                            ont = ont,
                            pvalueCutoff = pvalCutOff,
                            pAdjustMethod = padjustMethod, minGSSize = minGS, maxGSSize = maxGS)

      if(org!="sce"){go_enrich <-setReadable(x = go_enrich, OrgDb = orgDB, keyType = "UNIPROT")
      }else{
        go_enrich <- go_enrich
        print("Uniprot IDs for yeast cannot be translated to corresponding symbols")}

      go_enrich.list[[i+k]] <- go_enrich
      pathways.list[[i+k]] <- as.data.frame(go_enrich)
    }
    k = i-1

    names(pathways.list)[c(i+k,i+k+1)] <- paste(conditions, contrasts[i], sep = "_")
    names(go_enrich.list)[c(i+k,i+k+1)] <- paste(conditions, contrasts[i], sep = "_")
  }
  dir.create(paste(path1,"/",Fraction,"/Enrichment_analysis",sep = ""), showWarnings = FALSE)
  writexl::write_xlsx(x = pathways.list, path = paste(path1,"/",Fraction,"/Enrichment_analysis/",Fraction,"_ORA.xlsx", sep = ""), col_names = TRUE, format_headers = TRUE)
  return(go_enrich.list)
}else{
  stop("Accepted values for enrich are: gsea or ora")
}
}
