#' @title GSEAPlots
#'
#' @description Visualization of GSEA or ORA
#'
#' @param gseData Output from clusterProfiler (as list)
#' @param fraction Can be either "Proteome" or "Enriched". Indicates the type of input data used
#' @param plotType Can be one of "dotPlot", "cNetPlot", "heatPlot", "treePlot", "gseaPlot", "ridgePlot"
#' @param circular Logical, default is FALSE. If TRUE, circular visualization will be used for cNetPlot
#' @param colorEdge Logical, default is FALSE. Parameter to be used if circular = "TRUE" and plotType ' "cNetPlot
#' @param nodeLabel Can be one of "gene", "category", "all", "none". Parameter to be used if circular = "TRUE" and plotType ' "cNetPlot
#' @param cexLabelCategory Parameter to be used if circular = "TRUE" and plotType ' "cNetPlot (Check enrichplot package for further info). Defaults i 1.2
#' @param cexLabelGene Parameter to be used if circular = "TRUE" and plotType ' "cNetPlot (Check enrichplot package for further info). Default is 0.8
#' @param colorCcategory Parameter to be used if circular = "TRUE" and plotType ' "cNetPlot (Check enrichplot package for further info). Default is "black"
#' @param colorGene Parameter to be used if circular = "TRUE" and plotType ' "cNetPlot (Check enrichplot package for further info). Default is "black"
#' @param nCluster Number of clusters to be generated. Parameter used when plotType = "treePlot". Defualt is 5
#' @param showCategory Number of GO terms to be displayed in the plot. Default is 10
#' @param pvalCutOff P-value cut off for significant protein(s), site(s) and GO terms. Default is 0.05
#' @param org Database of the organism. Drosophila melanogaster = "dme", Mus muscuslus ' "mmu", Homo sapiens = "hsa", Saccharomyces cerevisae = "sce". Default is "dme"
#'
#' @return Generates plots from GO term analysis
#'
#' @export
#' @importFrom dplyr "%>%"
#'
################################################################################

# Various types of plot to represent the GO terms

GSEAPlots <- function(gseData, fraction = fraction, enrich = c('gsea', 'ora'),
                      plotType = c("dotPlot", "cNetPlot", "heatPlot", "treePlot", "gseaPlot", "ridgePlot"),
                      circular = FALSE, colorEdge = FALSE, nodeLabel = c("gene", "category", "all", "none"), cexLabelCategory = 1.2, cexLabelGene = 0.8, colorCcategory = "black", colorGene = "black",
                      nCluster = 5, showCategory = 10,
                      org = "dme", pvalCutOff = 0.05){
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
  dir.create(paste(getwd(),"/Results/Enrichment_analysis",sep = ""),showWarnings = FALSE)

  for(i in 1:length(gseData)){

    if(nrow(subset(gseData[[i]]@result, p.adjust < pvalCutOff)) == 0){next}

    # if(colnames(gseData[[i]]@result)[5] != "NES" & is.null(geneList) != TRUE){
    #   geneList <- geneList
    # }else if(colnames(gseData[[i]]@result)[5] != "NES" & is.null(geneList) == TRUE){
    #   geneList <- NULL
    # }else if(colnames(gseData[[i]]@result)[5] == "NES" & is.null(geneList) == TRUE){
    #   geneList <- gseData[[i]]@geneList
    # }

    if(plotType == "dotPlot"){
      pdf(file = paste(getwd(),"/Results/Enrichment_analysis/",fraction,"_",names(gseData[i]),"_enrichmentPlots-dotPlot.pdf",sep = ""), paper = "a4r")
      if(enrich == "gsea"){print(enrichplot::dotplot(gseData[[i]], showCategory=showCategory) + ggtitle(names(gseData[i])) +
                                   facet_grid(.~.sign))
      }else{print(enrichplot::dotplot(gseData[[i]], showCategory=showCategory) +
                    ggtitle(names(gseData[i])))
      }
      dev.off()
    }else if(plotType == "cNetPlot"){
      pdf(file = paste(getwd(),"/Results/Enrichment_analysis/",fraction,"_",names(gseData[i]),"_enrichmentPlots-cNetPlot.pdf",sep = ""), width = 5, height = 4, paper = "a4r")
      print(enrichplot::cnetplot(x = gseData[[i]], categorySize = "pvalue", foldChange=geneList, showCategory=showCategory,
                                 circular = circular, colorEdge = colorEdge, node_label = nodeLabel,
                                 cex_label_category = cexLabelCategory, cex_label_gene = cexLabelGene)
            + ggtitle(names(gseData[i])))
      dev.off()
    }else if(plotType == "heatPlot"){
      pdf(file = paste(getwd(),"/Results/Enrichment_analysis/",fraction,"_",names(gseData[i]),"_enrichmentPlots-heatPlot.pdf",sep = ""), width = 5, height = 4, paper = "a4r")
      print(enrichplot::heatplot(x = gseData[[i]], foldChange = geneList, showCategory = showCategory) +
              ggtitle(names(gseData[i])))
      dev.off()
    }else if(plotType == "treePlot"){
      pdf(file = paste(getwd(),"/Results/Enrichment_analysis/",fraction,"_",names(gseData[i]),"_enrichmentPlots-treePlot.pdf",sep = ""), paper = "a4r")
      edox <- pairwise_termsim(gseData[[i]])
      print(enrichplot::treeplot(edox, showCategory = showCategory, cluster.params = list(n = nCluster)) +
              ggtitle(names(gseData[i])))
      dev.off()
    }else if(plotType == "gseaPlot"){
      if(enrich == "ora"){
        print("GSEA plot cannot be used for Over-representation Analysis")
        stop()
      }
      pdf(file = paste(getwd(),"/Results/Enrichment_analysis/",fraction,"_",names(gseData[i]),"_enrichmentPlots-gseaPlot.pdf",sep = ""), paper = "a4r")
      if(is.numeric(showCategory) == TRUE){
        print(enrichplot::gseaplot2(gseData[[i]], geneSetID = 1:showCategory, pvalue_table = TRUE, title = names(gseData[i]), subplots = 1))
      }else{
        print(enrichplot::gseaplot2(gseData[[i]], geneSetID = showCategory, pvalue_table = TRUE, title = names(gseData[i]), subplots = 1))
      }
      dev.off()
    } else if(plotType == "ridgePlot"){
      if(enrich == "ora"){
        print("Ridge plot cannot be used for Over-representation Analysis")
        stop()
      }
      pdf(file = paste(getwd(),"/Results/Enrichment_analysis/",fraction,"_",names(gseData[i]),"_enrichmentPlots-ridgePlot.pdf",sep = ""), paper = "a4r")
      print(enrichplot::ridgeplot(x = gseData[[i]], showCategory = showCategory, fill = "p.adjust", decreasing = TRUE) +
              ggtitle(names(gseData[i])))
      dev.off()
    }
  }
}
