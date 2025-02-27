#' @title volcanoPlot
#'
#' @description Generates heatmap from pheatmap
#'
#' @param proteinList Vector with favourite proteins (or their corresponding sites) to be named in the volcano plot irrespective of its statistical significance
#' @param name.sigProteins Logical (Default is FALSE). If TRUE, names all significant protein(s) or site(s)
#' @param resData Data frame(s) (as list) from limma differential expression analysis
#'
#' @return Generates volcanoPlot
#'
#' @export
################################################################################
# Function for generating volcano plots

volcanoPlot <- function(proteinList, name.sigProteins = FALSE, resData){

  dir.create(paste(path1,"/",Fraction,"/VolcanoPlots",sep = ""), showWarnings = FALSE)
  pdf(paste(path1,"/",Fraction,"/VolcanoPlots/",Fraction,"_volcanoPlot.pdf",sep = ""), paper = "a4r")

  for(i in seq_along(1:length(contrasts))){
    conditions <- unlist(strsplit(gsub("_vs_|[()]", ",", contrasts[i]), ","))
    data.sub <- resData[[i]]
    rownames(data.sub) <- resData[[i]]$name
    proteinData <- data.sub[data.sub$symbol %in% proteinList,]

    # To set the corresponding column number for log-fc, qvalue and non-adj pvalue
    if(Fraction == "Proteome"){
      col_num = 14
      logFC_col_num = 6
      pval_col_num = 11
    } else if(Fraction == "Proteome" & filter.protein.type == "complete"){
      col_num = 12
      logFC_col_num = 4
      pval_col_num = 9
    } else if(Fraction == "Enriched"){
      col_num = 17
      logFC_col_num = 9
      pval_col_num = 14
    } else{
      col_num = 15
      logFC_col_num = 7
      pval_col_num = 12
    }

    signList <- rownames(data.sub[abs(data.sub[,logFC_col_num])>lfcCutOff & data.sub[,col_num] < pvalCutOff,])
    sigmaList <- rownames(data.sub[abs(data.sub[,logFC_col_num])>lfcCutOff & data.sub[,ncol(data.sub)] < sigmaCutOff,])
    signData <- data.sub[rownames(data.sub) %in% signList,]
    sigmaData <- data.sub[rownames(data.sub) %in% sigmaList,]
    sigmaData <- anti_join(sigmaData, signData, "name")
    allSigData <- rbind(sigmaData, signData)

    data.sub <- anti_join(data.sub, signData, "name")
    data.sub <- anti_join(data.sub, sigmaData, "name")

    volcano <- ggplot(mapping = aes(x = data.sub[,logFC_col_num], y = -log10(data.sub[,pval_col_num]))) +
      geom_point(size=2, shape =1, color = "grey") +
      theme_pubr() + # change theme
      labs(x = expression(log[2]~(FoldChange)), y = expression(-log[10]~(P-value)), subtitle = contrasts[i]) +
      {if(lfcCutOff!=0){geom_vline(xintercept = c(-lfcCutOff,lfcCutOff), colour = "#E10600FF", lty = "dashed")}} + # Add fold change cutoffs
      geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
      geom_text_repel(data=proteinData, aes(x = proteinData[,logFC_col_num], y = -log10(proteinData[,pval_col_num]), label = symbol), size = 2,
                      color = ifelse(proteinData[,col_num] < pvalCutOff,'#00FF00','#000000'),
                      fontface = "bold", point.padding = unit(0.5, 'lines'), box.padding = unit(.25, 'lines'), max.overlaps = Inf) + # Adjusted p value based on limma adjusted P-value
      geom_point(data = proteinData, aes(x = proteinData[,logFC_col_num], y = -log10(proteinData[,pval_col_num])), shape = 1, size = 2,
                 color = ifelse(proteinData[,col_num] < pvalCutOff,'#00FF00','#000000')) +
      geom_point(data = signData, mapping = aes(x = signData[,logFC_col_num], y = -log10(signData[,pval_col_num])), size = 2, color = ifelse(signData[,logFC_col_num] <0, "#00239CFF", "#E10600FF")) + # Adjusted p value based on limma adjusted P-value
      {if(name.sigProteins==TRUE){geom_text_repel(data=signData, aes(x = signData[,logFC_col_num], y = -log10(signData[,pval_col_num]), label = symbol), size = 2,
                                                  color = ifelse(signData[,logFC_col_num] < 0,'#00239CFF','#E10600FF'),
                                                  fontface = "bold", point.padding = unit(0.5, 'lines'), box.padding = unit(.25, 'lines'), max.overlaps = Inf)}} +

      geom_point(data = sigmaData, mapping = aes(x = sigmaData[,logFC_col_num], y = -log10(sigmaData[,pval_col_num])), shape = 1, size = 2, color = ifelse(sigmaData[,logFC_col_num] <0, "#00239CFF", "#E10600FF")) + # Adjusted p value based on limma adjusted P-value
      {if(name.sigProteins==TRUE){geom_text_repel(data=sigmaData, aes(x = sigmaData[,logFC_col_num], y = -log10(sigmaData[,pval_col_num]), label = symbol), size = 2,
                                                  color = ifelse(sigmaData[,logFC_col_num] < 0,'#00239CFF','#E10600FF'),
                                                  fontface = "italic", point.padding = unit(0.5, 'lines'), box.padding = unit(.25, 'lines'), max.overlaps = Inf)}} +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold.italic", colour = "#000000")) +
      theme(axis.text.y = element_text(size = 10, face = "bold.italic", colour = "#000000")) +
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0))) +
      theme(axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0))) +
      theme(axis.title.y = element_text(size = 12, face = "bold", angle = 90)) +
      theme(axis.title.x = element_text(size = 12, face = "bold")) +
      theme(plot.subtitle = element_text(size = 16, face = "bold")) +
      scale_x_continuous(breaks = seq(from = round(min(allSigData[,logFC_col_num]-0.5), 0), to = round(max(allSigData[,logFC_col_num]+0.5), 0), by = 1)) +
      scale_y_continuous(breaks = seq(from = 0, to = max(-log10(allSigData[,pval_col_num])) + 0.5, by = 1)) +
      geom_text(aes(x = min(resData[[i]][,logFC_col_num] + 1), y = 0, label = paste("Significant: ", (nrow(signData) + nrow(sigmaData))), fontface = "bold.italic"), colour = "darkgreen") +
      {if(length(conditions)==2){geom_text(aes(x = min(resData[[i]][,logFC_col_num] + 1), y = max(-log10(resData[[i]][,pval_col_num]) + 1), label = paste("Increased in", conditions[2], sep = " "), fontface = "bold"), size = 5, colour = "blue")}} +
      {if(length(conditions)==2){geom_text(aes(x = max(resData[[i]][,logFC_col_num] - 1), y = max(-log10(resData[[i]][,pval_col_num]) + 1), label = paste("Increased in", conditions[1], sep = " "), fontface = "bold"), size = 5, colour = "red")}}
    print(volcano)
  }
  dev.off()
}

