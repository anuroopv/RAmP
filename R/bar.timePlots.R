#' @title bar.timePlots
#'
#' @description Generates barplots or timeseries plots
#'
#' @param fav.proteins Vector with favourite proteins (or their corresponding sites) for which corresponding bar or timeseries plots are required
#' @param timeSeries Logical (Default is FALSE). If TRUE, line plots will be generated instead of barplots
#' @param imputed.data Imputed data obtained after QC.filter
#' @param lfq.data Output from editData function
#'
#' @return Generates bar plot (or timeseries) for the selected protein(s) or site(s)
#'
#' @export
#'
#' @importFrom Rmisc summarySE
################################################################################
# Barplot of for selected proteins/sites

bar.timePlots <- function(imputed.data, fav.proteins, timeSeries = FALSE, lfq.data){

  lfq.data.sub <- lfq.data[lfq.data$symbol %in% fav.proteins,]

  if(!is.null(fav.proteins) & nrow(lfq.data.sub != 0)){

    lfq.data.sub <- lfq.data[lfq.data$symbol %in% fav.proteins,]

    conditions <- unlist(strsplit(gsub("_vs_|[()]", ",", contrasts), ","))
    label_idx <- sampleTable[sampleTable$condition %in% conditions,]$label

    data.cont <- lfq.data.sub[,colnames(lfq.data.sub) %in% label_idx]
    data.cont$symbol <- lfq.data.sub$symbol
    data.cont$Uniprot <- lfq.data.sub$Uniprot

    if(Fraction == "Enriched"){data.cont$Sequence <- lfq.data.sub$Sequence[[1]]
    data.cont$Sequence <- make.unique(data.cont$Sequence)}

    if(Fraction == "Proteome"){names(data.cont)[1:(ncol(data.cont)-2)] <- sampleTable[grep(Fraction, sampleTable$fraction),]$condition[match(names(data.cont)[1:(ncol(data.cont)-2)], sampleTable[grep(Fraction, sampleTable$fraction),]$label)]
    names(data.cont)[1:(ncol(data.cont)-2)] <- make.unique(names(data.cont)[1:(ncol(data.cont)-2)])
    data.cont <- gather(data.cont, condition, measurement, 1:(length(colnames(data.cont))-2), factor_key=TRUE)
    }else{names(data.cont)[1:(ncol(data.cont)-3)] <- sampleTable[grep(Fraction, sampleTable$fraction),]$condition[match(names(data.cont)[1:(ncol(data.cont)-3)], sampleTable[grep(Fraction, sampleTable$fraction),]$label)]
    names(data.cont)[1:(ncol(data.cont)-3)] <- make.unique(names(data.cont)[1:(ncol(data.cont)-3)])
    data.cont <- gather(data.cont, condition, measurement, 1:(length(colnames(data.cont))-3), factor_key=TRUE)
    }

    data.cont$condition <- gsub("[.].*","",data.cont$condition)
    data.cont$logIntensity <- log2(data.cont$measurement)
    data.cont <- data.cont[complete.cases(data.cont$measurement),]

    summary.data <- if(Fraction == "Proteome"){summarySE(data = data.cont, measurevar = "measurement", groupvars = c("Uniprot", "symbol", "condition"), conf.interval = 0.95)}else{summarySE(data = data.cont, measurevar = "measurement", groupvars = c("Uniprot", "symbol", "condition", "Sequence"), conf.interval = 0.95)}
    summary.data$condition <- gsub(".*[.]", "", summary.data$condition)
    summary.data$genotype <- as.factor(gsub("[_].*", "", summary.data$condition))

    # if(nrow(summary.data) != 0){
      if(timeSeries == FALSE){
        summary.data$condition <- gsub(".*[_]", "", summary.data$condition)
        plot2 <- ggplot(data = summary.data, aes(x=condition, y=measurement, fill = genotype)) +
          geom_bar(stat="identity", position = position_dodge(.9)) +
          theme_pubr() +
          labs(x = "", y = "LFQ (or iBAQ)") +
          geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se), width = 0.2, position = position_dodge(0.9)) +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
          theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
          theme(strip.text = element_text(size=4)) +
          theme(legend.title = element_blank(), plot.subtitle = element_text(size = 10, face = "bold"), plot.caption = element_text(size = 8)) +
          theme(axis.title.x = element_text(size = 10, face = "bold")) +
          theme(axis.title.y = element_text(size = 10, face = "bold")) +
          {if(Fraction=="Proteome"){facet_wrap_paginate(~symbol, nrow = 1, ncol = 3, labeller = label_context, scales = "free_y")}else{facet_wrap_paginate(~symbol+Sequence,  nrow = 1, ncol = 3, labeller = label_context, scales = "free_y")}}
      }else{
        summary.data$condition <- as.numeric(gsub(".*[_]", "", summary.data$condition))
        plot2 <- ggplot(data = summary.data, aes(x=condition, y=measurement, color = genotype)) +
          geom_line(aes(group=1)) +
          geom_point() +
          theme_pubr() +
          labs(x = "Time", y = "LFQ (or iBAQ)") +
          geom_errorbar(aes(ymin=measurement-se, ymax=measurement+se), width = 0.2, position = position_dodge(0.9)) +
          scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
          theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
          theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
          theme(strip.text = element_text(size=4)) +
          theme(legend.title = element_blank(), plot.subtitle = element_text(size = 10, face = "bold"), plot.caption = element_text(size = 8)) +
          theme(axis.title.x = element_text(size = 10, face = "bold")) +
          theme(axis.title.y = element_text(size = 10, face = "bold")) +
          {if(Fraction=="Proteome"){facet_wrap_paginate(~symbol, nrow = 1, ncol = 3, labeller = label_context, scales = "free_y")}else{facet_wrap_paginate(~symbol+Sequence, nrow = 1, ncol = 3, labeller = label_context, scales = "free_y")}}
      }

      dir.create(paste(path1,"/",Fraction,"/BarPlots",sep = ""), showWarnings = FALSE)
      pdf(paste(path1,"/",Fraction,"/BarPlots/",Fraction,"_IntensityPlots_selectedFeatures.pdf",sep = ""), paper = "USr")
      for(i in 1:(n_pages(plot2))){
        if(Fraction=="Proteome"){print(plot2 + facet_wrap_paginate(~symbol, ncol = 3, nrow = 1, labeller = label_context, page = i, scales = "free_y"))
        }else{print(plot2 + facet_wrap_paginate(~symbol+Sequence, ncol = 3, nrow = 1, labeller = label_context, page = i, scales = "free_y"))
        }
      }
      dev.off()
    }else if(!is.null(fav.proteins) & nrow(lfq.data.sub == 0)){
      print("The given protein(s) (or site(s)) were NOT found. Check the data or the format of the protein name. Provided symbols should exactly match the symbols given in the database, including the case of each letter")
    } else{
      print("Favourite protein(s) or site(s) not provided")
    }
  }
