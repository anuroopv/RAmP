#' @title motif.analysis
#'
#' @description Generates motifs using the peptide information of the enriched data (both exclusive and non-exclusive significant ones) using rmotifx
#'
#' @param Ex.data Data frame(s) (as list) containing exclsive sites (generated only of filter.protein.type = "fraction")
#' @param nonEx.data Data frame(s) (as list) containing non-exclusive sites (Obtained after limma differential expression analysis)
#' @param fasta fasta file from uniprot (same fasta flies used for MaxQuant search)
#' @param aa Central amino acid. For example, "K" (for acetylome), "S", "T", "Y", "STY" (for phosphoproteome). Default is "K".
#' @param seq.width Width of the sequence for motif search. Default is 15
#' @param min.seqs This threshold refers to the minimum number of times you wish each of your extracted motifs to occur in the data set. Default is 20 and is usually appropriate, although this parameter may be adjusted to yield more specific or less specific motifs.
#' @param p.value.motif The p-value threshold for the binomial probability. This is used for the selection of significant residue/position pairs in the motif. Default is 1e-05 and is suggested to maintain a low false positive rate in standard protein motif analyses.
#' @param pvalCutOff P-value cut off for significant protein(s), site(s) and GO terms. Default is 0.05
#' @param sigmaCutOff PI-value cut off for significant protein(s) and site(s). Default is 0.05 (Refer Xiao et al, 2014 and Hostrup et al, 2022 for details on PI-value)
#' @param lfcCutOff Log-fold change cut off. Default is 0
#' @param contrasts Mentions the conditions to be compared. Ex. MUTANT_vs_WILDTYPE or (MUTANT-A_vs_WILDTYPE-A)_vs_(MUTANT-B_vs_WILDTYPE-B). This how the contrasts should be provided and the values should be the same as the one given in condition column of sampleTable. (MUTANT-A_vs_WILDTYPE-A)_vs_(MUTANT-B_vs_WILDTYPE-B): This type can be used for complex data when comparing the interaction between two conditions such genotype and time
#' @param sampleTable .xlsx file containing information about the samples. Three columns are mandatory (label, condition and replicate)
#'
#' @return Image of the most representated motif(s) and the corresponding information as a file (Check rmotifx package for more information)
#'
#' @export
#'
#' @import PTMphinder
#' @import rmotifx
#' @import ggseqlogo
################################################################################

# Motif analysis of significant enriched sites

motif.analysis <- function(raw.enrichData, Ex.data, nonEx.data, fasta, aa = "K", seq.width = 15, min.seqs = 20, p.value.motif = 1e-5,
                           contrasts, sampleTable,
                           pvalCutOff = 0.05, sigmaCutOff = 0.05, lfcCutOff = 0){

  enrichdata.sub <- raw.enrichData %>%
    select(ends_with("Probabilities"),
           starts_with("Number"),
           starts_with("Peptide.IDs"))

  colnames(enrichdata.sub) <- c("Sequence.probability", "Number.of.sites", "Peptide.IDs")

  enrich.sub <- list()
  for(i in 1:length(contrasts)){

    conditions <- unlist(strsplit(gsub("_vs_|[()]", ",", contrasts[i]), ","))

    j=1
    if(is.null(Ex.data) == FALSE & str_count(contrasts[i], "vs") == 1){
      enrich.excl <- Ex.data[[i]]
      enrich.excl.sub1 <- subset(enrich.excl, enrich.excl[,grep(conditions[j], colnames(enrich.excl))] <= (length(unique(sampleTable$replicate))-2) & enrich.excl[,grep(conditions[j+1], colnames(enrich.excl))] >= (length(unique(sampleTable$replicate))-1))
      enrich.excl.sub1 <- enrich.excl.sub1[,c(1,4:ncol(enrich.excl.sub1))]
      colnames(enrich.excl.sub1) <- c("name", "Uniprot", "symbol", "Sequence.probability")

      enrich.excl.sub2 <- subset(enrich.excl, enrich.excl[,grep(conditions[j+1], colnames(enrich.excl))] <= (length(unique(sampleTable$replicate))-2) & enrich.excl[,grep(conditions[j], colnames(enrich.excl))] >= (length(unique(sampleTable$replicate))-1))
      enrich.excl.sub2 <- enrich.excl.sub2[,c(1,4:ncol(enrich.excl.sub2))]
      colnames(enrich.excl.sub2) <- c("name", "Uniprot", "symbol", "Sequence.probability")
    }else{
      enrich.excl.sub1 <- NULL
      enrich.excl.sub2 <- NULL
    }

    enrich.logFC <- nonEx.data[[i]]
    enrich.logFC.sub <- enrich.logFC[,grep(contrasts[i], colnames(enrich.logFC))]

    enrich.logFC.sub$name <- enrich.logFC$name
    enrich.logFC.sub$Uniprot <- enrich.logFC$Uniprot
    enrich.logFC.sub$symbol <- enrich.logFC$symbol
    enrich.logFC.sub$Sequence.probability <- enrich.logFC[[5]]

    enrich.logFC.sub.1 <- subset(enrich.logFC.sub, (enrich.logFC.sub[,grep("q.val", colnames(enrich.logFC.sub))] < pvalCutOff | enrich.logFC.sub[,grep("SIGMA", colnames(enrich.logFC.sub))] < sigmaCutOff) & enrich.logFC.sub[,grep("logFC", colnames(enrich.logFC.sub))] > lfcCutOff)
    enrich.logFC.sub.1 <- enrich.logFC.sub.1[,c(12:ncol(enrich.logFC.sub.1))]

    enrich.logFC.sub.2 <- subset(enrich.logFC.sub, (enrich.logFC.sub[,grep("q.val", colnames(enrich.logFC.sub))] < pvalCutOff | enrich.logFC.sub[,grep("SIGMA", colnames(enrich.logFC.sub))] < sigmaCutOff) & enrich.logFC.sub[,grep("logFC", colnames(enrich.logFC.sub))] < lfcCutOff)
    enrich.logFC.sub.2 <- enrich.logFC.sub.2[,c(12:ncol(enrich.logFC.sub.2))]


    if(is.null(Ex.data) == FALSE){if(identical(x = colnames(enrich.excl.sub1), y = colnames(enrich.logFC.sub.1))==FALSE){stop()}}
    enrich.sub1 <- rbind(enrich.logFC.sub.1, enrich.excl.sub1)
    enrich.sub1 <- enrich.sub1[!duplicated(enrich.sub1),]

    enrich.sub1 <- merge(enrich.sub1, enrichdata.sub, by = "Sequence.probability")

    if(is.null(Ex.data) == FALSE){if(identical(x = colnames(enrich.excl.sub2), y = colnames(enrich.logFC.sub.2))==FALSE){stop()}}
    enrich.sub2 <- rbind(enrich.logFC.sub.2, enrich.excl.sub2)
    enrich.sub2 <- enrich.sub2[!duplicated(enrich.sub2),]

    enrich.sub2 <- merge(enrich.sub2, enrichdata.sub, by = "Sequence.probability")

    enrich.sub <- list(enrich.sub1, enrich.sub2)
    names(enrich.sub) <- paste(conditions, contrasts[i], sep="_")

    for(i in 1:length(conditions)){

      motif.data <- data.frame()
      motif.data <- as.data.frame(str_extract_all(enrich.sub[[i]]$Sequence.probability, "\\d+\\.*\\d*", simplify = TRUE))

      # Extract probabilities
      motif.data <- as.data.frame(sapply(motif.data, as.numeric)) # Make probabilities as numeric
      motif.data <- motif.data*100 # Probability to percentage
      motif.data <- unite(motif.data, "PTM_Score", sep = ";", remove = TRUE) # Join info from multiple probabilities for each site
      motif.data$PTM_Score <- gsub(";NA", "", motif.data$PTM_Score) # Remove NAs

      # Extract peptide sequence
      motif.data$Peptide_Seq <- enrich.sub[[i]]$Sequence.probability
      motif.data$Peptide_Seq <- gsub("([0-9]+)", "", motif.data$Peptide_Seq)
      motif.data$Peptide_Seq <- gsub("\\.", "", motif.data$Peptide_Seq)
      motif.data$Peptide_Seq <- gsub("\\(|\\)", "", motif.data$Peptide_Seq)

      # Obtain total number of sites possible for modification
      motif.data$Total_Sites <- enrich.sub[[i]]$Number.of.sites

      # Obtain peptide IDs
      motif.data$Identifier <- enrich.sub[[i]]$Peptide.IDs
      motif.data$Identifier <- gsub("[;].*", "", motif.data$Identifier)
      motif.data$Identifier <- paste("Pep", motif.data$Identifier ,sep="")

      # Obtain Uniprot IDs
      motif.data$Protein_ID <- enrich.sub[[i]]$Uniprot

      # Obtain the details of amino acid and corresponding sites in the peptide sequence

      a <- str_match_all(enrich.sub[[i]]$Sequence.probability, "(.*?)\\(\\d+\\.*\\d*\\)")
      names(a) <- motif.data$Peptide_Seq
      b <- do.call(rbind, lapply(a, data.frame))
      b <- tibble::rownames_to_column(b, "Peptide_Seq")
      b$Peptide_Seq <- gsub(".[0-9]", "", b$Peptide_Seq)
      b$X1 <- NULL
      colnames(b)[2] <- "Sequences"
      b$value <- b$Sequences
      b$count <- nchar(b$value)
      b$site <- substring(b$value, nchar(b$value))

      b <- setDT(b)[, Mod_pos := cumsum(count), Peptide_Seq][]
      #b %>% group_by(Peptide_Seq) %>% mutate(Mod_pos = cumsum(count))
      b$PTM_site <- paste(b$site, b$Mod_pos, sep = "")

      b.new <- as.data.frame(aggregate(PTM_site~ Peptide_Seq, data = b, paste, collapse = ";"))
      colnames(b.new) <- c("Peptide_Seq", "PTM_Loc")

      b.new$Peptide_Seq <- gsub("[0-9]+","", b.new$Peptide_Seq)

      data2 <- merge(motif.data, b.new, by = "Peptide_Seq")
      data2 <- data2[!duplicated(data2), ]
      data2 <- data2[,c(4,5,1,3,6,2)]

      # Obtain Uniprot fasta database

      ref_db <- fasta
      parseDB.data <- parseDB(ref = ref_db, db_source =  "UP", filt =  TRUE)

      s <- unlist(parseDB.data[,2])

      # Extract proteome-specific background motifs (may take a while depending on database size)
      extractBack <- extractBackground(s, aa, seq.width)

      phindPTMs.data <- phindPTMs(data = data2, reftab = parseDB.data)

      foreground_Seqs <- unlist(strsplit(phindPTMs.data[,"Flank_Seq"], split = "[.]"))
      foreground_Seqs_Filtered <- foreground_Seqs[which(lapply(foreground_Seqs, nchar)==seq.width)]

      # Run motif-x using foreground and background sequences from PTMphinder functions above
      motifx.data <- motifx(foreground_Seqs_Filtered, extractBack, central.res = aa, min.seqs = min.seqs, pval.cutoff = p.value.motif)

      # View motifx output
      dir.create(paste(getwd(),"/Results/Enriched/Motif_analysis",sep = ""), showWarnings = TRUE)
      pdf(paste(getwd(),"/Results/Enriched/Motif_analysis/",names(enrich.sub[i]),"_motif.pdf", sep = ""), width = 8, height = 3)

      if(!is.null(motifx.data)){print(ggseqlogo(motifx.data$motif, seq_type='aa', method = 'bits') + ggtitle(names(enrich.sub[i])))
        print(ggseqlogo(motifx.data$motif, seq_type='aa', method = 'prob')  + ggtitle(names(enrich.sub[i])))
      }else{
        print("No motifs matched")
      }
      dev.off()
      writexl::write_xlsx(path = paste(getwd(),"/Results/Enriched/Motif_analysis/",names(enrich.sub[i]),"_motif.xlsx", sep = ""), x = motifx.data, col_names = TRUE, format_headers = TRUE)
    }
  }
}
