#' Boxplots of PSI distributions by splicing type.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}
#' @return a ggplot object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' boxplot_PSI_levels(hypoxia_filt, type = "RI")
#' @export
#' @import ggplot2
#' @importFrom stats median
#' @import methods

boxplot_PSI_levels <- function(events, type = c("A3SS", "A5SS", "SE", "RI",
                                                "MXE")){
    
    Sample <- NULL
    
    if(!is(events, "Maser")){
      stop("Parameter events has to be a maser object.")
    }
    
    type <- match.arg(type)
    
    events <- as(events, "list")
    PSI <- events[[paste0(type,"_","PSI")]]

    PSI_long <- reshape2::melt(PSI)
    colnames(PSI_long) <- c("ID", "Sample", "PSI")

    Condition <- rep("NA",nrow(PSI_long))
    idx.cond1 <- grep(paste0("^", events$conditions[1]), x = PSI_long$Sample,
                      perl = TRUE)
    idx.cond2 <- grep(paste0("^", events$conditions[2]), x = PSI_long$Sample,
                      perl = TRUE)
    Condition[idx.cond1] <- events$conditions[1]
    Condition[idx.cond2] <- events$conditions[2]

    PSI_long <- cbind(PSI_long, Condition)

    ggplot(PSI_long, aes(x = Sample, y = PSI, fill = Condition, color = Condition)) +
        #geom_boxplot() +
        geom_violin(trim = FALSE, alpha = 0.8) +
        stat_summary(fun.y=median, geom="point", size=2, color="black") +
        theme_bw() +
        theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
              axis.text.y = element_text(size=12),
              axis.title.x = element_text(face="plain", colour="black", 
                                          size=12),
              axis.title.y = element_text(face="plain", colour="black", 
                                          size=12),
              legend.title=element_blank(),
              legend.text = element_text(face="plain", colour="black", 
                                         size=12)) +
        ylab(paste(type, "PSI")) +
        xlab("Sample") +
        scale_y_continuous(limits=c(-0.1, 1.05)) +
        scale_fill_manual(values = c("blue", "red") ) +
        scale_color_manual(values = c("blue", "red") )


}

#' Proportion of events by splicing type.
#' 
#' @param events a maser object.
#' @param fdr numeric, FDR (False Discovery Rate) cutoff.
#' @param deltaPSI numeric, absolute minimum PSI (Percent spliced-in) change
#' @return a ggplot object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' splicingDistribution(hypoxia_filt)
#' @export
#' @import ggplot2
#' @import methods
#' 
splicingDistribution <- function(events, fdr = 0.05, deltaPSI = 0.1){
    
    if(!is(events, "Maser")){
      stop("Parameter events has to be a maser object.")
    }
    
    events <- as(events, "list")
    
  # Plot distribution of splicing events per condition
    as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
    nevents_cond1 <- rep(0, length(as_types))
    nevents_cond2 <- rep(0, length(as_types))
    
    FDR <- NULL
    IncLevelDifference <- NULL
    Condition <- NULL
    Proportion <- NULL
    Type <- NULL

    for (i in 1:length(as_types)) {

        stats <- events[[paste0(as_types[i],"_","stats")]]
        cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)
        cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < 
                                 (-1*deltaPSI))

        nevents_cond1[i] <- length(cond1$ID)
        nevents_cond2[i] <- length(cond2$ID)
    }

    nevents_prop1 <- nevents_cond1/sum(nevents_cond1)
    nevents_prop2 <- nevents_cond2/sum(nevents_cond2)
    condition <- c(rep(events$conditions[1], length(as_types)),
                   rep(events$conditions[2], length(as_types)))
    condition <- factor(condition, levels = c(events$conditions[1],
                                              events$conditions[2]))

    df.plot <- data.frame(Condition = condition,
                          Type = c(as_types, as_types),
                          Proportion = c(nevents_prop1, nevents_prop2))

    ggplot(df.plot, aes(x = Condition, y = Proportion,
                        colour = Type, fill = Type)) +
        geom_bar(stat = "identity", alpha = 0.6) +
        theme_bw() +
        theme(axis.text.y = element_text(size=12, angle = 0, hjust = 0.5,
                                         face = "plain"),
              axis.text.x = element_text(size=12, angle = 0, hjust = 0.5,
                                         face = "plain"),
              axis.title.x = element_text(face="plain", colour="black",
                                          size=12),
              axis.title.y = element_text(face="plain", colour="black",
                                          size=12),
              legend.title=element_blank(),
              legend.text = element_text(face="plain", colour="black",
                                         size=12),
              panel.grid=element_blank()
        ) +
        ylab("Proportion of splicing events") +
        xlab("") +
        scale_fill_brewer(palette="Set2") +
        scale_color_brewer(palette="Set2") +
        coord_flip()


}

#' Volcano plot of splicing events.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}
#' @param fdr numeric, FDR (False Discovery Rate) cutoff.
#' @param deltaPSI numeric, absolute minimum PSI (Percent spliced-in) change
#' @return a ggplot object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' volcano(hypoxia_filt, type = "SE")
#' @export
#' @import ggplot2
#' @importFrom dplyr filter
#' @import methods
volcano <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE"),
                    fdr = 0.05, deltaPSI = 0.1){
    
    if(!is(events, "Maser")){
      stop("Parameter events has to be a maser object.")
    }
    
    type <- match.arg(type)
    
    events <- as(events, "list")
    
    IncLevelDifference <- NULL
    Status <- NULL

    stats <- events[[paste0(type,"_","stats")]]
    cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)
    cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < 
                             (-1*deltaPSI))

    status <- rep("Not significant", times = nrow(stats))
    status[stats$ID %in% cond1$ID] <- events$conditions[1]
    status[stats$ID %in% cond2$ID] <- events$conditions[2]

    FDR <- stats$FDR
    idx_zero <- which(stats$FDR == 0)
    idx_min_nonzero <- max(which(stats$FDR == 0))+1
    FDR[idx_zero] <- FDR[idx_min_nonzero]
    log10pval <- -1*log10(FDR)

    plot.df <- data.frame(ID = stats$ID,
                          deltaPSI = stats$IncLevelDifference,
                          log10pval = log10pval,
                          Status = factor(status,
                                          levels = c("Not significant",
                                                     events$conditions[1],
                                                     events$conditions[2])))
    if(length(unique(status)) < 3){
        colors <-  c("blue", "red")
    } else{
        colors <-  c("grey","blue", "red")
    }

    ggplot(plot.df, aes(x=deltaPSI, y=log10pval, colour=Status)) +
        geom_point(aes(colour = Status)) +
        scale_colour_manual(values = colors) +
        theme_bw() +
        theme(axis.text.x = element_text(size=12), axis.text.y = 
                element_text(size=12),
              axis.title.x = element_text(face="plain", colour="black", 
                                          size=12),
              axis.title.y = element_text(face="plain", colour="black", 
                                          size=12),
              panel.grid.minor = element_blank(),
              plot.background = element_blank()
        ) +
        labs(title="", x = "Log10 Adj. Pvalue",
             y = "Delta PSI" )
}

#' Dotplot representation of splicing events.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}
#' @param fdr numeric, FDR (False Discovery Rate) cutoff.
#' @param deltaPSI numeric, absolute minimum PSI (Percent spliced-in) change
#' @return a ggplot object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' dotplot(hypoxia_filt, type = "SE")
#' @export
#' @import ggplot2
#' @importFrom dplyr filter
#' @import methods

dotplot <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE"),
                    fdr = 0.05, deltaPSI = 0.1){
    
    if(!is(events, "Maser")){
      stop("Parameter events has to be a maser object.")
    }
    
    type <- match.arg(type)
    
    events <- as(events, "list")
    FDR <- NULL
    IncLevelDifference <- NULL
    Status <- NULL
    
    stats <- events[[paste0(type,"_","stats")]]
    cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)
    cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < 
                             (-1*deltaPSI))

    idx.cond1 <- seq(1, events$n_cond1, 1)
    idx.cond2 <- seq(events$n_cond1+1, events$n_cond1+events$n_cond2, 1)

    PSI <- events[[paste0(type,"_","PSI")]]
    psi1 <- rowMeans(PSI[, idx.cond1], na.rm = TRUE)
    psi2 <- rowMeans(PSI[, idx.cond2], na.rm = TRUE)

    status <- rep("Not significant", times = nrow(PSI))
    status[rownames(PSI) %in% cond1$ID] <- events$conditions[1]
    status[rownames(PSI) %in% cond2$ID] <- events$conditions[2]


    plot.df <- data.frame(ID = stats$ID, psi1 = psi1, psi2 = psi2,
                          Status = factor(status,
                                          levels = c("Not significant",
                                                     events$conditions[1],
                                                     events$conditions[2])))
    if(length(unique(status)) < 3){
        colors <-  c("blue", "red")
    } else{
        colors <-  c("grey","blue", "red")
    }

    ggplot2::ggplot(plot.df, ggplot2::aes(x=psi1, y=psi2, colour=Status)) +
        geom_point(aes(colour = Status)) +
        scale_colour_manual(values = colors) +
        theme_bw() +
        theme(axis.text.x = element_text(size=12), axis.text.y = 
                element_text(size=12),
              axis.title.x = element_text(face="plain", colour="black",
                                          size=12),
              axis.title.y = element_text(face="plain", colour="black",
                                          size=12),
              panel.grid.minor = element_blank(),
              plot.background = element_blank()
              ) +
        labs(title="", x = paste("Average", events$conditions[1]),
                       y = paste("Average", events$conditions[2]) )
}

#' Prinicipal component analysis of PSI distributions.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}
#' @return a ggplot object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' pca(hypoxia_filt, type = "RI")
#' @export
#' @import ggplot2
#' @importFrom stats prcomp
#' @import methods
#' 
pca <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")){
    
    if(!is(events, "Maser")){
      stop("Parameter events has to be a maser object.")
    }
    
    type <- match.arg(type)
    events <- as(events, "list")

    PC1 <- NULL
    PC2 <- NULL
    Condition <- NULL
    Samples <- NULL
    
    idx.cond1 <- seq(1, events$n_cond1, 1)
    idx.cond2 <- seq(events$n_cond1+1, events$n_cond1+events$n_cond2, 1)

    PSI <- events[[paste0(type,"_","PSI")]]
    pheno <- rep(NA, ncol(PSI))
    pheno[ idx.cond1 ] <- events$conditions[1]
    pheno[ idx.cond2 ] <- events$conditions[2]
    pheno <- factor(pheno, levels = c(events$conditions[1],
                                      events$conditions[2]))

    res <- rowSums(PSI)
    idx.na <- which(is.na(res))
    if(length(idx.na)>0){
        PSI_notna <- PSI[-1*idx.na, ]
    } else {
        PSI_notna <- PSI
    }

    # PCA analysis
    my.pc <- prcomp(PSI_notna, center = FALSE, scale = FALSE)
    my.rot <- my.pc$r
    my.sum <- summary(my.pc)
    my.imp <- my.sum$importance

    df.pca <- data.frame(PC1 = my.rot[,1],
                         PC2 = my.rot[,2],
                         Condition = pheno,
                         Samples = colnames(PSI))

    percentVar <- round(c(my.imp[2,"PC1"], my.imp[2,"PC2"]) * 100)

    ggplot(df.pca, aes(PC1, PC2, color=Condition, label = Samples ))+
        geom_point(size=5) +
        #geom_text(vjust = 1, hjust = 0) +
        scale_colour_manual(values = c("blue", "red")) +
        theme_bw() +
        theme(axis.text.x = element_text(size=12), axis.text.y = 
                element_text(size=12),
              axis.title.x = element_text(face="plain", colour="black", 
                                          size=12),
              axis.title.y = element_text(face="plain", colour="black", 
                                          size=12),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.background = element_blank(),
              legend.title = element_blank(),
              legend.position = "right" ) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance"))


}

#' @importFrom dplyr filter
#' @import ggplot2
#' @import methods
viewTopSplicedGenes <- function(events, types = c("A3SS", "A5SS", "SE", "RI",
                                                  "MXE"), n = 20){
  if(!is(events, "Maser")){
    stop("Parameter events has to be a maser object.")
  }
  
  types <- match.arg(types, several.ok = TRUE)
  events <- as(events, "list")
  
  type <- NULL
  count <- NULL
  desc <- NULL
  total <- NULL
  
  geneList <- lapply(types, function(atype){
    annot <- events[[paste0(atype,"_","events")]]
    return(annot$geneSymbol)
  })
  geneList <- unique(unlist(geneList))
  
  geneList_counts <- data.frame(gene = rep(geneList, each = length(types)), 
                            type = rep(types, length(geneList)), 
                            counts = rep(0, length(geneList)*length(types)))
  
  counts <- lapply(seq_along(geneList), function(i){
    gene_counts <- countGeneEvents(as(events, "Maser"), geneList[i])
    return(gene_counts$count)
  })
  geneList_counts$count <- unlist(counts)

  geneList_counts_filt <- dplyr::filter(geneList_counts, type %in% types)
  res <- dplyr::group_by(geneList_counts_filt, gene)
  res2 <- dplyr::summarise(res, total = sum(count)) 
  rankedGenes <- dplyr::arrange(res2, desc(total))
  
  genes_plot <- dplyr::filter(geneList_counts_filt, 
                              gene %in% rankedGenes$gene[1:n])
  
  genes_plot$gene <- factor(genes_plot$gene, levels = rankedGenes$gene[1:n])
  
  ggplot(genes_plot, aes(x=gene, y=count, colour = type, fill = type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1), 
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    scale_fill_brewer(palette="BrBG") +
    scale_color_brewer(palette="BrBG") +
    ylab("Splicing events") +
    xlab("Gene")
  
  
}

#' Visualization of splicing events annotation using an interactive data table.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}
#' @return a datatables object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' hypoxia_top <- topEvents(hypoxia_filt)
#' display(hypoxia_top, type = "SE")
#' @export
#' @importFrom DT datatable
#' @import methods
display <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")){
  
  if(!is(events, "Maser")){
    stop("Parameter events has to be a maser object.")
  }
  
  type <- match.arg(type)
  events <- as(events, "list")
  
  data <- create_stats(as(events, "Maser"), type)
  DT::datatable(data, options = list(
                        pageLength = 10,
                        filter = "none",
                        searchHighlight = TRUE,
                        rownames = FALSE,
                        style = "bootstrap"
                        ),
                escape = FALSE,
                rownames = FALSE,
                selection = "none",
                filter = 'top'
                )
  
}


  