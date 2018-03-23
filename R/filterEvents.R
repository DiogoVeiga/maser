
filterByCoverage <- function(events, avg_reads = 5){
    
    if(!is.maser(events)){
      stop("Parameter events has to be a maser object.")
    }
  
    as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
    # Re-create events list by coverage filtering
    events_filt <- events

    for (type in as_types){

        # Find event ids with avg_reads > threshold
        counts <- events[[paste0(type,"_","counts")]]
        res_id <- rownames(counts)[rowMeans(counts) > avg_reads]

        # filter read counts matrix
        events_filt[[paste0(type,"_","counts")]] <-
                            counts[ rownames(counts) %in% res_id, ]

        # filter PSI matrix
        PSI <- events[[paste0(type,"_","PSI")]]
        events_filt[[paste0(type,"_","PSI")]] <-
                            PSI[ rownames(PSI) %in% res_id, ]

        # filter rMATS stats
        stats <- events[[paste0(type,"_","stats")]]
        events_filt[[paste0(type,"_","stats")]] <- dplyr::filter(stats,
                                                    ID %in% res_id)

        # Filter Genomic ranges of alternative splicing events
        grl <- events[[paste0(type,"_","gr")]]
        grl_new <- grl
        for (exon in names(grl)) {

            exon.gr <- grl[[exon]]
            grl_new[[exon]] <- exon.gr[exon.gr$ID %in% res_id, ]

        }
        events_filt[[paste0(type,"_","gr")]] <- grl_new

        # Filter Event annotation
        annot <- events[[paste0(type,"_","events")]]
        events_filt[[paste0(type,"_","events")]] <-
                    dplyr::filter(annot, ID %in% res_id)

        #cat("Selecting  ", type, length(res_id), " events\n")

    } # for each event type

    return(events_filt)

}

topEvents <- function(events, fdr = 0.05, deltaPSI = 0.1){
    
    if(!is.maser(events)){
      stop("Parameter events has to be a maser object.")
    }
    
    as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
    events_top <- events

    for (type in as_types){

        stats <- events[[paste0(type,"_","stats")]]
        res <- dplyr::filter(stats, FDR < fdr,
                         abs(IncLevelDifference) > deltaPSI)

        # filter read counts matrix
        counts <- events[[paste0(type,"_","counts")]]
        events_top[[paste0(type,"_","counts")]] <- counts[ rownames(counts) %in% res$ID, ]

        # filter PSI matrix
        PSI <- events[[paste0(type,"_","PSI")]]
        events_top[[paste0(type,"_","PSI")]] <- PSI[ rownames(PSI) %in% res$ID, ]

        # filter rMATS stats
        events_top[[paste0(type,"_","stats")]] <- dplyr::filter(stats,
                                                                ID %in% res$ID)

        # Filter Genomic ranges of alternative splicing events
        grl <- events[[paste0(type,"_","gr")]]
        grl_new <- grl
        for (exon in names(grl)) {

            exon.gr <- grl[[exon]]
            grl_new[[exon]] <- exon.gr[exon.gr$ID %in% res$ID, ]

        }
        events_top[[paste0(type,"_","gr")]] <- grl_new

        # Filter Event annotation
        annot <- events[[paste0(type,"_","events")]]
        events_top[[paste0(type,"_","events")]] <- dplyr::filter(annot,
                                                                 ID %in% res$ID)

        #cat("Selecting  ", type, length(res$ID), " events\n")

    } #each event type

    # Number of samples condition 1 and 2
    events_top[["n_cond1"]] <- events[["n_cond1"]]
    events_top[["n_cond2"]] <- events[["n_cond2"]]

    return(events_top)

}

geneEvents <- function(events, geneS, fdr = 0.05, deltaPSI = 0.1){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  events_top <- events
  
  for (type in as_types){
    
    annot <- events[[paste0(type,"_","events")]]
    stats <- events[[paste0(type,"_","stats")]]
    res <- dplyr::filter(stats, FDR < fdr,
                         abs(IncLevelDifference) > deltaPSI)
    
    resAnnot <- dplyr::filter(annot, geneSymbol %in% geneS)
    
    keepIDs <- intersect(res$ID,resAnnot$ID)
    
    # filter read counts matrix
    counts <- events[[paste0(type,"_","counts")]]
    events_top[[paste0(type,"_","counts")]] <- counts[ rownames(counts) %in% keepIDs, ]
    
    # filter PSI matrix
    PSI <- events[[paste0(type,"_","PSI")]]
    events_top[[paste0(type,"_","PSI")]] <- PSI[ rownames(PSI) %in% keepIDs, ]
    
    # filter rMATS stats
    events_top[[paste0(type,"_","stats")]] <- dplyr::filter(stats,
                                                            ID %in% keepIDs)
    
    # Filter Genomic ranges of alternative splicing events
    grl <- events[[paste0(type,"_","gr")]]
    grl_new <- grl
    for (exon in names(grl)) {
      
      exon.gr <- grl[[exon]]
      grl_new[[exon]] <- exon.gr[exon.gr$ID %in% keepIDs, ]
      
    }
    events_top[[paste0(type,"_","gr")]] <- grl_new
    
    # Filter Event annotation
    events_top[[paste0(type,"_","events")]] <- dplyr::filter(annot,
                                                             ID %in% keepIDs)
    
    #cat("Selecting  ", type, length(keepIDs), " events\n")
    
  } #each event type
  
  # Number of samples condition 1 and 2
  events_top[["n_cond1"]] <- events[["n_cond1"]]
  events_top[["n_cond2"]] <- events[["n_cond2"]]
  
  return(events_top)
  
}


# Internal function
countGeneEvents <- function(events, geneS){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  event_counts <- c()

  for (type in as_types){
    
    annot <- events[[paste0(type,"_","events")]]
    resAnnot <- dplyr::filter(annot, geneSymbol %in% geneS)
    event_counts <- c(event_counts, length(resAnnot$ID))
    
  } #each event type
  
  res <- data.frame(gene = geneS, type = as_types, count = event_counts)
  return(res)
  
}


