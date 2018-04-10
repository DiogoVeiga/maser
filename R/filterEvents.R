filterByEventId <- function(events, event_id, type){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  event_filt <- list()
  
  annot <- events[[paste0(type,"_","events")]]
  idx.event <- grep(as.numeric(event_id), annot$ID)
  
  if(length(idx.event)==0){
    stop(cat("Event id not found."))
  }
  
  # filter read counts matrix
  counts <- events[[paste0(type,"_","counts")]]
  events_filt[[paste0(type,"_","counts")]] <-
    counts[ rownames(counts) %in% event_id, , drop = FALSE]
  
  # filter PSI matrix
  PSI <- events[[paste0(type,"_","PSI")]]
  events_filt[[paste0(type,"_","PSI")]] <-
    PSI[ rownames(PSI) %in% event_id, , drop = FALSE]
  
  # filter rMATS stats
  stats <- events[[paste0(type,"_","stats")]]
  events_filt[[paste0(type,"_","stats")]] <- dplyr::filter(stats,
                                                           ID %in% event_id)
  
  # Filter Genomic ranges of alternative splicing events
  grl <- events[[paste0(type,"_","gr")]]
  grl_new <- grl
  for (exon in names(grl)) {
    
    exon.gr <- grl[[exon]]
    grl_new[[exon]] <- exon.gr[exon.gr$ID %in% event_id, ]
    
  }
  events_filt[[paste0(type,"_","gr")]] <- grl_new
  
  # Filter Event annotation
  annot <- events[[paste0(type,"_","events")]]
  events_filt[[paste0(type,"_","events")]] <-
    dplyr::filter(annot, ID %in% event_id)
  
  
  events_filt[["n_cond1"]] <- events$n_cond1
  events_filt[["n_cond2"]] <- events$n_cond2
  events_filt[["conditions"]] <- events$conditions
  
  class(events_filt) <- "maser"
  return(events_filt)
  
}

#' Filter splice events based on coverage.
#' 
#' @param events a maser object.
#' @param avg_reads numeric, average number of reads covering the splice event.
#' @return a maser object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' @export

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
                  counts[ rownames(counts) %in% res_id, , drop = FALSE]

        # filter PSI matrix
        PSI <- events[[paste0(type,"_","PSI")]]
        events_filt[[paste0(type,"_","PSI")]] <-
                    PSI[ rownames(PSI) %in% res_id, , drop = FALSE]

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

#' Filter splice events based on false discovery rate and PSI change.
#' 
#' @param events a maser object.
#' @param fdr numeric, FDR (False Discovery Rate) cutoff.
#' @param deltaPSI numeric, absolute minimum PSI (Percent spliced-in) change
#' @return a maser object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' ## To select all events with minimum 10% change in PSI, and FDR lower than 0.01 
#' hypoxia_top <- topEvents(hypoxia, fdr = 0.01, deltaPSI = 0.1)
#' @export
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
        events_top[[paste0(type,"_","counts")]] <- 
            counts[ rownames(counts) %in% res$ID, , drop = FALSE]

        # filter PSI matrix
        PSI <- events[[paste0(type,"_","PSI")]]
        events_top[[paste0(type,"_","PSI")]] <- 
                PSI[ rownames(PSI) %in% res$ID, , drop = FALSE]

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

#' Retrieve splice events for a given gene.
#' 
#' @param events a maser object.
#' @param geneS a character indicating the gene symbol.
#' @param fdr numeric, FDR cutoff.
#' @param deltaPSI numeric, minimum PSI change. 
#' @return a maser object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_mib2 <- geneEvents(hypoxia_filt, geneS = "MIB2", fdr = 0.05, 
#' deltaPSI = 0.1)
#' @export
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
    events_top[[paste0(type,"_","counts")]] <- 
            counts[rownames(counts) %in% keepIDs, , drop = FALSE]
    
    # filter PSI matrix
    PSI <- events[[paste0(type,"_","PSI")]]
    events_top[[paste0(type,"_","PSI")]] <- 
            PSI[ rownames(PSI) %in% keepIDs, , drop = FALSE]
    
    # filter rMATS stats
    events_top[[paste0(type,"_","stats")]] <- 
                dplyr::filter(stats, ID %in% keepIDs)
    
    # Filter Genomic ranges of alternative splicing events
    grl <- events[[paste0(type,"_","gr")]]
    grl_new <- grl
    for (exon in names(grl)) {
      
      exon.gr <- grl[[exon]]
      grl_new[[exon]] <- exon.gr[exon.gr$ID %in% keepIDs, ]
      
    }
    events_top[[paste0(type,"_","gr")]] <- grl_new
    
    # Filter Event annotation
    events_top[[paste0(type,"_","events")]] <- 
                    dplyr::filter(annot, ID %in% keepIDs)
    #cat("Selecting  ", type, length(keepIDs), " events\n")
    
  } #each event type
  
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


