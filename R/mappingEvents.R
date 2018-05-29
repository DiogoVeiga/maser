#' Mapping of splice events to UniprotKB protein features.
#' 
#' @param events a maser object with transcript and protein identifiers.
#' @param tracks a character vector indicating valid UniprotKB features or 
#' categories.
#' @param by a character vector, possible values 
#' are \code{c("feature", "category")}.
#' @param ncores number of cores for multithreading (available only in OSX and Linux 
#' machines). If Windows, \code{ncores} will be set to 1 automatically.
#' @return a maser object with protein feature annotation.
#' @details This function performs mapping of splicing events to protein
#'  features available in the UniprotKB database. Annotation tracks of protein
#'  features mapped to the hg38 build of the human genome are retrieved from the
#'  public UniprotKB FTP. The function will overlap exons involved in the splice
#'  event with the feature genomic coordinates retrieved from UniprotKB.
#' 
#' Annotation can be executed either by feature or category. If categories are 
#' provided, all features within the category group will be included for 
#' annotation.
#' 
#' Thus, batch annotation is enabled either by using \code{by = category} or 
#' by providing mutilple features in the \code{tracks} argument.
#' 
#' Visualization of protein features can be done 
#' using \code{\link{plotUniprotKBFeatures}}.
#'  
#' @examples
#' ## Create the maser object
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' 
#' ## Ensembl GTF annotation for SRSF6
#' gtf_path <- system.file("extdata", file.path("GTF", "SRSF6_Ensembl85.gtf"),
#'  package = "maser")
#' ens_gtf <- rtracklayer::import.gff(gtf_path)
#' 
#' ## Retrieve gene specific splice events
#' srsf6_events <- geneEvents(hypoxia_filt, geneS = "SRSF6")
#' 
#' ## Map splicing events to transcripts
#' srsf6_mapped <- mapTranscriptsToEvents(srsf6_events, ens_gtf)
#' 
#' ## Annotate splice events with protein domains
#' srsf6_annot <- mapProteinFeaturesToEvents(srsf6_mapped, tracks = "domain")
#' head(annotation(srsf6_annot, "SE"))
#' 
#' @seealso \code{\link{plotUniprotKBFeatures}}
#' @export
#' @import GenomicRanges
#' @importFrom dplyr filter
#' @importFrom parallel mclapply
#' @import methods

mapProteinFeaturesToEvents <- function(events, tracks, by = c("feature", 
                                        "category"), ncores = 1){
  
  by <- match.arg(by)
  
  if(!is(events, "Maser")){
    stop("Parameter events has to be a maser object.")
  }
  
  if(.Platform$OS.type == "windows"){
    ncores = 1
  }
  
  df <- availableFeaturesUniprotKB()
  
  if(by == "feature"){
    
    if (!any(tracks %in% as.vector(df$Name))){
      stop(cat("\"tracks\" arg is invalid."))
    }
    features <- tracks
  }
  
  Category <- NULL
  
  if(by == "category"){
    
    if (!any(tracks %in% as.vector(df$Category))){
      stop(cat("\"tracks\" arg is invalid."))
    }
    
    df_filt <- dplyr::filter(df, Category %in% tracks)
    features <- as.vector(df_filt$Name)
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  events <- as(events, "list")
  
  # Add UniprotKB features annotation
  events_with_features <- events
  
  # Create GRanges for each UniprotKB feature in features argument
  features_Gr <- mclapply(features, createGRangesUniprotKBtrack, 
                          mc.cores = ncores)
  features_Gr <- GRangesList(features_Gr)
  names(features_Gr) <- features
  
  for (type in as_types){
    
    # Retrieve Events annotation
    annot <- events[[paste0(type,"_","events")]]
    grl <- events[[paste0(type,"_","gr")]]
    idx.cols <- grep("^list_ptn_", colnames(annot))
    
    if (nrow(annot) == 0){
      next
    }
    
    #Create matrix for storing results
    annot_uniprotkb <- matrix("NA", nrow = nrow(annot), ncol = length(features))
    colnames(annot_uniprotkb) <- features
    
    lapply(seq(1, nrow(annot)), function(i){ #for each event
      
      eventGr <- GRangesList()
      for (exon in names(grl)){
        eventGr[[paste0(exon)]] <- grl[[paste0(exon)]][i]
      }
      
      protein_ids <- unique(c(annot[i, idx.cols[1]], annot[i, idx.cols[2]]))
      
      lapply(seq_along(features), function(j){ #for each feature
        
        ovl_gr <- overlappingFeatures(features_Gr[[j]], eventGr)
        #ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ]
        ovl_gr_filt <- ovl_gr
        
        aux <- unique(as.character(ovl_gr_filt$Name))
        aux <- gsub(" ", "", aux)
        res <- paste0(aux, collapse = ",")
        if (!res == ""){
          annot_uniprotkb[i,j] <<- res  
        }
        
      })
      
    })
    
    # for (i in 1:nrow(annot)) {
    #   
    #   # Genomic ranges of alternative splicing event
    #   eventGr <- GRangesList()
    #   for (exon in names(grl)){
    #     eventGr[[paste0(exon)]] <- grl[[paste0(exon)]][i]
    #   }
    #   
    #   protein_ids <- unique(c(annot[i, idx.cols[1]], annot[i, idx.cols[2]]))
    #   
    #   
    #   for (j in 1:length(features)) {
    #     ovl_gr <- overlappingFeatures(features_Gr[[j]], eventGr)
    #     ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ] 
    #     
    #     aux <- unique(as.character(ovl_gr_filt$Name))
    #     aux <- gsub(" ", "", aux)
    #     res <- paste0(aux, collapse = ",")
    #     if (!res == ""){
    #       annot_uniprotkb[i,j] <- res  
    #     }
    #     
    #   } #all features
    #   
    #   
    # }#all events
    
    #write annotation to maser object
    events_with_features[[paste0(type,"_","events")]] <- 
      cbind(annot, annot_uniprotkb)

  }#all types
  
  return(as(events_with_features, "Maser"))
}


mapENSTtoUniprotKB <- function(enst_ids){
  
  if (enst_ids == ""){
    return(c(""))
  }
  
  aux <- strsplit(enst_ids,",")[[1]]
  tokens <- strsplit(aux, "\\.")
  
  enst_trans <- vapply(seq_along(aux),
                       function(i) tokens[[i]][[1]], character(1))
  
  idx.map <- match(enst_trans, UKB_ENST_map$ENST_ID)
  
  return(as.vector(UKB_ENST_map$UniprotKB_ID[idx.map]))
  
}

mapProteinsToEvents <- function(events){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  
  # Add UniprotKB ID
  events_with_ptn <- events
  
  for (type in as_types){
    
    # Retrieve Events annotation
    annot <- events[[paste0(type,"_","events")]]
    idx.cols <- grep("^txn_", colnames(annot))
    
    if (nrow(annot) == 0){
      next
    }
    
    list_ptn_a <- rep("", nrow(annot))
    list_ptn_b <- rep("", nrow(annot))
    
    for (i in 1:nrow(annot)) {
      
      
      list_ptn_a[i] <- paste(mapENSTtoUniprotKB(annot[i,idx.cols[1]]), 
                                        collapse = ",")
      list_ptn_b[i] <- paste(mapENSTtoUniprotKB(annot[i,idx.cols[2]]), 
                                        collapse = ",")
      
    } #for all events in annot
    
    annot[["list_ptn_a"]] <- list_ptn_a
    annot[["list_ptn_b"]] <- list_ptn_b
    
    events_with_ptn[[paste0(type,"_","events")]] <- annot
    
    
  } #for each AS type
  
  return(events_with_ptn)
  
}


#' Mapping of splice events to Ensembl transcripts.
#' 
#' @param events a maser object.
#' @param gtf a \code{GRanges} object obtained from an Ensembl or Gencode GTF
#'  file using the hg38 build of the human genome.
#' @param ncores number of cores for multithreading (available only in OSX and Linux 
#' machines). If Windows, \code{ncores} will be set to 1 automatically.
#' @return a maser object with transcript and protein identifiers.
#' @details This function performs mapping of splice events in the maser object
#'  to Ensembl transcripts by overlapping exons involved in the splice event to
#'  the transcript models provided in the GTF. 
#'  
#'  Each type of splice event requires a specific mapping procedure
#'   (described below).
#'  
#'  The mapping will also add Uniprot identifiers when the ENST transcript 
#'  encodes for a protein. 
#'    
#'  Visualization of affected transcripts can be done 
#'  using \code{\link{plotTranscripts}}.
#'  
#'   \describe{
#'     \item{\strong{Exon skipping}}{}
#'     \item{Inclusion transcript(s)}{Transcript(s) overlapping the cassette
#'      exon, as well both flanking exons (i.e upstream and downstream exons).}
#'     \item{Skipping transcript(s)}{Transcript(s) overlapping both flanking 
#'     exons but not the cassettte exon.}
#'   }
#'   
#'   \describe{
#'     \item{\strong{Intron retention}}{}
#'     \item{Retention transcript(s)}{Transcript(s) overlapping exactly the
#'      retained intron.}
#'     \item{Skipping transcript(s)}{Transcript(s) where intron is spliced out
#'      and overlapping both flanking exons.}
#'   }
#'   
#'   \describe{
#'   \item{\strong{Mutually exclusive exons}}{}
#'     \item{Exon1 transcript(s)}{Transcript(s) overlapping the first exon 
#'     and both flanking exons.}
#'     \item{Exon2 transcript(s)}{Transcript(s) overlapping the second exon and
#'      both flanking exons.}
#'   }
#'   
#'   \describe{
#'     \item{\strong{Alternative 3' and 5' splice sites}}{}
#'     \item{Short exon transcript(s)}{Transcript(s) overlapping both short and 
#'                        downstream exons.}
#'     \item{Long exon transcript(s)}{Transcript(s) overlapping both long and
#'      downstream exons.}
#'   }
#'   
#' @examples
#' ## Create the maser object
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' 
#' ## Ensembl GTF annotation for SRSF6
#' gtf_path <- system.file("extdata", file.path("GTF", 
#'  "Ensembl85_examples.gtf.gz"), package = "maser")
#' ens_gtf <- rtracklayer::import.gff(gtf_path)
#'  
#' ## Retrieve gene specific splice events
#' srsf6_events <- geneEvents(hypoxia_filt, geneS = "SRSF6")
#' 
#' ## Map splicing events to transcripts
#' srsf6_mapped <- mapTranscriptsToEvents(srsf6_events, ens_gtf)
#' head(annotation(srsf6_mapped, "SE"))
#' 
#' @seealso \code{\link{plotTranscripts}}
#' @export
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom parallel mclapply
#' @importFrom dplyr inner_join
#' @import methods
#' 
mapTranscriptsToEvents <- function(events, gtf, ncores = 1){
  
  is_strict = TRUE
  
  if(!is(events, "Maser")){
    stop("Parameter events has to be a maser object.")
  }
  
  if (!class(gtf) == "GRanges"){
    stop(cat("\"gtf\" should be a GRanges object."))
  }
  
  if(.Platform$OS.type == "windows"){
    ncores = 1
  }

  #Add chr to seqnames - necessary for Gviz plots and compatible with maser()
  if(any(!grepl("chr", GenomeInfoDb::seqlevels(gtf)))){
    GenomeInfoDb::seqlevels(gtf) <- paste0("chr", 
                                           GenomeInfoDb::seqlevels(gtf)) 
  }
  
  gtf_exons <- gtf[gtf$type=="exon",]
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  events <- as(events, "list")
  
  # Add transcripts to events using mapping functions based on sequence overlap 
  events_with_txn <- events
  
  for (type in as_types){
    
    # Retrieve Events annotation
    annot <- events[[paste0(type,"_","events")]]
    
    if (nrow(annot) == 0){
      next
    }
    
    # Retrieve Events gene ranges
    grl <- events[[paste0(type,"_","gr")]]
    
    list_txn <- mclapply(seq_along(annot$ID), function(i) {
      
      # Genomic ranges of alternative splicing events
      eventGr <- lapply(names(grl), function(exon){
        grl[[exon]][i]
      })
      eventGr <- GRangesList(eventGr)
      names(eventGr) <- names(grl)
      
      if(type == "SE") {
        tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
        return(data.frame(ID = annot$ID[i],
                  txn_3exons = paste(tx_ids$txn_3exons, collapse = ","),
                  txn_2exons = paste(tx_ids$txn_2exons, collapse = ","),
                  stringsAsFactors = FALSE)
               )
        
      }
      
      if(type == "MXE") {
        tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)
        return(data.frame(ID = annot$ID[i],
                  txn_mxe_exon1 = paste(tx_ids$txn_mxe_exon1, collapse = ","),
                  txn_mxe_exon2 = paste(tx_ids$txn_mxe_exon2, collapse = ","),
                  stringsAsFactors = FALSE)
        )
      }
      
      if(type == "RI") {
        tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
        return(data.frame(ID = annot$ID[i],
          txn_nonRetention = paste(tx_ids$txn_nonRetention, collapse = ","),
          txn_retention = paste(tx_ids$txn_retention, collapse = ","),
          stringsAsFactors = FALSE)
        )
      }
      
      if(type == "A5SS") {
        
        #reverse strand becomes A3SS
        if (as.character(strand(eventGr[1])) == "+"){
          tx_ids <<- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
        }else{
          tx_ids <<- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
        }
        
        return(data.frame(ID = annot$ID[i],
                txn_short = paste(tx_ids$txn_short, collapse = ","),
                txn_long = paste(tx_ids$txn_long, collapse = ","),
                stringsAsFactors = FALSE)
        )
        
      }
      
      if(type == "A3SS") {
        
        #reverse strand becomes A5SS
        if (as.character(strand(eventGr[1])) == "+"){
          tx_ids <<- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
        }else{
          tx_ids <<- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
        }
        
        return(data.frame(ID = annot$ID[i],
                  txn_short = paste(tx_ids$txn_short, collapse = ","),
                  txn_long = paste(tx_ids$txn_long, collapse = ","),
                  stringsAsFactors = FALSE)
        )
        
      }

    }, mc.cores = ncores) #for all events in annot
    
    df_txn <- do.call(rbind, list_txn)

    # for (i in 1:nrow(annot)) {
    #   
    #   # Genomic ranges of alternative splicing events
    #   eventGr <- lapply(names(grl), function(exon){
    #     grl[[exon]][i]
    #   })
    #   eventGr <- GRangesList(eventGr)
    #   names(eventGr) <- names(grl)
    #   
    #   if(type == "SE") {
    #     tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
    #     list_txn_a[i] <- paste(tx_ids$txn_3exons, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_2exons, collapse = ",")
    #   }
    #   
    #   if(type == "MXE") {
    #     tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)
    #     list_txn_a[i] <- paste(tx_ids$txn_mxe_exon1, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_mxe_exon2, collapse = ",")
    #   }
    #   
    #   if(type == "RI") {
    #     tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
    #     list_txn_a[i] <- paste(tx_ids$txn_nonRetention, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_retention, collapse = ",")
    #   }
    #   
    #   if(type == "A5SS") {
    #     
    #     #reverse strand becomes A3SS
    #     if (as.character(strand(eventGr[1])) == "+"){
    #       tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
    #     }else{
    #       tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
    #     }
    #     
    #     list_txn_a[i] <- paste(tx_ids$txn_short, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_long, collapse = ",")
    #   }
    #   
    #   if(type == "A3SS") {
    #     
    #     #reverse strand becomes A5SS
    #     if (as.character(strand(eventGr[1])) == "+"){
    #       tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
    #     }else{
    #       tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
    #     }
    # 
    #     list_txn_a[i] <- paste(tx_ids$txn_short, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_long, collapse = ",")
    #   }
    #   
    # } #for all events in annot
    
    
    # Update Events annotation
    annot_new <- dplyr::inner_join(annot, df_txn, by = "ID")
    events_with_txn[[paste0(type,"_","events")]] <- annot_new
      
  
  } #for each event type
  
  events_with_ptn <- mapProteinsToEvents(events_with_txn)
  
  return(as(events_with_ptn, "Maser"))
  
}

#' @import GenomicRanges
mapTranscriptsSEevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_upstream, 
                                        gtf_exons, type = "any")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_target, 
                                        gtf_exons, type = "any")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_downstream, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range for inclusion event and skipping event
  intron.skipping <- GenomicRanges::GRanges(
                                     seqnames = seqnames(eventGr$exon_target),
                                     ranges = IRanges::IRanges(
                                     start = end(eventGr$exon_upstream) + 1,  
                                     end = start(eventGr$exon_downstream) - 1),
                                     strand = strand(eventGr$exon_target)
  )
  
  intron.inclusion <- GenomicRanges::GRanges(
                                     seqnames = seqnames(eventGr$exon_target),
                                     ranges = IRanges::IRanges(
                                     start = c(end(eventGr$exon_upstream) + 1,
                                                end(eventGr$exon_target) + 1),  
                                     end = c(start(eventGr$exon_target) - 1,
                                             start(eventGr$exon_downstream) -1)
                                             ),
                                     strand = strand(eventGr$exon_target)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.inclusion <- GenomicRanges::findOverlaps(intron.inclusion, 
                                                      gtf_exons, type = "any")
  mytx.ids.intron.inclusion <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.inclusion)]
  
  ovl.intron.skipping <- GenomicRanges::findOverlaps(intron.skipping, 
                                                     gtf_exons, type = "any")
  mytx.ids.intron.skipping <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.skipping)]
  
  #decide wich transcripts to plot in inclusion and skipping tracks
  if (is_strict){
    mytx.ids.3exons <- intersect(mytx.ids.e1, mytx.ids.e3)#has both flank exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and target exon
    mytx.ids.2exons <- intersect(mytx.ids.e1, mytx.ids.e3)
    
  }else {
    mytx.ids.3exons <- union(mytx.ids.e1, mytx.ids.e3)#has either flaking exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and target exon  
    mytx.ids.2exons <- union(mytx.ids.e1, mytx.ids.e3)
  }
  
  mytx.ids.3exons <- setdiff(mytx.ids.3exons, mytx.ids.intron.inclusion)
  #mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.3exons)
  mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.intron.skipping)
  
  tx_ids <- list()
  tx_ids[["txn_3exons"]] <- mytx.ids.3exons
  tx_ids[["txn_2exons"]] <- mytx.ids.2exons
  
  return(tx_ids)  
  
}

#' @import GenomicRanges
mapTranscriptsRIevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_upstream, 
                                        gtf_exons, type = "any")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_ir, 
                                        gtf_exons, type = "equal")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_downstream, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range from the retention event
  intron <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_ir),
                                   ranges = IRanges::IRanges(
                                     start = end(eventGr$exon_upstream) + 1,  
                                     end = start(eventGr$exon_downstream) - 1),
                                   strand = strand(eventGr$exon_ir)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron <- GenomicRanges::findOverlaps(intron, gtf_exons, type = "any")
  mytx.ids.intron <- gtf_exons$transcript_id[subjectHits(ovl.intron)]
  
  #decide wich transcripts to plot in retention and non-retention tracks
  if (is_strict){
    #has both upstream and downstream exons
    tx.ids.nonRetention <- intersect(mytx.ids.e1, mytx.ids.e3) 
    
  }else {
    #has either upstream and downstream exons
    tx.ids.nonRetention <- union(mytx.ids.e1, mytx.ids.e3) 
  }
  
  tx.ids.nonRetention <- setdiff(tx.ids.nonRetention, mytx.ids.intron)
  tx.ids.Retention <- mytx.ids.e2
  
  tx_ids <- list()
  tx_ids[["txn_nonRetention"]] <- tx.ids.nonRetention
  tx_ids[["txn_retention"]] <- tx.ids.Retention
  return(tx_ids)
  
}

#' @import GenomicRanges
mapTranscriptsMXEevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_1, 
                                        gtf_exons, type = "any")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_2, 
                                        gtf_exons, type = "any")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_upstream,
                                        gtf_exons, type = "any")
  ovl.e4 <- GenomicRanges::findOverlaps(eventGr$exon_downstream,
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  mytx.ids.e4 <- gtf_exons$transcript_id[subjectHits(ovl.e4)]
  
  #obtain intron range for inclusion event and skipping event
  intron.mxe.exon1 <- GenomicRanges::GRanges(
                                    seqnames = seqnames(eventGr$exon_1),
                                    ranges = IRanges::IRanges(
                                    start = c(end(eventGr$exon_upstream) + 1,
                                              end(eventGr$exon_1) + 1),  
                                    end = c(start(eventGr$exon_1) - 1,
                                            start(eventGr$exon_downstream) -1)
                                     ),
                                    strand = strand(eventGr$exon_1)
  )
  
  intron.mxe.exon2 <- GenomicRanges::GRanges(
                                      seqnames = seqnames(eventGr$exon_2),
                                      ranges = IRanges::IRanges(
                                      start = c(end(eventGr$exon_upstream) + 1,
                                                 end(eventGr$exon_2) + 1),  
                                      end = c(start(eventGr$exon_2) - 1,
                                             start(eventGr$exon_downstream) -1)
                                             ),
                                      strand = strand(eventGr$exon_2)
                                      )
  
  #find transcripts with exons overlapping intronic regions
  ovl.mxe.exon1 <- GenomicRanges::findOverlaps(intron.mxe.exon1, 
                                               gtf_exons, type = "any")
  mytx.ids.intron1 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon1)]
  
  ovl.mxe.exon2 <- GenomicRanges::findOverlaps(intron.mxe.exon2, 
                                               gtf_exons, type = "any")
  mytx.ids.intron2 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon2)]
  
  
  #decide wich transcripts to plot in inclusion and skipping tracks
  if (is_strict){
    #has both flanking exons
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1) #and exon1
    
    #has both flanking exons
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2) #and exon2
    
  }else {
    #has either flanking exons
    mytx.ids.mxe.exon1 <- union(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1) #and exon 1
    
    #has both flanking exons
    mytx.ids.mxe.exon2 <- union(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2) #and exon2
  }
  
  #remove transcripts with exons in intronic regions
  mytx.ids.mxe.exon1 <- setdiff(mytx.ids.mxe.exon1, mytx.ids.intron1)
  mytx.ids.mxe.exon2 <- setdiff(mytx.ids.mxe.exon2, mytx.ids.intron2)
  
  tx_ids <- list()
  tx_ids[["txn_mxe_exon1"]] <- mytx.ids.mxe.exon1
  tx_ids[["txn_mxe_exon2"]] <- mytx.ids.mxe.exon2
  return(tx_ids)

}

#' @import GenomicRanges
mapTranscriptsA5SSevent <- function(eventGr, gtf_exons){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_short, 
                                        gtf_exons, type = "equal")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_long, 
                                        gtf_exons, type = "equal")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_flanking, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range for short event and long event
  intron.short <- GenomicRanges::GRanges(
                                      seqnames = seqnames(eventGr$exon_short),
                                      ranges = IRanges::IRanges(
                                      start = end(eventGr$exon_short) + 1,  
                                      end = start(eventGr$exon_flanking) - 1),
                                      strand = strand(eventGr$exon_short)
  )
  
  intron.long <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_long),
                                        ranges = IRanges::IRanges(
                                        start = end(eventGr$exon_long) + 1,  
                                        end = start(eventGr$exon_flanking) - 1),
                                        strand = strand(eventGr$exon_long)
  )
  
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.short <- GenomicRanges::findOverlaps(intron.short, 
                                                  gtf_exons, type = "any")
  mytx.ids.intron.short <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
  
  ovl.intron.long <- GenomicRanges::findOverlaps(intron.long, 
                                                 gtf_exons, type = "any")
  mytx.ids.intron.long <- gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
  
  #decide wich transcripts to plot in short and long tracks
  mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
  mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
  
  #remove transcripts with exons overlapping intronic regions
  mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
  mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
  
  tx_ids <- list()
  tx_ids[["txn_short"]] <- mytx.ids.short
  tx_ids[["txn_long"]] <- mytx.ids.long
  return(tx_ids)
  
}

#' @import GenomicRanges
mapTranscriptsA3SSevent <- function(eventGr, gtf_exons){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_short, 
                                        gtf_exons, type = "equal")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_long, 
                                        gtf_exons, type = "equal")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_flanking, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range for short event and long event
  intron.short <- GenomicRanges::GRanges(
                                      seqnames = seqnames(eventGr$exon_short),
                                      ranges = IRanges::IRanges(
                                      start = end(eventGr$exon_flanking) + 1,  
                                      end = start(eventGr$exon_short) - 1),
                                      strand = strand(eventGr$exon_short)
  )
  
  intron.long <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_long),
                                      ranges = IRanges::IRanges(
                                      start = end(eventGr$exon_flanking) + 1,  
                                      end = start(eventGr$exon_long) - 1),
                                      strand = strand(eventGr$exon_long)
  )  
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.short <- GenomicRanges::findOverlaps(intron.short, 
                                                  gtf_exons, type = "any")
  mytx.ids.intron.short <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
  
  ovl.intron.long <- GenomicRanges::findOverlaps(intron.long, 
                                                 gtf_exons, type = "any")
  mytx.ids.intron.long <- gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
  
  #decide wich transcripts to plot in short and long tracks
  mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
  mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
  
  #remove transcripts with exons overlapping intronic regions
  mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
  mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
  
  tx_ids <- list()
  tx_ids[["txn_short"]] <- mytx.ids.short
  tx_ids[["txn_long"]] <- mytx.ids.long
  return(tx_ids)
}
  