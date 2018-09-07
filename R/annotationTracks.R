createAnnotationTrack_event <- function(eventGr, type){
  
  if(type == "A3SS") {
    event_track <- createAnnotationTrackA3SS_event(eventGr)
  }
  
  if(type == "A5SS") {
    event_track <- createAnnotationTrackA5SS_event(eventGr)
  }
  
  if (type == "SE"){
    event_track <- createAnnotationTrackSE_event(eventGr)
  }
  
  if (type == "RI"){
    event_track <- createAnnotationTrackRI_event(eventGr)
  }
  
  if (type == "MXE"){
    event_track <- createAnnotationTrackMXE_event(eventGr)
  }
  
  return(event_track)
  
}

createAnnotationTrack_transcripts <- function(eventGr, gtf_exons, 
                                              type, is_strict){
  
  if(type == "A3SS") {
    txn_tracks <- createAnnotationTrackA3SS_transcripts(eventGr, gtf_exons)
  }
  
  if(type == "A5SS") {
    txn_tracks <- createAnnotationTrackA5SS_transcripts(eventGr, gtf_exons)
  }
  
  if (type == "SE"){
    txn_tracks <- createAnnotationTrackSE_transcripts(eventGr, gtf_exons,
                                                      is_strict)
  }
  
  if (type == "RI"){
    txn_tracks <- createAnnotationTrackRI_transcripts(eventGr, gtf_exons,
                                                      is_strict)
  }
  
  if (type == "MXE"){
    txn_tracks <- createAnnotationTrackMXE_transcripts(eventGr, gtf_exons,
                                                      is_strict)
  }
  
  return(txn_tracks)
  
}
#' @importFrom dplyr filter
createExonTable <- function(gtf_exons, ids){
  
  transcript_id <- NULL
  # Recover exons of transcripts for the inclusion track using transcript
  #IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% ids)
  
  # Create data frame for transcript track - follow the model 
  #from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon", 
                        "transcript")
  
  return(res.df)
}

#' @importFrom Gviz GeneRegionTrack
#' @importFrom GenomicRanges GRanges
createTxnTrack <- function(res.df, trackLabel, featureName){
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- featureName
    txnTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                            name = trackLabel, 
                                            transcriptAnnotation = "transcript",
                                            col = NULL,
                                            col.line = NULL)  
  }else {
    txnTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                            name = trackLabel, 
                                            transcriptAnnotation = "transcript",
                                            col = NULL,
                                            col.line = NULL)  
  }
  return(txnTrack)
  
}

#' @importFrom Gviz GeneRegionTrack
createAnnotationTrackSE_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
  
  # Inclusion track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_3exons)
  inclusionTrack <- createTxnTrack(res.df, "Inclusion", "Inclusion")
  
  # Skipping track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_2exons)
  skippingTrack <- createTxnTrack(res.df, "Skipping", "Skipping")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackRI_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
  
  # Retention track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_retention)
  retention_Track <- createTxnTrack(res.df, "Retention", "Retention")
  
  # Non-retention track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_nonRetention)
  nonRetention_Track <- createTxnTrack(res.df, "Non-retention", "Non_Retention")
  
  txn_tracks <- list("inclusionTrack" = retention_Track,
                     "skippingTrack" = nonRetention_Track)
  return(txn_tracks)
  
}

createAnnotationTrackMXE_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)

  # MXE Exon 1 track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_mxe_exon1)
  inclusionTrack <- createTxnTrack(res.df, "MXE Exon 1", "MXE_Exon1")
  
  # MXE Exon 2 track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_mxe_exon2)
  skippingTrack <- createTxnTrack(res.df, "MXE Exon 2", "MXE_Exon2")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackA5SS_transcripts <- function(eventGr, gtf_exons){
  
  #reverse strand becomes A3SS
  if (as.character(strand(eventGr[1])) == "+"){
    tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
  }else{
    tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
  }
  
  # Short exon track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_short)
  inclusionTrack <- createTxnTrack(res.df, "A5SS Short", "A5SS_Short")
  
  # Long track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_long)
  skippingTrack <- createTxnTrack(res.df, "A5SS Long", "A5SS_Long")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackA3SS_transcripts <- function(eventGr, gtf_exons){
  
  #reverse strand becomes A5SS
  if (as.character(strand(eventGr[1])) == "+"){
    tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
  }else{
    tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
  }
  
  # Short exon track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_short)
  inclusionTrack <- createTxnTrack(res.df, "A3SS Short", "A3SS_Short")
  
  # Long exon track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_long)
  skippingTrack <- createTxnTrack(res.df, "A3SS Long", "A3SS_Long")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackSE_event <- function(eventGr){
  
  transcript_id <- NULL
  trackGr <- c(unlist(eventGr), unlist(eventGr[2:3]))
  trackGr$group <- rep(c("Inclusion", "Skipping"), c(3, 2))
  trackGr$type <- rep("Exon skipping", 5)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", 
                                       shape = "box",
                                       stacking = "squish", 
                                       id = "Exon skipping",
                                       col = NULL,
                                       col.line = NULL)
  Gviz::feature(event_track) <- rep(c("Inclusion", "Skipping"), c(3, 2))

  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackRI_event <- function(eventGr){
  
  transcript_id <- NULL
  trackGr <- c(eventGr$exon_ir, eventGr$exon_upstream, 
               eventGr$exon_downstream)
  trackGr$group <- rep(c("Retention", "Non-retention"), c(1, 2))
  trackGr$type <- rep("Intron retention", 3)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                          groupAnnotation = "group", 
                          shape = "box",
                          stacking = "squish", 
                          id = "Intron retention",
                          col = NULL,
                          col.line = NULL)
  
  Gviz::feature(event_track) <- rep(c("Retention", "Non_Retention"),
                                    c(1, 2))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
#' 
createAnnotationTrackMXE_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_upstream, eventGr$exon_1, 
               eventGr$exon_downstream, eventGr$exon_upstream,
               eventGr$exon_2, eventGr$exon_downstream)
  
  trackGr$group <- rep(c("MXE_Exon1", "MXE_Exon2"), c(3, 3))
  trackGr$type <- rep("Mutually Exclusive Exons", 3)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                      groupAnnotation = "group", shape = "box",
                      stacking = "squish", 
                      id = "Mutually Exclusive Exons",
                      col = NULL,
                      col.line = NULL)
  
  Gviz::feature(event_track) <- rep(c("MXE_Exon1", "MXE_Exon2"), c(3, 3))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackA5SS_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_short, eventGr$exon_flanking,
               eventGr$exon_long, eventGr$exon_flanking)
  trackGr$group <- rep(c("A5SS Short", "A5SS Long"), c(2, 2))
  trackGr$type <- rep("A5SS", 4)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                groupAnnotation = "group", 
                                shape = "box",
                                stacking = "squish", id = "A5SS",
                                col = NULL,
                                col.line = NULL)
  Gviz::feature(event_track) <- rep(c("A5SS_Short", "A5SS_Long"), c(2, 2))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackA3SS_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_flanking, eventGr$exon_short,
               eventGr$exon_flanking, eventGr$exon_long)
  trackGr$group <- rep(c("A3SS Short", "A3SS Long"), c(2, 2))
  trackGr$type <- rep("A3SS", 4)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                groupAnnotation = "group", 
                                shape = "box",
                                stacking = "squish", 
                                id = "A3SS",
                                col = NULL,
                                col.line = NULL)
  Gviz::feature(event_track) <- rep(c("A3SS_Short", "A3SS_Long"), c(2, 2))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom GenomicRanges GRanges
#' @importFrom parallel mclapply
createUniprotKBtracks <- function(eventGr, features, protein_ids, ncores = 1){
  
  options(ucscChromosomeNames=FALSE)
  
  if(.Platform$OS.type == "windows"){
    ncores = 1
  }
  
  uniprotTracks <- mclapply(seq_along(features), function(i){
    
    feature_gr <- createGRangesUniprotKBtrack(features[i])
    ovl_gr <- overlappingFeatures(feature_gr, eventGr)
    #ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ] 
    ovl_gr_filt <- ovl_gr
    
    uniq_features <- match(unique(as.vector(ovl_gr_filt$Name)), 
                           as.vector(ovl_gr_filt$Name))
    
    ovl_gr_filt_uniq <- ovl_gr_filt[uniq_features, ]
    
    if (length(ovl_gr_filt_uniq) > 0 ){
      
      track <- Gviz::AnnotationTrack(range = ovl_gr_filt_uniq, 
                                     name = features[i], 
                                     id = ovl_gr_filt_uniq$Name, 
                                     showFeatureId = TRUE,
                                     fontcolor.feature = "darkblue",
                                     #fill = "aliceblue", 
                                     shape = "box",
                                     col = NULL,
                                     col.line = NULL)  
      
    }else {
      track <- Gviz::AnnotationTrack(range = GRanges(), name = features[i])
    }
    
    return(track)
    
  }, mc.cores = ncores)

  return(uniprotTracks)  
  
}



createPSITrack_event <- function(eventGr, PSI_event, groups, type, zoom){
  
  if(type == "A3SS" || type == "A5SS" ) {
    psi_track <- createPSIDataTrack(eventGr, PSI_event, groups, zoom, 
                                    eventGr$exon_long)
  }
  
  if (type == "SE"){
    psi_track <- createPSIDataTrack(eventGr, PSI_event, groups, zoom, 
                                    eventGr$exon_target)
  }
  
  if (type == "RI"){
    psi_track <- createPSIDataTrack(eventGr, PSI_event, groups, zoom, 
                                    eventGr$exon_ir)
  }
  
  if (type == "MXE"){
    psi_track <- createPSITrackMXE_event(eventGr, PSI_event, groups, zoom)
  }
  
  return(psi_track)
  
}

#' @importFrom Gviz DataTrack
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges values
createPSIDataTrack <- function(eventGr, PSI_event, groups, zoom, exonGr){

  if(zoom){
    trackGr <- exonGr  
  }else{
    #create space for boxplot plotting
    trackGr <- range(unlist(eventGr))
    start(trackGr) <- start(trackGr) - 200
    end(trackGr) <- end(trackGr) + 200
  }
  
  values(trackGr) <- PSI_event 
  psi_track <- Gviz::DataTrack(trackGr, 
                               name = "PSI", 
                               groups = groups,
                               feature = groups,
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("#0000ff", "#ff0000"),
                               col = c("#0000ff", "#ff0000"),
                               col.line = c("#0000ff", "#ff0000"))
  
  return(psi_track)
  
}


#' @importFrom Gviz DataTrack
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges values
createPSITrackMXE_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- c(eventGr$exon_1, eventGr$exon_2)
    start(trackGr) <- start(trackGr) - 50
    end(trackGr) <- end(trackGr) + 50
  }else{
    #create space for boxplot plotting
    trackGr <- c( range(c(range(eventGr$exon_upstream), 
                          range(eventGr$exon_1))),
                  
                  range(c(range(eventGr$exon_2), 
                          range(eventGr$exon_downstream)))
    )
    start(trackGr) <- start(trackGr) - 200
    end(trackGr) <- end(trackGr) + 200
  }
  #rMATS PSI levels are for Exon 1
  data <- rbind(PSI_event, 1-PSI_event )
  values(trackGr) <- data
  psi_track <- Gviz::DataTrack(trackGr, 
                               name = "PSI", 
                               groups = groups,
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("blue", "red"),
                               col = c("blue", "red"))
  return(psi_track)
  
  
}
