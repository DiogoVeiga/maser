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


createAnnotationTrackSE_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  transcript_id <- NULL
  tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
  
  # Inclusion track
  # Recover exons of transcripts for the inclusion track using transcript
  #IDs
  # AnnotationDbi::keytypes(gtf_txdb)
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_3exons)
  
  # Create data frame for inclusion track - follow the model 
  #from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon", 
                        "transcript")
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- "Inclusion"
    inclusionTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                    name = "Inclusion", 
                                    transcriptAnnotation = "transcript")  
  }else {
    inclusionTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                      name = "Inclusion", 
                                      transcriptAnnotation = "transcript")  
  }
  
  # Skipping track
  # Recover exons of transcripts for the skipping track using 
  # transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_2exons)
  
  # Create data frame for inclusion track - follow the model 
  # from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon", 
                        "transcript")
  
  if (nrow(res.df) > 0){
    res.df$feature <- "Skipping"
    skippingTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                      name = "Skipping", 
                                      transcriptAnnotation = "transcript")  
  }else {
    skippingTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                      name = "Skipping", 
                                      transcriptAnnotation = "transcript")  
  }
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackRI_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  transcript_id <- NULL
  tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
  
  # Retention track
  # Recover exons of transcripts for the retention track using 
  # transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_retention)
  
  # Create data frame for inclusion track - follow the model from 
  # data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon",
                        "transcript")
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- "Retention"
    retention_Track <- Gviz::GeneRegionTrack(range = res.df,
                                      name = "Retention", 
                                      transcriptAnnotation = "transcript")  
  }else {
    retention_Track <- Gviz::GeneRegionTrack(range = GRanges(),
                                      name = "Retention", 
                                      transcriptAnnotation = "transcript")  
  }
  
  # Non-retention track
  # Recover exons of transcripts for the non-retention track using
  # transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_nonRetention)
  
  # Create data frame for inclusion track - follow the model
  # from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon",
                        "transcript")
  
  if (nrow(res.df) > 0){
    res.df$feature <- "Non_retention"
    nonRetention_Track <- Gviz::GeneRegionTrack(range = res.df, 
                                      name = "Non-retention", 
                                      transcriptAnnotation = "transcript",
                                      feature = "Non_Retention")  
  }else {
    nonRetention_Track <- Gviz::GeneRegionTrack(range = GRanges(),
                                      name = "Non-retention", 
                                      transcriptAnnotation = "transcript",
                                      feature = "Non_Retention")  
  }
  
  txn_tracks <- list("inclusionTrack" = retention_Track,
                     "skippingTrack" = nonRetention_Track)
  return(txn_tracks)
  
}

createAnnotationTrackMXE_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  transcript_id <- NULL
  tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)

  # MXE Exon 1 track
  # Recover exons of transcripts for the MXE Exon 1 track using 
  # transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_mxe_exon1)
  
  # Create data frame for inclusion track - follow the model
  # from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon", 
                        "transcript")
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- "MXE_Exon1"
    inclusionTrack <- Gviz::GeneRegionTrack(range = res.df,
                                    name = "MXE Exon 1", 
                                    transcriptAnnotation = "transcript")  
  }else {
    inclusionTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                    name = "MXE Exon 1", 
                                    transcriptAnnotation = "transcript")  
  }
  
  # MXE Exon 2 track
  # Recover exons of transcripts for the MXE Exon 2 track using
  #transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_mxe_exon2)
  
  # Create data frame for inclusion track - follow the model from
  # data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon", 
                        "transcript")
  
  if (nrow(res.df) > 0){
    res.df$feature <- "MXE_Exon2"
    skippingTrack <- Gviz::GeneRegionTrack(range = res.df,
                                    name = "MXE Exon 2", 
                                    transcriptAnnotation = "transcript")  
  }else {
    skippingTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                    name = "MXE Exon 2", 
                                    transcriptAnnotation = "transcript")  
  }
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackA5SS_transcripts <- function(eventGr, gtf_exons){
  
  transcript_id <- NULL
  #reverse strand becomes A3SS
  if (as.character(strand(eventGr[1])) == "+"){
    tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
  }else{
    tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
  }
  
  # Short track
  # Recover exons of transcripts for the short track using transcript IDs
  # AnnotationDbi::keytypes(gtf_txdb)
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_short)
  
  # Create data frame for inclusion track - follow the model from
  # data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon",
                        "transcript")
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- "A5SS_Short"
    inclusionTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                      name = "A5SS Short", 
                                      transcriptAnnotation = "transcript")  
  }else {
    inclusionTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                    name = "A5SS Short", 
                                    transcriptAnnotation = "transcript")  
  }
  
  # Long track
  # Recover exons of transcripts for the long track using transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_long)
  
  # Create data frame for inclusion track - follow the model
  # from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon",
                        "transcript")
  
  if (nrow(res.df) > 0){
    res.df$feature <- "A5SS_Long"
    skippingTrack <- Gviz::GeneRegionTrack(range = res.df,
                                      name = "A5SS Long", 
                                      transcriptAnnotation = "transcript")  
  }else {
    skippingTrack <- Gviz::GeneRegionTrack(range = GRanges(), 
                                      name = "A5SS Long", 
                                      transcriptAnnotation = "transcript")  
  }
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackA3SS_transcripts <- function(eventGr, gtf_exons){
  
  transcript_id <- NULL
  #reverse strand becomes A5SS
  if (as.character(strand(eventGr[1])) == "+"){
    tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
  }else{
    tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
  }
  
  # Short track
  # Recover exons of transcripts for the short track using transcript IDs
  # AnnotationDbi::keytypes(gtf_txdb)
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_short)
  
  # Create data frame for inclusion track - follow the model 
  # from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon",
                        "transcript")
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- "A3SS_Short"
    inclusionTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                      name = "A3SS Short", 
                                      transcriptAnnotation = "transcript")  
  }else {
    inclusionTrack <- Gviz::GeneRegionTrack(range = GRanges(), 
                                      name = "A3SS Short", 
                                      transcriptAnnotation = "transcript")  
  }
  
  # Long track
  # Recover exons of transcripts for the long track using transcript IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% tx_ids$txn_long)
  
  # Create data frame for inclusion track - follow the model 
  # from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", 
                    "exon_id", "transcript_name")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon",
                        "transcript")
  
  if (nrow(res.df) > 0){
    res.df$feature <- "A3SS_Long"
    skippingTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                    name = "A3SS Long", 
                                    transcriptAnnotation = "transcript")  
  }else {
    skippingTrack <- Gviz::GeneRegionTrack(range = GRanges(), 
                                    name = "A3SS Long", 
                                    transcriptAnnotation = "transcript")  
  }
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackSE_event <- function(eventGr){
  
  transcript_id <- NULL
  trackGr <- c(unlist(eventGr), unlist(eventGr[2:3]))
  trackGr$group <- rep(c("Inclusion", "Skipping"), c(3, 2))
  trackGr$type <- rep("Exon skipping", 5)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", 
                                       shape = "box",
                                       stacking = "squish", 
                                       id = "Exon skipping")
  Gviz::feature(event_track) <- rep(c("Inclusion", "Skipping"), c(3, 2))

  return(event_track)
  
}

createAnnotationTrackRI_event <- function(eventGr){
  
  transcript_id <- NULL
  trackGr <- c(eventGr$exon_ir, eventGr$exon_upstream, 
               eventGr$exon_downstream)
  trackGr$group <- rep(c("Retention", "Non-retention"), c(1, 2))
  trackGr$type <- rep("Intron retention", 3)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                          groupAnnotation = "group", shape = "box",
                          stacking = "squish", id = "Intron retention")
  
  Gviz::feature(event_track) <- rep(c("Retention", "Non_Retention"),
                                    c(1, 2))
  
  
  return(event_track)
  
}

createAnnotationTrackMXE_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_upstream, eventGr$exon_1, 
               eventGr$exon_downstream, eventGr$exon_upstream,
               eventGr$exon_2, eventGr$exon_downstream)
  
  trackGr$group <- rep(c("MXE_Exon1", "MXE_Exon2"), c(3, 3))
  trackGr$type <- rep("Mutually Exclusive Exons", 3)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                      groupAnnotation = "group", shape = "box",
                      stacking = "squish", id = "Mutually Exclusive Exons")
  
  Gviz::feature(event_track) <- rep(c("MXE_Exon1", "MXE_Exon2"), c(3, 3))
  
  
  return(event_track)
  
}

createAnnotationTrackA5SS_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_short, eventGr$exon_flanking,
               eventGr$exon_long, eventGr$exon_flanking)
  trackGr$group <- rep(c("A5SS Short", "A5SS Long"), c(2, 2))
  trackGr$type <- rep("A5SS", 4)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                groupAnnotation = "group", shape = "box",
                                stacking = "squish", id = "A5SS")
  Gviz::feature(event_track) <- rep(c("A5SS_Short", "A5SS_Long"), c(2, 2))
  
  
  return(event_track)
  
}

createAnnotationTrackA3SS_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_flanking, eventGr$exon_short,
               eventGr$exon_flanking, eventGr$exon_long)
  trackGr$group <- rep(c("A3SS Short", "A3SS Long"), c(2, 2))
  trackGr$type <- rep("A3SS", 4)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                groupAnnotation = "group", shape = "box",
                                stacking = "squish", id = "A3SS")
  Gviz::feature(event_track) <- rep(c("A3SS_Short", "A3SS_Long"), c(2, 2))
  
  
  return(event_track)
  
}

#deprecated rtracklayer too slow - modify this function 
# to use createGRangesUniprotKB
createUniprotUCSCtrack_localization <- function(eventGr, genome){
  
  uniprotTracks <- list()
  
  session <- rtracklayer::browserSession("UCSC")
  rtracklayer::genome(session) <- genome
  
  #tableNames(ucscTableQuery(mySession, track="uniprot"))
  # Tables to query in the Uniprot track
  # Creates a track for each table related to cell localization
  tables_uniprot <- c("unipLocExtra", "unipLocTransMemb", 
                      "unipLocCytopl")
  names_uniprot <- c("Extra", "TransMemb", 
                      "Cytop")
  
  #Define region around splicing event
  region <- range(unlist(eventGr))
  start(region) <- start(region) - 10
  end(region) <- end(region) + 10
  GenomeInfoDb::genome(region) <- "hg38"
  
  for(i in 1:length(tables_uniprot)){
    
    query <- rtracklayer::ucscTableQuery(session, track = "uniprot", 
                            table = tables_uniprot[i],
                            range = region)
    query_gr <- rtracklayer::track(query)
    query_table <- rtracklayer::getTable(query)
    
    if (length(query_gr) > 0 ){
      
      track <- Gviz::AnnotationTrack(range = query_gr, genome = genome, 
                                     name = names_uniprot[i], 
                                     id = query_gr$name, 
                                     showFeatureId = TRUE,
                                     fill = query_gr$itemRgb,
                                     shape = "arrow")  
        
    }else {
      track <- Gviz::AnnotationTrack(range = GRanges(), 
                                     name = names_uniprot[i])
    }
    
    uniprotTracks[[tables_uniprot[i]]] <- track  
  }
  
  return(uniprotTracks)  

}

createUniprotKBtracks <- function(eventGr, features, protein_ids){
  
  uniprotTracks <- list()
  
  for(i in 1:length(features)){
    
    feature_gr <- createGRangesUniprotKBtrack(features[i])
    ovl_gr <- overlappingFeatures(feature_gr, eventGr)
    ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ] 
    
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
                                     shape = "box")  
      
    }else {
      track <- Gviz::AnnotationTrack(range = GRanges(), name = features[i])
    }
    
    uniprotTracks[[features[i]]] <- track  
  }
  
  return(uniprotTracks)  
  
}



createPSITrack_event <- function(eventGr, PSI_event, groups, type, zoom){
  
  if(type == "A3SS") {
    psi_track <- createPSITrackA3SS_event(eventGr, PSI_event, groups, zoom)
  }
  
  if(type == "A5SS") {
    psi_track <- createPSITrackA5SS_event(eventGr, PSI_event, groups, zoom)
  }
  
  if (type == "SE"){
    psi_track <- createPSITrackSE_event(eventGr, PSI_event, groups, zoom)
  }
  
  if (type == "RI"){
    psi_track <- createPSITrackRI_event(eventGr, PSI_event, groups, zoom)
  }
  
  if (type == "MXE"){
    psi_track <- createPSITrackMXE_event(eventGr, PSI_event, groups, zoom)
  }
  
  return(psi_track)
  
}

createPSITrackSE_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- eventGr$exon_target  
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
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("blue", "red"),
                               col = c("blue", "red"))
  return(psi_track)
  
}

createPSITrackRI_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- eventGr$exon_ir  
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
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("blue", "red"),
                               col = c("blue", "red"))
  return(psi_track)
  
  
}

createPSITrackMXE_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- c(eventGr$exon_1, eventGr$exon_2)
    start(trackGr) <- start(trackGr) - 50
    end(trackGr) <- end(trackGr) + 50
  }else{
    #create space for boxplot plotting
    trackGr <- c( range(unlist(eventGr$exon_upstream,eventGr$exon_1)),
                  range(unlist(eventGr$exon_2,eventGr$exon_downstream))
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

createPSITrackA5SS_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- eventGr$exon_long  
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
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("blue", "red"),
                               col = c("blue", "red"))
  return(psi_track)
  
}

createPSITrackA3SS_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- eventGr$exon_long  
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
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("blue", "red"),
                               col = c("blue", "red"))
  return(psi_track)
  
}