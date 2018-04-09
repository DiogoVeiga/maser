#' Boxplots of Percent spliced-in levels for gene events.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}
#' @param show_replicates logical, add data points for individual replicates     
#' @return a ggplot object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' hypoxia_mib2 <- geneEvents(hypoxia_filt, geneS = "MIB2")
#' plotGenePSI(hypoxia_mib2, type = "SE", show_replicates = T)
#' @export

plotGenePSI <- function(events, type, show_replicates = TRUE){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select AS events."))
  }
  
  PSI <- events[[paste0(type,"_","PSI")]]
  PSI_long <- reshape2::melt(PSI)
  colnames(PSI_long) <- c("ID", "Sample", "PSI")
  
  Condition <- rep("NA",nrow(PSI_long))
  idx.cond1 <- grep(paste0("^", events$conditions[1]), x = PSI_long$Sample,
                    perl = T )
  idx.cond2 <- grep(paste0("^", events$conditions[2]), x = PSI_long$Sample,
                    perl = T )
  Condition[idx.cond1] <- events$conditions[1]
  Condition[idx.cond2] <- events$conditions[2]
  
  PSI_long <- cbind(PSI_long, Condition)
  
  if (show_replicates){
    
    ggplot(PSI_long, aes(x = Condition, y = PSI, fill = Condition, color = Condition)) +
      #geom_boxplot() +
      geom_violin(trim = F, alpha = 0.6) +
      geom_jitter(position=position_jitter(0.05), size = 2) +
      theme_bw() +
      theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
            axis.text.y = element_text(size=12),
            axis.title.x = element_text(face="plain", colour="black", size=12),
            axis.title.y = element_text(face="plain", colour="black", size=12),
            legend.title=element_blank(),
            legend.text = element_text(face="plain", colour="black", size=12)) +
      ylab(paste(type, "PSI")) +
      #xlab("Sample") +
      scale_y_continuous(limits=c(-0.1, 1.05)) +
      scale_fill_manual(values = c("blue", "red") ) +
      scale_color_manual(values = c("blue", "red") ) +
      facet_grid(. ~ ID)

  } else{  
  
  ggplot(PSI_long, aes(x = Condition, y = PSI, fill = Condition, color = Condition)) +
    #geom_boxplot() +
    geom_violin(trim = F) +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    ylab(paste(type, "PSI")) +
    #xlab("Sample") +
    scale_y_continuous(limits=c(-0.1, 1.05)) +
    scale_fill_manual(values = c("blue", "red") ) +
    scale_color_manual(values = c("blue", "red") ) +  
    facet_grid(. ~ ID)
    
  }
  
}


plotTranscripts_old_txdb <- function(events, type, event_id, gtf_txdb, 
                            is_strict = FALSE){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select gene-specific AS events."))
  }
  
  # Genomic ranges of alternative splicing events
  grl <- events[[paste0(type,"_","gr")]]
  idx.event <- grep(as.numeric(event_id), grl[[1]]$ID)
  
  eventGr <- GRangesList()
  for (feature in names(grl)){
    eventGr[[paste0(feature)]] <- grl[[paste0(feature)]][idx.event]
  }
  
  eventTrack <- createAnnotationTrack(eventGr, type)
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicFeatures::transcriptsByOverlaps(gtf_txdb,eventGr[1])
  ovl.e2 <- GenomicFeatures::transcriptsByOverlaps(gtf_txdb,eventGr[2])
  ovl.e3 <- GenomicFeatures::transcriptsByOverlaps(gtf_txdb,eventGr[3])
  
  mytx.ids.e1 <- ovl.e1$tx_name
  mytx.ids.e2 <- ovl.e2$tx_name
  mytx.ids.e3 <- ovl.e3$tx_name
  
  if (is_strict){
    mytx.ids.3exons <- intersect(mytx.ids.e1, mytx.ids.e3) #has both flaking exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and the target exon
    mytx.ids.2exons <- intersect(mytx.ids.e1, mytx.ids.e3)
    
  }else {
    mytx.ids.3exons <- union(mytx.ids.e1, mytx.ids.e3) #has either flaking exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and the target exon  
    mytx.ids.2exons <- union(mytx.ids.e1, mytx.ids.e3)
  }
  
  mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.3exons)
  
  if(length(mytx.ids.3exons)==0){
    stop(cat("No transcripts matching the splicing event were found."))
  }
  
  # Inclusion track
  # Recover exons of transcripts for the inclusion track using transcript IDs
  # AnnotationDbi::keytypes(gtf_txdb)
  res <- AnnotationDbi::select(gtf_txdb, mytx.ids.3exons, 
                               columns = AnnotationDbi::columns(gtf_txdb),
                               keytype ="TXNAME")
  
  # Create data frame for inclusion track - follow the model from data(geneModels)
  res.df <- res[, c("EXONCHROM", "EXONSTART", "EXONEND", "EXONSTRAND", "EXONNAME", "TXNAME")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon","transcript")
  inclusionTrack <- GeneRegionTrack(range = res.df, name = "Inclusion", 
                            transcriptAnnotation = "transcript")
  
  # Skipping track
  # Recover exons of transcripts for the skipping track using transcript IDs
  res <- AnnotationDbi::select(gtf_txdb, mytx.ids.2exons, 
                               columns = AnnotationDbi::columns(gtf_txdb),
                               keytype ="TXNAME")
  
  # Create data frame for inclusion track - follow the model from data(geneModels)
  res.df <- res[, c("EXONCHROM", "EXONSTART", "EXONEND", "EXONSTRAND", "EXONNAME", "TXNAME")]
  colnames(res.df) <- c("chromosome","start","end","strand","exon","transcript")
  skippingTrack <- GeneRegionTrack(range = res.df, name = "Skipping", 
                                    transcriptAnnotation = "transcript")
  
  
  plotTracks(list(eventTrack, inclusionTrack, skippingTrack), 
             col.line = NULL, col = NULL,
             Inclusion = "darkred", Skipping = "darkblue")
  
  
}


#' Mapping and visualization of Ensembl transcripts affected by splicing.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}.
#' @param event_id numeric, event identifier.
#' @param gtf a \code{GRanges}, Ensembl or Gencode GTF using the hg38 build of the human genome.
#' @param zoom logical, zoom to the genomic coordinates of the splice event.
#' @param show_PSI logical, display the PSI track.     
#' @return a Gviz object.
#' @details This is a wrapper function for performing both mapping and visualization
#'  of Ensembl transcripts that are compatible with the splice event. This function
#'  calls \code{\link{mapTranscriptsToEvents}} for transcript mapping, which in 
#'  turn uses \code{\link[GenomicRanges:findOverlaps]{GenomicRanges}} for transcript
#'  overlapping. The \code{\link[Gviz:plotTracks]{GViz}} package is used for creating
#'  annotation tracks for genomic visualization of splice events. 
#'  
#'  Each type of splice event requires a specific overlapping rule (described below),
#'  and a customized \code{Gviz} plot is created for each splicing type.
#'   
#'   \describe{
#'     \item{\strong{Exon skipping}}{}
#'     \item{Inclusion track}{Transcript(s) overlapping the cassette exon, as well 
#'                            both flanking exons (i.e upstream and downstream exons).}
#'     \item{Skipping track}{Transcript(s) overlapping both flanking exons but not the
#'                            cassettte exon.}
#'   }
#'   
#'   \describe{
#'     \item{\strong{Intron retention}}{}
#'     \item{Retention track}{Transcript(s) overlapping exactly the retained intron.}
#'     \item{Skipping track}{Transcript(s) where intron is spliced out and overlapping
#'                              both flanking exons.}
#'   }
#'   
#'   \describe{
#'   \item{\strong{Mutually exclusive exons}}{}
#'     \item{Exon1 track}{Transcript(s) overlapping the first exon and both 
#'                       flanking exons.}
#'     \item{Exon2 track}{Transcript(s) overlapping the second exon and both
#'                        flanking exons.}
#'   }
#'   
#'   \describe{
#'     \item{\strong{Alternative 3' and 5' splice sites}}{}
#'     \item{Short exon track}{Transcript(s) overlapping both short and 
#'                        downstream exons.}
#'     \item{Long exon track}{Transcript(s) overlapping both long and downstream 
#'                           exons.}
#'   }
#'   
#' @examples
#' ## Create the maser object
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' 
#' ## Retrieve Ensembl GTF annotation
#' ah <- AnnotationHub::AnnotationHub() 
#' ens_gtf <- qhs[["AH51014"]] #Homo_sapiens.GRCh38.85.gtf 
#' 
#' ## Retrieve gene specific splice events
#' srsf6_events <- geneEvents(hypoxia_filt, geneS = "SRSF6")
#' 
#' ## Plot exon skipping event
#' plotTranscripts(srsf6_events, type = "SE", event_id = 33209, gtf = ens_gtf)
#' 
#' @seealso \code{\link{mapTranscriptsToEvents}}
#' @export

plotTranscripts <- function(events, type, event_id, gtf, 
                            zoom = FALSE, show_PSI = TRUE){
  
  is_strict = TRUE #affects exon skipping, MXE, RI
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  if (!class(gtf) == "GRanges"){
    stop(cat("\"gtf\" should be a GRanges class."))
  }
  gtf <- GenomeInfoDb::keepSeqlevels(gtf, c(paste0(seq(1:22)), "X", "Y"),
                                     pruning.mode = "coarse")
  seqlevels(gtf, pruning.mode="coarse") <- c(paste0("chr", seq(1:22)),
                                             "chrX", "chrY")

  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select gene-specific AS events."))
  }
  
  # Genomic ranges of alternative splicing events
  grl <- events[[paste0(type,"_","gr")]]
  idx.event <- grep(as.numeric(event_id), grl[[1]]$ID)
  
  eventGr <- GRangesList()
  for (feature in names(grl)){
    eventGr[[paste0(feature)]] <- grl[[paste0(feature)]][idx.event]
  }
  eventTrack <- createAnnotationTrack_event(eventGr, type)
  
  gtf_exons <- gtf[gtf$type=="exon",]
  txnTracks <- createAnnotationTrack_transcripts(eventGr, gtf_exons,
                                            type, is_strict)
  
  if (show_PSI){
    PSI <- events[[paste0(type,"_","PSI")]]
    PSI_event <- PSI[idx.event, , drop = FALSE]
    groups <- factor(c(rep(events$conditions[1], events$n_cond1),
                       rep(events$conditions[2], events$n_cond2)),
                     levels = events$conditions)
    psiTrack <- createPSITrack_event(eventGr, PSI_event, groups, type, zoom)  
    trackList <- list(psiTrack, eventTrack, txnTracks$inclusionTrack, 
                      txnTracks$skippingTrack)  
  }else{
    trackList <- list(eventTrack, txnTracks$inclusionTrack, 
                      txnTracks$skippingTrack)  
  }
  
  
  if (zoom){
    Gviz::plotTracks(trackList, 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple",
                     from = start(range(unlist(eventGr))) - 500,
                     to = end(range(unlist(eventGr))) + 500)  
  }else {
    Gviz::plotTracks(trackList,
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple")  
  }

}


#' Mapping and visualization of UniprotKB protein features affected by splicing.
#' 
#' @param events a maser object.
#' @param type character indicating splice type. Possible values are
#'    \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}.
#' @param event_id numeric, event identifier.
#' @param gtf a \code{GRanges}, Ensembl or Gencode GTF using the hg38 build of the human genome.
#' @param features a character vector indicating valid UniprotKB features.
#' @param zoom logical, zoom to the genomic coordinates of the splice event.
#' @param show_transcripts logical, display transcripts track.
#' @return a Gviz object.
#' @details This is a wrapper function for performing both mapping and visualization
#'  of protein features affected by the splice event. This function
#'  calls \code{\link{mapProteinFeaturesToEvents}} for mapping of protein features to 
#'  splice events. 
#'  
#' The \code{\link[Gviz:plotTracks]{GViz}} package is used for creating
#'  annotation tracks for genomic visualization. 
#' 
#' Multiple protein annotation tracks can be created using the \code{features}
#' argument.
#' 
#' @examples
#' ## Create the maser object
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' 
#' ## Retrieve Ensembl GTF annotation
#' ah <- AnnotationHub::AnnotationHub() 
#' ens_gtf <- qhs[["AH51014"]] #Homo_sapiens.GRCh38.85.gtf 
#' 
#' ## Retrieve gene specific splice events
#' srsf6_events <- geneEvents(hypoxia_filt, geneS = "SRSF6")
#' 
#' ## Map splicing events to transcripts
#' srsf6_mapped <- mapTranscriptsToEvents(srsf6_events, ens_gtf)
#' 
#' ## Annotate splice events with protein features
#' categories <- c("Domain_and_Sites", "Topology")
#' srsf6_annot <- mapProteinFeaturesToEvents(srsf6_mapped, categories, by="category")
#' head(annot(srsf6_annot, "SE"))
#' 
#' ## Plot splice event, transcripts and protein features
#' plotUniprotKBFeatures(srsf6_annot, "SE", event_id = 33209, gtf = ens_gtf, 
#'   features = c("domain", "chain", "mod-res"), show_transcripts = T)
#' 
#' @seealso \code{\link{mapProteinFeaturesToEvents}}
#' @export


plotUniprotKBFeatures <- function(events, type, event_id, gtf,
                                  features, zoom = FALSE,
                                  show_transcripts = FALSE){
  
  is_strict = FALSE
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  if (!class(gtf) == "GRanges"){
    stop(cat("\"gtf\" should be a GRanges class."))
  }
  gtf <- GenomeInfoDb::keepSeqlevels(gtf, c(paste0(seq(1:22)), "X", "Y"),
                                     pruning.mode = "coarse")
  seqlevels(gtf, pruning.mode="coarse") <- c(paste0("chr", seq(1:22)),
                                             "chrX", "chrY")
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select gene-specific AS events."))
  }
  

  # Genomic ranges of alternative splicing events
  grl <- events[[paste0(type,"_","gr")]]
  idx.event <- grep(as.numeric(event_id), annot$ID)
  
  eventGr <- GRangesList()
  for (feature in names(grl)){
    eventGr[[paste0(feature)]] <- grl[[paste0(feature)]][idx.event]
  }
  
  eventTrack <- createAnnotationTrack_event(eventGr, type)
  
  gtf_exons <- gtf[gtf$type=="exon",]
  idx.cols <- grep("^list_ptn_", colnames(annot))
  protein_ids <- unique(c(annot[idx.event, idx.cols[1]], annot[idx.event, idx.cols[2]]))
  uniprotTracks <- createUniprotKBtracks(eventGr, features, protein_ids)
  
  if (show_transcripts){
    txnTracks <- createAnnotationTrack_transcripts(eventGr, gtf_exons,
                                                   type, is_strict)
    trackList <- c(list(eventTrack, txnTracks$inclusionTrack, 
                      txnTracks$skippingTrack),
                   uniprotTracks)
  }else {
    trackList <- c(list(eventTrack), uniprotTracks)
  }
  
  if (zoom){
    Gviz::plotTracks(trackList, 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple",
                     from = start(range(unlist(eventGr))) - 500,
                     to = end(range(unlist(eventGr))) + 500)  
  }else {
    Gviz::plotTracks(trackList, 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple")  
  }
  
  
}

# deprecated - uses UCSC rtracklayer but too slow
plotUniprotUCSCFeatures <- function(events, type, event_id, gtf, 
                            is_strict = FALSE, zoom = FALSE, 
                            show_transcripts = FALSE){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select gene-specific AS events."))
  }
  
  # Genomic ranges of alternative splicing events
  grl <- events[[paste0(type,"_","gr")]]
  idx.event <- grep(as.numeric(event_id), grl[[1]]$ID)
  
  eventGr <- GRangesList()
  for (feature in names(grl)){
    eventGr[[paste0(feature)]] <- grl[[paste0(feature)]][idx.event]
  }
  
  eventTrack <- createAnnotationTrack_event(eventGr, type)
  
  gtf_exons <- gtf[gtf$type=="exon",]
  uniprotTracks <- createUniprotUCSCtrack_localization(eventGr, "hg38")
  
  if (show_transcripts){
    txnTracks <- createAnnotationTrack_transcripts(eventGr, gtf_exons,
                                                   type, is_strict)
    trackList <- list(eventTrack, txnTracks$inclusionTrack, 
                      txnTracks$skippingTrack,
                      uniprotTracks$unipLocExtra, 
                      uniprotTracks$unipLocTransMemb,
                      uniprotTracks$unipLocCytopl)
  }else {
    trackList <- list(eventTrack, uniprotTracks$unipLocExtra, 
                      uniprotTracks$unipLocTransMemb,
                      uniprotTracks$unipLocCytopl)
  }
  
  if (zoom){
    Gviz::plotTracks(trackList, 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple",
                     from = start(range(unlist(eventGr))) - 500,
                     to = end(range(unlist(eventGr))) + 500)  
  }else {
    Gviz::plotTracks(trackList, 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple")  
  }
  
}


