
PSI_geneEvents <- function(events, type, show_replicates = TRUE){
  
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


plotTranscripts_old_txdb <- function(gene_events, type, event_id, gtf_txdb, 
                            is_strict = FALSE){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- gene_events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select gene-specific AS events."))
  }
  
  # Genomic ranges of alternative splicing events
  grl <- gene_events[[paste0(type,"_","gr")]]
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


plotTranscripts <- function(gene_events, type, event_id, gtf, 
                            is_strict = FALSE, zoom = FALSE){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- gene_events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select gene-specific AS events."))
  }
  
  # Genomic ranges of alternative splicing events
  grl <- gene_events[[paste0(type,"_","gr")]]
  idx.event <- grep(as.numeric(event_id), grl[[1]]$ID)
  
  eventGr <- GRangesList()
  for (feature in names(grl)){
    eventGr[[paste0(feature)]] <- grl[[paste0(feature)]][idx.event]
  }
  
  eventTrack <- createAnnotationTrack_event(eventGr, type)
  
  gtf_exons <- gtf[gtf$type=="exon",]
  txnTracks <- createAnnotationTrack_transcripts(eventGr, gtf_exons,
                                            type, is_strict)
  if (zoom){
    Gviz::plotTracks(list(eventTrack, txnTracks$inclusionTrack, 
                          txnTracks$skippingTrack), 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple",
                     from = start(range(unlist(eventGr))) - 500,
                     to = end(range(unlist(eventGr))) + 500)  
  }else {
    Gviz::plotTracks(list(eventTrack, txnTracks$inclusionTrack, 
                          txnTracks$skippingTrack), 
                     col.line = NULL, col = NULL,
                     Inclusion = "orange", Skipping = "purple",
                     Retention = "orange", Non_Retention = "purple",
                     MXE_Exon1 = "orange", MXE_Exon2 = "purple",
                     A5SS_Short = "orange", A5SS_Long = "purple",
                     A3SS_Short = "orange", A3SS_Long = "purple")  
  }

}

