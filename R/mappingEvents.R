#internal function? called inside mapTranscriptsEvents
mapTranscriptsToProtein <- function(){
  
}


# User accessible function
# Will add transcript IDs to all events
mapTranscriptsToEvents <- function(events, gtf, is_strict = TRUE){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  gtf_exons <- gtf[gtf$type=="exon",]
  
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
    
    list_txn_a <- c()
    list_txn_b <- c()
    
    for (i in 1:nrow(annot)) {
      
      # Genomic ranges of alternative splicing events
      eventGr <- GRangesList()
      for (feature in names(grl)){
        eventGr[[paste0(feature)]] <- grl[[paste0(feature)]][i]
      }
      
      if(type == "SE") {
        tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
        list_txn_a <- c(list_txn_a, paste(tx_ids$txn_3exons, collapse = ",") )
        list_txn_b <- c(list_txn_b, paste(tx_ids$txn_2exons, collapse = ",") )
      }
      
      if(type == "MXE") {
        tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)
        list_txn_a <- c(list_txn_a, paste(tx_ids$txn_mxe_exon1, collapse = ",") )
        list_txn_b <- c(list_txn_b, paste(tx_ids$txn_mxe_exon2, collapse = ",") )
      }
      
      if(type == "RI") {
        tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
        list_txn_a <- c(list_txn_a, paste(tx_ids$txn_nonRetention, collapse = ",") )
        list_txn_b <- c(list_txn_b, paste(tx_ids$txn_retention, collapse = ",") )
      }
      
      if(type == "A5SS") {
        
        #reverse strand becomes A3SS
        if (as.character(strand(eventGr[1])) == "+"){
          tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
        }else{
          tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
        }
        
        list_txn_a <- c(list_txn_a, paste(tx_ids$txn_short, collapse = ",") )
        list_txn_b <- c(list_txn_b, paste(tx_ids$txn_long, collapse = ",") )
      }
      
      if(type == "A3SS") {
        
        #reverse strand becomes A5SS
        if (as.character(strand(eventGr[1])) == "+"){
          tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
        }else{
          tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
        }

        list_txn_a <- c(list_txn_a, paste(tx_ids$txn_short, collapse = ",") )
        list_txn_b <- c(list_txn_b, paste(tx_ids$txn_long, collapse = ",") )
      }
    } #for all events in annot
    
    annot[[names(tx_ids)[1]]] <- list_txn_a
    annot[[names(tx_ids)[2]]] <- list_txn_b
    
    events_with_txn[[paste0(type,"_","events")]] <- annot
      
  
  } #for each event type
  
  return(events_with_txn)
  
}

mapTranscriptsSEevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  tx_ids <- list()
  
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
  intron.skipping <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_target),
                                            ranges = IRanges(
                                              start = end(eventGr$exon_upstream) + 1,  
                                              end = start(eventGr$exon_downstream) - 1),
                                            strand = strand(eventGr$exon_target)
  )
  
  intron.inclusion <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_target),
                                             ranges = IRanges(
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
  mytx.ids.intron.inclusion <- gtf_exons$transcript_id[subjectHits(ovl.intron.inclusion)]
  
  ovl.intron.skipping <- GenomicRanges::findOverlaps(intron.skipping, 
                                                     gtf_exons, type = "any")
  mytx.ids.intron.skipping <- gtf_exons$transcript_id[subjectHits(ovl.intron.skipping)]
  
  #decide wich transcripts to plot in inclusion and skipping tracks
  if (is_strict){
    mytx.ids.3exons <- intersect(mytx.ids.e1, mytx.ids.e3) #has both flaking exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and the target exon
    mytx.ids.2exons <- intersect(mytx.ids.e1, mytx.ids.e3)
    
  }else {
    mytx.ids.3exons <- union(mytx.ids.e1, mytx.ids.e3) #has either flaking exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and the target exon  
    mytx.ids.2exons <- union(mytx.ids.e1, mytx.ids.e3)
  }
  
  mytx.ids.3exons <- setdiff(mytx.ids.3exons, mytx.ids.intron.inclusion)
  #mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.3exons)
  mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.intron.skipping)
  
  tx_ids[["txn_3exons"]] <- mytx.ids.3exons
  tx_ids[["txn_2exons"]] <- mytx.ids.2exons
  
  return(tx_ids)  
  
}

mapTranscriptsRIevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  tx_ids <- list()
  
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
                                   ranges = IRanges(
                                     start = end(eventGr$exon_upstream) + 1,  
                                     end = start(eventGr$exon_downstream) - 1),
                                   strand = strand(eventGr$exon_ir)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron <- GenomicRanges::findOverlaps(intron, gtf_exons, type = "any")
  mytx.ids.intron <- gtf_exons$transcript_id[subjectHits(ovl.intron)]
  
  #decide wich transcripts to plot in retention and non-retention tracks
  if (is_strict){
    tx.ids.nonRetention <- intersect(mytx.ids.e1, mytx.ids.e3) #has both upstream and downstream exons
    
  }else {
    tx.ids.nonRetention <- union(mytx.ids.e1, mytx.ids.e3) #has either upstream and downstream exons
  }
  
  tx.ids.nonRetention <- setdiff(tx.ids.nonRetention, mytx.ids.intron)
  tx.ids.Retention <- mytx.ids.e2
  
  tx_ids[["txn_nonRetention"]] <- tx.ids.nonRetention
  tx_ids[["txn_retention"]] <- tx.ids.Retention
  
  return(tx_ids)
  
}

mapTranscriptsMXEevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  tx_ids <- list()
  
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
  intron.mxe.exon1 <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_1),
                                             ranges = IRanges(
                                               start = c(end(eventGr$exon_upstream) + 1,
                                                         end(eventGr$exon_1) + 1),  
                                               end = c(start(eventGr$exon_1) - 1,
                                                       start(eventGr$exon_downstream) -1)
                                             ),
                                             strand = strand(eventGr$exon_1)
  )
  
  intron.mxe.exon2 <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_2),
                                             ranges = IRanges(
                                               start = c(end(eventGr$exon_upstream) + 1,
                                                         end(eventGr$exon_2) + 1),  
                                               end = c(start(eventGr$exon_2) - 1,
                                                       start(eventGr$exon_downstream) -1)
                                             ),
                                             strand = strand(eventGr$exon_2)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.mxe.exon1 <- GenomicRanges::findOverlaps(intron.mxe.exon1, gtf_exons, type = "any")
  mytx.ids.intron1 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon1)]
  
  ovl.mxe.exon2 <- GenomicRanges::findOverlaps(intron.mxe.exon2, gtf_exons, type = "any")
  mytx.ids.intron2 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon2)]
  
  
  #decide wich transcripts to plot in inclusion and skipping tracks
  if (is_strict){
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.e3, mytx.ids.e4) #has both flanking exons
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1) #and exon1
    
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.e3, mytx.ids.e4) #has both flanking exons
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2) #and exon2
    
  }else {
    mytx.ids.mxe.exon1 <- union(mytx.ids.e3, mytx.ids.e4) #has either flanking exons
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1) #and exon 1
    
    mytx.ids.mxe.exon2 <- union(mytx.ids.e3, mytx.ids.e4) #has both flanking exons
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2) #and exon2
  }
  
  #remove transcripts with exons in intronic regions
  mytx.ids.mxe.exon1 <- setdiff(mytx.ids.mxe.exon1, mytx.ids.intron1)
  mytx.ids.mxe.exon2 <- setdiff(mytx.ids.mxe.exon2, mytx.ids.intron2)
  
  
  tx_ids[["txn_mxe_exon1"]] <- mytx.ids.mxe.exon1
  tx_ids[["txn_mxe_exon2"]] <- mytx.ids.mxe.exon2
  
  return(tx_ids)

}

mapTranscriptsA5SSevent <- function(eventGr, gtf_exons){
  
  tx_ids <- list()
  
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
  intron.short <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_short),
                                         ranges = IRanges(
                                           start = end(eventGr$exon_short) + 1,  
                                           end = start(eventGr$exon_flanking) - 1),
                                         strand = strand(eventGr$exon_short)
  )
  
  intron.long <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_long),
                                        ranges = IRanges(
                                          start = end(eventGr$exon_long) + 1,  
                                          end = start(eventGr$exon_flanking) - 1),
                                        strand = strand(eventGr$exon_long)
  )
  
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.short <- GenomicRanges::findOverlaps(intron.short, gtf_exons, type = "any")
  mytx.ids.intron.short <- gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
  
  ovl.intron.long <- GenomicRanges::findOverlaps(intron.long, gtf_exons, type = "any")
  mytx.ids.intron.long <- gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
  
  #decide wich transcripts to plot in short and long tracks
  mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
  mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
  
  #remove transcripts with exons overlapping intronic regions
  mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
  mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
  
  tx_ids[["txn_short"]] <- mytx.ids.short
  tx_ids[["txn_long"]] <- mytx.ids.long
  
  return(tx_ids)
  
}

mapTranscriptsA3SSevent <- function(eventGr, gtf_exons){
  
  tx_ids <- list()
  
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
  intron.short <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_short),
                                         ranges = IRanges(
                                           start = end(eventGr$exon_flanking) + 1,  
                                           end = start(eventGr$exon_short) - 1),
                                         strand = strand(eventGr$exon_short)
  )
  
  intron.long <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_long),
                                        ranges = IRanges(
                                          start = end(eventGr$exon_flanking) + 1,  
                                          end = start(eventGr$exon_long) - 1),
                                        strand = strand(eventGr$exon_long)
  )  
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.short <- GenomicRanges::findOverlaps(intron.short, gtf_exons, type = "any")
  mytx.ids.intron.short <- gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
  
  ovl.intron.long <- GenomicRanges::findOverlaps(intron.long, gtf_exons, type = "any")
  mytx.ids.intron.long <- gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
  
  #decide wich transcripts to plot in short and long tracks
  mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
  mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
  
  #remove transcripts with exons overlapping intronic regions
  mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
  mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
  
  tx_ids[["txn_short"]] <- mytx.ids.short
  tx_ids[["txn_long"]] <- mytx.ids.long
  
  return(tx_ids)
}
  