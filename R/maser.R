
create_GRanges_ASS <- function(events){
  
  exon_long <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$longExonStart_0base+1,
                     end = events$longExonEnd),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  exon_short <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$shortES+1,
                     end = events$shortEE),
    strand = events$strand,
    ID = as.vector(events$ID),
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
    
  )
  exon_flanking <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$flankingES+1,
                     end = events$flankingEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  grl <- GRangesList("exon_long" = exon_long, "exon_short" = exon_short,
                     "exon_flanking" = exon_flanking)
  return(grl)
}

create_GRanges_ES <- function(events){
  
  exon_target <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$exonStart_0base+1,
                     end = events$exonEnd),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon_upstream <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$upstreamES+1,
                     end = events$upstreamEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon_downstream <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$downstreamES+1,
                     end = events$downstreamEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  grl <- GRangesList("exon_target" = exon_target,
                     "exon_upstream" = exon_upstream,
                     "exon_downstream" = exon_downstream)
  return(grl)
}

create_GRanges_IR <- function(events){
  
  exon_ir <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$riExonStart_0base+1,
                     end = events$riExonEnd),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon_upstream <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$upstreamES+1,
                     end = events$upstreamEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon_downstream <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$downstreamES+1,
                     end = events$downstreamEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  grl <- GRangesList("exon_ir" = exon_ir,
                     "exon_upstream" = exon_upstream,
                     "exon_downstream" = exon_downstream)
  return(grl)
}

create_GRanges_MXE <- function(events){
  
  exon1 <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$X1stExonStart_0base+1,
                     end = events$X1stExonEnd),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon2 <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$X2ndExonStart_0base+1,
                     end = events$X2ndExonEnd),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon_upstream <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$upstreamES+1,
                     end = events$upstreamEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  exon_downstream <- GRanges(
    seqnames = events$chr,
    ranges = IRanges(start = events$downstreamES+1,
                     end = events$downstreamEE),
    strand = events$strand,
    ID = events$ID,
    GeneID = events$GeneID,
    geneSymbol = events$geneSymbol
  )
  
  grl <- GRangesList("exon_1" = exon1,
                     "exon_2" = exon2,
                     "exon_upstream" = exon_upstream,
                     "exon_downstream" = exon_downstream)
  return(grl)
}

create_GRanges <- function(events, type){
  
  if(type == "A3SS" || type == "A5SS") {
    grl <- create_GRanges_ASS(events)
  }
  
  if (type == "SE"){
    grl <- create_GRanges_ES(events)
  }
  
  if (type == "RI"){
    grl <- create_GRanges_IR(events)
  }
  
  if (type == "MXE"){
    grl <- create_GRanges_MXE(events)
  }
  
  return(grl)
  
}

#' Create a maser object by importing rMATS splice events.
#' 
#' @param path a character specifiying the folder containing rMATS output files.
#' @param cond_labels a character vector of length 2 describing labels for experimental conditions.
#' @param rtype a character indicating the read type. Possible values are \code{c("ReadsOnTargetAndJunction", "JunctionCountOnly")}. 
#' @return A maser object.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' @export
maser <- function(path, cond_labels,
                  rtype = "ReadsOnTargetAndJunction"){
  
  rmats_out <- list.files(path, pattern = rtype, full.names = F)
  mats <- list()
  
  if(rtype == "JunctionCountOnly"){
    counts.col <- c("IJC_SAMPLE_1", "IJC_SAMPLE_2")
  }else{
    counts.col <- c("IC_SAMPLE_1", "IC_SAMPLE_2")
  }
  
  # For each AS type
  for (f in rmats_out) {
    
    events <- read.table(paste0(path,"/",f), sep = "\t",
                         stringsAsFactors = F, header = T)
    type <-  unlist(strsplit(f, ".", fixed = T))[1]
    
    # prepare read counts matrix
    inc1 <- strsplit(events[ , counts.col[1]], ",")
    inc2 <- strsplit(events[ , counts.col[2]], ",")
    
    reads.mat <- matrix(0, nrow = length(inc1),
                        ncol = length(inc1[[1]]) +
                          length(inc2[[1]]))
    for (i in 1:length(inc1)) {
      reads.mat[i, ] <- suppressWarnings(
                      c(as.numeric(inc1[[i]]), as.numeric(inc2[[i]]) )
                      )
    }
    rownames(reads.mat) <- events$ID
    col_names <- c(paste0(cond_labels[1], "_",
                          seq(1, length(inc1[[1]]), 1)),
                   paste0(cond_labels[2], "_",
                          seq(1, length(inc2[[1]]), 1)) )
    
    colnames(reads.mat) <- col_names
    mats[[paste0(type,"_","counts")]] <- reads.mat
    
    # prepare PSI matrix
    inc1 <- strsplit(events[ , "IncLevel1"], ",")
    inc2 <- strsplit(events[ , "IncLevel2"], ",")
    
    reads.mat <- matrix(0, nrow = length(inc1),
                        ncol = length(inc1[[1]]) +
                          length(inc2[[1]]))
    
    for (i in 1:length(inc1)) {
      reads.mat[i, ] <- suppressWarnings(
                          c(as.numeric(inc1[[i]]), as.numeric(inc2[[i]]) )
                        )
    }
    rownames(reads.mat) <- events$ID
    colnames(reads.mat) <- col_names
    mats[[paste0(type,"_","PSI")]] <- reads.mat
    
    # Number of samples condition 1 and 2
    mats[["n_cond1"]] <- length(inc1[[1]])
    mats[["n_cond2"]] <- length(inc2[[1]])
    
    mats[["conditions"]] <- cond_labels
    
    
    # rMATS stats
    mats[[paste0(type,"_","stats")]] <-
      events[ , c("ID", "PValue", "FDR", "IncLevelDifference")]
    
    # Genomic ranges of alternative splicing events
    grl <- create_GRanges(events, type)
    mats[[paste0(type,"_","gr")]] <- grl
    
    # Event annotation
    mats[[paste0(type,"_","events")]] <-
      events[ , c("ID", "GeneID", "geneSymbol")]
    
    #cat("Importing ", f, " \n")
    
    
  }
  
  class(mats) <- "maser"
  return(mats)
  
}


#summary.maser <- function(x) "maser"
#print.maser <- function(x) "maser"

#' Print a maser object.
#' 
#' @param x a maser object.
#' @param ... further arguments.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' print(hypoxia)
#' @export
print.maser <- function(x, ...){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  counts <- rep(0, length(as_types))
  total <- 0
  line_events <- ""
  
  for (type in as_types) {
    PSI <- x[[paste0(type,"_","PSI")]]
    nevents <- nrow(PSI)
    line_events <- paste0(line_events, 
                          type, ".......... ", nevents, " events\n")
    total <- total + nevents
  }
  
  line1 <- paste0("A maser object with ", total, " splicing events.\n\n")
  line2 <- paste0("Samples description: \n", "Label=", x$conditions[1], "     n=", x$n_cond1, 
                  " replicates\n")
  line3 <- paste0( "Label=", x$conditions[2], "     n=", x$n_cond1, 
                   " replicates\n\n")
  line4 <- paste0("Splicing events: \n")
  
  cat(paste0(line1, line2, line3, line4, line_events))
  
}

summary.maser <- function(x){
  
  print.maser(x)
}

is.maser <- function(x){
  
  attributes <- c("A3SS_counts", "A3SS_PSI", "n_cond1", "n_cond2", "conditions",
                  "A3SS_stats", "A3SS_gr", "A3SS_events", "A5SS_counts", 
                  "A5SS_PSI",  "A5SS_stats", "A5SS_gr", "A5SS_events",
                  "MXE_counts", "MXE_PSI", "MXE_stats", "MXE_gr", "MXE_events",
                  "RI_counts", "RI_PSI", "RI_stats", "RI_gr", "RI_events",
                  "SE_counts", "SE_PSI", "SE_stats", "SE_gr", "SE_events")
  
  #check if all elements are allowed attributes, and that all attributes are present
  # and has been created using the constructor maser()
  if (any(!names(x) %in% attributes)  || any(!attributes %in% names(x)) || 
      !class(x) == "maser"){
    return(FALSE)
  }else{
    return(TRUE)    
  }
  
}

asDataFrame <- function(events, type){
  
  annot <- events[[paste0(type,"_","events")]]
  stats <- events[[paste0(type,"_","stats")]]
  PSI <- events[[paste0(type,"_","PSI")]]
  grl <- events[[paste0(type,"_","gr")]]
  
  df <- annot
  df <- cbind(df, stats[,2:4])
  
  PSI_1 <- c("NA", nrow(df))
  PSI_2 <- c("NA", nrow(df))
  
  idx.cond1 <- seq(1, events$n_cond1, 1)
  idx.cond2 <- seq(events$n_cond1+1, events$n_cond1+events$n_cond2, 1)
  
  for (i in 1:nrow(df)) {
    PSI_1[i] <- paste(PSI[i, idx.cond1], collapse = ",")
    PSI_2[i] <- paste(PSI[i, idx.cond2], collapse = ",")
  } 
  
  df <- cbind(df, PSI_1 = PSI_1, PSI_2 = PSI_2)
  
  exon_df <- data.frame(Chr = as.character(seqnames(grl[[1]])),
                        Strand = as.character(strand(grl[[1]])))
  
  for (i in 1:length(names(grl))) {
    res <- as.data.frame(grl[[i]])
    coord <- paste0(res$start, "-", res$end)  
    label <- names(grl[i])
    exon_df <- cbind(exon_df, Exon = coord)
  }
  
  colnames(exon_df) <- c("Chr", "Strand", names(grl))
  
  df <- cbind(df, exon_df)
  return(df)
  
}

#' Retrieve PSI (percent spliced in) values from a maser object.
#' 
#' @param events a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a data.frame.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(PSI(hypoxia, "SE"))
#' @export
PSI <- function(events, type){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  return(events[[paste0(type,"_","PSI")]])
  
}

#' Retrieve raw read counts values from a maser object.
#' 
#' @param events a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a data.frame.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(counts(hypoxia, "SE"))
#' @export
counts <- function(events, type){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  return(events[[paste0(type,"_","counts")]])
  
}

#' Retrieve annotation of splice events from a maser object.
#' 
#' @param events a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a data.frame.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(annot(hypoxia, "SE"))
#' @export
annot <- function(events, type){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  return(events[[paste0(type,"_","events")]])
  
}

#' Retrieve rMATS stats of differential splicing from a maser object.
#' 
#' @param events a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a data.frame.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(stats(hypoxia, "SE"))
#' @export
stats <- function(events, type){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  return(asDataFrame(events, type))
  
}
