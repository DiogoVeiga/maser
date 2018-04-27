#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
create_GRanges_exon <- function(chr, start, end, strand, id, gene_id, 
                                gene_symbol){
  
  if(any(!grepl("chr", chr))){
    echr <- paste0("chr", chr)
  }else{
    echr <- chr
  }
  
  exon <- GRanges(seqnames = echr, 
    ranges = IRanges::IRanges(start = start, end = end), strand = strand,
    ID = id, GeneID = gene_id, geneSymbol = gene_symbol)
  
  return(exon)
  
}

create_GRanges_ASS <- function(events){
  
  exon_long <- create_GRanges_exon(events$chr, events$longExonStart_0base+1,
    events$longExonEnd, events$strand, events$ID, events$GeneID, 
    events$geneSymbol)
  
  exon_short <- create_GRanges_exon(events$chr, events$shortES+1, events$shortEE,
              events$strand, events$ID, events$GeneID, events$geneSymbol)
  
  exon_flanking <- create_GRanges_exon(events$chr, events$flankingES+1, 
              events$flankingEE, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  grl <- GRangesList("exon_long" = exon_long, "exon_short" = exon_short,
                     "exon_flanking" = exon_flanking)
  return(grl)

}

create_GRanges_ES <- function(events){
  
 
  exon_target <- create_GRanges_exon(events$chr, events$exonStart_0base+1,
              events$exonEnd, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  exon_upstream <- create_GRanges_exon(events$chr, events$upstreamES+1,
              events$upstreamEE, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  exon_downstream <- create_GRanges_exon(events$chr, events$downstreamES+1,
              events$downstreamEE, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  grl <- GRangesList("exon_target" = exon_target,
                     "exon_upstream" = exon_upstream,
                     "exon_downstream" = exon_downstream)
  return(grl)
}

create_GRanges_IR <- function(events){
  
  exon_ir <- create_GRanges_exon(events$chr, events$riExonStart_0base+1,
              events$riExonEnd, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  exon_upstream <- create_GRanges_exon(events$chr, events$upstreamES+1,
              events$upstreamEE, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  exon_downstream <- create_GRanges_exon(events$chr, events$downstreamES+1,
              events$downstreamEE, events$strand, events$ID, events$GeneID, 
              events$geneSymbol)
  
  grl <- GRangesList("exon_ir" = exon_ir,
                     "exon_upstream" = exon_upstream,
                     "exon_downstream" = exon_downstream)
  return(grl)
}

create_GRanges_MXE <- function(events){
  
  exon1 <- create_GRanges_exon(events$chr, events$X1stExonStart_0base+1,
                  events$X1stExonEnd, events$strand, events$ID, events$GeneID, 
                  events$geneSymbol)
  
  exon2 <- create_GRanges_exon(events$chr, events$X2ndExonStart_0base+1,
                  events$X2ndExonEnd, events$strand, events$ID, events$GeneID, 
                  events$geneSymbol)
  
  exon_upstream <- create_GRanges_exon(events$chr, events$upstreamES+1,
                  events$upstreamEE, events$strand, events$ID, events$GeneID, 
                  events$geneSymbol)
  
  exon_downstream <- create_GRanges_exon(events$chr, events$downstreamES+1,
                  events$downstreamEE, events$strand, events$ID, events$GeneID, 
                  events$geneSymbol)

  grl <- GRangesList("exon_1" = exon1, "exon_2" = exon2,
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

#' Create a maser object by importing rMATS splicing events.
#' 
#' @param path a character specifiying the folder containing rMATS output files.
#' @param cond_labels a character vector of length 2 describing labels for 
#' experimental conditions.
#' @param ftype a character indicating the rMATS file type.
#' Possible values are \code{c("ReadsOnTargetAndJunctionCounts",
#' "JunctionCountOnly", "JCEC", "JC")}. 
#' @return A maser object.
#' @details This function creates a maser object by importing rMATS output. 
#' \code{ftype} indicates which rMATS files to import. 
#' \code{ReadsOnTargetandJunction or JunctionCountOnly} are used in rMATS 3.2.5 
#' or lower. Newer versions (>4.0.1) use \code{"JCEC" or "JC"} nomenclature.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' @export
maser <- function(path, cond_labels,
                  ftype = c("ReadsOnTargetAndJunctionCounts", 
                            "JunctionCountOnly",
                            "JCEC", "JC")){
  
  ftype <- match.arg(ftype)
  
  rmats_out <- list.files(path, pattern = paste0(ftype, ".txt"),
                          full.names = FALSE)
  mats <- maser_empty()
  
  if(ftype == "ReadsOnTargetAndJunctionCounts"){
    counts.col <- c("IC_SAMPLE_1", "IC_SAMPLE_2")
  }else{
    counts.col <- c("IJC_SAMPLE_1", "IJC_SAMPLE_2")    
  }
  
  if(!grepl("/$", path)){
    path <- paste0(path, "/")
  }
  
  # For each AS type
  for (f in rmats_out) {
    
    events <- read.table(paste0(path, f), sep = "\t",
                         stringsAsFactors = FALSE, header = TRUE)
    type <-  unlist(strsplit(f, ".", fixed = TRUE))[1]
    
    # prepare read counts matrix
    inc1 <- strsplit(events[ , counts.col[1]], ",")
    inc2 <- strsplit(events[ , counts.col[2]], ",")
    
    reads.inc1 <- suppressWarnings(matrix(as.numeric(unlist(inc1)), 
                 nrow = length(inc1), ncol = length(inc1[[1]]), byrow = TRUE))
    
    reads.inc2 <- suppressWarnings(matrix(as.numeric(unlist(inc2)), 
                 nrow = length(inc2), ncol = length(inc2[[1]]), byrow = TRUE))
    
    reads.mat <- cbind(reads.inc1, reads.inc2)
    rownames(reads.mat) <- events$ID
    col_names <- c(paste0(cond_labels[1], "_", seq(1, length(inc1[[1]]), 1)),
                   paste0(cond_labels[2], "_", seq(1, length(inc2[[1]]), 1)) )
    
    colnames(reads.mat) <- col_names
    mats[[paste0(type,"_","counts")]] <- reads.mat
    
    # prepare PSI matrix
    inc1 <- strsplit(events[ , "IncLevel1"], ",")
    inc2 <- strsplit(events[ , "IncLevel2"], ",")
    
    reads.inc1 <- suppressWarnings(matrix(as.numeric(unlist(inc1)), 
                   nrow = length(inc1), ncol = length(inc1[[1]]), byrow = TRUE))
    
    reads.inc2 <- suppressWarnings(matrix(as.numeric(unlist(inc2)), 
                   nrow = length(inc2), ncol = length(inc2[[1]]), byrow = TRUE))
    
    reads.mat <- cbind(reads.inc1, reads.inc2)
    
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
#' @return a character.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' print(hypoxia)
#' @export
print.maser <- function(x, ...){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  total <- 0
  line_events <- ""
  nevents <- 0
  
  lapply(as_types, function(type){
    PSI <- x[[paste0(type,"_","PSI")]]
    nevents <- nrow(PSI)
    line_events <<- paste0(line_events,
                          type, ".......... ", nevents, " events\n")
    total <<- total + nevents
  })
  
  line1 <- paste0("A maser object with ", total, " splicing events.\n\n")
  line2 <- paste0("Samples description: \n", "Label=", x$conditions[1],
                  "     n=", x$n_cond1, 
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
  
  #check if all elements are allowed attributes, and that all attributes are
  # present and has been created using the constructor maser()
  if (any(!names(x) %in% attributes)  || any(!attributes %in% names(x)) || 
      !class(x) == "maser"){
    return(FALSE)
  }else{
    return(TRUE)    
  }
  
}

#' @importFrom GenomicRanges GRangesList
maser_empty <- function(){
  
  as_types = c("A3SS", "A5SS", "SE", "RI", "MXE")
  x <- list()
  
  for (type in as_types) {
    x[[paste0(type,"_","counts")]] <- matrix(nrow = 0, ncol = 0)
    x[[paste0(type,"_","PSI")]] <- matrix(nrow = 0, ncol = 0)
    x[[paste0(type,"_","stats")]] <- data.frame()
    x[[paste0(type,"_","events")]] <- data.frame()
    x[[paste0(type,"_","gr")]] <- GRangesList()
  }
  x$n_cond1 <- NA
  x$n_cond2 <- NA
  x$conditions <- "NA"
  
  class(x) <- "maser"
  return(x)

}

asDataFrame <- function(events, type){
  
  annot <- events[[paste0(type,"_","events")]]
  stats <- events[[paste0(type,"_","stats")]]
  PSI <- events[[paste0(type,"_","PSI")]]
  grl <- events[[paste0(type,"_","gr")]]
  
  df <- annot
  df <- cbind(df, stats[,2:4])
  
  idx.cond1 <- seq(1, events$n_cond1, 1)
  idx.cond2 <- seq(events$n_cond1+1, events$n_cond1+events$n_cond2, 1)
  
  PSI_1 <- vapply(seq_along(PSI[,1]), function(i){
    paste(PSI[i, idx.cond1], collapse = ",")
    }, character(1)
  )
  
  PSI_2 <- vapply(seq_along(PSI[,1]), function(i){
    paste(PSI[i, idx.cond2], collapse = ",")
  }, character(1)
  )
  
  df <- cbind(df, PSI_1 = PSI_1, PSI_2 = PSI_2)
  
  exon_df <- data.frame(Chr = as.character(seqnames(grl[[1]])),
                        Strand = as.character(strand(grl[[1]])))
  
  lapply(seq_along(grl), function(i){
    res <- as.data.frame(grl[[i]])
    coord <- paste0(res$start, "-", res$end)  
    exon_df <<- cbind(exon_df, Exon = coord)
  })
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
PSI <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  type <- match.arg(type)
  
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
counts <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  type <- match.arg(type)
  
  return(events[[paste0(type,"_","counts")]])
  
}

#' Retrieve annotation of splicing events from a maser object.
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
annot <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  type <- match.arg(type)
  
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
stats <- function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")){
  
  if(!is.maser(events)){
    stop("Parameter events has to be a maser object.")
  }
  
  type <- match.arg(type)
  
  return(asDataFrame(events, type))
  
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to maser")
}
