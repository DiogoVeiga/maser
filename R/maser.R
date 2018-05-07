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
#' @import methods
maser <- function(path, cond_labels,
                  ftype = c("ReadsOnTargetAndJunctionCounts", 
                            "JunctionCountOnly",
                            "JCEC", "JC")){
  
  ftype <- match.arg(ftype)
  
  rmats_out <- list.files(path, pattern = paste0(ftype, ".txt"),
                          full.names = FALSE)
  mats <- new("Maser")
  
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
    slot(mats, paste0(type,"_","counts")) <- reads.mat
    
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
    slot(mats, paste0(type,"_","PSI")) <- reads.mat
    
    # Number of samples condition 1 and 2
    mats@n_cond1 <- length(inc1[[1]])
    mats@n_cond2 <- length(inc2[[1]])
    mats@conditions <- cond_labels

    # rMATS stats
    slot(mats, paste0(type,"_","stats")) <-
      events[ , c("ID", "PValue", "FDR", "IncLevelDifference")]
    
    # Genomic ranges of alternative splicing events
    grl <- create_GRanges(events, type)
    slot(mats, paste0(type,"_","gr")) <- grl
    
    # Event annotation
    slot(mats, paste0(type,"_","events")) <-
      events[ , c("ID", "GeneID", "geneSymbol")]
    

  }

  return(mats)
  
}

#' S4 class to represent splicing events imported from rMATS.
#'
setClass("Maser",
         slots = list(A3SS_counts = "matrix", A3SS_PSI = "matrix",
                      A3SS_stats = "data.frame", A3SS_events = "data.frame",
                      A3SS_gr = "GRangesList",
                      A5SS_counts = "matrix", A5SS_PSI = "matrix",
                      A5SS_stats = "data.frame", A5SS_events = "data.frame",
                      A5SS_gr = "GRangesList",
                      SE_counts = "matrix", SE_PSI = "matrix",
                      SE_stats = "data.frame", SE_events = "data.frame",
                      SE_gr = "GRangesList",
                      RI_counts = "matrix", RI_PSI = "matrix",
                      RI_stats = "data.frame", RI_events = "data.frame",
                      RI_gr = "GRangesList",
                      MXE_counts = "matrix", MXE_PSI = "matrix",
                      MXE_stats = "data.frame", MXE_events = "data.frame",
                      MXE_gr = "GRangesList",
                      n_cond1 = "numeric",
                      n_cond2 = "numeric",
                      conditions = "character"
                      )
  )

setValidity2("Maser", function(object){
  if(!length(object@conditions) == 2){
    return("'conditions' slot must be a character vector of size 2.")
  }
  TRUE
})

setMethod("show", "Maser", function(object){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  total <- 0
  line_events <- ""
  nevents <- 0
  
  lapply(as_types, function(type){
    
    PSI <- slot(object, paste0(type,"_","PSI"))
    nevents <- nrow(PSI)
    line_events <<- paste0(line_events,
                           type, ".......... ", nevents, " events\n")
    total <<- total + nevents
  })
  
  conditions <- slot(object, "conditions")
  line1 <- paste0("A Maser object with ", total, " splicing events.\n\n")
  line2 <- paste0("Samples description: \n", "Label=", conditions[1],
                  "     n=", slot(object, "n_cond1"), 
                  " replicates\n")
  line3 <- paste0( "Label=", conditions[2], "     n=", slot(object, "n_cond2"), 
                   " replicates\n\n")
  line4 <- paste0("Splicing events: \n")
  
  cat(paste0(line1, line2, line3, line4, line_events))
  
  
})

setAs("Maser", "list", function(from){
  
  mapply(function(y) {
    #apply as.list if the slot is again an user-defined object
    #therefore, as.list gets applied recursively
    if (inherits(slot(from,y),"Maser")) {
      as.list(slot(from,y))
    } else {
      #otherwise just return the slot
      slot(from,y)
    }
  },
  slotNames(class(from)),
  SIMPLIFY=FALSE)
  
})

setAs("list", "Maser", function(from){
  y <- new("Maser")
  if (!any(names(from) %in% slotNames(y))){
    return("Invalid slot names.")
  }
  
  lapply(names(from), function(aslot){
    slot(y, paste0(aslot)) <<- from[[paste0(aslot)]]
  })
  return(y)
})

#' @import methods
create_stats <- function(events, type){
  
  annot <- slot(events, paste0(type,"_","events"))
  stats <- slot(events, paste0(type,"_","stats"))
  PSI <- slot(events, paste0(type,"_","PSI"))
  grl <- slot(events, paste0(type,"_","gr"))
  
  df <- annot
  df <- cbind(df, stats[,2:4])
  
  idx.cond1 <- seq(1, events@n_cond1, 1)
  idx.cond2 <- seq(events@n_cond1+1, events@n_cond1+events@n_cond2, 1)
  
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
#' @return a matrix.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(PSI(hypoxia, "SE"))
#' @export
#' @import methods
setGeneric("PSI", function(events, type) standardGeneric("PSI"))

#' Retrieve PSI (percent spliced in) values from a maser object.
#' 
#' @param events a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a matrix.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(PSI(hypoxia, "SE"))
#' @export
#' @import methods
setMethod("PSI", signature(events="Maser", type="character"),
    function(events, type = c("A3SS", "A5SS", "SE", "RI", "MXE")) {
    
      type <- match.arg(type)
      return(slot(events, paste0(type,"_","PSI")))            
      
})


#' Retrieve raw read counts values from a maser object.
#' 
#' @param object a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a matrix.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(counts(hypoxia, "SE"))
#' @export
#' @import methods
#' @importMethodsFrom BiocGenerics counts

setMethod("counts", "Maser", 
      function(object, type = c("A3SS", "A5SS", "SE", "RI", "MXE"))  {
        
        type <- match.arg(type)
        return(slot(object, paste0(type,"_","counts")))
})


#' Retrieve annotation of splicing events from a maser object.
#' 
#' @param object a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a data.frame.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(annotation(hypoxia, "SE"))
#' @export
#' @import methods
#' @importMethodsFrom BiocGenerics annotation

setMethod("annotation", "Maser", 
      function(object, type = c("A3SS", "A5SS", "SE", "RI", "MXE"))  {

        type <- match.arg(type)
        return(slot(object, paste0(type,"_","events")))
})


#' Retrieve rMATS stats of differential splicing from a maser object.
#' 
#' @param object a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @return a data.frame.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(summary(hypoxia, "SE"))
#' @export
#' @import methods
setMethod("summary", "Maser", 
    function(object, type = c("A3SS", "A5SS", "SE", "RI", "MXE"))  {
            
      type <- match.arg(type)
      return(create_stats(object, type))
})

#' Retrieve genomic ranges of splicing events from a maser object.
#' 
#' @param x a maser object.
#' @param type a character indicating the splice type. Possible values 
#' are  \code{c("A3SS", "A5SS", "SE", "RI", "MXE")}. 
#' @param ... additional arguments.
#' @return a GRangesList.
#' @examples
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' head(granges(hypoxia, type = "SE"))
#' @import methods
setMethod("granges", "Maser", 
    function(x, type = c("A3SS", "A5SS", "SE", "RI", "MXE"), ...)  {
      
      type <- match.arg(type)
      return(slot(x, paste0(type,"_","gr")))
})


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to maser")
}
