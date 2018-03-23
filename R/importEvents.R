
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

#create maser object by importing rMATS events
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
            reads.mat[i, ] <- c(as.numeric(inc1[[i]]), as.numeric(inc2[[i]]) )
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
            reads.mat[i, ] <- c(as.numeric(inc1[[i]]), as.numeric(inc2[[i]]) )
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

        cat("Importing ", f, " \n")


    }

    class(mats) <- "maser"
    return(mats)

}


#summary.maser <- function(x) "maser"
#print.maser <- function(x) "maser"

print.maser <- function(x){
  
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

