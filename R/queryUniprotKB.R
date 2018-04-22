
#' @importFrom utils read.table
#' @importFrom methods as
createGRangesUniprotKBtrack <- function(track_name){
  
  Name <- NULL
  track_df <- urlTracksUniprotKB()
  
  if(!track_name %in% track_df$Name){
    stop(cat("Unknown track name."))
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  bed <- read.table(as.character(track$URL), header = FALSE, sep = "\t", 
                    quote = NULL,
                    stringsAsFactors = FALSE)
  
  colnames(bed)[1:6] <- c("chr", "start", "end", "Uniprot_ID", "V5", "strand") 
  res <- strsplit(bed$V14, ";")
  
  name <- vapply(seq_along(res), function(i){

      if (length(unlist(res[i])) > 1){
        return(paste0(bed$Uniprot_ID[i], ":", res[[i]][[2]]))
      }else{
        return(paste0(bed$Uniprot_ID[i], ":", "NA"))
      }

    }, character(1)
  )
  
  bed <- cbind(bed, Name = name)
  bed.gr <- methods::as(bed, "GRanges")
  
  GenomeInfoDb::genome(bed.gr) <- "hg38"

  
  return(bed.gr)
  
}

#' Query available human protein features in UniprotKB.
#' 
#' @return a data.frame.
#' @examples
#' head(availableFeaturesUniprotKB(), 10)
#' @export
availableFeaturesUniprotKB <- function(){
  
  Name <- NULL
  trackMetadata <- paste0("ftp://ftp.uniprot.org/pub/databases/uniprot/",
                    "current_release/knowledgebase/genome_annotation_tracks/",
                    "UP000005640_9606_tracks.txt")
  
  data <- readLines(trackMetadata)
  
  track_df <- data.frame()
  
  track_meta <- lapply(seq_along(data), function(i){
    aux <- gsub("\"", "", data[i])
    tokens <- strsplit(aux, split = "=", fixed = FALSE)
    values <- tokens[[1]]
    
    trackName <- gsub(" description", "", values[2])
    trackName <- gsub("UniProtKB ", "", trackName)
    trackDesc <- gsub(" type", "", values[3])
    return(data.frame(Name = trackName, Description = trackDesc))
    
  })
  track_df <- do.call(rbind, track_meta)
  
  
  track_df_filt <- dplyr::filter(track_df, 
                                 !Name %in% c("non_std_aa", "peptide",
                                            "UP000005640_9606_proteome",
                                            "UP000005640_9606_variants"))
  
  dm_sites <- c("act-site", "binding", "Ca-binding", "coiled", "DNA-bind",
                "domain", "metal", "motif", "NP bind", "region", "repeat",
                "site", "Zn-fing")
  mp <- c("chain", "signal", "transit", "propep", "init_meth")
  mut <- c("mutagen")
  ptm <- c("carbohyd", "crosslnk", "disulfide", "lipid", "mod-res")
  sf <- c("helix", "turn", "strand")
  topo <- c("intramem", "topo-dom", "transmem")
  
  category <- rep("NA", nrow(track_df_filt))
  
  category[track_df_filt$Name %in% dm_sites] <- "Domain_and_Sites"
  category[track_df_filt$Name %in% mp] <- "Molecule_Processing"
  category[track_df_filt$Name %in% mut] <- "Mutagenesis"
  category[track_df_filt$Name %in% ptm] <- "PTM"
  category[track_df_filt$Name %in% sf] <- "Structural_Features"
  category[track_df_filt$Name %in% topo] <- "Topology"
  
  track_df_filt <- cbind(track_df_filt, Category = category)
  
  return(track_df_filt)
}

urlTracksUniprotKB <- function(){
  
  trackMetadata <- paste0("ftp://ftp.uniprot.org/pub/databases/uniprot/",
    "current_release/knowledgebase/genome_annotation_tracks/",
    "UP000005640_9606_tracks.txt")
  
  data <- readLines(trackMetadata)
  
  track_df <- data.frame()
  
  track_meta <- lapply(seq_along(data), function(i){
    
    aux <- gsub("\"", "", data[i])
    tokens <- strsplit(aux, split = "=", fixed = FALSE)
    
    values <- tokens[[1]]
    
    trackName <- gsub(" description", "", values[2])
    trackName <- gsub("UniProtKB ", "", trackName)
    
    trackDesc <- gsub(" type", "", values[3])
    
    trackUrl <- gsub(" url", "", values[8])
    trackUrl <- gsub("_hub", "_beds", trackUrl)
    trackUrl <- gsub(".bb", ".bed", trackUrl)
    trackUrl <- gsub("/hg38", "", trackUrl)
    
    return(data.frame(Name = trackName, Description = trackDesc, 
                      URL = as.character(trackUrl)))
    
  })
  track_df <- do.call(rbind, track_meta)
  
  return(track_df)
}

overlappingFeatures <- function(feature_gr, eventGr){
  
  #Define region around splicing event
  region <- range(unlist(eventGr))
  start(region) <- start(region) - 10
  end(region) <- end(region) + 10
  GenomeInfoDb::genome(region) <- "hg38"
  
  ov <- findOverlaps(eventGr, feature_gr)
  ov_features <- feature_gr[subjectHits(ov)]
  
  return(ov_features)
  
}