
createGRangesUniprotKBtrack <- function(track_name){
  
  track_df <- urlTracksUniprotKB()
  
  if(!track_name %in% track_df$Name){
    stop(cat("Unknown track name."))
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  bed <- read.table(as.character(track$URL), header = F, sep = "\t", quote = NULL,
                    stringsAsFactors = FALSE)
  
  colnames(bed)[1:6] <- c("chr", "start", "end", "Uniprot_ID", "V5", "strand") 
  res <- strsplit(bed$V14, ";")
  
  name <- rep("NA", length(res))
  for (i in 1:length(res)) {
    if (length(unlist(res[i])) > 1){
      name[i] <- paste0(bed$Uniprot_ID[i], ":", res[[i]][[2]])
    }else{
      name[i] <- paste0(bed$Uniprot_ID[i], ":", "NA")
    }
  }
  bed <- cbind(bed, Name = name)
  bed.gr <- as(bed, "GRanges")
  genome(bed.gr) <- "hg38"
  seqlevels(bed.gr, pruning.mode="coarse") <- c(paste0("chr", seq(1:22)), "chrX", "chrY")
  
  
  return(bed.gr)
  
}

# possibly remove - use createGRangesUniprotKBrack instead
downloadUniprotKBtrack <- function(track_name, destfolder = getwd()){
  
  
  track_df <- urlTracksUniprotKB()
  
  if(!track_name %in% track_df$Name){
    stop(cat("Unknown track name."))
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  
  dir_path <- file.path(destfolder, "UP000005640_9606_beds")
  
  if (!dir.exists(dir_path)){
    dir.create(dir_path)  
  }
  
  bedFile <- basename(as.character(track$URL))
  bedPath <- file.path(dir_path, bedFile)
  download.file(as.character(track$URL), destfile = bedPath) #will trigger error if URL does not work
  
  bed <- read.table(as.character(track$URL), header = F, sep = "\t", quote = NULL,
                    stringsAsFactors = FALSE)
  
  colnames(bed)[1:4] <- c("chr", "start", "end", "Uniprot_ID") 
  res <- strsplit(bed$V14, ";")
  
  name <- rep("NA", length(res))
  for (i in 1:length(res)) {
    #name <- c(name, res[[1]][[2]])
    name[i] <- res[[i]][[2]]
  }
  bed <- cbind(bed, Name = name)
  bed.gr <- as(bed, "GRanges")
}

#' Query available protein features in UniprotKB.
#' 
#' @return a data.frame.
#' @examples
#' head(availableFeaturesUniprotKB(), 10)
#' @export
availableFeaturesUniprotKB <- function(){
  
  trackMetadata <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_tracks.txt"
  data <- readLines(trackMetadata)
  
  track_df <- data.frame()
  
  for (i in 1:length(data)){
    
    aux <- gsub("\"", "", data[i])
    tokens <- strsplit(aux, split = "=", fixed = F)
    
    values <- tokens[[1]]
    
    trackName <- gsub(" description", "", values[2])
    trackName <- gsub("UniProtKB ", "", trackName)
    
    trackDesc <- gsub(" type", "", values[3])
    
    track_df <- rbind(track_df, 
                      data.frame(Name = trackName, Description = trackDesc))

  }
  
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

# internal function
urlTracksUniprotKB <- function(){
  
  trackMetadata <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/genome_annotation_tracks/UP000005640_9606_tracks.txt"
  data <- readLines(trackMetadata)
  
  track_df <- data.frame()
  
  for (i in 1:length(data)){
    
    aux <- gsub("\"", "", data[i])
    tokens <- strsplit(aux, split = "=", fixed = F)
    
    values <- tokens[[1]]
    
    trackName <- gsub(" description", "", values[2])
    trackName <- gsub("UniProtKB ", "", trackName)
    
    trackDesc <- gsub(" type", "", values[3])
    
    trackUrl <- gsub(" url", "", values[8])
    trackUrl <- gsub("_hub", "_beds", trackUrl)
    trackUrl <- gsub(".bb", ".bed", trackUrl)
    trackUrl <- gsub("/hg38", "", trackUrl)
    
    track_df <- rbind(track_df, 
                      data.frame(Name = trackName, Description = trackDesc, 
                                 URL = as.character(trackUrl)))
    
  }
  
  return(track_df)
}

# internal function
overlappingFeatures <- function(feature_gr, eventGr){
  
  #Define region around splicing event
  region <- range(unlist(eventGr))
  start(region) <- start(region) - 10
  end(region) <- end(region) + 10
  genome(region) <- "hg38"
  
  ov <- findOverlaps(eventGr, feature_gr)
  ov_features <- feature_gr[subjectHits(ov)]
  
  return(ov_features)
  
}