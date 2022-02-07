
#' @importFrom data.table fread
#' @importFrom methods as
#' @importFrom dplyr filter
createGRangesUniprotKBtrack <- function(track_name){
  
  Name <- NULL
  track_df <- urlTracksUniprotKB()
  
  if(!track_name %in% track_df$Name){
    stop(cat("Unknown track name."))
  }
  
  track <- dplyr::filter(track_df, Name %in% track_name)
  # bed <- read.table(as.character(track$URL), header = FALSE, sep = "\t", 
  #                   quote = NULL,
  #                   stringsAsFactors = FALSE)
  
  bed <- data.table::fread(as.character(track$URL), header = FALSE, sep = "\t", 
                           quote = "", stringsAsFactors = FALSE,
                           data.table = FALSE, showProgress = FALSE)
  
  colnames(bed)[1:6] <- c("chr", "start", "end", "Uniprot_ID", "V5", "strand") 
  res <- strsplit(bed$V14, ";")
  
  name <- vapply(seq_along(res), function(i){

      if (length(unlist(res[i])) > 1){
        #return(paste0(bed$Uniprot_ID[i], ":", res[[i]][[2]]))
        return(paste0(bed$Uniprot_ID[i], ":", bed$V14[i]))
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
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr arrange
availableFeaturesUniprotKB <- function(){
  
  track_df <- read.csv("https://raw.githubusercontent.com/DiogoVeiga/maser_aux/main/FeaturesUniprotKB.csv")
  track_df <- dplyr::arrange(track_df, Category, Name)
  track_df <- dplyr::select(track_df, c("Name", "Description", "Category"))
  
  return(track_df)
}

#' @importFrom dplyr filter
urlTracksUniprotKB <- function(){
  
  track_df <- read.csv("https://raw.githubusercontent.com/DiogoVeiga/maser_aux/main/FeaturesUniprotKB.csv")
  track_df <- dplyr::select(track_df, c("Name", "Description", "URL"))
  
  return(track_df)
}

#' @import GenomicRanges
#' @import GenomeInfoDb
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