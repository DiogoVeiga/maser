#Parsing Uniprot ID mapping file
#Dowloaded in March 2018

path <- file.path(system.file("data-raw", package = "maser"),
                "HUMAN_9606_idmapping_EnsemblTRS.dat")

dat <- read.table(path, header = F, sep = "\t", quote = NULL,
                  stringsAsFactors = FALSE)
colnames(dat) <- c("UniprotKB_isoform_ID", "Mapping", "ENST_ID")

unip_id <- rep("NA", nrow(dat))

tokens <- strsplit(dat$UniprotKB_isoform_ID, "-")

for (i in 1:nrow(dat)) {
  unip_id[i] <- tokens[[i]][[1]]  
}

UKB_ENST_map <- cbind(UniprotKB_ID = unip_id, dat)
devtools::use_data(UKB_ENST_map, internal = TRUE)
