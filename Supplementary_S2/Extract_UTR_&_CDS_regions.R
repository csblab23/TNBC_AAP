## This script contains the code to extract sequence region of Untranslated regions (5'UTR; 3'UTR) & Coding sequence (CDS) in BED files to find the shared and unique regions ##

# Load required packages
library(GenomicFeatures)
library(rtracklayer)

# Specify path to GTF file
gtf_path <- "gencode.v44.basic.annotation.gtf"  

# Define transcript IDs of interest
normal_tx <- "ENST00000707077.1"   # pr1077 transcript
cancer_tx <- "ENST00000456174.6"   # pr1079 transcript

# Build TxDb object
txdb <- makeTxDbFromGFF(gtf_path, format = "gtf")

## Extract regions by transcript
# CDS
cds_list <- cdsBy(txdb, by = "tx", use.names = TRUE)
# 5′ UTR
utr5_list <- fiveUTRsByTranscript(txdb, use.names = TRUE)
# 3′ UTR
utr3_list <- threeUTRsByTranscript(txdb, use.names = TRUE)

# Function to export BED files
export_if_exists <- function(gr_list, tx_id, filename) {
  if (tx_id %in% names(gr_list)) {
    export(gr_list[[tx_id]], con = filename, format = "BED")
    message("File exported: ", filename)
  } else {
    message("Transcript not found: ", tx_id, " in ", filename)
  }
}

# --- For pr1077 transcript ---
export_if_exists(utr5_list, normal_tx, "normal_5utr.bed")
export_if_exists(utr3_list, normal_tx, "normal_3utr.bed")
export_if_exists(cds_list,  normal_tx, "normal_cds.bed")

# --- For pr1079 transcript ---
export_if_exists(utr5_list, cancer_tx, "cancer_5utr.bed")
export_if_exists(utr3_list, cancer_tx, "cancer_3utr.bed")
export_if_exists(cds_list,  cancer_tx, "cancer_cds.bed")

#############################################################
#############################################################
