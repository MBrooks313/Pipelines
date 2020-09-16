########################################
# This script takes the FIMO output and converts it to BED format. It uses the
# TRANSFAC matrix ID to get the gene symbol replacing the FIMO output gene symbol
# from a two column tab-separated text file, ex. V$MYOD_01       MYOD1
#
# Written Sept 10th, 2020 by Matthew J Brooks
# Written for Fimo output from meme/5.1.0
# Tested in R/4.0


##############################
# USAGE: Rscript Fimo2bed_v1_1.R <sample> <matrix_txt>
#
# samp = full path to the fimo output including the file name, i.e. /path/to/fimo/sample1/fimo.tsv
# matrix = full path to the matrix2gene text file for TRANSFAC, i.e. /path/o/transfac/matrixMAPgene_Mm_v2.txt
#
##############################


#------------------------------------------------------------------------------#
# Fimo output to bed file function
fimo2bed <- function(samp, matrix){

  #Load libraries
  library("data.table")
  library(GenomicRanges)
  library(rtracklayer)

  #Load TRANSFAC matrix to symbol data
  transfac <- read.delim(matrix, sep = "\t", header = T, stringsAsFactors = F)

  # ----------------------------------------#
  # Process the fimo file

  #Load FIMO output
  fimo <- read.table(samp, sep = "\t", header = T, stringsAsFactors = F, blank.lines.skip = T, comment.char = "#")

  #Split the fimo$sequence_name column to get chr and start information
  chr_split <- strsplit(fimo$sequence_name, split = ":", )
  chr <- unlist(lapply(chr_split, function(x){x[[1]]}))
  pos_split <- strsplit(unlist(lapply(chr_split, function(x){x[[2]]})), "-")
  start <- as.numeric(unlist(lapply(pos_split, function(x){x[[1]]})))


  # ----------------------------------------#
  # Generate the BED file

  #Reorder FIMO column output to BED format
  fimo.bed <- cbind(chr = chr,
                    start = start + fimo$start,
                    end = start + fimo$stop,
                    fimo[,c(2,7,6,1,8:10)])

  #Sort the BED file based on chromosome and location
  fimo.bed <- fimo.bed[order(fimo.bed$chr, fimo.bed$start, fimo.bed$end),]

  #Replace the FIMO output gene symbol with the official gene symbol
  fimo.bed$motif_alt_id <- transfac[match(gsub("V_", "V\\$", fimo.bed$motif_id), transfac$Matrix), 3]

  # Make a GRanges object from the fimo.bed
  gr <- makeGRangesFromDataFrame(fimo.bed, keep.extra.columns = T)

  # Loop through each TF and reduce, combine and sort to get final non-redundant bed file
  final.bed <- GRanges()
  for (i in unique(fimo.bed$motif_alt_id)){
    tmp_gr <- gr[gr$motif_alt_id == i]
    tmp_red <- reduce(tmp_gr)
    mcols(tmp_red, level="within")$name <- i
    final.bed <- append(final.bed, tmp_red)
  }
  final.bed <- sort(final.bed)


  # ----------------------------------------#
  #Output the rearranged fimo output and BED file
  samp.dir <- gsub("fimo.tsv", "", samp)
  samp.base <- gsub(".+\\/(.+)\\/", "\\1", samp.dir)

  write.table(fimo.bed, file = paste0(samp.dir, samp.base, "_fimo.tsv"), row.names = F, quote = F, sep = "\t")
  export(object = final.bed, con = paste0(samp.dir, samp.base, "_fimo.bed"), format = "bed")
}
#------------------------------------------------------------------------------#


# Get the arguments from cml or swarm script
args <- commandArgs(trailingOnly=TRUE)

# Run function with the arguments provided by cml
fimo2bed(args[1], args[2])


#----------------------#
#Test data
#samp <- "/Volumes/data/CutNRun/Mouse/Analysis/20200309_16-29/Meme/fimo/Analysis_200910/Nrl_P28_200910/fimo.tsv"
#matrix <- "/Volumes/data/Index/TRANSFAC_Pro_data/TFP_2017.3_data/dat/matrixMAPgene_Mm_v2.txt"
#----------------------#
