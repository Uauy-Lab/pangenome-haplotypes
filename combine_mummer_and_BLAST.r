## Script written by Jemima Brinton 2019
## Aim - combine haplotype blocks from nucmer and BLAST alignments. First use mummer/nucmer blocks, then check BLAST blocks - any not represented by the mummer blocks and >= bin size should be added to the final set 

library(GenomicRanges)
##functions

check_block <- function(data, EI_vars, mummer_blocks){
  #print(data)
  if ((data["ref"] %in% EI_vars) & (data["query"] %in% EI_vars)){
    count_overlaps <- 0
  } else {
    if ((data["ref"] %in% EI_vars)){
      mummer_gdf <- mummer_blocks[mummer_blocks$query == data["ref"] & mummer_blocks$ref == data["query"], c("chrom", "block_start", "block_end")]
    } else {
      mummer_gdf <- mummer_blocks[mummer_blocks$query == data["query"] & mummer_blocks$ref == data["ref"], c("chrom", "block_start", "block_end")]
    }
    if(nrow(mummer_gdf)>0){
      colnames(mummer_gdf) <- c("chr", "start", "end")
      #print(mummer_gdf)
      mummer_gr <- makeGRangesFromDataFrame(mummer_gdf)
    
      to_check_gdf <- data.frame(chr = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)
      to_check_gdf[1,"chr"] <- data["chrom"]
      to_check_gdf[1,"start"] <- data[("ref_start")]
      to_check_gdf[1,"end"] <- data[("ref_end")]
      #print(to_check_gdf)
    
      to_check_gr <- makeGRangesFromDataFrame(to_check_gdf)
    
      overlaps <- findOverlaps(to_check_gr, mummer_gr)
      count_overlaps <- length(unique(subjectHits(overlaps)))
    } else {
      count_overlaps <- 0
    }
  }
  return(count_overlaps)
  #print("finish")
}

###

mummer_dir <- "X:/brintonj/haplotype/whole_genome_mummer/blocks/"
BLAST_dir <- "X:/brintonj/haplotype/whole_genome_blast/blocks/"

chromosomes <- list.files(mummer_dir)

bin_size <- 5000000

combined_blocks <- data.frame(ref = character(),
                              query = character(),
                              chrom = character(),
                              ref_start = numeric(),
                              ref_end = numeric(),
                              ref_start_transcript = character(),
                              ref_end_transcript = character(),
                              cs_start = numeric(),
                              cs_end = numeric(),
                              cs_start_transcript = character(),
                              cs_end_transcript = character(),
                              source = character())

for (chr in chromosomes){
  mummer_blocks <- read.table(paste0(mummer_dir, chr, "/mummer_blocks_chr", chr, ".min20000.5Mb_bins.txt"),
                              header = TRUE, 
                              sep = "\t",
                              stringsAsFactors = FALSE)
  
  BLAST_blocks <- read.table(paste0(BLAST_dir, "chr", chr, "/haplotype_blocks_BLAST_25_gene_window_2000bp_ref_coords.txt"),
                             header = TRUE, 
                             sep = "\t",
                             stringsAsFactors = FALSE)
  
  #filter the BLAST for blocks only above the mummer bin size
  #remove spelta
  BLAST_blocks$chrom <- paste0("chr", chr)
  BLAST_blocks <- BLAST_blocks[!(BLAST_blocks$ref == "spelta" | BLAST_blocks$query == "spelta"),]
  BLAST_blocks$block_size <- BLAST_blocks$ref_end - BLAST_blocks$ref_start
  BLAST_bin <- BLAST_blocks[BLAST_blocks$block_size >= bin_size,]
  BLAST_bin <- BLAST_bin[complete.cases(BLAST_bin$block_start),]

  EI_vars <- c("cadenza", "claire", "paragon", "robigus", "weebil")
  BLAST_bin$overlap <- apply(BLAST_bin, 1, check_block, EI_vars = EI_vars, mummer_blocks = mummer_blocks)
  BLAST_bin$source <- "BLAST"
  BLAST_bin$in_reciprocal <- NA
  
  mummer_blocks$source <- "mummer"
  
  BLAST_to_add <- BLAST_bin[BLAST_bin$overlap == 0, c("ref", 
                                                      "query", 
                                                      "chrom", 
                                                      "ref_start", 
                                                      "ref_end", 
                                                      "ref_start_transcript", 
                                                      "ref_end_transcript",
                                                      "block_start", 
                                                      "block_end",
                                                      "start_transcript",
                                                      "end_transcript",
                                                      "source")]
  
  colnames(BLAST_to_add)[c(8:11)] <- c("cs_start", "cs_end", "cs_start_transcript", "cs_end_transcript")
  
  mummer_blocks$cs_start <- NA
  mummer_blocks$cs_end <- NA
  mummer_blocks$cs_start_transcript <- NA
  mummer_blocks$cs_end_transcript <- NA
  mummer_blocks$ref_start_transcript <- NA
  mummer_blocks$ref_end_transcript <- NA
  
  mummer_to_add <- mummer_blocks[, c("ref", 
                                     "query", 
                                     "chrom", 
                                     "block_start", 
                                     "block_end", 
                                     "ref_start_transcript", 
                                     "ref_end_transcript",
                                     "cs_start", 
                                     "cs_end",
                                     "cs_start_transcript",
                                     "cs_end_transcript",
                                     "source")]
colnames(mummer_to_add)[c(4,5)] <- c("ref_start", "ref_end")

mummer_BLAST <- rbind(mummer_to_add, BLAST_to_add)

combined_blocks <- rbind(combined_blocks, mummer_BLAST)
}

write.table(combined_blocks, file = "X:/brintonj/haplotype/whole_genome_mummer_BLAST_5mbp_blocks_combined.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

ref_coords_only <- combined_blocks[,c(1:5, 12)]

write.table(ref_coords_only, file = "X:/brintonj/haplotype/whole_genome_mummer_BLAST_5mbp_blocks_combined_ref_coords.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## we also want to format the dataframe with the blocks being given a unique identifier
## columns we want: assembly chromosome start end block_no assembly_coords

EI_vars <- c("cadenza", "claire", "paragon", "robigus", "weebil")

ref_coords_only_copy <- ref_coords_only

ref_coords_only_copy$block_no <- NA
ref_coords_only_copy$assembly_coords <- NA

ref_coords_only_copy$block_size <- ref_coords_only_copy$ref_end - ref_coords_only_copy$ref_start
ref_coords_only_copy <- ref_coords_only_copy[order(ref_coords_only_copy$block_size, decreasing = TRUE),]

block_no <- 1


for(i in seq(1, nrow(ref_coords_only_copy))){
  data <- ref_coords_only_copy[i,]
  
  if(is.na(ref_coords_only_copy[i, "block_no"])){
    ref <- ref_coords_only_copy[i,"ref"]
    query <- ref_coords_only_copy[i,"query"]
    chrom <- ref_coords_only_copy[i,"chrom"]
    
    ref_coords_only_copy[i,"block_no"] <- block_no
    
    matching_rows <- which(((ref_coords_only_copy$ref == query) &
                              (ref_coords_only_copy$query == ref) &
                             (ref_coords_only_copy$chrom == chrom)) & (((ref_coords_only_copy$ref_start >= data$ref_start & ref_coords_only_copy$ref_start <= data$ref_end) |
                            (ref_coords_only_copy$ref_end >= data$ref_start & ref_coords_only_copy$ref_end <= data$ref_end)) | 
                            ((data$ref_start >= ref_coords_only_copy$ref_start & data$ref_start <= ref_coords_only_copy$ref_end) |
                            (data$ref_end >= ref_coords_only_copy$ref_start & data$ref_end <= ref_coords_only_copy$ref_end))))

    
    if(length(matching_rows) > 0){
      ref_coords_only_copy[matching_rows, "block_no"] <- block_no
    } else if (query %in% EI_vars) {
      last_line <- nrow(ref_coords_only_copy)
      ref_coords_only_copy <- rbind(ref_coords_only_copy, ref_coords_only_copy[i,])
      ref_coords_only_copy[last_line+1,"block_no"] <- block_no
      ref_coords_only_copy[last_line+1,"ref"] <- query
      ref_coords_only_copy[last_line+1,"query"] <- ref
      ref_coords_only_copy[last_line+1,"assembly_coords"] <- ref
    }
    if(query %in% EI_vars & ref %in% EI_vars){
      ref_coords_only_copy[i,"assembly_coords"] <- "IWGSCv1.1"
      ref_coords_only_copy[matching_rows,"assembly_coords"] <- "IWGSCv1.1"
    } else if (query %in% EI_vars) {
      ref_coords_only_copy[i,"assembly_coords"] <- ref
      ref_coords_only_copy[matching_rows,"assembly_coords"] <- ref
    } else if (ref %in% EI_vars) {
      ref_coords_only_copy[i,"assembly_coords"] <- query
      ref_coords_only_copy[matching_rows,"assembly_coords"] <- query
    } else {
      ref_coords_only_copy[i,"assembly_coords"] <- ref
      ref_coords_only_copy[matching_rows,"assembly_coords"] <- query
    }
    block_no <- block_no + 1
  } else {
    print("sorted")
  }
}

ref_coords_only_copy <- ref_coords_only_copy[order(ref_coords_only_copy$chrom, ref_coords_only_copy$ref, ref_coords_only_copy$ref_start),]

#combine with the name of the ref chroms

ref_suffix <- read.table("X:/brintonj/haplotype/assemblies_suffix.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ref_coords_only_copy$assembly_chrom <- ref_coords_only_copy$chrom

for (i in seq(1, nrow(ref_coords_only_copy))){
  assembly <- ref_coords_only_copy[i, "assembly_coords"]
  chrom <- ref_coords_only_copy[i, "chrom"]
  if (assembly == "IWGSCv1.1"){
  } else {
    ref_coords_only_copy[i, "assembly_chrom"] <- paste0(chrom, "__", ref_suffix[ref_suffix$name == assembly, "suffix"])
  }
}

colnames(ref_coords_only_copy)[8] <- "ref_assembly"
colnames(ref_coords_only_copy)[10] <- "ref_chrom"

write.table(ref_coords_only_copy, file = "X:/brintonj/haplotype/whole_genome_mummer_BLAST_5mbp_blocks_combined_ref_coords_block_numbers.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

