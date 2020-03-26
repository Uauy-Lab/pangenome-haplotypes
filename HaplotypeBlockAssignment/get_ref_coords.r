## Get the position of the BLAST blocks with respect to the specific reference assemblies, rather than chinese spring RefSeq
library(stringr)

## define functions here
get_ref_start <- function(data, ref){
  ref_start_data <- ref_gff_gene_only[ref_gff_gene_only$cs_id == data["start_transcript"], ]
  if (nrow(ref_start_data) > 1){
    ref_start <- ref_start_data[ref_start_data$score == max(ref_start_data$score), "start"]
    if (length(ref_start) > 1){
      temp <- ref_start_data[ref_start_data$score == max(ref_start_data$score), ]
      ref_start <- temp[temp$start == min(temp$start), "start"]
    }
  } else {
    ref_start <- ref_start_data[, "start"]
  }
  ref_start
}

get_ref_start_transcript <- function(data, ref){
  ref_start_transcript_data <- ref_gff_gene_only[ref_gff_gene_only$cs_id == data["start_transcript"],]
  if (nrow(ref_start_transcript_data) > 1){
    ref_start_transcript <- ref_start_transcript_data[ref_start_transcript_data$score == max(ref_start_transcript_data$score), "var_id"]
    if (length(ref_start_transcript) > 1){
      temp <- ref_start_transcript_data[ref_start_transcript_data$score == max(ref_start_transcript_data$score), ]
      ref_start_transcript <- temp[temp$start == min(temp$start), "var_id"]
    }
  } else {
    ref_start_transcript <- ref_start_transcript_data[, "var_id"]
  }
}

get_ref_end <- function(data, ref){
  ref_end_data <- ref_gff_gene_only[ref_gff_gene_only$cs_id == data["end_transcript"], ]
  if (nrow(ref_end_data) > 1){
    ref_end <- ref_end_data[ref_end_data$score == max(ref_end_data$score), "end"]
    if (length(ref_end) > 1){
      temp <- ref_end_data[ref_end_data$score == max(ref_end_data$score), ]
      ref_end <- temp[temp$start == min(temp$start), "end"]
    }
  } else {
    ref_end <- ref_end_data[, "end"]
  }
  ref_end
}

get_ref_end_transcript <- function(data, ref){
  ref_end_transcript_data <- ref_gff_gene_only[ref_gff_gene_only$cs_id == data["end_transcript"], ]
  if (nrow(ref_end_transcript_data) > 1){
    ref_end_transcript <- ref_end_transcript_data[ref_end_transcript_data$score == max(ref_end_transcript_data$score), "var_id"]
    if (length(ref_end_transcript) > 1){
      temp <- ref_end_transcript_data[ref_end_transcript_data$score == max(ref_end_transcript_data$score), ]
      ref_end_transcript <- temp[temp$start == min(temp$start), "var_id"]
    }
  } else {
    ref_end_transcript <- ref_end_transcript_data[, "var_id"]
  }
  ref_end_transcript
}

## end of functions
data_dir <- "X:/brintonj/haplotype/whole_genome_blast/blocks/"

chromosomes <- list.files(data_dir)

for (chr in chromosomes){
  BLAST_path <- paste0(data_dir, chr, "/haplotype_blocks_BLAST_25_gene_window_2000bp.txt")

  BLAST <- read.table(BLAST_path,
                    sep = "\t", 
                    header = TRUE,
                    stringsAsFactors = FALSE)

  head(BLAST)

  BLAST$ref <- str_split_fixed(BLAST$aln_type, "->", 2)[,1]
  BLAST$query <- str_split_fixed(BLAST$aln_type, "->", 2)[,2]

  #make the opposite ref query scenario for the blast
  BLAST_copy <- BLAST
  BLAST_copy$ref <- BLAST$query
  BLAST_copy$query <- BLAST$ref

  BLAST_recip <- rbind(BLAST, BLAST_copy)

  EI_vars <- c("cadenza", "claire", "paragon", "robigus", "weebil")
  EI_blocks <- subset(BLAST_recip, ref %in% EI_vars)

  references <- unique(BLAST_recip$ref)
  references <- references[!(references %in% EI_vars)]

  ref_gff_dir <- "W:/assemblies/releasePGSBv2.0/gff/"

  BLAST_coords <- data.frame(block_no=numeric(),
                           block_start = numeric(),
                           block_end = numeric(),
                           aln_type = character(),
                           start_transcript = character(),
                           end_transcript = character(),
                           window = numeric(),
                           ref = character(),
                           query = character(),
                           ref_start = numeric(),
                           ref_start_transcript = character(),
                           ref_end = numeric(),
                           ref_end_transcript = character())
  
  for (ref in references) {
  
    ref_gff_path <- paste0(ref_gff_dir, ref, ".gff")

    ref_gff <- read.table(ref_gff_path, 
                      sep = "\t", 
                      header = FALSE,
                      stringsAsFactors = FALSE)


    ref_gff_gene_only <- ref_gff[ref_gff$V3 == "gene",]

    ref_gff_gene_only <- cbind(ref_gff_gene_only, str_split_fixed(ref_gff_gene_only$V9, ";", 4))
    colnames(ref_gff_gene_only) <- c("chr", "annotation", "biotype", "start", "end", "score", "strand", "missing", "info", "ref_id", "cs_id", "frame", "note")

    ref_gff_gene_only$var_id <- gsub("ID=", "", ref_gff_gene_only$ref_id)
    ref_gff_gene_only$cs_id <- gsub("srcmodel=", "", ref_gff_gene_only$cs_id)
    ref_gff_gene_only$cs_id <- str_split_fixed(ref_gff_gene_only$cs_id, "\\.", 2)[,1]


    ## this bit gets the reference blocks
  ## we also need to make sure that we incorportate the EI varieties because we can't get their position (scaffolds) so we will use the query position for this too
    BLAST_ref <- BLAST_recip[BLAST_recip$ref == ref,]
    BLAST_ref <- rbind(BLAST_ref, EI_blocks[EI_blocks$query == ref,])
    BLAST_ref <- BLAST_ref[complete.cases(BLAST_ref),]

    BLAST_ref$ref_start <- apply(BLAST_ref, 1, FUN=get_ref_start, ref = ref)
    BLAST_ref$ref_start_transcript <- apply(BLAST_ref, 1, FUN=get_ref_start_transcript, ref = ref)
    BLAST_ref$ref_end <- apply(BLAST_ref, 1, FUN=get_ref_end, ref = ref)
    BLAST_ref$ref_end_transcript <- apply(BLAST_ref, 1, FUN=get_ref_end_transcript, ref = ref)

    BLAST_coords <- rbind(BLAST_coords, BLAST_ref)
  }

  BLAST_coords$block_no <- as.character(BLAST_coords$block_no)
  BLAST_coords$block_start <- as.numeric(as.character(BLAST_coords$block_start))
  BLAST_coords$block_end <- as.numeric(as.character(BLAST_coords$block_end))
  BLAST_coords$aln_type <- as.character(BLAST_coords$aln_type)
  BLAST_coords$start_transcript <- as.character(BLAST_coords$start_transcript)
  BLAST_coords$end_transcript <- as.character(BLAST_coords$end_transcript)
  BLAST_coords$window <- as.numeric(as.character(BLAST_coords$window))
  BLAST_coords$ref <- as.character(BLAST_coords$ref)
  BLAST_coords$query <- as.character(BLAST_coords$query)
  BLAST_coords$ref_start <- as.numeric(as.character(BLAST_coords$ref_start))
  BLAST_coords$ref_start_transcript <- as.character(BLAST_coords$ref_start_transcript)
  BLAST_coords$ref_end <- as.numeric(as.character(BLAST_coords$ref_end))
  BLAST_coords$ref_end_transcript <- as.character(BLAST_coords$ref_end_transcript)
  
  write.table(BLAST_coords,
              file = paste0(data_dir, chr, "/haplotype_blocks_BLAST_25_gene_window_2000bp_ref_coords.txt"),
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
}

## Unfortunately this misses out EI EI comparisons, so we need to add them back
for (chr in chromosomes){
  BLAST_path <- paste0(data_dir, chr, "/haplotype_blocks_BLAST_25_gene_window_2000bp.txt")
  
  BLAST <- read.table(BLAST_path,
                      sep = "\t", 
                      header = TRUE,
                      stringsAsFactors = FALSE)

  BLAST_coords_path <- paste0(data_dir, chr, "/haplotype_blocks_BLAST_25_gene_window_2000bp_ref_coords.txt")
  
  BLAST_coords <- read.table(BLAST_coords_path,
                      sep = "\t", 
                      header = TRUE,
                      stringsAsFactors = FALSE)
  
  BLAST$ref <- str_split_fixed(BLAST$aln_type, "->", 2)[,1]
  BLAST$query <- str_split_fixed(BLAST$aln_type, "->", 2)[,2]

  EI_BLAST <- BLAST[BLAST$ref %in% EI_vars & BLAST$query %in% EI_vars,]  
  
  EI_BLAST$ref_start <- EI_BLAST$block_start
  EI_BLAST$ref_start_transcript <- EI_BLAST$start_transcript
  EI_BLAST$ref_end <- EI_BLAST$block_end
  EI_BLAST$ref_end_transcript <- EI_BLAST$end_transcript
  
  EI_recip <- EI_BLAST
  EI_recip$ref <- EI_BLAST$query
  EI_recip$query <- EI_BLAST$ref
  EI_recip <- rbind(EI_BLAST, EI_recip)
  
  BLAST_coords <- rbind(BLAST_coords, EI_recip)
  
  write.table(BLAST_coords,
              file = paste0(data_dir, chr, "/haplotype_blocks_BLAST_25_gene_window_2000bp_ref_coords.txt"),
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
}
