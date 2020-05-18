library(stringr)
library(plyr)
library(reshape2)

## functions

read_pairwise_position <- function(blast_path_gz, gtf, outfile){
  all_comp <- read.table(gzfile(blast_path_gz), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  head(all_comp)
  
  #get the refseq position of the genes
  gtf$transcript_id <- str_split_fixed(gtf$V9, ";", 2)[,1]
  gtf$transcript_id <- gsub("transcript_id ", "", gtf$transcript_id)
  
  gene_only <- gtf[gtf$V3 == "gene",]
  
  gene_positions <- gene_only[,c("transcript_id", "V1", "V4",  "V5", "V7")]
  colnames(gene_positions) <- c("transcript", "chr", "start", "end", "strand")
  gene_positions$transcript <- gsub("ID=", "", gene_positions$transcript)
  
  #now add positions to the pairwise comparison file
  
  all_comp_positions <- merge(all_comp, gene_positions, all.x = TRUE, all.y = FALSE)
  head(all_comp_positions)
  
  write.table(all_comp_positions, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(all_comp_positions)
}

separate_output_identities <- function(comp_data, outdir, suffix){
  varieties <- unique(comp_data$var_query)
  for (i in seq(1, length(varieties))){
    variety <- varieties[i]
    print(variety)
    variety_data <- comp_data[grep(variety, comp_data$aln_type),]
    write.table(variety_data, file = paste0(outdir, variety, suffix, ".tab"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}

read_pairwise_position_no_header <- function(blast_path_gz, gtf, outfile, col_names){
  all_comp <- read.table(gzfile(blast_path_gz), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(all_comp) <- col_names
  head(all_comp)
  
  #get the refseq position of the genes
  gtf$transcript_id <- str_split_fixed(gtf$V9, ";", 2)[,1]
  gtf$transcript_id <- gsub("transcript_id ", "", gtf$transcript_id)
  
  gene_only <- gtf[gtf$V3 == "gene",]
  
  gene_positions <- gene_only[,c("transcript_id", "V1", "V4",  "V5", "V7")]
  colnames(gene_positions) <- c("transcript", "chr", "start", "end", "strand")
  gene_positions$transcript <- gsub("ID=", "", gene_positions$transcript)
  
  #now add positions to the pairwise comparison file
  
  all_comp_positions <- merge(all_comp, gene_positions, all.x = TRUE, all.y = FALSE)
  head(all_comp_positions)
  
  write.table(all_comp_positions, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  return(all_comp_positions)
}

separate_output_identities <- function(comp_data, outdir, suffix){
  varieties <- unique(comp_data$var_query)
  for (i in seq(1, length(varieties))){
    variety <- varieties[i]
    print(variety)
    variety_data <- comp_data[grep(variety, comp_data$aln_type),]
    write.table(variety_data, file = paste0(outdir, variety, suffix, ".tab"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}
calculate_pid_windows <- function(aln_data, window_size = 20, drop_no = 2){
  aln_data <- aln_data[order(aln_data$start),]
  aln <- unique(aln_data$aln_type)
  mean_pidents <- data.frame(start = numeric(), end = numeric(), start_position = numeric(), end_position = numeric(), pident_mean = numeric(), aln_type = character()) 
  
  for(j in seq(1, nrow(aln_data)-window_size)){
    window_start = j
    start_position = aln_data[window_start, "start"]
    window_end = (window_start + window_size) - 1
    end_position = aln_data[window_end, "end"]
    pidents <- aln_data[c(window_start:window_end), "pident"]
    pidents <- pidents[order(pidents, decreasing = TRUE)]
    pidents_sub <- pidents[c(1:(length(pidents)-drop_no))]
    pident_mean <- mean(pidents_sub)
    aln_to_add = data.frame(start = window_start, 
                            end = window_end, 
                            start_position = start_position, 
                            end_position = end_position, 
                            pident_mean = pident_mean, 
                            aln_type = aln)
    mean_pidents <- rbind(mean_pidents, aln_to_add)
  }
  return(mean_pidents)
}


assign_blocks <-function(mean_pidents){
  block_positions <- data.frame(block_no = numeric(), block_start = numeric(), block_end = numeric(), aln_type = character())
  
  mean_pidents_copy <- mean_pidents
  mean_pidents_copy$block_no <- NA
  block_no = 1
  new_block <- "no"
  
  for (i in seq(1, nrow(mean_pidents_copy))){
    print(i)
    if(mean_pidents_copy[i, "pident_mean"] < 100){
      print("less than 100%")
      mean_pidents_copy[i, "block_no"] <- NA
    } else if (mean_pidents_copy[i, "pident_mean"] == 100){
      if(new_block == "no"){
        mean_pidents_copy[i, "block_no"] <- block_no
      } else if (new_block == "yes"){
        if(mean_pidents_copy[i, "start_position"] < end_prev_block){
          block_no <- block_no-1
          mean_pidents_copy[i, "block_no"] <- block_no
        } else {
          mean_pidents_copy[i, "block_no"] <- block_no
        }
      }
      if(i == nrow(mean_pidents_copy)){
        print("end")
      } else if (mean_pidents_copy[i+1, "pident_mean"] == 100){
        block_no <- block_no
        new_block <- "no"
        print("same block")
      } else if (mean_pidents_copy[i+1, "pident_mean"] < 100){
        block_no <- block_no+1
        new_block <- "yes"
        print("new_block")
        end_prev_block <- mean_pidents_copy[i, "end_position"]
      }
    }
  }
  return(mean_pidents_copy)
}

summarise_blocks <- function(mean_pidents_copy){
  #now get the start and end of each block
  blocks <- unique(mean_pidents_copy$block_no)
  blocks <- blocks[complete.cases(blocks)]
  
  block_positions <- data.frame(block_no = numeric(), block_start = numeric(), block_end = numeric(), aln_type = character())
  
  for(block in blocks){
    block_data <- subset(mean_pidents_copy, block_no == block)
    block_start <- min(block_data$start_position)
    block_end <- max(block_data$end_position)
    aln_type <- unique(block_data$aln_type)
    to_add <- data.frame(block_no = block, block_start = block_start, block_end = block_end, aln_type = aln_type)
    block_positions <- rbind(block_positions, to_add)
  }
  return(block_positions)
}

calculate_blocks <- function (data_dir, window_sizes, gene_no_window, drop_no, chrom = "6A"){
  
  overall_blocks <- data.frame(block_no = numeric(), block_start = numeric(), block_end = numeric(), aln_type = character(), window = character())
  
  for(window in window_sizes){
    print(window)
    window_dir <- paste0(data_dir, window)
    id_data_path <- paste0(window_dir, "/varieties_", chrom, "_identities_", window, "_refseq_positions.tab")
    id_data <- read.table(id_data_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    id_data_noN <- id_data[id_data$Ns_total == 0,]
    
    window_blocks <- data.frame(block_no = numeric(), block_start = numeric(), block_end = numeric(), aln_type = character(), window_size = character())
    
    alns <- unique(id_data_noN$aln_type)
    for(aln in alns){
      print(aln)
      aln_data <- id_data_noN[id_data_noN$aln_type == aln,]
      #block_info <- calculate_hap_blocks(aln_data = aln_data)
      
      aln_mean_pidents <- calculate_pid_windows(aln_data = aln_data, window_size = gene_no_window, drop_no = drop_no)
      aln_mean_pidents <- aln_mean_pidents[complete.cases(aln_mean_pidents),]
      if(max(aln_mean_pidents$pident_mean) < 100){
        block_info <- data.frame(block_no = NA, block_start = NA, block_end = NA, aln_type = aln, window_size = window)
      } else {
        aln_mean_pidents_copy <- assign_blocks(mean_pidents = aln_mean_pidents)
        block_info <- summarise_blocks(mean_pidents_copy = aln_mean_pidents_copy)
        block_info$window_size <- window
      }
      window_blocks <- rbind(window_blocks, block_info)
    }
    overall_blocks <- rbind(overall_blocks, window_blocks)
  }
  outfile <- paste0(data_dir, "haplotype_blocks_aln_", gene_no_window, "_gene_window_drop", drop_no, ".tsv")
  write.table(overall_blocks, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  return(overall_blocks)
}

## End of functions

## calculate the blocks for 1D, 2B, 3B, 4D, 6A, 7A 

raw_data_dir <- "W:/ramirezr/SM1/final_pairwise_blast/"
out_dir <- "X:/brintonj/haplotype/whole_genome_blast/precision_recall/varying_window_blocks/"

pairwise_header <- c("transcript", "query", "subject", "var_query", "var_subject", "aln_type", "length", "pident", "Ns_query", "Ns_subject", "Ns_total", "Flanking")
#same gtf for all (using RefSeq order for now)
HC_gtf <- read.table("W://WGAv1.0//annotation//IWGSC_v1.1_HC_20170706.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
LC_gtf <- read.table("W://WGAv1.0//annotation//IWGSC_v1.1_LC_20170706.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ALL_gtf <- rbind(HC_gtf, LC_gtf)

chromosomes_to_calculate <- c("1D", "2B", "3B", "4D", "6A", "7A")

for (chrom in chromosomes_to_calculate){
  files_to_run <- list.files(raw_data_dir)[grep(paste0("varieties_", chrom, "_identities"), list.files(raw_data_dir))]
  
  for (i in seq(1, length(files_to_run))){
    file <- files_to_run[i]
    window <- gsub(paste0("varieties_", chrom, "_identities_"), "", file)
    window <- gsub("\\.tab.gz", "", window)
    window_dir <- paste0(out_dir, window, "/")
    
    dir.create(window_dir)
   
    all_comp_positions <- read_pairwise_position_no_header(blast_path_gz = paste0(raw_data_dir, file), gtf = ALL_gtf, outfile = paste0(window_dir, "/varieties_", chrom, "_identities_", window, "_refseq_positions.tab"), col_names = pairwise_header)
    
    #separate_output_identities(comp_data = all_comp_positions, outdir = window_dir, suffix = paste0( "_", chrom, "_identities_refseq_positions_", window))
  }
  
  ## ok now we need to call the blocks for all the different window sizes and criteria
  
  window_sizes <- c("cds", "0bp", "1000bp", "2000bp", "5000bp")
  
  #OK now we want to calculate the blocks for all the different comparisons for the different window sizes
  
  gene_no_window <- 10
  
  window_10_blocks <- calculate_blocks(data_dir = out_dir, 
                                       window_sizes = window_sizes, 
                                       gene_no_window = gene_no_window, 
                                       drop_no = ceiling(gene_no_window*0.1),
                                       chrom = chrom)
  
  window_10_blocks$block_size <- window_10_blocks$block_end-window_10_blocks$block_start
  window_10_blocks$chr <- paste0("chr", chrom)
  window_10_blocks$gene_block <- gene_no_window
  
  gene_no_window <- 15
  
  window_15_blocks <- calculate_blocks(data_dir = out_dir, 
                                       window_sizes = window_sizes, 
                                       gene_no_window = gene_no_window, 
                                       drop_no = ceiling(gene_no_window*0.1),
                                       chrom = chrom)
  
  window_15_blocks$block_size <- window_15_blocks$block_end-window_15_blocks$block_start
  window_15_blocks$chr <- paste0("chr", chrom)
  window_15_blocks$gene_block <- gene_no_window
  
  gene_no_window <- 25
  
  window_25_blocks <- calculate_blocks(data_dir = out_dir, 
                                       window_sizes = window_sizes, 
                                       gene_no_window = gene_no_window, 
                                       drop_no = ceiling(gene_no_window*0.1),
                                       chrom = chrom)
  
  window_25_blocks$block_size <- window_25_blocks$block_end-window_25_blocks$block_start
  window_25_blocks$chr <- paste0("chr", chrom)
  window_25_blocks$gene_block <- gene_no_window
  
  gene_no_window <- 30
  
  window_30_blocks <- calculate_blocks(data_dir = out_dir, 
                                       window_sizes = window_sizes, 
                                       gene_no_window = gene_no_window, 
                                       drop_no = ceiling(gene_no_window*0.1),
                                       chrom = chrom)
  
  window_30_blocks$block_size <- window_30_blocks$block_end-window_30_blocks$block_start
  window_30_blocks$chr <- paste0("chr", chrom)
  window_30_blocks$gene_block <- gene_no_window
  
  gene_no_window <- 20
  
  window_20_blocks <- calculate_blocks(data_dir = out_dir, 
                                       window_sizes = window_sizes, 
                                       gene_no_window = gene_no_window, 
                                       drop_no = ceiling(gene_no_window*0.1),
                                       chrom = chrom)
  
  window_20_blocks$block_size <- window_20_blocks$block_end-window_20_blocks$block_start
  window_20_blocks$chr <- paste0("chr", chrom)
  window_20_blocks$gene_block <- gene_no_window
  
  all_window_blocks <- rbind(window_10_blocks, window_15_blocks, window_20_blocks, window_25_blocks, window_30_blocks)
  
  write.table(all_window_blocks, file = paste0(out_dir, "combined_", chrom, "_all_flanking_all_windown_sizes_BLAST_blocks.tsv"),
              sep = "\t", 
              col.names = TRUE, 
              row.names = FALSE, 
              quote = FALSE)
}

#####
