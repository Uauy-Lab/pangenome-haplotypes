## Script written by Jemima Brinton 2019
## Aim - define haplotype blocks across whole wheat genome based on pairwise nucmer alignments between genome assemblies

###FUNCTIONS
library(dplyr)
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(viridis)
library(stringr)


#Read Delta functions and plotting functions (adapted) used are from https://jmonlong.github.io/Hippocamplus/2017/09/19/mummerplots-with-ggplot2/
#Reads in delta file output from nucmer into R and summarises into dataframe
readDelta <- function(deltafile){
  lines = scan(deltafile, 'a', sep='\n', quiet=TRUE)
  lines = lines[-1]
  lines.l = strsplit(lines, ' ')
  lines.len = lapply(lines.l, length) %>% as.numeric
  lines.l = lines.l[lines.len != 1]
  lines.len = lines.len[lines.len != 1]
  head.pos = which(lines.len == 4)
  head.id = rep(head.pos, c(head.pos[-1], length(lines.l)+1)-head.pos)
  mat = matrix(as.numeric(unlist(lines.l[lines.len==7])), 7)
  res = as.data.frame(t(mat[1:5,]))
  colnames(res) = c('rs','re','qs','qe','error')
  res$qid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 2))
  res$rid = unlist(lapply(lines.l[head.id[lines.len==7]], '[', 1)) %>% gsub('^>', '', .)
  res$strand = ifelse(res$qe-res$qs > 0, '+', '-')
  res
}

#calculates percentage identity and mid point for each mummer alignment
calculate_perc_id_mid_points <- function(data){
  data$r_length <- (data$re - data$rs)
  data$perc_id <- ((data$r_length - data$error)/data$r_length)*100
  data$perc_id_factor <- data$perc_id
  data[data$perc_id < 100, "perc_id_factor"] <- "<100"
  data$perc_id_factor <- as.factor(data$perc_id_factor)
  data$r_mid <- (data$rs + data$re)/2
  data$q_mid <- (data$qs + data$qe)/2
  return(data)
}

#reads in delta file, calculates percentage ID and filters for a minimum alignment size
pre_plot_analysis <- function(delta_path, min_size = 20000){
  data = readDelta(delta_path)
  data <- calculate_perc_id_mid_points(data)
  data_filt <- data[data$r_length >= min_size,]
  return(data_filt)
}

# Plots diagonal scatter plot of a nucmer pairwise alignment with points coloured by percentage identity of alignment - x axis is reference, y axis is query
plot_by_perc_id_cap <- function(data, xmin = 0, xmax = max(data$re), cap_lower, cap_upper){
  ggplot(data[data$rs > xmin & data$re < xmax,], aes(x=rs, xend=re, y=qs, yend=qe, colour=perc_id, shape = strand)) + geom_segment() +
    geom_point(alpha=.5) +
    theme_bw() + 
    theme(strip.text.y=element_text(angle=180, size=5),
          strip.background=element_blank()) +
    xlab(unique(data$rid)) + ylab(unique(data$qid)) +
    scale_colour_viridis(limits = c(cap_lower, cap_upper))
}

#scatter plot against reference position, y axis is percentage identity colour is size of alignment
plot_perc_id_v_ref <- function(data, xmin = 0, xmax = max(data$re), ymin = 97, ymax = 100){
  ggplot(data[data$r_mid > xmin & data$r_mid < xmax,], aes(x=r_mid, y=perc_id, colour=r_length)) +
    theme_bw() + xlab(data$rid) + ylab(paste0('percentage ID v ', data$qid)) +
    geom_point(alpha=.5) +
    ylim(ymin, ymax) +
    scale_colour_viridis()
}

#bins alignment dataframe by reference position
bin_data <- function(data, bin_size, max_chrom_size){
  bins <- seq(0, max_chrom_size, by = bin_size)
  data$bin <- NA
  for (i in bins){
    data$bin <- ifelse(((data$r_mid > (i-bin_size)) & (data$r_mid < (i-1))), i, data$bin)
  }
  return(data)
}

#plots median percentage identity of each bin against reference position, colour of point is based on if the median percentage identity passes the identity threshold
plot_medians_line_colour_cut_off <- function(data, bin_size = 10000000, cut_off = 99.99, ymin = 97, ymax = 100, max_chrom_size){
  comparison_filt_bin <- bin_data(data, bin_size = bin_size, max_chrom_size = max_chrom_size)
  comparison_medians <- aggregate(perc_id ~ bin, data = comparison_filt_bin, FUN=median)
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- NA
  comparison_medians$cut_off <- ifelse(comparison_medians$perc_id >= cut_off, paste0(">=",cut_off), paste0("<",cut_off))
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)
  
  ggplot(comparison_medians, aes(x=bin, y = perc_id_median, colour = cut_off)) +
    geom_line(colour = "grey") + 
    geom_point(size = 1)  + 
    ylim(ymin,ymax) +
    scale_colour_manual(values = c("#73D055FF", "#440154FF")) +
    labs(colour = "% id") +
    #scale_colour_viridis() +
    xlab(data$rid) +
    
    ggtitle(paste0(data$rid, " v ", data$qid, " BinSize: ", bin_size, " CutOff: ", cut_off))
}

#plots boxplots of percentage identity of each bin against reference position, colour of boxplot is based on if the median percentage identity passes the identity threshold
plot_boxplots_median_colour_cut_off <- function(data, bin_size = 10000000, cut_off = 99.99, ymin = 97, ymax = 100, max_chrom_size){
  comparison_filt_bin <- bin_data(data, bin_size = bin_size, max_chrom_size = max_chrom_size)
  comparison_medians <- aggregate(perc_id ~ bin, data = comparison_filt_bin, FUN=median)
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- NA
  comparison_medians$cut_off <- ifelse(comparison_medians$perc_id >= cut_off, paste0(">=",cut_off), paste0("<",cut_off))
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)
  
  comparison_to_plot <- merge(comparison_filt_bin, comparison_medians)
  
  ggplot(comparison_to_plot, aes(x=bin, y = perc_id, group = bin, fill = cut_off)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(ymin,ymax) +
    scale_fill_manual(values = c("#73D055FF", "#440154FF")) +
    labs(fill = "% id") +
    #scale_colour_viridis() +
    xlab(data$rid) +
    ggtitle(paste0(data$rid, " v ", data$qid, " BinSize: ", bin_size, " CutOff: ", cut_off))
}

#assigns haplotype blocks based on bins of alignments exceeding the 99.99 % id threshold
assign_blocks_mummer <-function(median_cutoffs){
  median_cutoffs_copy <- median_cutoffs
  median_cutoffs_copy$block_no <- NA
  block_no = 1
  
  for (i in seq(1, nrow(median_cutoffs_copy))){
    print(i)
    if(median_cutoffs_copy[i, "perc_id_median"] < 99.99){
      median_cutoffs_copy[i, "block_no"] <- NA
    } else if (median_cutoffs_copy[i, "perc_id_median"] >= 99.99){
      median_cutoffs_copy[i, "block_no"] <- block_no
      if (i > (nrow(median_cutoffs_copy)-3)){
        print("coming to end")
      } else if ((median_cutoffs_copy[i+1, "perc_id_median"] < 99.99) & (median_cutoffs_copy[i+2, "perc_id_median"] < 99.99) & (median_cutoffs_copy[i+3, "perc_id_median"] < 99.99)){
        block_no <- block_no + 1
      }
    }
  }
  return(median_cutoffs_copy)
}

#summarises adjacent bins into haplotype blocks
summarise_blocks <- function(median_cutoffs_copy){
  #now get the start and end of each block
  blocks <- unique(median_cutoffs_copy$block_no)
  blocks <- blocks[complete.cases(blocks)]
  
  block_positions <- data.frame(block_no = numeric(), block_start = numeric(), block_end = numeric(), ref = character(), query = character())
  
  for(block in blocks){
    block_data <- subset(median_cutoffs_copy, block_no == block)
    block_start <- min(block_data$bin) - 5000000
    block_end <- max(block_data$bin)
    ref <- unique(block_data$ref)
    query <- unique(block_data$query)
    to_add <- data.frame(block_no = block, block_start = block_start, block_end = block_end, ref = ref, query = query)
    block_positions <- rbind(block_positions, to_add)
  }
  return(block_positions)
}

#checks if haplotype blocks are also called in the reciprocal mummer alignment
check_reciprocal <- function(data, all_ref_query_coords, chrom){
  ref <- data["ref"]
  query <- data["query"]
  
  ref_vars <- as.character(unique(all_ref_query_coords$ref))
  
  if(query %in% ref_vars){
    recip_table <- all_ref_query_coords[(all_ref_query_coords$ref == query) & (all_ref_query_coords$query == ref), ]
    if (nrow(recip_table) == 0){
      in_reciprocal <- "N"
    } else {
      
      ref_gdf <- data.frame(chr=character(), start=numeric(), end=numeric())
      to_add <- c(data["chrom"], data["block_start"], data["block_end"])
      ref_gdf <- rbind(ref_gdf, to_add)
      colnames(ref_gdf) <- c("chr", "start", "end")
      
      
      query_gdf <- recip_table[, c("chrom", "block_start", "block_end")]
      colnames(query_gdf) <- c("chr", "start", "end")
      
      # make genomic ranges objects for each
      ref_gr <- makeGRangesFromDataFrame(ref_gdf)
      query_gr <- makeGRangesFromDataFrame(query_gdf)
      
      overlaps <- findOverlaps(ref_gr, query_gr)
      #precision first - how many BLAST blocks found in the  mummer haplotype?
      count_overlaps <- length(unique(subjectHits(overlaps)))
      if (count_overlaps > 0){
        in_reciprocal <- "Y"
        
      } else {
        in_reciprocal <- "N"
      }
    }
  } else {
    in_reciprocal <- NA
  }
  return(in_reciprocal)
}

##END OF FUNCTIONS

## plots the median cut offs and call the blocks for each chromosome
base_dir <- "X:/brintonj/haplotype/whole_genome_mummer/"

plot_base <- paste0(base_dir, "plots")
blocks_base <- paste0(base_dir, "blocks")

dir.create(plot_base)
dir.create(blocks_base)

#read in the fasta indexes so we can set the max chromosome length
chrom_lengths <- read.table("W:/assemblies/releasePGSBv2.0/genome/combined_fai.txt", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
#remove spelta and synthetic from this table as we are not including them
chrom_lengths <- chrom_lengths[!grepl("ash", chrom_lengths[,1]),]
chrom_lengths <- chrom_lengths[!grepl("tsp", chrom_lengths[,1]),]

#get a list of all the varieties to plot as references (chromosome level as references only)
varieties_to_plot <- list.files(paste0(base_dir, "aln"))
chromosomes_to_plot <- list.files(paste0(base_dir, "aln/", varieties_to_plot[1]))

#we will only use the 20kb filter here (this is already filtered using mummer)
min_size <- 20000
cut_off <- 99.99
chromosomes_to_plot <- chromosomes_to_plot[!(chromosomes_to_plot == "1A")]

for (i in seq(1, length(chromosomes_to_plot))){
  chrom <- paste0("chr", chromosomes_to_plot[i])
  print(chrom)
  max_chrom_size <- max(chrom_lengths[grepl(chrom, chrom_lengths[,1]),2])
  
  blocks_chrom_dir <- paste0(blocks_base, "/", chromosomes_to_plot[i])
  plots_chrom_dir <- paste0(plot_base, "/", chromosomes_to_plot[i])
  
  dir.create(blocks_chrom_dir)
  dir.create(plots_chrom_dir)
  
  chrom_blocks_coords <- data.frame(block_no = numeric(), 
                                    block_start = numeric(), 
                                    block_end = numeric(), 
                                    ref = character(), 
                                    query = character(),
                                    chrom = character(),
                                    ref_chrom = character())

  
  for (variety in varieties_to_plot){
    print(variety)
    data_dir <- paste0(base_dir, "aln/", variety, "/", chromosomes_to_plot[i], "/")
    plot_dir <- paste0(plots_chrom_dir, "/", variety, "/")
    dir.create(plot_dir)
    
    all_files <- list.files(data_dir)
    #This gets just the filtered deltas 
    filtered_delta <- all_files[grep("_L20Kb_rq.delta", all_files)]
    
    for (comparison in filtered_delta){
      comparison_delta_path <- paste0(data_dir, "/", comparison)
      comparison_id <- str_split_fixed(comparison, "\\.", 3)[1]
      
      ref <- variety
      if (ref == "sy_mattis"){
        query <- str_split_fixed(comparison_id, "_", 4)[4]
      } else {
        query <- str_split_fixed(comparison_id, "_", 3)[3]
      }
      
      if (query == "sy"){
        query <- "sy_mattis"
      }
    
      if(file.info(comparison_delta_path)$size == 0){
        print("file empty")
      } else {
        #read and filter for the minimum size
        comparison_filt <- pre_plot_analysis(delta_path = comparison_delta_path, min_size = min_size)
        print("data read in and filtered")
        
        #plot the diagonal dot plots
        diagonal_dot_plot <- plot_by_perc_id_cap(comparison_filt, cap_lower = 97, cap_upper = 100)
        ggsave(diagonal_dot_plot, file = paste0(plot_dir, comparison_id, ".", chrom, ".diagonal_dot.min", min_size, ".png"), dpi = 300, height = 5, width = 9)
        
        #plot the percentage dotplot capped at 97%
        perc_dot_plot_97 <- plot_perc_id_v_ref(comparison_filt, ymin = 97)
        ggsave(perc_dot_plot_97, file = paste0(plot_dir, comparison_id, ".", chrom, ".percentage_dot.min", min_size, "_cap97.png"), dpi = 300, height = 5, width = 9)
       
        #now plot the median cut offs
        perc_line_5M <- plot_medians_line_colour_cut_off(comparison_filt, bin_size = 5000000, max_chrom_size = max_chrom_size)
        ggsave(perc_line_5M, file = paste0(plot_dir, comparison_id, ".", chrom, ".percentage_line.min", min_size, "_5Mb_bin_99.99median.png"), dpi = 300, height = 5, width = 9)
        
        perc_boxplot_5M <- plot_boxplots_median_colour_cut_off(comparison_filt, bin_size = 5000000, max_chrom_size = max_chrom_size)
        ggsave(perc_boxplot_5M, file = paste0(plot_dir, comparison_id, ".", chrom, ".percentage_boxplot.min", min_size, "_5Mb_bin_99.99median.png"), dpi = 300, height = 5, width = 9)
        
        #also output table with the medians of the bins
        
        comparison_filt_bin <- bin_data(comparison_filt, bin_size = 5000000, max_chrom_size = max_chrom_size)
        comparison_medians <- aggregate(perc_id ~ bin, data = comparison_filt_bin, FUN=median)
        colnames(comparison_medians)[2] <- "perc_id_median"
        comparison_medians$cut_off <- NA
        comparison_medians$cut_off <- ifelse(comparison_medians$perc_id >= cut_off, paste0(">=",cut_off), paste0("<",cut_off))
        comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)
        
        comparison_medians$ref <- ref
        comparison_medians$query <- query
        
        write.table(comparison_medians, file = paste0(blocks_chrom_dir, "/", comparison_id, ".", chrom, ".percentage_medians.min", min_size, "_5Mb_bin.txt"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
        
        #now assign the blocks
        ref_query_blocks <- assign_blocks_mummer(comparison_medians)
        ref_query_coords <- summarise_blocks(ref_query_blocks)
        if (nrow(ref_query_coords) > 0){
          ref_query_coords$chrom <- chrom
          ref_query_coords$ref_chrom <- unique(comparison_filt$rid)
        
          chrom_blocks_coords <- rbind(chrom_blocks_coords, ref_query_coords)
        }
      }
    }
  }
  #now check if blocks are in reciprocals
  chrom_blocks_coords$in_reciprocal <- apply(chrom_blocks_coords, 1, check_reciprocal, all_ref_query_coords = chrom_blocks_coords)
  
  write.table(chrom_blocks_coords, 
              file = paste0(blocks_chrom_dir, "/mummer_blocks_", chrom, ".min", min_size, ".5Mb_bins.txt"), 
              sep = "\t",
              col.names = TRUE, 
              row.names = FALSE, 
              quote = FALSE)
}



