## Script written by Jemima Brinton 2019
## Aim - define haplotype blocks across whole wheat genome based on pairwise BLAST alignments 


##DEFINE FUNCTIONS
library(stringr)
library(ggplot2)
library("ggdendro")
library(plyr)
library(reshape2)
library("grid")

#read in BLAST alignments and get positions 
read_pairwise_position <- function(blast_path_gz, gtf, outfile, write_table = TRUE){
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
  
  if(write_table == TRUE){
  write.table(all_comp_positions, file = outfile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
  return(all_comp_positions)
}

#cap percentage identity at a certain value
cap_data <- function(all_comp_positions, cap_value){
  all_comp_positions$capped_ident <- all_comp_positions$pident
  all_comp_positions[all_comp_positions$capped_ident < cap_value, "capped_ident"] <- 0
  return(all_comp_positions)
}

#filter for the number of missing alignments
filter_aln_no <- function(comp_data, missing_limit){
  #total_number of alignments
  total_aln <- length(unique(comp_data$aln_type))
  aln_limit <- total_aln - missing_limit
  
  transcript_counts<- data.frame(table(comp_data$transcript))
  transcripts_to_keep <- transcript_counts[transcript_counts$Freq >= aln_limit, "Var1"]
  
  comp_data_filt <- subset(comp_data, transcript %in% transcripts_to_keep)
  
  nrow(comp_data)
  nrow(comp_data_filt)
  length(unique(comp_data$transcript))
  length(unique(comp_data_filt$transcript))
  
  return(comp_data_filt)
}


#this zooms in on an region and plots the 100% identity and orders alignments via heirachical clustering
plot_zoom_100_clust <- function(data, start, end, outfile){
  data <- data[order(data$start),]
  data_zoom <- data[(data$end > start) & (data$start < end),]
  #print(table(data_zoom$aln_type))
  zoom_mat <- dcast(data_zoom[,c("aln_type", "transcript", "capped_ident")], aln_type ~ transcript)
  zoom.matrix <- as.matrix(zoom_mat[,c(2:(ncol(zoom_mat)-1))])
  dim(zoom.matrix)
  rownames(zoom.matrix) <- zoom_mat$aln_type
  zoom.dendro <- as.dendrogram(hclust(d=dist(x = zoom.matrix)))
  dendro.plot <- ggdendrogram(data = zoom.dendro, rotate = TRUE)
  zoom.order <- order.dendrogram(zoom.dendro)
  data_zoom$aln_type <- factor(x = data_zoom$aln_type, levels = zoom_mat$aln_type[zoom.order], ordered = TRUE)
  plot <- ggplot(data_zoom , aes(x=factor(transcript),y=aln_type)) +
    geom_tile(aes(fill = data_zoom$capped_ident)) +
    scale_fill_distiller(limits=c(0,100), type='div', palette="YlGnBu", na.value = "gray94", direction = 1, name = "% ID") +
    theme(axis.text.x=element_text(angle=90, hjust=1)) 
  ggsave(plot = plot, file =outfile, dpi = 600, height = 8, width = 18)
  return(plot)
}

#calculate sliding windows of percentage identity based on a certain number of consective transcripts (window_size) and dropping a specific numer with the lowest identity before calculating the mean % id
calculate_pid_windows <- function(aln_data, window_size = 20, drop_no = 2){
  aln_data <- aln_data[order(aln_data$start),]
  aln <- unique(aln_data$aln_type)
  mean_pidents <- data.frame(start = numeric(), 
                             end = numeric(), 
                             start_position = numeric(), 
                             end_position = numeric(), 
                             pident_mean = numeric(), 
                             aln_type = character(),
                             start_transcript = character(),
                             end_transcript = character()) 
  
  for(j in seq(1, nrow(aln_data)-window_size)){
    window_start = j
    start_position = aln_data[window_start, "start"]
    start_transcript = aln_data[window_start, "transcript"]
    window_end = (window_start + window_size) - 1
    end_position = aln_data[window_end, "end"]
    end_transcript = aln_data[window_end, "transcript"]
    pidents <- aln_data[c(window_start:window_end), "pident"]
    pidents <- pidents[order(pidents, decreasing = TRUE)]
    pidents_sub <- pidents[c(1:(length(pidents)-drop_no))]
    pident_mean <- mean(pidents_sub)
    aln_to_add = data.frame(start = window_start, 
                            end = window_end, 
                            start_position = start_position, 
                            end_position = end_position, 
                            pident_mean = pident_mean, 
                            aln_type = aln, 
                            start_transcript = start_transcript,
                            end_transcript = end_transcript)
    mean_pidents <- rbind(mean_pidents, aln_to_add)
  }
  return(mean_pidents)
}

#assign haplotype blocks based on the sliding window mean of percentage identity
assign_blocks <-function(mean_pidents){
  mean_pidents_copy <- mean_pidents
  mean_pidents_copy$block_no <- NA
  block_no = 1
  new_block <- "no"
  
  for (i in seq(1, nrow(mean_pidents_copy))){
    #print(i)
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

#summarise adjacted haplotype blocks calculated from sliding window mean
summarise_blocks <- function(mean_pidents_copy){
  #now get the start and end of each block
  blocks <- unique(mean_pidents_copy$block_no)
  blocks <- blocks[complete.cases(blocks)]
  
  block_positions <- data.frame(block_no = numeric(), 
                                block_start = numeric(), 
                                block_end = numeric(), 
                                aln_type = character(),
                                start_transcript = character(),
                                end_transcript = character())
  
  for(block in blocks){
    block_data <- subset(mean_pidents_copy, block_no == block)
    block_start <- min(block_data$start_position)
    start_transcript <- unique(block_data[block_data$start_position == block_start, "start_transcript"])
    block_end <- max(block_data$end_position)
    end_transcript <- unique(block_data[block_data$end_position == block_end, "end_transcript"])
    aln_type <- unique(block_data$aln_type)
    to_add <- data.frame(block_no = block, 
                         block_start = block_start, 
                         block_end = block_end, 
                         aln_type = aln_type,
                         start_transcript = start_transcript,
                         end_transcript = end_transcript)
    
    block_positions <- rbind(block_positions, to_add)
  }
  return(block_positions)
}

calculate_blocks_single_window <- function (data, gene_no_window, drop_no, window, out_dir, chrom){
  
  overall_blocks <- data.frame(block_no = numeric(), 
                               block_start = numeric(), 
                               block_end = numeric(), 
                               aln_type = character(),
                               start_transcript = character(),
                               end_transcript = character(),
                               window = character())
  
  alns <- unique(data$aln_type)
    for(aln in alns){
      print(aln)
      aln_data <- data[data$aln_type == aln,]
      #block_info <- calculate_hap_blocks(aln_data = aln_data)
      
      aln_mean_pidents <- calculate_pid_windows(aln_data = aln_data, window_size = gene_no_window, drop_no = drop_no)
      aln_mean_pidents <- aln_mean_pidents[complete.cases(aln_mean_pidents),]
      if(max(aln_mean_pidents$pident_mean) < 100){
        block_info <- data.frame(block_no = NA, 
                                 block_start = NA, 
                                 block_end = NA, 
                                 aln_type = aln, 
                                 start_transcript = NA,
                                 end_transcript = NA,
                                 window = window)
      } else {
        aln_mean_pidents_copy <- assign_blocks(mean_pidents = aln_mean_pidents)
        block_info <- summarise_blocks(mean_pidents_copy = aln_mean_pidents_copy)
        block_info$window <- window
      }
    overall_blocks <- rbind(overall_blocks, block_info)
    }
  
  #outfile <- paste0(out_dir, "haplotype_blocks_aln_", gene_no_window, "_gene_window_drop", drop_no, ".tsv")
  #write.table(overall_blocks, file = outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  return(overall_blocks)
}

## END OF FUNCTIONS

raw_data_path <- "W:/ramirezr/SM1/pairwise_blast_nov_2019/varieties_all_identities_2000bp.tab.gz"
out_dir <- "X:/brintonj/haplotype/whole_genome_blast/"
plot_dir <- paste0(out_dir, "plots/")
blocks_dir <- paste0(out_dir, "blocks/")

dir.create(plot_dir)
dir.create(blocks_dir)

#same gtf for all (using RefSeq order for now)
HC_gtf <- read.table("W://WGAv1.0//annotation//IWGSC_v1.1_HC_20170706.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
LC_gtf <- read.table("W://WGAv1.0//annotation//IWGSC_v1.1_LC_20170706.gff3", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
ALL_gtf <- rbind(HC_gtf, LC_gtf)

raw_data <- read_pairwise_position(blast_path_gz = raw_data_path, 
                                   gtf = ALL_gtf, 
                                   write_table = FALSE)

#make a capped column in the all positions data_frame
all_comp_positions <- cap_data(all_comp_positions = raw_data, cap_value = 100)

#filter for Ns
all_comp_positions_noN <- all_comp_positions[all_comp_positions$Ns_total == 0,]

#plot how many Ns we lose for each chromosome
gene_counts_all <- data.frame(table(all_comp_positions[,c("aln_type", "chr")]))
gene_counts_all$filter <- "all"
gene_counts_noN <- data.frame(table(all_comp_positions_noN[,c("aln_type", "chr")]))
gene_counts_noN$filter <- "noN"

gene_counts_filter <- rbind(gene_counts_all, gene_counts_noN)

ggplot(gene_counts_filter, aes(x = chr, y = Freq, fill = filter)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("number of genes")

ggsave(file = paste0(out_dir, "gene_counts_pairwise_blast_whole_genome_all_v_noN.png"), height = 5, width = 6, dpi = 300)

chromosomes <- unique(all_comp_positions_noN$chr)
#get rid of chrUn - we need physical order
chromosomes <- chromosomes[!(chromosomes == "chrUn")]

window <- unique(all_comp_positions_noN$Flanking)

#use a sliding window size of 25 consecutive transcripts
gene_no_window <- 25

for (chr in chromosomes){
  chr_data <- all_comp_positions_noN[all_comp_positions_noN$chr == chr,]
  chrom_plot_dir <- paste0(plot_dir, chr, "/")
  dir.create(chrom_plot_dir, recursive = TRUE)
  
  chrom_blocks_dir <- paste0(blocks_dir, chr, "/")
  dir.create(chrom_blocks_dir, recursive = TRUE)
  
  varieties <- unique(chr_data$var_query)

    for (i in seq(1, length(varieties))){
      variety <- varieties[i]
    
      variety_data <- chr_data[grep(variety, chr_data$aln_type),]
    
      #now filter for those in all comparisons
      variety_data_all <- filter_aln_no(comp_data = variety_data, missing_limit = 0)
    
      #exclude those transcripts missing more than two comparisons
      variety_data_2 <- filter_aln_no(comp_data = variety_data, missing_limit = 2)
    
      #plot heatmap with all data
      all_plot <- plot_zoom_100_clust(variety_data, 
                                    start = 0, 
                                    end = ceiling(max(variety_data$end)), 
                                    outfile = paste(chrom_plot_dir, variety, "_", window, "_noN_100cap_clustered.png", sep = ""))
    
      #plot heatmap with no missing data
      all_plot <- plot_zoom_100_clust(variety_data_all, 
                                    start = 0, 
                                    end = ceiling(max(variety_data$end)), 
                                    outfile = paste(chrom_plot_dir, variety, "_", window, "_0missing_noN_100cap_clustered.png", sep = ""))
    
    #plot heatmap with up to 2 missing data points
    all_plot <- plot_zoom_100_clust(variety_data_2, 
                                    start = 0, 
                                    end = ceiling(max(variety_data$end)), 
                                    outfile = paste(chrom_plot_dir, variety, "_", window, "_max2missing_noN_100cap_clustered.png", sep = ""))
    }

  #now do the block calculation
  chrom_blocks <- calculate_blocks_single_window(data = chr_data,
                                                 gene_no_window = gene_no_window, 
                                                 drop_no = ceiling(gene_no_window*0.1), 
                                                 window = window, 
                                                 out_dir = chrom_blocks_dir, 
                                                 chrom = chr)
  
  write.table(chrom_blocks, file = paste0(chrom_blocks_dir, "haplotype_blocks_BLAST_", gene_no_window, "_gene_window_", window, "bp.txt"),
              sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)

  
}
