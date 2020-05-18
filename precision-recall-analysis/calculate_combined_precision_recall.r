#Jemima Brinton 2020
#Script to carry out the precision recall calculations for the mummer and blast blocks

library(GenomicRanges)
library(ggplot2)
## Functions

read_BLAST_remove_EI <- function(blast_blocks_path, chrom, EI_vars, gene_window_sizes){
  blast_blocks <- data.frame(assembly = character(),
                             reference = character(),
                             chromosome	= character(),
                             start = numeric(),
                             end = numeric(),
                             block_no = numeric(),
                             chr_length = numeric(),
                             window_size = numeric(),
                             gene_window = numeric())
  
  for (gene_window in gene_window_sizes){
    temp <- read.table(paste0(blast_blocks_path, "coord_converted.chr", chrom, "_haplotype_blocks_aln_", gene_window, "_gene_window_drop", ceiling(gene_window*0.1), ".tsv"),
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE)
    
    temp$gene_window <- gene_window
    temp$window_block <- paste0(temp$block_no, "_", gene_window, "gene")
    blast_blocks <- rbind(blast_blocks, temp)
  }
  
  #remove EI v EI comparisons as we can't compare with mummer
  
  block_numbers <- unique(blast_blocks$window_block)
  
  blocks_to_remove <- character()
  
  for (block_number in block_numbers){
    assemblies <- unique(blast_blocks[blast_blocks$window_block == block_number, "assembly"])
    if (length(assemblies) > 2){
      print("ERROR")
      print(block_number)
      print(assemblies)
    } else if((assemblies[1] %in% EI_vars) & (assemblies[2] %in% EI_vars)){
      blocks_to_remove <- c(blocks_to_remove, block_number)
      #print(block_number)
      #print(assemblies)
      
    }
  }
  
  blast_blocks_to_use <- blast_blocks[!(blast_blocks$window_block %in% blocks_to_remove),]
  blast_blocks_to_use[blast_blocks_to_use$window_size == "cds", "window_size"] <- "cdsbp"
  return(blast_blocks_to_use)
}

precision_recall_flanking_window_sizes <- function(overall_blocks_to_use, mummer){
  precision_recall <- data.frame(aln_type = character(), window_size = character(), BLAST_blocks = numeric(), BLAST_found = numeric(), mummer_blocks = numeric(), mummer_found = numeric(), precision = numeric(), recall = numeric())
  
  window_sizes <- unique(overall_blocks_to_use$window_size)
  alns <- unique(mummer$aln_type)
  
  for (window_size in window_sizes){
    window_overall_blocks <- overall_blocks_to_use[overall_blocks_to_use$window_size == window_size,]
    block_numbers_window <- unique(window_overall_blocks$block_no)
    for (aln in alns){
      mummer_aln <- mummer[mummer$aln_type == aln,]
      
      ref <- unique(mummer_aln$ref)
      query <- unique(mummer_aln$query)
      
      aln_blocks <- numeric()
      #get the corresponding blast block IDs
      for (block_number in block_numbers_window){
        assemblies <- unique(window_overall_blocks[window_overall_blocks$block_no == block_number, "assembly"])
        if((ref %in% assemblies) & (query %in% assemblies)){
          aln_blocks <- c(aln_blocks, block_number)
        }
      }
      
      window_aln <- window_overall_blocks[(window_overall_blocks$block_no %in% aln_blocks) & (window_overall_blocks$assembly == ref),]
      
      
      mummer_blocks <- nrow(mummer_aln)
      BLAST_blocks <- nrow(window_aln)
      
      if(mummer_blocks == 0 | BLAST_blocks == 0){
        BLAST_found <- 0
        mummer_found <- 0
      } else {
        mummer_gdf <- mummer_aln[,c("ref_chrom", "block_start", "block_end")]
        colnames(mummer_gdf) <- c("chr", "start", "end")
        
        window_gdf <- window_aln[,c("chromosome", "start", "end")]
        colnames(window_gdf) <- c("chr", "start", "end")
        
        # make genomic ranges objects for each
        mummer_gr <- makeGRangesFromDataFrame(mummer_gdf)
        window_gr <- makeGRangesFromDataFrame(window_gdf)
        
        overlaps <- findOverlaps(mummer_gr, window_gr)
        #precision first - how many BLAST blocks found in the  mummer haplotype?
        BLAST_found <- length(unique(subjectHits(overlaps)))
        #recall - how many mummer blocks found in the  BLAST haplotype?
        mummer_found <- length(unique(queryHits(overlaps)))
      }
      
      precision <- BLAST_found/BLAST_blocks
      recall <- mummer_found/mummer_blocks
      precision_recall_aln <- data.frame(aln_type = aln, window_size, BLAST_blocks, BLAST_found, mummer_blocks, mummer_found, precision, recall)
      precision_recall <- rbind(precision_recall, precision_recall_aln)
    }
  }
  return(precision_recall)
}

f1 <- function(precision, recall){
  f1_score <- 2*((precision*recall)/(precision+recall))
  return(f1_score)
}

precision_recall_all_windows <- function(blast_blocks_to_use, gene_window_sizes){
  precision_recall_all_blocks <- data.frame(aln_type = character(), 
                                            window_size = character(), 
                                            BLAST_blocks = numeric(), 
                                            BLAST_found = numeric(), 
                                            mummer_blocks = numeric(), 
                                            mummer_found = numeric(), 
                                            precision = numeric(), 
                                            recall = numeric(),
                                            gene_block = character())
  
  
  
  for (gene_window in gene_window_sizes){
    gene_window_data <- blast_blocks_to_use[blast_blocks_to_use$gene_window == gene_window,]
    
    
    precision_recall_block <- precision_recall_flanking_window_sizes(overall_blocks_to_use = gene_window_data, mummer = mummer_blocks)
    
    precision_recall_block$gene_block <- gene_window
    precision_recall_all_blocks <- rbind(precision_recall_all_blocks, precision_recall_block)
  }
  
  
  
  precision_recall_all_blocks$window_size <- factor(precision_recall_all_blocks$window_size, levels = c("cdsbp", "0bp", "1000bp", "2000bp", "5000bp"))
  precision_recall_all_blocks$gene_block <- as.character(precision_recall_all_blocks$gene_block)
  precision_recall_all_blocks$gene_block <- factor(precision_recall_all_blocks$gene_block, levels = c("10", "15", "20", "25", "30"))
  
  precision_recall_all_blocks$f1_score <- f1(precision = precision_recall_all_blocks$precision, recall = precision_recall_all_blocks$recall)
  return(precision_recall_all_blocks)
}

plots_precision_recall <- function(precision_recall_all_blocks, chrom, plot_dir){
  precision <- ggplot(precision_recall_all_blocks, aes(x = window_size, y = precision, fill = gene_block)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1")
  
  ggsave(precision, file = paste0(plot_dir, "chr", chrom, "precision_window_sizes_varying_gene_blocks.png"), 
         height = 3, width = 6)
  
  recall <- ggplot(precision_recall_all_blocks, aes(x = window_size, y = recall, fill = gene_block)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1")
  
  ggsave(recall, file = paste0(plot_dir, "chr", chrom, "recall_window_sizes_varying_gene_blocks.png"), 
         height = 3, width = 6)
  
  f1 <- ggplot(precision_recall_all_blocks, aes(x = window_size, y = f1_score, fill = gene_block)) +
    geom_boxplot() +
    scale_fill_brewer(palette = "Set1")
  
  ggsave(f1, file = paste0(plot_dir, "chr", chrom, "f1_score_window_sizes_varying_gene_blocks.png"),
         height = 3, width = 6)
  
}
## end of functions

base_dir <- "C:/Users/brintonj/Documents/2020_03_WFH/haplotype_manuscript/precision_recall/"
data_dir <- paste0(base_dir, "all_chrom_data_to_use/")
plot_dir <- paste0(base_dir, "plots/")

dir.create(plot_dir)

mummer_blocks_path <- "X:/brintonj/haplotype/whole_genome_mummer/blocks/"

EI_vars <- c("cadenza", "claire", "paragon", "robigus", "weebil")

mummer_bin_size <- 5000000

gene_window_sizes <- c(10, 15, 20, 25,30)

chromosomes <- c("1D", "2B", "3B", "4D", "6A", "7A")

precision_recall_all_blocks_all_chroms <- data.frame(aln_type = character(), 
                                          window_size = character(), 
                                          BLAST_blocks = numeric(), 
                                          BLAST_found = numeric(), 
                                          mummer_blocks = numeric(), 
                                          mummer_found = numeric(), 
                                          precision = numeric(), 
                                          recall = numeric(),
                                          gene_block = character(),
                                          f1 = numeric(),
                                          chrom = character())

for (chrom in chromosomes){
  print(chrom)
  #read in mummer blocks
  mummer_blocks <- read.table(paste0(mummer_blocks_path, chrom, "/mummer_blocks_chr", chrom, ".min20000.5Mb_bins.txt"),
                              header = TRUE, 
                              stringsAsFactors = FALSE)
  
  mummer_blocks$aln_type <- paste0(mummer_blocks$ref, "->", mummer_blocks$query)
  
  #read in blast blocks
  blast_blocks_to_use <- read_BLAST_remove_EI(blast_blocks_path = data_dir,
                                              chrom = chrom,
                                              EI_vars = EI_vars,
                                              gene_window_sizes = gene_window_sizes)
  
  precision_recall_all_blocks <- precision_recall_all_windows(blast_blocks_to_use = blast_blocks_to_use,
                                                              gene_window_sizes = gene_window_sizes)
  
  write.table(precision_recall_all_blocks, file = paste0(plot_dir, "chr", chrom, "precision_recall_f1_table.tsv"),
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE)
  
  plots_precision_recall(precision_recall_all_blocks = precision_recall_all_blocks,
                         chrom = chrom,
                         plot_dir = plot_dir)
  
  precision_recall_all_blocks$chrom <- chrom
  
  precision_recall_all_blocks_all_chroms <- rbind(precision_recall_all_blocks_all_chroms, precision_recall_all_blocks)
  
}

write.table(precision_recall_all_blocks_all_chroms, file = paste0(plot_dir, "combined_precision_recall_f1_table.tsv"),
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE)

plots_precision_recall(precision_recall_all_blocks = precision_recall_all_blocks_all_chroms,
                       chrom = "combined",
                       plot_dir = plot_dir)

