

library(stringr)
library(ggplot2)


#Compare all the snps against chinese spring on 6A
#Takes output from mummer show-snps


sequence_dir <- "X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/6A/"
output_dir <- "X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/"

ref <- "chinese"

sequence_files <- list.files(sequence_dir, pattern = paste0(ref, "_v_"))
sequence_files

snps_positions_combined <- integer()

for (file in sequence_files){
  
  ref_query_temp <- str_split_fixed(file, "\\.", 2)[1]
  
  query <- str_split_fixed(ref_query_temp, "_", 3)[3]
  

  if(query == "sy"){
    query <- "sy_mattis"
  }
  
  snp_flanks <- read.table(file = paste0(sequence_dir, file),
                          sep = "\t",
                          skip = 4,
                          header = FALSE,
                          stringsAsFactors = FALSE)

  colnames(snp_flanks) <- c("P1", "SUB1", "SUB2", "P2", "BUFF", "DIST", "REF_SEQ", "QUERY_SEQ", "FRM", "TAGS", "REF", "QUERY")
  
  #remove any snps that are Ns and also get rid of all irrelevant info - i.e. only keep chinese spring position and alleles
  snp_flanks_noN <- snp_flanks[!((snp_flanks$SUB1 == "N") | (snp_flanks$SUB2 == "N")), c("P1", "SUB1", "SUB2")]
  nrow(snp_flanks_noN)
  
  #also remove any where the alleles are "." i.e. an indel
  snps_only <- snp_flanks_noN[!((snp_flanks_noN$SUB1 == ".") | (snp_flanks_noN$SUB2 == ".")),]
  nrow(snps_only)
  
  snp_positions <- unique(snps_only$P1)
  length(snp_positions)
  
  snp_positions_to_keep <- snp_positions[!(snp_positions %in% snps_positions_combined)]
  
  snps_positions_combined <- c(snps_positions_combined, snp_positions_to_keep)
  print(file)
}

head(snps_positions_combined)
length(snps_positions_combined)

snps_positions_combined <- snps_positions_combined[order(snps_positions_combined)]
hist(snps_positions_combined)

snps_positions_combined_df <- data.frame(chrom = "chr6A", pos = snps_positions_combined)

write.table(snps_positions_combined_df, 
            file = "Y:/Publications/Haplotypes/Figures/figure_3/pangenome_snp_dist_6A_CS_ref.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
