library(dplyr)
library(ggplot2)
##calculate for 6A
base_dir <- "X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/"

seq_files <- list.files(base_dir, pattern = "tsv")

length(seq_files)
sample_no <- 10000
hap_blocks <- data.frame(matrix( nrow = length(seq_files)*sample_no, ncol = 6))
colnames(hap_blocks) <- c("P1", "REF_SEQ", "REF", "QUERY", "hap_block", "REF_CPLX")
not_block <- data.frame(matrix(nrow = length(seq_files)*sample_no, ncol = 6))
colnames(not_block) <- c("P1", "REF_SEQ", "REF", "QUERY", "hap_block", "REF_CPLX")

stats <- data.frame(matrix(nrow = length(seq_files), ncol = 13))
colnames(stats) <- c("FILE", 
                     "N_HAP", 
                     "N_NOT", 
                     "MEDIAN_HAP", 
                     "MEDIAN_NOT", 
                     "MEAN_HAP", 
                     "MEAN_NOT", 
                     "P_VALUE", 
                     "MEDIAN_HAP_SUB", 
                     "MEDIAN_NOT_SUB", 
                     "MEAN_HAP_SUB", 
                     "MEAN_NOT_SUB", 
                     "P_VALUE_SUB")

row_start <- 1

for (i in seq(1, length(seq_files))){
  file <- seq_files[i]
  name <- gsub("_snps_ref_complexity.tsv", "", file)
  seqs <- read.table(file = paste0(base_dir, file),
                     sep = "\t",
                     header = TRUE, 
                     stringsAsFactors = FALSE)
  
  seqs$hap_block <- factor(seqs$hap_block, levels = c("Y", "N"))
  
  distr <- ggplot(seqs, aes(x = REF_CPLX, group = hap_block, fill = hap_block)) +
            geom_density(alpha = 0.5) +
            ggtitle(paste0(file, "\nN_hap: ", nrow(seqs[seqs$hap_block == "Y",]), "; N_not: ", nrow(seqs[seqs$hap_block == "N",]) ))
  

  pdf(paste0("X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/combined/", name, "_complexity_distribution.pdf"), height = 3, width = 3.5)
  print(distr)
  dev.off()
  
  head(seqs)
  nrow(seqs)
  counts <- data.frame(table(seqs$hap_block))
  
  medians <- data.frame(aggregate(REF_CPLX ~ hap_block, data = seqs, FUN = median))
  means <- data.frame(aggregate(REF_CPLX ~ hap_block, data = seqs, FUN = mean))
  
  stats[i,"FILE"] <- file
  stats[i,"N_HAP"] <- counts[counts$Var1 == "Y", "Freq"]
  stats[i,"N_NOT"] <- counts[counts$Var1 == "N", "Freq"]
  
  if (nrow(medians[medians$hap_block == "Y",]) > 0){
  stats[i,"MEDIAN_HAP"] <- medians[medians$hap_block == "Y","REF_CPLX"]
  stats[i,"MEAN_HAP"] <- means[means$hap_block == "Y","REF_CPLX"]
  }
  if (nrow(medians[medians$hap_block == "N",]) > 0){
  stats[i,"MEDIAN_NOT"] <- medians[medians$hap_block == "N","REF_CPLX"]
  stats[i,"MEAN_NOT"] <- means[means$hap_block == "N","REF_CPLX"]
  }
  
  if ((nrow(medians[medians$hap_block == "N",]) > 0) & (nrow(medians[medians$hap_block == "Y",]) > 0)){
  stats[i,"P_VALUE"] <- (pairwise.wilcox.test(seqs$REF_CPLX,seqs$hap_block,  p.adjust.method = "BH"))$p.value
  }
  
  row_end <- sample_no*i
  if((nrow(seqs[seqs$hap_block == "Y",]) < sample_no) | (nrow(seqs[seqs$hap_block == "N",]) < sample_no)){
    print("not enough samples")
    hap_blocks[c(row_start:row_end),"REF"] <- file
    not_block[c(row_start:row_end),"REF"] <- file
  } else {
    seqs$hap_block <- as.character(seqs$hap_block)
    hap_samples <- sample_n(seqs[seqs$hap_block == "Y",], sample_no)
    hap_blocks[c(row_start:row_end),] <- hap_samples
    not_samples <- sample_n(seqs[seqs$hap_block == "N",], sample_no)
    not_block[c(row_start:row_end),] <- not_samples
    
    #also want to calculate the subsampled stats per alignment
    subsamples <- rbind(hap_samples, not_samples)
    
    medians_sub <- data.frame(aggregate(REF_CPLX ~ hap_block, data = subsamples, FUN = median))
    means_sub <- data.frame(aggregate(REF_CPLX ~ hap_block, data = subsamples, FUN = mean))
    
    stats[i,"MEDIAN_HAP_SUB"] <- medians_sub[medians_sub$hap_block == "Y","REF_CPLX"]
    stats[i,"MEAN_HAP_SUB"] <- means_sub[means_sub$hap_block == "Y","REF_CPLX"]
    
    stats[i,"MEDIAN_NOT_SUB"] <- medians_sub[medians_sub$hap_block == "N","REF_CPLX"]
    stats[i,"MEAN_NOT_SUB"] <- means_sub[means_sub$hap_block == "N","REF_CPLX"]
    
    stats[i,"P_VALUE_SUB"] <- (pairwise.wilcox.test(subsamples$REF_CPLX, subsamples$hap_block,  p.adjust.method = "BH"))$p.value
  }
  row_start <- row_end + 1
}

write.table(stats, file = paste0("X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/combined/stats_seq_complexity_per_chrom", sample_no, "random.txt"),
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE,
            quote = FALSE)

head(hap_blocks)
head(not_block)

haps_noNA <- hap_blocks[complete.cases(hap_blocks$P1),]
haps_noNA$hap_block <- "Y"
not_noNA <- not_block[complete.cases(not_block$P1),]
not_noNA$hap_block <- "N"

combined <- rbind(haps_noNA, not_noNA)
sink(paste0("X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/combined/stats_seq_overall_", sample_no, "random.txt"))
median(haps_noNA$REF_CPLX)
median(not_noNA$REF_CPLX)

pairwise.wilcox.test(combined$REF_CPLX, combined$hap_block,  p.adjust.method = "BH")

sink()

ggplot(combined, aes(x = hap_block, y = REF_CPLX, group = hap_block)) +
  geom_boxplot()

ggsave(paste0("X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/combined/stats_seq_overall_", sample_no, "random_boxplot.png"))

#plot density distribution
combined$hap_block <- factor(combined$hap_block, levels = c("Y", "N"))

distr_dens <- ggplot(combined, aes(x = REF_CPLX, group = hap_block, fill = hap_block)) +
geom_density(alpha = 0.5) +
  ggtitle(paste0("6A random ", sample_no, "\nN_hap: ", nrow(combined[combined$hap_block == "Y",]), "; N_not: ", nrow(combined[combined$hap_block == "N",]) ))

pdf(paste0("X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/combined/stats_seq_overall_", sample_no, "random_density.pdf"), height = 3, width = 3.5)
print(distr_dens)
dev.off()

##now calculate for 2B
base_dir <- "X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/2B/complexity/"
out_dir <- paste0(base_dir, "stats/")
dir.create(out_dir)
chrom <- "2B"
seq_files <- list.files(base_dir, pattern = "tsv")

length(seq_files)
sample_no <- 10000
hap_blocks <- data.frame(matrix( nrow = length(seq_files)*sample_no, ncol = 6))
colnames(hap_blocks) <- c("P1", "REF_SEQ", "REF", "QUERY", "hap_block", "REF_CPLX")
not_block <- data.frame(matrix(nrow = length(seq_files)*sample_no, ncol = 6))
colnames(not_block) <- c("P1", "REF_SEQ", "REF", "QUERY", "hap_block", "REF_CPLX")

stats <- data.frame(matrix(nrow = length(seq_files), ncol = 13))
colnames(stats) <- c("FILE", 
                     "N_HAP", 
                     "N_NOT", 
                     "MEDIAN_HAP", 
                     "MEDIAN_NOT", 
                     "MEAN_HAP", 
                     "MEAN_NOT", 
                     "P_VALUE", 
                     "MEDIAN_HAP_SUB", 
                     "MEDIAN_NOT_SUB", 
                     "MEAN_HAP_SUB", 
                     "MEAN_NOT_SUB", 
                     "P_VALUE_SUB")

row_start <- 1

for (i in seq(1, length(seq_files))){
  file <- seq_files[i]
  name <- gsub("_snps_ref_complexity.tsv", "", file)
  seqs <- read.table(file = paste0(base_dir, file),
                     sep = "\t",
                     header = TRUE, 
                     stringsAsFactors = FALSE)
  
    seqs$hap_block <- factor(seqs$hap_block, levels = c("Y", "N"))
  
  distr <- ggplot(seqs, aes(x = REF_CPLX, group = hap_block, fill = hap_block)) +
    geom_density(alpha = 0.5) +
    ggtitle(paste0(file, "\nN_hap: ", nrow(seqs[seqs$hap_block == "Y",]), "; N_not: ", nrow(seqs[seqs$hap_block == "N",]) ))
  
  pdf(paste0(out_dir, name, "_", chrom, "_complexity_distribution.pdf"), height = 3, width = 3.5)
  print(distr)
  dev.off()
  
  #ggsave(distr, file = paste0(out_dir, name, "_", chrom, "_complexity_distribution.png"), height = 5, width = 6)
  
  head(seqs)
  nrow(seqs)
  counts <- data.frame(table(seqs$hap_block))
  
  medians <- data.frame(aggregate(REF_CPLX ~ hap_block, data = seqs, FUN = median))
  means <- data.frame(aggregate(REF_CPLX ~ hap_block, data = seqs, FUN = mean))
  
  stats[i,"FILE"] <- file
  stats[i,"N_HAP"] <- counts[counts$Var1 == "Y", "Freq"]
  stats[i,"N_NOT"] <- counts[counts$Var1 == "N", "Freq"]
  
  if (nrow(medians[medians$hap_block == "Y",]) > 0){
    stats[i,"MEDIAN_HAP"] <- medians[medians$hap_block == "Y","REF_CPLX"]
    stats[i,"MEAN_HAP"] <- means[means$hap_block == "Y","REF_CPLX"]
  }
  if (nrow(medians[medians$hap_block == "N",]) > 0){
    stats[i,"MEDIAN_NOT"] <- medians[medians$hap_block == "N","REF_CPLX"]
    stats[i,"MEAN_NOT"] <- means[means$hap_block == "N","REF_CPLX"]
  }
  
  if ((nrow(medians[medians$hap_block == "N",]) > 0) & (nrow(medians[medians$hap_block == "Y",]) > 0)){
    stats[i,"P_VALUE"] <- (pairwise.wilcox.test(seqs$REF_CPLX,seqs$hap_block,  p.adjust.method = "BH"))$p.value
  }
  
  row_end <- sample_no*i
  if((nrow(seqs[seqs$hap_block == "Y",]) < sample_no) | (nrow(seqs[seqs$hap_block == "N",]) < sample_no)){
    print("not enough samples")
    hap_blocks[c(row_start:row_end),"REF"] <- file
    not_block[c(row_start:row_end),"REF"] <- file
  } else {
    seqs$hap_block <- as.character(seqs$hap_block)
    hap_samples <- sample_n(seqs[seqs$hap_block == "Y",], sample_no)
    hap_blocks[c(row_start:row_end),] <- hap_samples
    not_samples <- sample_n(seqs[seqs$hap_block == "N",], sample_no)
    not_block[c(row_start:row_end),] <- not_samples
    
    #also want to calculate the subsampled stats per alignment
    subsamples <- rbind(hap_samples, not_samples)
    
    medians_sub <- data.frame(aggregate(REF_CPLX ~ hap_block, data = subsamples, FUN = median))
    means_sub <- data.frame(aggregate(REF_CPLX ~ hap_block, data = subsamples, FUN = mean))
    
    stats[i,"MEDIAN_HAP_SUB"] <- medians_sub[medians_sub$hap_block == "Y","REF_CPLX"]
    stats[i,"MEAN_HAP_SUB"] <- means_sub[means_sub$hap_block == "Y","REF_CPLX"]
    
    stats[i,"MEDIAN_NOT_SUB"] <- medians_sub[medians_sub$hap_block == "N","REF_CPLX"]
    stats[i,"MEAN_NOT_SUB"] <- means_sub[means_sub$hap_block == "N","REF_CPLX"]
    
    stats[i,"P_VALUE_SUB"] <- (pairwise.wilcox.test(subsamples$REF_CPLX, subsamples$hap_block,  p.adjust.method = "BH"))$p.value
  }
  row_start <- row_end + 1
}

write.table(stats, file = paste0(out_dir, "stats_seq_complexity_per_chrom_", chrom, sample_no, "random.txt"),
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE,
            quote = FALSE)

head(hap_blocks)
head(not_block)

haps_noNA <- hap_blocks[complete.cases(hap_blocks$P1),]
haps_noNA$hap_block <- "Y"
not_noNA <- not_block[complete.cases(not_block$P1),]
not_noNA$hap_block <- "N"

combined <- rbind(haps_noNA, not_noNA)
sink(paste0(out_dir, "stats_seq_overall_", chrom, sample_no, "random.txt"))
median(haps_noNA$REF_CPLX)
median(not_noNA$REF_CPLX)

pairwise.wilcox.test(combined$REF_CPLX, combined$hap_block,  p.adjust.method = "BH")

sink()

ggplot(combined, aes(x = hap_block, y = REF_CPLX, group = hap_block)) +
  geom_boxplot()

ggsave(paste0(out_dir, "stats_seq_overall_", chrom, sample_no, "random_boxplot.png"))
combined$hap_block <- factor(combined$hap_block, levels = c("Y", "N"))

distr_dens <- ggplot(combined, aes(x = REF_CPLX, group = hap_block, fill = hap_block)) +
  geom_density(alpha = 0.5) +
  ggtitle(paste0(file, "\nN_hap: ", nrow(combined[combined$hap_block == "Y",]), "; N_not: ", nrow(combined[combined$hap_block == "N",]) ))

pdf(paste0(out_dir, "stats_seq_overall_", chrom, sample_no, "random_density.pdf"), height = 3, width = 3.5)
print(distr_dens)
dev.off()



#Also want to remove any sequences that have any Ns in them from the calculations
base_dir <- "X:/brintonj/haplotype/whole_genome_mummer/surrounding_seq/2B/complexity/"
out_dir <- paste0(base_dir, "stats/noN_unique/")
dir.create(out_dir)
chrom <- "2B"
seq_files <- list.files(base_dir, pattern = "tsv")

length(seq_files)
sample_no <- 10000
hap_blocks <- data.frame(matrix( nrow = length(seq_files)*sample_no, ncol = 6))
colnames(hap_blocks) <- c("P1", "REF_SEQ", "REF", "QUERY", "hap_block", "REF_CPLX")
not_block <- data.frame(matrix(nrow = length(seq_files)*sample_no, ncol = 6))
colnames(not_block) <- c("P1", "REF_SEQ", "REF", "QUERY", "hap_block", "REF_CPLX")

stats <- data.frame(matrix(nrow = length(seq_files), ncol = 13))
colnames(stats) <- c("FILE", 
                     "N_HAP", 
                     "N_NOT", 
                     "MEDIAN_HAP", 
                     "MEDIAN_NOT", 
                     "MEAN_HAP", 
                     "MEAN_NOT", 
                     "P_VALUE", 
                     "MEDIAN_HAP_SUB", 
                     "MEDIAN_NOT_SUB", 
                     "MEAN_HAP_SUB", 
                     "MEAN_NOT_SUB", 
                     "P_VALUE_SUB")

row_start <- 1

for (i in seq(1, length(seq_files))){
  file <- seq_files[i]
  name <- gsub("_snps_ref_complexity.tsv", "", file)
  seqs_raw <- read.table(file = paste0(base_dir, file),
                     sep = "\t",
                     header = TRUE, 
                     stringsAsFactors = FALSE)
  
  seqs_noN <- seqs_raw[grep("N", seqs_raw$REF_SEQ),]

  seqs <- unique(seqs_noN)
  
  seqs$hap_block <- factor(seqs$hap_block, levels = c("Y", "N"))
  
  distr <- ggplot(seqs, aes(x = REF_CPLX, group = hap_block, fill = hap_block)) +
    geom_density(alpha = 0.5) +
    ggtitle(paste0(file, "\nN_hap: ", nrow(seqs[seqs$hap_block == "Y",]), "; N_not: ", nrow(seqs[seqs$hap_block == "N",]) ))
  
  
  ggsave(distr, file = paste0(out_dir, name, "_", chrom, "_complexity_distribution_noN_unique.png"), height = 5, width = 6)
  
  head(seqs)
  nrow(seqs)
  counts <- data.frame(table(seqs$hap_block))
  
  medians <- data.frame(aggregate(REF_CPLX ~ hap_block, data = seqs, FUN = median))
  means <- data.frame(aggregate(REF_CPLX ~ hap_block, data = seqs, FUN = mean))
  
  stats[i,"FILE"] <- file
  stats[i,"N_HAP"] <- counts[counts$Var1 == "Y", "Freq"]
  stats[i,"N_NOT"] <- counts[counts$Var1 == "N", "Freq"]
  
  if (nrow(medians[medians$hap_block == "Y",]) > 0){
    stats[i,"MEDIAN_HAP"] <- medians[medians$hap_block == "Y","REF_CPLX"]
    stats[i,"MEAN_HAP"] <- means[means$hap_block == "Y","REF_CPLX"]
  }
  if (nrow(medians[medians$hap_block == "N",]) > 0){
    stats[i,"MEDIAN_NOT"] <- medians[medians$hap_block == "N","REF_CPLX"]
    stats[i,"MEAN_NOT"] <- means[means$hap_block == "N","REF_CPLX"]
  }
  
  if ((nrow(medians[medians$hap_block == "N",]) > 0) & (nrow(medians[medians$hap_block == "Y",]) > 0)){
    stats[i,"P_VALUE"] <- (pairwise.wilcox.test(seqs$REF_CPLX,seqs$hap_block,  p.adjust.method = "BH"))$p.value
  }
  
  row_end <- sample_no*i
  if((nrow(seqs[seqs$hap_block == "Y",]) < sample_no) | (nrow(seqs[seqs$hap_block == "N",]) < sample_no)){
    print("not enough samples")
    hap_blocks[c(row_start:row_end),"REF"] <- file
    not_block[c(row_start:row_end),"REF"] <- file
  } else {
    seqs$hap_block <- as.character(seqs$hap_block)
    hap_samples <- sample_n(seqs[seqs$hap_block == "Y",], sample_no)
    hap_blocks[c(row_start:row_end),] <- hap_samples
    not_samples <- sample_n(seqs[seqs$hap_block == "N",], sample_no)
    not_block[c(row_start:row_end),] <- not_samples
    
    #also want to calculate the subsampled stats per alignment
    subsamples <- rbind(hap_samples, not_samples)
    
    medians_sub <- data.frame(aggregate(REF_CPLX ~ hap_block, data = subsamples, FUN = median))
    means_sub <- data.frame(aggregate(REF_CPLX ~ hap_block, data = subsamples, FUN = mean))
    
    stats[i,"MEDIAN_HAP_SUB"] <- medians_sub[medians_sub$hap_block == "Y","REF_CPLX"]
    stats[i,"MEAN_HAP_SUB"] <- means_sub[means_sub$hap_block == "Y","REF_CPLX"]
    
    stats[i,"MEDIAN_NOT_SUB"] <- medians_sub[medians_sub$hap_block == "N","REF_CPLX"]
    stats[i,"MEAN_NOT_SUB"] <- means_sub[means_sub$hap_block == "N","REF_CPLX"]
    
    stats[i,"P_VALUE_SUB"] <- (pairwise.wilcox.test(subsamples$REF_CPLX, subsamples$hap_block,  p.adjust.method = "BH"))$p.value
  }
  row_start <- row_end + 1
}

write.table(stats, file = paste0(out_dir, "stats_seq_complexity_per_chrom_", chrom, sample_no, "random_noN_unique.txt"),
            sep = "\t", 
            col.names = TRUE, 
            row.names = FALSE,
            quote = FALSE)

head(hap_blocks)
head(not_block)

haps_noNA <- hap_blocks[complete.cases(hap_blocks$P1),]
haps_noNA$hap_block <- "Y"
not_noNA <- not_block[complete.cases(not_block$P1),]
not_noNA$hap_block <- "N"

combined <- rbind(haps_noNA, not_noNA)
sink(paste0(out_dir, "stats_seq_overall_", chrom, sample_no, "random_noN_unique.txt"))
median(haps_noNA$REF_CPLX)
median(not_noNA$REF_CPLX)

pairwise.wilcox.test(combined$REF_CPLX, combined$hap_block,  p.adjust.method = "BH")

sink()


ggplot(combined, aes(x = hap_block, y = REF_CPLX, group = hap_block)) +
  geom_boxplot()

ggsave(paste0(out_dir, "stats_seq_overall_", chrom, sample_no, "random_boxplot_noN_unique.png"))
