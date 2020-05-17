#!/bin/bash

cd pangenome-web

input_folder=blocks

for chr in `ls $input_folder | grep chr` ; do

	input="${input_folder}/${chr}/haplotype_blocks_BLAST_25_gene_window_2000bp.txt"
	output="${input_folder}/${chr}/haplotype_blocks_BLAST_25_gene_window_2000bp_converted_${chr}_10genesBlocks.csv"
	echo $chr
	rake  "haplotypes:convert_gene_coordinates[${input},${output},Wheat]"

done

