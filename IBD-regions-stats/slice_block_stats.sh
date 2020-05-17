#!/bin/bash -e

cd pangenome-web

output_folder=block_stats_by_slice_500kbp/
mkdir -p $output_folder
folder=regions
regions_path=$folder/cs_regions.bed
block_sizes='Tae_5000kbp Tae_2500kbp Tae_1000kbp    '
slice_size=500000

for f in $block_sizes ;do
	echo $f
	output="${output_folder}/Block_slice500k_stats_${f}.tsv"
	echo $output
	rake "haplotypes:export_haplotype_block_stats_in_points[${output},${f},Wheat,$slice_size]" &
done

wait

for f in $block_sizes ;do
	echo $f
	output="${output_folder}/Block_slice500k_stats_${f}.tsv"
	echo $output
	gzip  -f $output &
	gzip  -f $output.missing &
done 