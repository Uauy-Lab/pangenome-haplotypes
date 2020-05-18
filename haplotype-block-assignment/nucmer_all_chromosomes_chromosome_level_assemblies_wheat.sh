## this will be just to write and submit the scripts for each chromosome

cd /nbi/group-data/ifs/NBI/Cristobal-Uauy/Jemima/haplotype/mummer_whole_genome

chr_names=("1A" "1B" "1D" "2A" "2B" "2D" "3A" "3B" "3D" "4A" "4B" "4D" "5A" "5B" "5D" "6A" "6B" "6D" "7A" "7B" "7D")


for CHR_NAME in ${chr_names[*]}
	do
	echo $CHR_NAME
	echo \
'#!/bin/bash -e
#SBATCH --mem=30G
#SBATCH --mail-type=END,FAIL 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p RG-Cristobal-Uauy,jic-long
#SBATCH --array=0-99
#SBATCH -o ./run_logs/run_nucmer.%a.%N.%j.out
#SBATCH -e ./run_logs/run_nucmer.%a.%N.%j.err



t_id=$SLURM_ARRAY_TASK_ID

x=$(($t_id / 10))
q_x=$(($t_id % 10))


declare -a names=("arinalrfor" "chinese" "jagger" "julius" "lancer" "landmark" "mace" "norin61" "stanley" "sy_mattis")

source mummer-3.23


echo $t_id
echo $x
echo $q_x
echo $c_x


REF_NAME=${names[$x]}
QUERY_NAME=${names[$q_x]}


echo $REF_NAME
echo $QUERY_NAME
echo '$CHR_NAME'


ASSEMBLIES_PATH="/nbi/group-data/ifs/NBI/Cristobal-Uauy/ramirezr/Pangenomes/Bn"
OUTPUT_DIR="/jic/scratch/groups/Cristobal-Uauy/brintonj/haplotype/whole_genome_mummer/aln/$REF_NAME/'$CHR_NAME'"
CHROM_DIR="/jic/scratch/groups/Cristobal-Uauy/brintonj/haplotype/whole_genome_mummer/fasta/'$CHR_NAME'"

mkdir -p $OUTPUT_DIR

if [ "$REF_NAME" == "$QUERY_NAME" ]
then
	echo "Reference and Query are the same"
else
	echo "Reference and Query are different, continuing with mummer alignments"
	#run nucmer
	nucmer --mum --delta $CHROM_DIR/$REF_NAME.'$CHR_NAME'.fa $CHROM_DIR/$QUERY_NAME.'$CHR_NAME'.fa --prefix $OUTPUT_DIR/$REF_NAME.v.$QUERY_NAME.'$CHR_NAME'
	#filter for alignments of at least 20kb and rq
	delta-filter -l 20000 -r -q $OUTPUT_DIR/$REF_NAME.v.$QUERY_NAME.'$CHR_NAME'.delta > $OUTPUT_DIR/$REF_NAME.v.$QUERY_NAME.'$CHR_NAME'.filtered_L20Kb_rq.delta
fi
' > /nbi/group-data/ifs/NBI/Cristobal-Uauy/Jemima/haplotype/mummer_whole_genome/scripts/chrom_aln_scripts/nucmer_aln_$CHR_NAME.sh

sbatch /nbi/group-data/ifs/NBI/Cristobal-Uauy/Jemima/haplotype/mummer_whole_genome/scripts/chrom_aln_scripts/nucmer_aln_$CHR_NAME.sh
done