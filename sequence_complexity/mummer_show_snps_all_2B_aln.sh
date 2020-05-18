#!/bin/bash -e
#SBATCH --mem=10000
#SBATCH --mail-type=END,FAIL 
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -p jic-short,RG-Cristobal-Uauy
#SBATCH --array=0-255
#SBATCH -o ./run_logs/show_snps_2B.%a.%N.%j.out
#SBATCH -e ./run_logs/show_snps_2B.%a.%N.%j.err

CHR_NAME="2B"

t_id=$SLURM_ARRAY_TASK_ID

x=$(($t_id / 15))
q_x=$(($t_id % 15))


declare -a names=("arinalrfor" "chinese" "jagger" "julius" "lancer" "landmark" "mace" "norin61" "stanley" "sy_mattis" "cadenza" "claire" "paragon" "robigus" "weebil")

source mummer-3.23


echo $t_id
echo $x
echo $q_x
echo $c_x


REF_NAME=${names[$x]}
QUERY_NAME=${names[$q_x]}


echo $REF_NAME
echo $QUERY_NAME
echo $CHR_NAME


ALN_DIR='/jic/scratch/groups/Cristobal-Uauy/brintonj/6A_region/mummer/whole_2B_aln/'$REF_NAME
ALN_ID=$REF_NAME'_v_'$QUERY_NAME'.all_'$CHR_NAME'_filtered_L20Kb_rq'

OUTPUT_DIR='/jic/scratch/groups/Cristobal-Uauy/brintonj/haplotype/whole_genome_mummer/surrounding_seq/'$CHR_NAME

mkdir -p $OUTPUT_DIR

if [ "$REF_NAME" == "$QUERY_NAME" ]
then
	echo "Reference and Query are the same"
else
	echo "Reference and Query are different, continuing with mummer alignments"
	#extract the 10bp up and downstream of each snp
	show-snps -C -T -x 10 $ALN_DIR/$ALN_ID.delta > $OUTPUT_DIR/$ALN_ID'_10bpFlank.snps'
fi

