# script for splitting a bam file per chromosome
# Usage: bash split_bam_per_chromosome.sh $sample_name $bamfile
SAMPLE_NAME=$1
BAM=$2

for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT; do

	samtools view $BAM $chr -b >  "output/split_bams/"$SAMPLE_NAME"_bqsr_chr"$chr".bam"
	samtools index "output/split_bams/"$SAMPLE_NAME"_bqsr_chr"$chr".bam"

done




