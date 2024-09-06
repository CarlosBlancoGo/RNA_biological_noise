#!/bin/bash
#SBATCH --job-name=tama_merge_samples
#SBATCH --time 1-00:00:00
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10gb
#SBATCH --output=tama_merge_samples.out
#SBATCH --error=tama_merge_samples.err


source activate tama

filtered_gtfs=$(find "/storage/gge/Carlos/aging/biological_noise/SQANTI3_filter_results/no_fl_counts_filter" -type f -name "*filtered.gtf" | sort -V)

# First create the fofn with the filtered tama_collapse models in gtf format. 
fofn="/storage/gge/Carlos/aging/biological_noise/SQANTI3_filter_results/no_fl_counts_filter/filtered_tama_models_gtfs.fofn"
quant_fofn="/storage/gge/Carlos/aging/data/read_counts/filtered_tama_models_quant.fofn"

for sample in $filtered_gtfs;
do
	echo $sample >> "$fofn"
done


# GTF to BED12
while read -r line; do
	tama_bed=$(dirname ${line})/$(basename $line ".gtf").bed
	python /storage/gge/Carlos/tama/tama_go/format_converter/tama_format_gff_to_bed12_cupcake.py \
		$line \
		$tama_bed
	echo -e "${tama_bed}\tcapped\t1,1,1\t$(basename $line ".gtf")" >> /storage/gge/Carlos/aging/biological_noise/SQANTI3_filter_results/no_fl_counts_filter/filtered_tama_models_beds_merge.fofn

	cuant_file=$(dirname ${line})/$(basename $line ".gtf").counts.tsv
	sample_name=$(basename ${line} | sed 's/.filtered.gtf//g')
	class_file=$(dirname ${line})/*classification.txt
	echo -e "ids\t$(basename $line ".gtf")" > $cuant_file
	grep -v Artifact $class_file | tail -n +2 | cut -f1,22 >> $cuant_file
	echo $cuant_file >> ${quant_fofn}

done < $fofn


# Run tama_merge
python /storage/gge/Carlos/tama/tama_merge.py -f /storage/gge/Carlos/aging/biological_noise/SQANTI3_filter_results/no_fl_counts_filter/filtered_tama_models_beds_merge.fofn -p /storage/gge/Carlos/aging/data/transcript_models/merged_filtered_tama_isoforms -a 50 -m 0 -z 50 -d merge_dup

# Convertir el bed con las isoformas de las merged samples de formato bed12 a gtf
cd /storage/gge/Carlos/aging/data/transcript_models
merged_bed="/storage/gge/Carlos/aging/data/transcript_models/merged_filtered_tama_isoforms.bed"

# Convertir a GTF usando TAMA GO
python /storage/gge/Carlos/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py $merged_bed brain_filtered_tama_isoforms.gtf

conda deactivate
source activate my_r_env #conda environment with R and tidyverse installed

Rscript /home/apadepe/documenting_NIH/jorge/scripts/flair/merge_flair_quatification_files_after_tama_merge.R \
	$quant_fofn \
	/storage/gge/Carlos/aging/data/transcript_models/merged_filtered_tama_isoforms_merge.txt \
	/storage/gge/Carlos/aging/data/read_counts/merged_filtered_samples.counts.tsv \
	"isoseq"