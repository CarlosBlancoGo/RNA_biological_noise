#!/bin/bash

#SBATCH --job-name=tama_collapse_and_process_nanopore
#SBATCH --time=7-00:00:00
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=130gb
#SBATCH --output=tama_collapse_and_process_nanopore_%A_%a.out
#SBATCH --error=tama_collapse_and_process_nanopore_%A_%a.err
#SBATCH --array=0-8

# Load the needed modules and conda env
source activate tama
module load samtools

# Establish the variables and routes for the needed files
WD="/storage/gge/Carlos/aging/data_nanopore/subsampled_bams_seed123/"
if [ ! -d "$WD" ]; then
    mkdir -p "$WD"
fi
bam_files=("${WD}subsampled_bams"/B{1,3}*subsampled.bam)
reference_genome="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"

# Select the BAM file for the current array job
bam_file="${bam_files[$SLURM_ARRAY_TASK_ID]}"
name=$(basename ${bam_file} .bam)

cd $WD

# Create necessary directories
directories=("split_files" "transcript_models" "read_counts")
for dir in "${directories[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir -p "$dir"
  fi
done

# Split the bam files into different sam chunks
cd split_files
if [ ! -d "${name}_split" ]; then
    mkdir -p "${name}_split"
fi
cd ${WD}split_files/${name}_split

# Split the bam file
python /storage/gge/Carlos/tama/tama_go/split_files/tama_mapped_sam_splitter.py ${bam_file} 25 ${name}_split

# Create the directory for the collapsed files
if [ ! -d "${name}_split_collapsed" ]; then
    mkdir -p "${name}_split_collapsed"
fi
cd ./${name}_split_collapsed

# Collapse the sam_files in parallel using GNU Parallel
sam_files=$(find ../ -name "*.sam")
export reference_genome
export name
parallel -j 3 python /storage/gge/Carlos/tama/tama_collapse.py -s {} -f ${reference_genome} -p {/.} -c 75 -log log_off -icm ident_map -a 50 -m 0 -z 50 ::: $sam_files

# Concatenate the bed12 files and the .txt files with the read count of the transcript models
cd ${WD}transcript_models
mkdir -p bed_files
cd bed_files

# Get the collapsed bed files corresponding to each bam file
bed_split_files=$(find "${WD}split_files/${name}_split/${name}_split_collapsed" -type f -name "*_split_*.bed" ! -name "*_trans_read.bed" | sort -V)

# Extract the sample name
sample_name=$(basename ${bam_file} _primary_aln_sorted_subsampled.bam)

# Concatenate the bed files
cat $bed_split_files > "tmp_${sample_name}_tama_models.bed"
echo "Collapsing and concatenation for ${sample_name} completed successfully."

# Concatenate the .txt files with the counts
cd ${WD}read_counts
counts_split_files=$(find "$WD/split_files/${name}_split/${name}_split_collapsed" -type f -name "*_trans_report.txt" | sort -V)
cat $counts_split_files > "tmp_${sample_name}_counts.txt"

# Update the transcript models IDs
cd ${WD}transcript_models/bed_files

# Update bed files
INPUT_FILE="tmp_${sample_name}_tama_models.bed"
OUTPUT_FILE="${sample_name}_tama_models.bed"
echo "Updating ${sample_name} transcript models IDs..."

global_transcript_counter=0
last_transcript_number=0

while IFS= read -r line || [[ -n "$line" ]]; do
    transcript_id=$(awk '{print $4}' <<< "$line" | cut -d ';' -f1 | grep -o -E 'G[0-9]+')
    isoform_id=$(awk '{print $4}' <<< "$line" | cut -d ';' -f2)
    transcript_number=$(grep -o -E '[0-9]+' <<< "$transcript_id")

    if [[ "$transcript_number" -ne "$last_transcript_number" ]]; then
        global_transcript_counter=$((global_transcript_counter + 1))
        last_transcript_number="$transcript_number"
    fi

    new_transcript_id="G${global_transcript_counter}"
    new_id="${new_transcript_id};${new_transcript_id}.${isoform_id#*.}"
    new_line=$(sed "s/$transcript_id;$isoform_id/$new_id/" <<< "$line")
    echo "$new_line"
done < "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Updating of ${sample_name} transcript models IDs completed successfully. Results are in ${OUTPUT_FILE}."
rm tmp_${sample_name}_tama_models.bed

# Update counts files
cd $WD/read_counts
INPUT_FILE="tmp_${sample_name}_counts.txt"
OUTPUT_FILE="${sample_name}_counts.txt"
echo -e "pbid\tcount_fl\tnorm_fl" > "$OUTPUT_FILE"

global_transcript_counter=0
last_transcript_number=0

while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" =~ ^transcript_id ]]; then
        continue
    fi
    transcript_id=$(awk '{print $1}' <<< "$line" | cut -d '.' -f1 | grep -o -E 'G[0-9]+')
    isoform_id=$(awk '{print $1}' <<< "$line")
    transcript_number=$(grep -o -E '[0-9]+' <<< "$transcript_id")

    if [[ "$transcript_number" -ne "$last_transcript_number" ]]; then
        global_transcript_counter=$((global_transcript_counter + 1))
        last_transcript_number="$transcript_number"
    fi

    new_transcript_id="G${global_transcript_counter}"
    new_id="${new_transcript_id}.${isoform_id#*.}"
    new_line=$(sed "s/$isoform_id/$new_id/" <<< "$line")
    new_line=$(echo "$new_line" | awk 'BEGIN{OFS="\t"} {print $1, $2, $2}')
    echo "$new_line" >> "$OUTPUT_FILE"
done < "$INPUT_FILE"

echo "Updating of ${sample_name} transcript models counts IDs completed successfully. Results are in ${OUTPUT_FILE}."
rm tmp_${sample_name}_counts.txt

# Convert the bed12 files to gtf
cd ${WD}transcript_models
mkdir -p gtf_files
cd gtf_files

INPUT_FILE="${WD}transcript_models/bed_files/${sample_name}_tama_models.bed"
OUTPUT_FILE="${sample_name}_tama_models.gtf"

python /storage/gge/Carlos/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py $INPUT_FILE $OUTPUT_FILE

pwd; date
