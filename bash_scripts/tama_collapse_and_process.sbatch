#!/bin/bash

#SBATCH --job-name=tama_collapse_and_process
#SBATCH --time 3-00:00:00
#SBATCH --qos=medium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=77gb
#SBATCH --output=tama_collapse_and_process_%A_%a.out
#SBATCH --error=tama_collapse_and_process_%A_%a.err
#SBATCH --array=0-8

# This script generates collapsed transcripts in gtf format using tama_collapse.py from the aligned and sorted bam files from young and old mice Iso-seq data samples (B31-B35 and B151-B154, respectively). 
# Because tama_collapse is very memory greedy, the bam files are first splitted into a minimum of 25 sam files with tama_mapped_sam_splitter.py. Then, the sam chunks are collapsed separately and the collapsed transcript models bed12 output files are concatenated to rebuild the whole transcript models files. The flnc reads are collapsed by their TSS and TTS with a threshold of 50 bps and by exon-exon junctions with a threshold of 0 bps. Because the IDs for the transcript models are reset with each chunk, the IDs must be updated in the final concatenated bed12 files. Finally, the bed12 output files for each sample are converted to Ensembl GTF without CDS represented using tama_convert_bed_gtf_ensembl_no_cds.py

# The same concatenation, processing and update of the transcript models IDs are performed for the *_trans_report.txt files of each sample, obtained with tama_collapse.py. These contain the number of reads collapsed in each transcript model, and are recommended when running SQANTI3 later as a mean to add quantification information.

pwd;date

# Load the needed modules and conda env
source activate tama
module load samtools

# And establish the variables and routes for the needed files
WD="/storage/gge/Carlos/aging"
bam_files=($(find "/storage/gge/nih/PacBio_IsoSeq/merged_reads/mouse_B100K0" "/storage/gge/nih/PacBio_IsoSeq/merged_reads/old_mouse" -type f -name "*.aln.sorted.bam"))
reference_genome="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"


# Select the BAM file for the current array job
bam_file="${bam_files[$SLURM_ARRAY_TASK_ID]}"
name=`basename ${bam_file} | sed 's/.bam//g'`

cd $WD

# Create the main directory where everything will be stored. 
main_directory="data"
if [ ! -d "$main_directory" ]; then
    mkdir "$main_directory"
fi
cd $main_directory

# Then create inside the different directories that will be needed to store all the files.
directories=("split_files" "transcript_models" "read_counts")

for dir in "${directories[@]}"; do
  if [ ! -d "$dir" ]; then
    mkdir "$dir"
  fi
done

# First split the bam files with the reads into different sam chunks. 
cd split_files

# Create the directories
if [ ! -d "${name}_split" ]; then
    mkdir "${name}_split"
fi

# Get inside its respective directory
cd /storage/gge/Carlos/aging/data/split_files/${name}_split

# Split the bam file
python /storage/gge/Carlos/tama/tama_go/split_files/tama_mapped_sam_splitter.py ${bam_file} 25 ${name}_split

# Create the directory for the collapsed files
if [ ! -d "${name}_split_collapsed" ]; then
    mkdir "${name}_split_collapsed"
fi
cd ./${name}_split_collapsed

# And collapse the sam_files
sam_files=$(find ../ -name "*.sam")
for sam_file in $sam_files;
do
    sam_name=`basename ${sam_file} | sed 's/.sam//g'`
	python /storage/gge/Carlos/tama/tama_collapse.py -s ${sam_file} -f ${reference_genome} -p ${sam_name} -log log_off -a 50 -m 0 -z 50 
done

# Once each chunk of each sample is collapsed, concatenate the bed12 files with the transcript models for each sample and the .txt files with the read count of the transcript models
cd $WD/$main_directory/transcript_models
if [ ! -d "bed_files" ]; then
    mkdir "bed_files"
fi
cd bed_files

# Get the collapsed bed files corresponding to each bam file
bed_split_files=$(find "$WD/${main_directory}/split_files/${name}_split/${name}_split_collapsed" -type f -name "*_split_*.bed" ! -name "*_trans_read.bed" | sort -V)

# Extract the sample name
sample_name=$(basename ${bam_file} | sed 's/.aln.sorted.bam//g')

# Concatenate the bed files
cat $bed_split_files > "tmp_${sample_name}_tama_models.bed"

echo "Collapsing and concatenation for ${sample_name} completed succesfully."

# Do the same for the .txt files with the counts
cd $WD/$main_directory/read_counts
counts_split_files=$(find "$WD/${main_directory}/split_files/${name}_split/${name}_split_collapsed" -type f -name "*_trans_report.txt" | sort -V)
cat $counts_split_files > "tmp_${sample_name}_counts.txt"

# Once we have the concatenated files for the transcript models and counts for each sample, we need to update the transcript models IDs
cd $WD/$main_directory/transcript_models/bed_files

# First for the bed files
# Initiate the processing of the file
INPUT_FILE="tmp_${sample_name}_tama_models.bed"
OUTPUT_FILE="${sample_name}_tama_models.bed"
echo "Updating ${sample_name} transcript models IDs..."

# Initialize ID number variables
global_transcript_counter=0
last_transcript_number=0

# Read each bed file
while IFS= read -r line || [[ -n "$line" ]]; do
	# Get the current ID of the transcript model
  	transcript_id=$(awk '{print $4}' <<< "$line" | cut -d ';' -f1 | grep -o -E 'G[0-9]+')
  	isoform_id=$(awk '{print $4}' <<< "$line" | cut -d ';' -f2)

	# Get the number of the transcript
	transcript_number=$(grep -o -E '[0-9]+' <<< "$transcript_id")
	
	# Verify if the transcript number has changed since the last ID
	if [[ "$transcript_number" -ne "$last_transcript_number" ]]; then
		# Increase the global gene counter and update the las number
		global_transcript_counter=$((global_transcript_counter + 1))
		last_transcript_number="$transcript_number" 
	fi
	
	# Build the new transcript and isoform ID
	new_transcript_id="G${global_transcript_counter}"
	new_id="${new_transcript_id};${new_transcript_id}.${isoform_id#*.}"
	
	# Change the transcript and isoform ID in the line and add it to the output
	new_line=$(sed "s/$transcript_id;$isoform_id/$new_id/" <<< "$line")
	echo "$new_line"
done < "$INPUT_FILE" > "$OUTPUT_FILE"

echo "Updating of ${sample_name} transcript models IDs completed succesfully. Results are in ${OUTPUT_FILE}."

# Now we can remove the temporal non-updated transcript models files
rm tmp_${sample_name}_tama_models.bed


# Then for the counts files
cd $WD/$main_directory/read_counts
# Initiate the processing of the file
INPUT_FILE="tmp_${sample_name}_counts.txt"
OUTPUT_FILE="${sample_name}_counts.txt"

# Add the needed header to the output fl abundance file
echo -e "pbid\tcount_fl\tnorm_fl" > "$OUTPUT_FILE"

echo "Updating ${sample_name} transcript models counts IDs..."

# Initialize ID number variables
global_transcript_counter=0
last_transcript_number=0

# Read each counts txt file
while IFS= read -r line || [[ -n "$line" ]]; do
       	if [[ "$line" =~ ^transcript_id ]]; then
		continue  # This ommits the header
	fi	    
    # Get the current ID of the transcript model
	transcript_id=$(awk '{print $1}' <<< "$line" | cut -d '.' -f1 | grep -o -E 'G[0-9]+')
	isoform_id=$(awk '{print $1}' <<< "$line")

    # Get the number of the transcript
	transcript_number=$(grep -o -E '[0-9]+' <<< "$transcript_id")

	# Verify if the transcript number has changed since the last ID
	if [[ "$transcript_number" -ne "$last_transcript_number" ]]; then
		# Increase the global gene counter and update the las number
		global_transcript_counter=$((global_transcript_counter + 1))
		last_transcript_number="$transcript_number"
	fi

	# Build the new transcript and isoform ID
	new_transcript_id="G${global_transcript_counter}"
	new_id="${new_transcript_id}.${isoform_id#*.}"

	# Change the transcript and isoform ID in the line and add it to the output
	new_line=$(sed "s/$isoform_id/$new_id/" <<< "$line")
	new_line=$(echo "$new_line" | awk 'BEGIN{OFS="\t"} {print $1, $2, $2}')
	echo "$new_line" >> "$OUTPUT_FILE"
done < "$INPUT_FILE"

echo "Updating of ${sample_name} transcript models counts IDs completed succesfully. Results are in ${OUTPUT_FILE}."

# And remove the temporal counts files
rm tmp_${sample_name}_counts.txt

# Lastly, convert the bed12 files with the transcript models to gtf and store them inside a different directory
cd $WD/$main_directory/transcript_models
if [ ! -d "gtf_files" ]; then
    mkdir "gtf_files"
fi
cd gtf_files

INPUT_FILE="$WD/$main_directory/transcript_models/bed_files/${sample_name}_tama_models.bed"
OUTPUT_FILE="${sample_name}_tama_models.gtf"

python /storage/gge/Carlos/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py $INPUT_FILE $OUTPUT_FILE


pwd;date
