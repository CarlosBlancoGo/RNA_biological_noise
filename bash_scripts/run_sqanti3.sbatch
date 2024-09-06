#!/bin/bash
#SBATCH --job-name=run_sqanti3
#SBATCH --time 1-00:00:00
#SBATCH --qos=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mail-type=ALL # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mem=30gb
#SBATCH --output=run_sqanti3_%A_%a.out
#SBATCH --error=run_sqanti3_%A_%a.err
#SBATCH --array=0-8

pwd;date

# This script fully runs SQANTI3 modules for the characterization of the Iso-seq mice brain data.  

source activate SQANTI3.env
export PYTHONPATH=$PYTHONPATH:/home/carblan/software/SQANTI3/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/home/carblan/software/SQANTI3/cDNA_Cupcake


# First we stablish the directory with the gtf files that will be processed, reference genome, reference annotation, polyA motifs, short reads and full length abundancy files.
WD="/storage/gge/Carlos/aging/"
cd $WD
if [ ! -d "sqanti3_results" ]; then
    mkdir "sqanti3_results"
fi
data="/storage/gge/Carlos/aging/data/"
reference_genome="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39_SIRV.fa"
reference_annotation="/storage/gge/genomes/mouse_ref_NIH/reference_genome/mm39.ncbiRefSeq_SIRV.gtf"
polyA_motifs="/home/carblan/software/SQANTI3/data/polyA_motifs/mouse_and_human.polyA_motif.txt"
cage_data="/storage/gge/genomes/mouse_ref_NIH/reference_genome/lft_mm39_CAGE.bed"
sr_dir="/storage/gge/nih/Illumina_short_reads/short_reads/NOVOGENE_stranded/mapped/"
fl_counts_dir="/storage/gge/Carlos/aging/data/read_counts/"


#########################################
# 1. Run SQANTI3 Quality Control module #
#########################################

# Create the array of files
mapfile -t gtf_files < <(find "$data/transcript_models/gtf_files/" -type f -name "*.gtf" | sort -V)

# Read the file corresponding to the task
gtf=${gtf_files[$SLURM_ARRAY_TASK_ID]}

# Get the transcript models file name without extension
gtf_file_name=$(basename $gtf .gtf)

# And the sample name
sample_name=$(echo $gtf_file_name | cut -d'_' -f1)

# Now find the sr BAM file corresponding to the sample name
sj=$(echo ${sr_dir}${sample_name}/*SJ.out.tab)
bam=$(echo ${sr_dir}${sample_name}/*.sortedByCoord.out.bam)

# Also find the fl abundance file contaning the number of reads collapsed in each transcript model
counts_file=$(echo ${fl_counts_dir}${sample_name}_counts.txt)

# Write sr BAM dir to file
echo $bam > tmp_$SLURM_ARRAY_TASK_ID.fofn

# Run SQANTI3 QC for each transcript models file
python /home/carblan/software/SQANTI3/sqanti3_qc.py \
	--dir "${WD}sqanti3_results/qc_results/${sample_name}" \
	--force_id_ignore \
	--output "${gtf_file_name}" -c $sj --min_ref_len "0" \
	--SR_bam tmp_$SLURM_ARRAY_TASK_ID.fofn \
	--CAGE_peak "$cage_data" --polyA_motif_list "$polyA_motifs" \
	-fl $counts_file \
	--report "both" --saturation \
	$gtf $reference_annotation $reference_genome


echo "Classification for ${gtf_file_name} finished."


######################################################################################################################
# 2. Run SQANTI3 Quality Control module on the reference transcriptome using the short-reads for each of the samples #
######################################################################################################################

# Run SQANTI3 QC on the reference transcriptome
python /home/carblan/software/SQANTI3/sqanti3_qc.py \
        --dir "${WD}sqanti3_results/qc_results/reference_transcriptome/${sample_name}" \
        --force_id_ignore \
        --output "${sample_name}_reference_transcriptome" -c $sj --min_ref_len "0" \
        --SR_bam tmp_$SLURM_ARRAY_TASK_ID.fofn \
        --CAGE_peak "$cage_data" --polyA_motif_list "$polyA_motifs" \
        -fl $counts_file \
        --report "both" --saturation \
        $reference_annotation $reference_annotation $reference_genome

rm tmp_$SLURM_ARRAY_TASK_ID.fofn

echo "Classification of the reference transcriptome for ${sample_name} short reads finished."


################################
# 3. Run SQANTI3 Filter module #
################################

# Get the classification file of the sample
classification_file=$(echo $WD/sqanti3_results/qc_results/${sample_name}/*classification.txt)

# Get the classification file name without extension
classification_file_name=$(basename "$classification_file" | sed 's/_classification\.txt$//')

# Run SQANTI3 filter module for each classification file
python /home/carblan/software/SQANTI3/sqanti3_filter.py rules \
       --output "${gtf_file_name}" \
       --dir "${WD}sqanti3_results/filter_results/${sample_name}" \
       --gtf $corrected_gtf --faa $corrected_faa \
       -j "${WD}scripts/technical_noise_filter.json" \
       $classification_file 
echo "Isoform filtering for ${gtf_file_name} finished."


################################
# 4. Run SQANTI3 Rescue module #
################################

# Get the filtered classification file of the sample
filter_classification_file=$(echo $WD/sqanti3_results/filter_results/${sample_name}/*classification.txt)

# Get the filtered classification file name without extension
filter_classification_file_name=$(basename "$filter_classification_file" | sed 's/_classification\.txt$//')

# Get the FASTA file output by SQANTI QC, which contains the full long read transcriptome
isoforms_fasta=$(echo ${WD}sqanti3_results/qc_results/${sample_name}/*_corrected.fasta)

# Get the long read-defined transcriptome filtered annotation
filtered_annotation=$(echo ${WD}sqanti3_results/filter_results/${sample_name}/*.filtered.gtf)

# Get the reference transcriptome classification file
refclassif=$(echo ${WD}sqanti3_results/qc_results/reference_transcriptome/${sample_name}/*.classification.txt)

# Run SQANTI3 rescue module
python /home/carblan/software/SQANTI3/sqanti3_rescue.py rules \
       --isoforms "${isoforms_fasta}" \
       --gtf "${filtered_annotation}" \
       --refGTF "${reference_annotation}" \
       --refGenome "${reference_genome}" \
       --refClassif "${refclassif}" \
	--output "${gtf_file_name}" \
       --dir "${WD}sqanti3_results/rescue_results/${sample_name}" \
       $filter_classification_file

echo "Isoform rescue for ${gtf_file_name} finished. SQANTI3 pipeline completed."

pwd;date