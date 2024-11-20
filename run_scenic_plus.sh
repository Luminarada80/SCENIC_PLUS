#!/bin/bash -l

#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem-per-cpu=16G
#SBATCH -o LOGS/scenic_plus_test.log
#SBATCH -e LOGS/scenic_plus_tes.err
#srun source /gpfs/Home/esm5360/miniconda3/envs/scenicplus

conda activate scenicplus
cd /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS

echo $(which python)

# Add the local pycisTopic to the python path so it's recognized as a module
export PYTHONPATH=$PYTHONPATH:/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/pycisTopic/src

LOG_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/LOGS"
SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS"
OUT_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/scenicplus/mESC_new_scenicplus/outs/"
QC_DIR="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/outs/qc"
DATA_DIR="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.SC_MO_TRN_BENCHMARKING.MIRA/SCENIC_PLUS.HABIBA/data"

# Paths to create a custom cistarget database
CISTARGET_SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/create_cisTarget_databases/"
GENOME_FASTA="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/mm10.fa"
CHROMSIZES="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/mm10.chrom.sizes"
REGION_BED="/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/scenicplus/mESC_new_scenicplus/outs/consensus_peak_calling/consensus_regions.bed"

# Function to run each Python step with timing and memory tracking
run_python_step() {
  step_name=$1
  script_path=$2
  shift 2  # Remove the first two arguments
  echo "Running $step_name..."
  
  # Use /usr/bin/time with verbose output to capture memory and timing
  /usr/bin/time -v /gpfs/Home/esm5360/miniconda3/envs/scenicplus/bin/python "$script_path" "$@" 2>> "${LOG_DIR}/${step_name}_time_mem.log"
}

run_bash_step() {
  step_name=$1
  script_path=$2
  shift 2  # Remove the first two arguments
  echo "Running $step_name..."
  
  # Use /usr/bin/time with verbose output to capture memory and timing
  /usr/bin/time -v "$script_path" "$@" 2>> "${LOG_DIR}/${step_name}_time_mem.log"
}

# run_python_step "step00_RNA_preprocessing" "${SCRIPT_DIR}/step00_RNA_preprocessing.py"

# run_python_step "step01_preprocessing_pseudobulk_profiles.py" "${SCRIPT_DIR}/step01_ATAC_preprocessing_pseudobulk_profiles.py"

# run_python_step "step02_ATAC_infering_consensus_peaks.py" "${SCRIPT_DIR}/step02_ATAC_infering_consensus_peaks.py"

# run_python_step "step01_ATAC_preprocessing" "${SCRIPT_DIR}/ATAC_preprocessing.py"

# Get the mm10 TSS data
# pycistopic tss get_tss \
#   --output "${QC_DIR}/tss.bed" \
#   --name "mmusculus_gene_ensembl" \
#   --to-chrom-source ucsc \
#   --ucsc mm10

# # Run pycistopic qc
# pycistopic qc \
#   --fragments "${DATA_DIR}/GSM6205427_E7.5_rep1_ATAC_fragments.tsv.gz" \
#   --regions "${OUT_DIR}/consensus_peak_calling/consensus_regions.bed" \
#   --tss "${QC_DIR}/tss.bed" \
#   --output "${QC_DIR}/mESC"

# run_python_step "step03_ATAC_qc.py" "${SCRIPT_DIR}/step03_ATAC_qc.py"

# run_python_step "step04_ATAC_creating_cistopic_object.py" "${SCRIPT_DIR}/step04_ATAC_creating_cistopic_object.py"

# run_python_step "step05_ATAC_run_models.py" "${SCRIPT_DIR}/step05_ATAC_run_models.py"

# run_python_step "step06_ATAC_model_selection.py" "${SCRIPT_DIR}/step06_ATAC_model_selection.py"

# run_python_step "step07_ATAC_topic_binarization_to_save_region_sets.py" "${SCRIPT_DIR}/step07_ATAC_topic_binarization_to_save_region_sets.py"

# Make sure that you have downloaded the create_cisTarget_databases directory from GitHub
# git clone https://github.com/aertslab/create_cisTarget_databases

# Download cluster-buster
# wget https://resources.aertslab.org/cistarget/programs/cbust > chmod a+x cbust
# export PATH=$(pwd):$PATH

# Download the motif collection
# mkdir -p "${SCRIPT_DIR}/aertslab_motif_colleciton"

# Download the motif collection to the folder
# wget -O "${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public.zip https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip"
# unzip -q "${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public.zip"

FASTA_FILE="${OUT_DIR}/mm10.mESC.with_1kb_bg_padding.fa"

# Ensure the directory exists
mkdir -p "$(dirname "${FASTA_FILE}")"

# Create an empty file
touch "${FASTA_FILE}"

module load bedtools/2.31.0

#Creates the custom cistarget database
run_bash_step "test_step08_ATAC_preprocessing" "${CISTARGET_SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
    ${GENOME_FASTA} \
    ${CHROMSIZES} \
    ${REGION_BED} \
    ${FASTA_FILE} \
    1000 \
    yes \

# Creating the ranking and score databases
CBDIR="${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
MOTIF_LIST="/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/motifs.txt"

export PATH=${SCRIPT_DIR}:$PATH

run_python_step "test_step09_creating_cistarget_databases" "${CISTARGET_SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    --min 5 \
    --max 1000 \
    -o ${OUT_DIR} \
    --bgpadding 1000 \
    -t 16