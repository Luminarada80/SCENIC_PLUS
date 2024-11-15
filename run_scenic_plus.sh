#!/bin/bash -l

#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem-per-cpu=16G
#SBATCH -o LOGS/scenic_plus.log
#SBATCH -e LOGS/scenic_plus.err

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

# Function to run each Python step with timing and memory tracking
run_step() {
  step_name=$1
  script_path=$2
  shift 2  # Remove the first two arguments
  echo "Running $step_name..."
  
  # Use /usr/bin/time with verbose output to capture memory and timing
  /usr/bin/time -v python3 "$script_path" "$@" 2>> "${LOG_DIR}/${step_name}_time_mem.log"
}

# run_step "step00_RNA_preprocessing.py" "${SCRIPT_DIR}/step00_RNA_preprocessing.py"

# run_step "step01_preprocessing_pseudobulk_profiles.py" "${SCRIPT_DIR}/step01_ATAC_preprocessing_pseudobulk_profiles.py"

# run_step "step02_ATAC_infering_consensus_peaks.py" "${SCRIPT_DIR}/step02_ATAC_infering_consensus_peaks.py"

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

# run_step "step03_ATAC_qc.py" "${SCRIPT_DIR}/step03_ATAC_qc.py"

# run_step "step04_ATAC_creating_cistopic_object.py" "${SCRIPT_DIR}/step04_ATAC_creating_cistopic_object.py"

run_step "step05_ATAC_run_models.py" "${SCRIPT_DIR}/step05_ATAC_run_models.py"

run_step "step06_ATAC_model_selection.py" "${SCRIPT_DIR}/step06_ATAC_model_selection.py"

run_step "step07_ATAC_topic_binarization_to_save_region_sets.py" "${SCRIPT_DIR}/step07_ATAC_topic_binarization_to_save_region_sets.py"
