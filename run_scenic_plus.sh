#!/bin/bash -l

###############################################################################
# SLURM DIRECTIVES
###############################################################################
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=16G
#SBATCH -o LOGS/scenic_plus_test.log
#SBATCH -e LOGS/scenic_plus_test.err

###############################################################################
# ENVIRONMENT SETUP
###############################################################################
set -e

conda activate scenicplus
cd /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS

echo "Python executable: $(which python)"

# Add the local pycisTopic to the python path so it is recognized as a module
export PYTHONPATH="${PYTHONPATH}:/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/pycisTopic/src"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib

###############################################################################
# DECIDE WHICH STEPS TO RUN
###############################################################################
STEP_01_RNA_PREPROCESSING=false
STEP_02_ATAC_PREPROCESSING=true
STEP_03_GET_TSS_DATA=true
STEP_04_CREATE_FASTA=true
STEP_05_CREATE_CISTARGET_MOTIF_DATABASES=TRUE
STEP_06_RUN_SNAKEMAKE_PIPELINE=TRUE

###############################################################################
# PATH SETUP
###############################################################################
LOG_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/LOGS"
SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS"
CISTARGET_SCRIPT_DIR="${SCRIPT_DIR}/create_cisTarget_databases"
TEMP_DIR="${SCRIPT_DIR}/tmp"
INPUT_DIR="${SCRIPT_DIR}/input"
OUTPUT_DIR="${SCRIPT_DIR}/outs"
QC_DIR="${OUTPUT_DIR}/qc"

###############################################################################
# INPUT FILES & DIRECTORIES
###############################################################################
CHROMSIZES="${INPUT_DIR}/mm10.chrom.sizes"
GENOME_FASTA="${INPUT_DIR}/mm10.fa"
RNA_FILE_NAME="multiomic_data_1000_cells_E7.5_rep1_RNA.csv"
ATAC_FILE_NAME="multiomic_data_1000_cells_E7.5_rep1_ATAC.csv"
MM10_BLACKLIST="${SCRIPT_DIR}/pycisTopic/blacklist/mm10-blacklist.v2.bed"
REGION_BED="${OUTPUT_DIR}/consensus_peak_calling/consensus_regions.bed"
FASTA_FILE="${INPUT_DIR}/mm10.mESC.with_1kb_bg_padding.fa"

echo "RNA Data File: $RNA_FILE_NAME"
echo "ATAC Data File: $ATAC_FILE_NAME"
echo "Genome FASTA File: $GENOME_FASTA"

###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
run_python_step() {
    step_name=$1
    script_path=$2
    shift 2  # Remove the first two arguments
    echo "Running ${step_name}..."
    /usr/bin/time -v /gpfs/Home/esm5360/miniconda3/envs/scenicplus/bin/python \
        "${script_path}" "$@" \
        2>> "${LOG_DIR}/${step_name}_time_mem.log"
}

run_bash_step() {
    step_name=$1
    script_path=$2
    shift 2  # Remove the first two arguments
    echo "Running ${step_name}..."
    /usr/bin/time -v "${script_path}" "$@" \
        2>> "${LOG_DIR}/${step_name}_time_mem.log"
}

# Function to check if a directory exists, and create it if it doesn't
check_or_create_dir() {
    local dir_path="$1"
    if [ ! -d "$dir_path" ]; then
        echo "Directory '$dir_path' does not exist. Creating it now..."
        mkdir -p "$dir_path"
    fi
}

# Function to check if a file exists, and exit with an error if it doesn't
check_file_exists() {
    local file_path="$1"
    if [ ! -f "$file_path" ]; then
        echo "Error: File '$file_path' does not exist!"
        exit 1
    fi
}

# Function to check if a directory exists, and exit with an error if it doesn't
check_dir_exists() {
    local dir_path="$1"
    if [ ! -d "$dir_path" ]; then
        echo "Error: Directory '$dir_path' does not exist!"
        exit 1
    fi
}

###############################################################################
# ADDITIONAL NOTES
###############################################################################
# 1. Download cluster-buster (cbust):
#    wget https://resources.aertslab.org/cistarget/programs/cbust
#    chmod a+x cbust
#    export PATH=$(pwd):$PATH
#
# 2. Download the motif collection:
#    mkdir -p "${SCRIPT_DIR}/aertslab_motif_colleciton"
#    wget -O "${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public.zip" \
#         "https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip"
#    unzip -q "${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public.zip"
#
# 3. For region_to_gene, ensure dask is a compatible version:
#    pip install "dask==2024.5.0"  (Ignore pip’s possible upgrade error)
#
# 4. Make sure bedtools is loaded if needed:
#    module load bedtools/2.31.0

###############################################################################
# CHECK PATHS
###############################################################################

# Check required directories
check_dir_exists "$CISTARGET_SCRIPT_DIR"
check_dir_exists "$SCRIPT_DIR"

# Check required files
check_file_exists "$ATAC_FILE_NAME"
check_file_exists "$RNA_FILE_NAME"
check_file_exists "$MM10_BLACKLIST"
check_file_exists "$GENOME_FASTA"

# Check to see if SCENIC+ generated directories exist or create them
check_or_create_dir "$LOG_DIR"
check_or_create_dir "$OUTPUT_DIR"
check_or_create_dir "$TEMP_DIR"
check_or_create_dir "$QC_DIR"


###############################################################################
# STEP 1: RNA PREPROCESSING
###############################################################################
if [ "$STEP_01_RNA_PREPROCESSING" = true ]; then
    run_python_step "Step 1: RNA preprocessing" "${SCRIPT_DIR}/Step01.RNA_preprocessing.py" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --rna_file_name "${RNA_FILE_NAME}";
fi

###############################################################################
# STEP 2: ATAC PREPROCESSING
###############################################################################
if [ "$STEP_02_ATAC_PREPROCESSING" = true ]; then
    run_python_step "Step 2: ATAC preprocessing" "${SCRIPT_DIR}/Step02.ATAC_preprocessing.py" \
        --input_dir "${INPUT_DIR}" \
        --output_dir "${OUTPUT_DIR}" \
        --tmp_dir "${TEMP_DIR}" \
        --atac_file_name "${ATAC_FILE_NAME}" \
        --mm10_blacklist "${MM10_BLACKLIST}";
fi

###############################################################################
# STEP 3: GET TSS DATA
###############################################################################
if [ "$STEP_03_GET_TSS_DATA" = true ]; then
    echo "Step 3: Getting Transcription Start Site data"
    pycistopic tss get_tss \
        --output "${QC_DIR}/tss.bed" \
        --name "mmusculus_gene_ensembl" \
        --to-chrom-source ucsc \
        --ucsc mm10;
fi

###############################################################################
# STEP 4: CREATE FASTA WITH PADDED BACKGROUND
###############################################################################
if [ "$STEP_04_CREATE_FASTA" = true ]; then
    module load bedtools/2.31.0
    run_bash_step "Step 4: Prepare fasta from consensus regions" \
        "${CISTARGET_SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh" \
        "${GENOME_FASTA}" \
        "${CHROMSIZES}" \
        "${REGION_BED}" \
        "${FASTA_FILE}" \
        1000 \
        yes;
fi

###############################################################################
# STEP 5: CREATE CISTARGET MOTIF DATABASES
###############################################################################
if [ "$STEP_05_CREATE_CISTARGET_MOTIF_DATABASES" = true ]; then
    CBDIR="${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
    MOTIF_LIST="${INPUT_DIR}/motifs.txt"

    run_python_step "Step 5: Create cistarget databases" \
        "${CISTARGET_SCRIPT_DIR}/create_cistarget_motif_databases.py" \
        -f "${FASTA_FILE}" \
        -M "${CBDIR}" \
        -m "${MOTIF_LIST}" \
        --min 5 \
        --max 1000 \
        -o "${OUTPUT_DIR}" \
        --bgpadding 1000 \
        -t 32;
fi

###############################################################################
# STEP 6: RUN SNAKEMAKE PIPELINE
###############################################################################
if [ "$STEP_06_RUN_SNAKEMAKE_PIPELINE" = true ]; then
    cd "${SCRIPT_DIR}/scplus_pipeline/Snakemake"
    echo "Step 6: Run SCENIC+ snakemake"
    snakemake --cores 32 --latency-wait 600 >> "${LOG_DIR}/Step06.snakemake_time_mem.log"
fi