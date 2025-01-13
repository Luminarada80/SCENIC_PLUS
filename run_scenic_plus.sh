#!/bin/bash -l

###############################################################################
# SLURM DIRECTIVES
###############################################################################
#SBATCH --job-name=SCENIC+
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=50G
#SBATCH -o LOGS/scenic_plus_test.log
#SBATCH -e LOGS/scenic_plus_test.err

###############################################################################
# ENVIRONMENT SETUP
###############################################################################
set -e

conda activate scenicplus
cd /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS

echo "Python executable: $(which python)"
echo ""

# Add the local pycisTopic to the python path so it is recognized as a module
export PYTHONPATH="${PYTHONPATH}:/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/pycisTopic/src"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/miniconda3/lib

NUM_CPU=${SLURM_CPUS_PER_TASK}
echo "Number of CPUs allocated: ${NUM_CPU}"
echo ""

###############################################################################
# DECIDE WHICH STEPS TO RUN
###############################################################################
STEP_01_RNA_PREPROCESSING=false
STEP_02_ATAC_PREPROCESSING=false
STEP_03_GET_TSS_DATA=false
STEP_04_CREATE_FASTA=false
STEP_05_CREATE_CISTARGET_MOTIF_DATABASES=false
STEP_06_RUN_SNAKEMAKE_PIPELINE=true

# Optional: Use precomputed cisTarget database
USE_PRECOMPUTED_CISTARGET_DB=true

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
CHROMSIZES="${INPUT_DIR}/hg38.chrom.sizes"
GENOME_FASTA="${INPUT_DIR}/hg38.fa"
RNA_FILE_NAME="K562_human_filtered_RNA.csv"
ATAC_FILE_NAME="K562_human_filtered_ATAC.csv"
MM10_BLACKLIST="${SCRIPT_DIR}/pycisTopic/blacklist/hg38-blacklist.v2.bed"
REGION_BED="${OUTPUT_DIR}/consensus_peak_calling/consensus_regions.bed"
FASTA_FILE="${INPUT_DIR}/hg38.K562.with_1kb_bg_padding.fa"

echo "Input files:"
echo "    RNA Data File: $RNA_FILE_NAME"
echo "    ATAC Data File: $ATAC_FILE_NAME"
echo "    Genome FASTA File: $GENOME_FASTA"
echo ""

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
        2> "${LOG_DIR}/${step_name}.log"
}

run_bash_step() {
    step_name=$1
    script_path=$2
    shift 2  # Remove the first two arguments
    echo "Running ${step_name}..."
    /usr/bin/time -v "${script_path}" "$@" \
        2> "${LOG_DIR}/${step_name}.log"
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
# CHECK PATHS AND DEPENDENCIES
###############################################################################

echo "Checking for all required directories and files"
# Check required directories
check_dir_exists "$CISTARGET_SCRIPT_DIR"
check_dir_exists "$SCRIPT_DIR"

# Check required files
check_file_exists "$INPUT_DIR/$ATAC_FILE_NAME"
check_file_exists "$INPUT_DIR/$RNA_FILE_NAME"
check_file_exists "$MM10_BLACKLIST"
check_file_exists "$GENOME_FASTA"

# Check to see if SCENIC+ generated directories exist or create them
check_or_create_dir "$LOG_DIR"
check_or_create_dir "$OUTPUT_DIR"
check_or_create_dir "$TEMP_DIR"
check_or_create_dir "$QC_DIR"

echo "    All required files and directories found"
echo ""

echo "Checking to see if Cluster-Buster is in the PATH"
# Check if 'cbust' is in the PATH
if ! command -v cbust &> /dev/null; then
    echo "    'cbust' not found in PATH. Setting it up..."

    # Download the cbust if its not in the script path
    if [ ! -f "${SCRIPT_DIR}/cbust" ]; then
        echo "    cbust not downloaded, downloading..."
        # Download cluster-buster (cbust)
        wget https://resources.aertslab.org/cistarget/programs/cbust -O cbust

    else
        echo "    cbust file found, adding to PATH"
    fi

    # Make it executable
    chmod a+x cbust

    # Add it to the PATH
    export PATH=$(pwd):$PATH
    echo "    'cbust' has been added to PATH."

else
    echo "    'cbust' is already in PATH."
fi
echo ""

echo "Checking if the Aertslab motif collection is downloaded"
# Check for the motif collection directory and file or download
MOTIF_DIR="${SCRIPT_DIR}/aertslab_motif_colleciton"
MOTIF_ZIP="${MOTIF_DIR}/v10nr_clust_public.zip"
MOTIF_URL="https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip"

# Check if the motif collection already exists
if [ ! -d "${MOTIF_DIR}/v10nr_clust_public" ]; then
    echo "    Motif collection not found. Downloading and extracting it now..."

    # Create the motif collection directory if it doesn't exist
    mkdir -p "${MOTIF_DIR}"

    # Download the motif collection zip file
    wget -O "${MOTIF_ZIP}" "${MOTIF_URL}"

    # Extract the zip file
    unzip -q "${MOTIF_ZIP}" -d "${MOTIF_DIR}"

    echo "    Motif collection downloaded and extracted to ${MOTIF_DIR}/v10nr_clust_public."
else
    echo "    Motif collection exists at ${MOTIF_DIR}/v10nr_clust_public."
fi
echo ""

# Update the n_cpu value in the config.yaml file
CONFIG_FILE="${SCRIPT_DIR}/scplus_pipeline/Snakemake/config/config.yaml"
sed -i "s/^\(\s*n_cpu:\s*\).*/\1${NUM_CPU}/" "${CONFIG_FILE}"
echo "Updated config.yaml with n_cpu: ${NUM_CPU}"
echo ""

# Download the precomputed cisTarget database
if [ "$USE_PRECOMPUTED_CISTARGET_DB" = true ] && [ "$STEP_05_CREATE_CISTARGET_MOTIF_DATABASES" = false ]; then

    echo "Using pre-computed cisTarget database"
    if [ -f "$INPUT_DIR/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather" ]; then
        echo "    Precomputed cisTarget hg38 regions_vs_motifs.rankings.feather file exists"
    else
        # Optional: Download the precomputed cisTarget database
        echo "    Downloading precomputed cisTarget database: rankings.feather file"
        curl -o -s -q  "$INPUT_DIR/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather" \
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
        echo "        Done!" 
    fi

    if [ -f "$INPUT_DIR/hg38_screen_v10_clust.regions_vs_motifs.scores.feather" ]; then
        echo "    Precomputed cisTarget hg38 regions_vs_motifs.scores.feather file exists"
    else
        # Optional: Download the precomputed cisTarget database
        echo "    Downloading precomputed cisTarget database: scores.feather file"
        curl -o -s -q "$INPUT_DIR/hg38_screen_v10_clust.regions_vs_motifs.scores.feather" \
            "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
        echo "        Done!"
    fi
    echo ""
fi

echo "Checks complete, starting pipeline"

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
        --name "hsapiens_gene_ensembl" \
        --to-chrom-source ucsc \
        --ucsc hg38 > "${LOG_DIR}/Step 3: Getting Transcription Start Site data.log" 2>&1;
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
        -t "${NUM_CPU}"
fi



###############################################################################
# STEP 6: RUN SNAKEMAKE PIPELINE
###############################################################################
if [ "$STEP_06_RUN_SNAKEMAKE_PIPELINE" = true ]; then
    echo "Step 6: Run SCENIC+ snakemake"
    

    # Check for Snakemake lock files
    LOCK_FILE="$SCRIPT_DIR/scplus_pipeline/Snakemake/.snakemake/locks"
    if [ -d "$LOCK_FILE" ]; then
        echo "    Lock files detected at $LOCK_FILE"

        # Use the SLURM job name for comparison
        JOB_NAME="${SLURM_JOB_NAME}"  # Dynamically retrieve the job name from SLURM

        # Check for running jobs with the same name, excluding the current job
        RUNNING_COUNT=$(squeue --name="$JOB_NAME" --noheader | wc -l)

        # If other SCENIC+ jobs are running and there are lock files, exit
        if [ "$RUNNING_COUNT" -gt 1 ]; then
            echo "    A job with the name '"$JOB_NAME"' is already running:"
            echo "    Exiting to avoid conflicts."
            exit 1
        
        # If no other SCENIC+ jobs are running and there are lock files, remove them
        else
            echo "    No other jobs with the name '"$JOB_NAME"', removing locks"
            rm -rf "$LOCK_FILE"
        fi
        
    else
        echo "    No lock files detected."
    fi

    echo "    Running snakemake"

    cd "${SCRIPT_DIR}/scplus_pipeline/Snakemake"
    snakemake --cores ${NUM_CPU} --latency-wait 600 > "${LOG_DIR}/Step 6: Snakemake.log" 2>&1;
fi