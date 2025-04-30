#!/bin/bash -l
#SBATCH --job-name="SCENIC+_${CELL_TYPE}_${SAMPLE_NAME}"
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 5
#SBATCH --mem=64G
LOG_DIR="${SCRIPT_DIR}/LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs/"

# DECIDE WHICH STEPS TO RUN
STEP_01_RNA_PREPROCESSING=false
STEP_02_ATAC_PREPROCESSING=false
STEP_03_GET_TSS_DATA=false
STEP_04_CREATE_FASTA=false
USE_PRECOMPUTED_CISTARGET_DB=false
STEP_06_RUN_SNAKEMAKE_PIPELINE=false
STEP_07_FORMAT_INFERRED_GRN=true

CONDA_ENV_NAME="test_scenicplus"

###############################################################################
# ENVIRONMENT SETUP
###############################################################################
set -euo pipefail
trap "echo 'An error occurred. Exiting...'; exit 1;" ERR

###############################################################################
# INPUT FILES & DIRECTORIES
###############################################################################


OUTPUT_DIR="${SCRIPT_DIR}/output/${CELL_TYPE}_${SAMPLE_NAME}_outs"
REGION_BED="${OUTPUT_DIR}/consensus_peak_calling/consensus_regions.bed"
CISTARGET_SCRIPT_DIR="${SCRIPT_DIR}/create_cisTarget_databases"
TEMP_DIR="${SCRIPT_DIR}/tmp/${CELL_TYPE}_${SAMPLE_NAME}_tmp"
QC_DIR="${OUTPUT_DIR}/qc"
MOTIF_DATABASE_DIR="${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/snapshots"

if [ $SPECIES == "human" ]; then
    ORGANISM_DIR="${SCRIPT_DIR}/organism_genome_files/human"
    ENSEMBL_SPECIES="hsapiens"
    MOTIF_ENRICHMENT_SPECIES="homo_sapiens"
    PYCISTOPIC_SPECIES="hsapiens_gene_ensembl"
    PYCISTOPIC_SPECIES_CODE="hg38"
    ANNOTATION_VERSION="v10nr_clust"
    CHROMSIZES="${ORGANISM_DIR}/hg38.chrom.sizes"
    GENOME_FASTA="${ORGANISM_DIR}/hg38.fa"
    FASTA_FILE="${INPUT_DIR}/hg38.${CELL_TYPE}.with_1kb_bg_padding.fa"
    MOTIF_ANNOT_FILE="${MOTIF_DATABASE_DIR}/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl"
    CISTARGET_RANKINGS_PRECOMP="hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
    CISTARGET_SCORES_PRECOMP="hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
    BLACKLIST="${SCRIPT_DIR}/pycisTopic/blacklist/hg38-blacklist.v2.bed"
fi

if [ $SPECIES == "mouse" ]; then
    ORGANISM_DIR="${SCRIPT_DIR}/organism_genome_files/mouse"
    ENSEMBL_SPECIES="mmusculus"
    MOTIF_ENRICHMENT_SPECIES="mus_musculus"
    PYCISTOPIC_SPECIES="mmusculus_gene_ensembl"
    PYCISTOPIC_SPECIES_CODE="mm10"
    ANNOTATION_VERSION="v10nr_clust"
    CHROMSIZES="${ORGANISM_DIR}/mm10.chrom.sizes"
    GENOME_FASTA="${ORGANISM_DIR}/mm10.fa"
    FASTA_FILE="${INPUT_DIR}/mm10.${CELL_TYPE}.with_1kb_bg_padding.fa"
    MOTIF_ANNOT_FILE="${MOTIF_DATABASE_DIR}/motifs-v10-nr.mgi-m0.00001-o0.0.tbl"
    CISTARGET_RANKINGS_PRECOMP="mm10_screen_v10_clust.regions_vs_motifs.rankings.feather"
    CISTARGET_SCORES_PRECOMP="mm10_screen_v10_clust.regions_vs_motifs.scores.feather"
    BLACKLIST="${SCRIPT_DIR}/pycisTopic/blacklist/mm10-blacklist.v2.bed"
fi

echo "Input files:"
echo "    RNA Data File: $RNA_FILE_NAME"
echo "    ATAC Data File: $ATAC_FILE_NAME"
echo "    Genome FASTA File: $GENOME_FASTA"
echo "    Cell Type: $CELL_TYPE"
echo "    Sample: $SAMPLE_NAME"
echo "    Species: $SPECIES"

###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
check_if_running() {
    echo ""
    echo "[INFO] Checking for SCENIC+ jobs running for the current sample..."
    # Use the SLURM job name for comparison
    JOB_NAME="${SLURM_JOB_NAME}"  # Dynamically retrieve the job name from SLURM

    # Check for running jobs with the same name, excluding the current job
    RUNNING_COUNT=$(squeue --name="$JOB_NAME" --noheader | wc -l)

    # If other SCENIC+ jobs are running and there are lock files, exit
    if [ "$RUNNING_COUNT" -gt 1 ]; then
        echo ""
        echo "[WARNING] A job with the name '"$JOB_NAME"' is already running:"
        echo "    Exiting to avoid conflicts."
        exit 1
    else
        echo "    - No jobs with the same name running, continuing"
    fi
}

determine_num_cpus() {
    echo ""
    echo "[INFO] Checking the number of CPUs available for parallel processing"
    if [ -z "${SLURM_CPUS_PER_TASK:-}" ]; then
        if command -v nproc &> /dev/null; then
            TOTAL_CPUS=$(nproc --all)
            case $TOTAL_CPUS in
                [1-15]) IGNORED_CPUS=1 ;;  # Reserve 1 CPU for <=15 cores
                [16-31]) IGNORED_CPUS=2 ;; # Reserve 2 CPUs for <=31 cores
                *) IGNORED_CPUS=4 ;;       # Reserve 4 CPUs for >=32 cores
            esac
            NUM_CPU=$((TOTAL_CPUS - IGNORED_CPUS))
            echo "    - Running locally. Detected $TOTAL_CPUS CPUs, reserving $IGNORED_CPUS for system tasks. Using $NUM_CPU CPUs."
        else
            NUM_CPU=1  # Fallback
            echo "    - Running locally. Unable to detect CPUs, defaulting to $NUM_CPU CPU."
        fi
    else
        NUM_CPU=${SLURM_CPUS_PER_TASK}
        echo "    - Running on SLURM. Number of CPUs allocated: ${NUM_CPU}"
    fi
}

run_python_step() {
    step_name=$1
    script_path=$2
    shift 2  # Remove the first two arguments
    echo "Running ${step_name}..."
    /usr/bin/time -v python3 \
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

activate_conda_env() {
    echo ""
    echo "[INFO] Attempting to load the specified Conda module"
    CONDA_BASE=$(conda info --base)
    if [ -z "$CONDA_BASE" ]; then
        echo ""
        echo "[ERROR] Conda base could not be determined. Is Conda installed and in your PATH?"
        exit 1
    fi

    source "$CONDA_BASE/bin/activate"
    if ! conda env list | grep -q "^$CONDA_ENV_NAME "; then
        echo ""
        echo "[ERROR] Conda environment '$CONDA_ENV_NAME' does not exist."
        echo "   - Attempting to create $CONDA_ENV_NAME environment..."
        # redirect both stdout and stderr into create_env.err
        {
        conda create --name "$CONDA_ENV_NAME" python=3.11.8 -y 
        conda activate "$CONDA_ENV_NAME"
        pip install "pip<24.1" 
        pip install ${SCRIPT_DIR}
        } > "${LOG_DIR}/create_conda_env.out" 2> "${LOG_DIR}/create_conda_env.out"
        if [[ $? -ne 0 ]]; then
            echo "[ERROR] Failed to create Conda environment, see ${LOG_DIR}/create_env.err for details."
            exit 1
        fi
    fi

    conda activate "$CONDA_ENV_NAME" || { echo "Error: Failed to activate Conda environment '$CONDA_ENV_NAME'."; exit 1; }
    echo "   - Successfully activated Conda environment: $CONDA_ENV_NAME"
    echo "   - Python executable: $(which python)"
}

install_scenic_plus() {
    echo ""
    echo "[INFO] Checking that scenicplus is installed"
    local repo_dir="${SCRIPT_DIR}/src/scenicplus"
    local logf="${LOG_DIR}/install_scenic_plus.log"

    python3 -c "import scenicplus" 2>/dev/null
    local is_installed=$?

    if [[ $is_installed -ne 0 || ! -d "$repo_dir" ]]; then
        echo "    - scenicplus not found in Python or directory missing, installing..."
        mkdir -p "$(dirname "$logf")"

        {
            echo "---- Cloning scenicplus ----"
            git clone https://github.com/Luminarada80/scenicplus_src.git "$repo_dir"
            echo "---- Installing scenicplus ----"
            pip install "${SCRIPT_DIR}"
        } > "$logf" 2>&1

        if [[ $? -ne 0 ]]; then
            echo "[ERROR] scenicplus install failed, see $logf"
            exit 1
        else
            echo "    - scenicplus installed; logs in $logf"
        fi
    else
        echo "    - scenicplus is already installed and directory exists"
    fi
}


install_pycistopic() {
    echo ""
    echo "[INFO] Checking that pycisTopic is installed"
    local repo_dir="${SCRIPT_DIR}/pycisTopic"
    local logf="${LOG_DIR}/install_pycistopic.log"

    if [[ ! -d "$repo_dir" ]]; then
        echo "    - pycisTopic directory not found, installing..."
        # make sure log directory exists
        mkdir -p "$(dirname "$logf")"

        # run both clone and install, capturing all output
        {
            echo "---- Cloning pycisTopic ----"
            git clone https://github.com/aertslab/pycisTopic.git "$repo_dir"
            echo "---- Installing pycisTopic ----"
            pip install -e "$repo_dir"
        } > "$logf" 2>&1

        if [[ $? -ne 0 ]]; then
            echo "[ERROR] pycisTopic install failed, see $logf"
            exit 1
        else
            echo "    - pycisTopic installed; log in $logf"
        fi
    else
        echo "    - pycisTopic directory exists"
    fi
}

add_pycistopic_to_path(){
    echo ""
    echo "[INFO] Adding pycisTopic to PATH"
    # Add the local pycisTopic to the python path so it is recognized as a module
    if [ -z "${PYTHONPATH+x}" ]; then
        export PYTHONPATH="${SCRIPT_DIR}/pycisTopic/src"
    else
        export PYTHONPATH="${PYTHONPATH}:${SCRIPT_DIR}/pycisTopic/src"
    fi

    if [ -z "${LD_LIBRARY_PATH+x}" ]; then
        export LD_LIBRARY_PATH="$HOME/miniconda3/lib"
    else
        export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$HOME/miniconda3/lib"
    fi
}

ensure_python_pkg() {
    local pkg="$1"
    if ! python3 - <<EOF 2>/dev/null
try:
    __import__("$pkg")
except ImportError:
    raise
EOF
    then
        echo "    - Python package '$pkg' not found; attempting to install…"
        # try conda first
        if conda install -y "$pkg" >> "$LOG_DIR/python_deps.log" 2>&1; then
            echo "        - Installed '$pkg' via conda"
        else
            echo "        - [WARN] Conda install failed for '$pkg'; falling back to pip"
            if pip install "$pkg" >> "$LOG_DIR/python_deps.log" 2>&1 ; then
                echo "            - Installed '$pkg' via pip"
            else
                echo ""
                echo "[ERROR] Could not install python package '$pkg'"
                exit 1
            fi
        fi
    else
        echo "    - Python package '$pkg' is already installed"
    fi
}

check_python_deps() {
    echo ""
    echo "[INFO] Checking python package requirements"
    local pkgs=(ruamel.yaml requests numpy pandas scanpy mudata)
    for p in "${pkgs[@]}"; do
        ensure_python_pkg "$p"
    done
} 

# Function to check if a directory exists, and create it if it doesn't
check_or_create_dir() {
    local dir_path="$1"
    if [ ! -d "$dir_path" ]; then
        echo "    - Directory '$dir_path' does not exist. Creating it now..."
        mkdir -p "$dir_path"
    fi
}

# Function to check if a file exists, and exit with an error if it doesn't
check_file_exists() {
    local file_path="$1"
    if [ ! -f "$file_path" ]; then
        echo "    - Error: File '$file_path' does not exist!"
        exit 1
    fi
}

# Function to check if a directory exists, and exit with an error if it doesn't
check_dir_exists() {
    local dir_path="$1"
    if [ ! -d "$dir_path" ]; then
        echo "    - Error: Directory '$dir_path' does not exist!"
        exit 1
    fi
}

generate_config() {
    python3 update_config_yaml.py \
    --cisTopic_obj_fname "${OUTPUT_DIR}/cistopic_obj.pkl" \
    --GEX_anndata_fname "${OUTPUT_DIR}/adata_final.h5ad" \
    --region_set_folder "${OUTPUT_DIR}/region_sets" \
    --ctx_db_fname "${ORGANISM_DIR}/${CISTARGET_RANKINGS_PRECOMP}" \
    --dem_db_fname "${ORGANISM_DIR}/${CISTARGET_SCORES_PRECOMP}" \
    --path_to_motif_annotations "${MOTIF_ANNOT_FILE}" \
    --combined_GEX_ACC_mudata "${OUTPUT_DIR}/ACC_GEX.h5mu" \
    --dem_result_fname "${OUTPUT_DIR}/dem_results.hdf5" \
    --ctx_result_fname "${OUTPUT_DIR}/ctx_results.hdf5" \
    --output_fname_dem_html "${OUTPUT_DIR}/dem_results.html" \
    --output_fname_ctx_html "${OUTPUT_DIR}/ctx_results.html" \
    --cistromes_direct "${OUTPUT_DIR}/cistromes_direct.h5ad" \
    --cistromes_extended "${OUTPUT_DIR}/cistromes_extended.h5ad" \
    --tf_names "${OUTPUT_DIR}/tf_names.txt" \
    --genome_annotation "${OUTPUT_DIR}/genome_annotation.tsv" \
    --chromsizes "${OUTPUT_DIR}/chromsizes.tsv" \
    --search_space "${OUTPUT_DIR}/search_space.tsv" \
    --tf_to_gene_adjacencies "${OUTPUT_DIR}/tf_to_gene_adj.tsv" \
    --region_to_gene_adjacencies "${OUTPUT_DIR}/region_to_gene_adj.tsv" \
    --eRegulons_direct "${OUTPUT_DIR}/eRegulons_direct.tsv" \
    --eRegulons_extended "${OUTPUT_DIR}/eRegulons_extended.tsv" \
    --AUCell_direct "${OUTPUT_DIR}/AUCell_direct.h5mu" \
    --AUCell_extended "${OUTPUT_DIR}/AUCell_extended.h5mu" \
    --scplus_mdata "${OUTPUT_DIR}/scplusmdata.h5mu" \
    --temp_dir "${SCRIPT_DIR}/tmp/${CELL_TYPE}/${SAMPLE_NAME}" \
    --n_cpu "${NUM_CPU}" \
    --seed 666 \
    --ensembl_species "${ENSEMBL_SPECIES}" \
    --motif_enrichment_species "${MOTIF_ENRICHMENT_SPECIES}" \
    --annotation_version "${ANNOTATION_VERSION}" \
    --output_config_path "${SCRIPT_DIR}/scplus_pipeline/Snakemake/config/${CELL_TYPE}_${SAMPLE_NAME}_config.yaml"
}

check_clusterbuster(){
    echo "[INFO] Checking to see if Cluster-Buster is in the PATH"
    # Check if 'cbust' is in the PATH
    if ! command -v cbust &> /dev/null; then
        echo "    - 'cbust' not found in PATH. Setting it up..."

        # Download the cbust if its not in the script path
        if [ ! -f "${SCRIPT_DIR}/cbust" ]; then
            echo "    cbust not downloaded, downloading..."
            # Download cluster-buster (cbust)
            wget -q https://resources.aertslab.org/cistarget/programs/cbust -O "${SCRIPT_DIR}/cbust"

        else
            echo "        - 'cbust' file found, adding to PATH..."
        fi

        # Make it executable
        chmod +x "${SCRIPT_DIR}/cbust"

        # Add it to the PATH
        export PATH="${SCRIPT_DIR}:$PATH"
        echo "    - 'cbust' has been added to PATH."

    else
        echo "    - 'cbust' is already in PATH."
    fi
    echo ""
}

check_aertslab_motif_collection() {
    local logf="${LOG_DIR}/check_aertslab_motif_collection.log"

    # Ensure the log directory exists
    mkdir -p "$(dirname "$logf")"

    {
        echo "[INFO] Checking if the Aertslab motif collection is downloaded"

        MOTIF_DIR="${SCRIPT_DIR}/aertslab_motif_colleciton"
        MOTIF_ZIP="${MOTIF_DIR}/v10nr_clust_public.zip"
        MOTIF_URL="https://resources.aertslab.org/cistarget/motif_collections/v10nr_clust_public/v10nr_clust_public.zip"

        if [[ ! -d "${MOTIF_DIR}/v10nr_clust_public" ]]; then
            echo "    - Not found, creating ${MOTIF_DIR}"
            mkdir -p "${MOTIF_DIR}"

            echo "    - Downloading to ${MOTIF_ZIP}"
            wget -q -O "${MOTIF_ZIP}" "${MOTIF_URL}"

            echo "    - Extracting archive"
            unzip -q "${MOTIF_ZIP}" -d "${MOTIF_DIR}"

            echo "[INFO] Downloaded & extracted to ${MOTIF_DIR}/v10nr_clust_public"
        else
            echo "[INFO] Motif collection already at ${MOTIF_DIR}/v10nr_clust_public"
        fi

        echo ""  # blank line for readability
    } >"$logf" 2>&1

    if [[ $? -ne 0 ]]; then
        echo "[ERROR] check_aertslab_motif_collection failed; see $logf"
        exit 1
    else
        echo "[INFO] check_aertslab_motif_collection succeeded; details in $logf"
    fi
}


check_organism_genome_files(){
    echo "[INFO] Checking to see if the organism genome directory contains the correct files"
    if [ ! -d "${ORGANISM_DIR}" ]; then
        mkdir -p "$ORGANISM_DIR"
    fi

    if [ "$SPECIES" == "human" ]; then
        if [ ! -f "${ORGANISM_DIR}/hg38.chrom.sizes" ]; then
            echo "    - hg38.chrom.sizes does not exist, downloading..."
            ORGANISM_CHROM_SIZE_LINK="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
            curl -L -o "${ORGANISM_DIR}/hg38.chrom.sizes" "${ORGANISM_CHROM_SIZE_LINK}"
            echo "        Done!"
        fi
        if [ ! -f "${ORGANISM_DIR}/hg38.fa" ]; then
            echo "    - hg38.fa does not exist, downloading..."
            ORGANISM_FASTA_LINK="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
            curl -L -o "${ORGANISM_DIR}/hg38.fa.gz" "${ORGANISM_FASTA_LINK}" | gunzip "${ORGANISM_DIR}/hg38.fa.gz"
            echo "        Done!"
        fi
    fi

    if [ "$SPECIES" == "mouse" ]; then
        if [ ! -f "${ORGANISM_DIR}/mm10.chrom.sizes" ]; then
            echo "    - mm10.chrom.sizes does not exist, downloading..."
            ORGANISM_CHROM_SIZE_LINK="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"
            curl -L -o "${ORGANISM_DIR}/mm10.chrom.sizes" "${ORGANISM_CHROM_SIZE_LINK}"
            echo "        Done!"
        fi
        if [ ! -f "${ORGANISM_DIR}/mm10.fa" ]; then
            echo "    - mm10.fa does not exist, downloading..."
            ORGANISM_FASTA_LINK="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
            curl -L -o "${ORGANISM_DIR}/mm10.fa.gz" "${ORGANISM_FASTA_LINK}" | gunzip "${ORGANISM_DIR}/mm10.fa.gz"
            echo "        Done!"
        fi
    fi

    # Download the precomputed cisTarget database to the organism genome file
    if [ "$USE_PRECOMPUTED_CISTARGET_DB" = true ]; then

        echo "Using pre-computed cisTarget database"

        # Ensure destination directory exists
        mkdir -p "$INPUT_DIR"

        # File: rankings.feather
        if [ -f "${ORGANISM_DIR}/${CISTARGET_RANKINGS_PRECOMP}" ]; then
            echo "    Precomputed cisTarget ${PYCISTOPIC_SPECIES_CODE} regions_vs_motifs.rankings.feather file exists"
        else
            echo "    Downloading rankings.feather file..."
            if [ "$SPECIES" == "human" ]; then
                RANKINGS_FEATHER_LINK="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
            
            elif [ "$SPECIES" == "mouse" ]; then
                RANKINGS_FEATHER_LINK="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather"

            fi
            curl -L -o "${ORGANISM_DIR}/${CISTARGET_RANKINGS_PRECOMP}" \
                "${RANKINGS_FEATHER_LINK}"
            if [ $? -eq 0 ]; then
                echo "        Done!"
            else
                echo "        Error: Failed to download rankings.feather file."
            fi
        fi

        # File: scores.feather
        if [ -f "${ORGANISM_DIR}/${CISTARGET_SCORES_PRECOMP}" ]; then
            echo "    Precomputed cisTarget ${PYCISTOPIC_SPECIES_CODE} regions_vs_motifs.scores.feather file exists"
        else
            echo "    Downloading scores.feather file..."
            if [ "$SPECIES" == "human" ]; then
                SCORES_FEATHER_LINK="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
            
            elif [ "$SPECIES" == "mouse" ]; then
                SCORES_FEATHER_LINK="https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/screen/mc_v10_clust/region_based/mm10_screen_v10_clust.regions_vs_motifs.scores.feather"
            fi

            curl -L -o "${ORGANISM_DIR}/${CISTARGET_SCORES_PRECOMP}" \
                "${SCORES_FEATHER_LINK}"
            if [ $? -eq 0 ]; then
                echo "        Done!"
            else
                echo "        Error: Failed to download scores.feather file."
            fi
        fi

        echo ""
    fi
}


###############################################################################
# CHECK PATHS AND DEPENDENCIES
###############################################################################

check_or_create_dir "$LOG_DIR"

cd "${SCRIPT_DIR}"

check_if_running
determine_num_cpus
activate_conda_env

install_scenic_plus
install_pycistopic
add_pycistopic_to_path
check_python_deps 

echo ""
echo "[INFO] Checking required directories and files"
# Check required directories
check_dir_exists "$CISTARGET_SCRIPT_DIR"
check_dir_exists "$SCRIPT_DIR"

# Check required files
check_file_exists "$INPUT_DIR/$ATAC_FILE_NAME"
check_file_exists "$INPUT_DIR/$RNA_FILE_NAME"

# Check to see if SCENIC+ generated directories exist or create them
check_or_create_dir "$OUTPUT_DIR"
check_or_create_dir "$TEMP_DIR"
check_or_create_dir "$QC_DIR"
check_or_create_dir "${SCRIPT_DIR}/formatted_inferred_GRNs"
check_or_create_dir "${SCRIPT_DIR}/scplus_pipeline/Snakemake/config"
echo "    - Done!"

# Generate the config file for the cell type and sample
generate_config

echo "    All required files and directories found"
echo ""

check_clusterbuster
check_aertslab_motif_collection
check_organism_genome_files

check_file_exists "$BLACKLIST"
check_file_exists "$GENOME_FASTA"

echo ""
echo "===== CHECKS COMPLETE ====="

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
        --blacklist "${BLACKLIST}" \
        --chromsize_file_path "${CHROMSIZES}";
fi

###############################################################################
# STEP 3: GET TSS DATA
###############################################################################
if [ "$STEP_03_GET_TSS_DATA" = true ]; then
    echo "Step 3: Getting Transcription Start Site data"
    /usr/bin/time -v pycistopic tss get_tss \
        --output "${QC_DIR}/tss.bed" \
        --name "${PYCISTOPIC_SPECIES}" \
        --to-chrom-source ucsc \
        --ucsc "${PYCISTOPIC_SPECIES_CODE}" > "${LOG_DIR}/Step 3: Getting Transcription Start Site data.log" 2>&1;
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

# ###############################################################################
# # STEP 5: CREATE CISTARGET MOTIF DATABASES
# ###############################################################################
# if [ "$STEP_05_CREATE_CISTARGET_MOTIF_DATABASES" = true ]; then
#     CBDIR="${SCRIPT_DIR}/aertslab_motif_colleciton/v10nr_clust_public/singletons"
#     MOTIF_LIST="${INPUT_DIR}/motifs.txt"

#     run_python_step "Step 5: Create cistarget databases" \
#         "${CISTARGET_SCRIPT_DIR}/create_cistarget_motif_databases.py" \
#         -f "${FASTA_FILE}" \
#         -M "${CBDIR}" \
#         -m "${MOTIF_LIST}" \
#         --min 5 \
#         --max 1000 \
#         -o "${OUTPUT_DIR}" \
#         --bgpadding 1000 \
#         -t "${NUM_CPU}"
# fi

###############################################################################
# STEP 6: RUN SNAKEMAKE PIPELINE
###############################################################################
if [ "$STEP_06_RUN_SNAKEMAKE_PIPELINE" = true ]; then
    echo "Step 6: Run SCENIC+ snakemake"

    # Define the Snakefile and target config path
    SNAKEFILE="${SCRIPT_DIR}/scplus_pipeline/Snakemake/workflow/Snakefile"
    NEW_CONFIG_PATH="config/${CELL_TYPE}_${SAMPLE_NAME}_config.yaml"

    echo "    Running snakemake"

    cd "${SCRIPT_DIR}/scplus_pipeline/Snakemake"
    /usr/bin/time -v snakemake \
        --nolock \
        --cores ${NUM_CPU} \
        --snakefile $SNAKEFILE \
        --latency-wait 600 \
        --configfile $NEW_CONFIG_PATH \
        > "${LOG_DIR}/Step 6: Snakemake.log" 2>&1;

    echo "Done!"
    echo ""
fi

if [ "$STEP_07_FORMAT_INFERRED_GRN" = true ]; then
    echo "Step 7: Format SCENIC+ inferred GRN (scplusmdata.h5mu)"
    if [ -f "$OUTPUT_DIR/scplusmdata.h5mu" ]; then
        echo "    File Found! Formatting..."
        run_python_step "Step 7: Format Inferred GRN" \
            "${SCRIPT_DIR}/Step07.format_inferred_grn.py" \
                --output_dir "${SCRIPT_DIR}/formatted_inferred_GRNs" \
                --inferred_grn_file "${OUTPUT_DIR}/scplusmdata.h5mu" \
                --cell_type "$CELL_TYPE" \
                --sample_name "$SAMPLE_NAME"
        echo "    DONE! Formatted GRN saved as 'scenic_plus_inferred_grn_${CELL_TYPE}.tsv'"
    else
        echo "    ERROR! formatting inferred GRN 'scplusmdata.h5mu': File not found in the output directory"
    fi
fi