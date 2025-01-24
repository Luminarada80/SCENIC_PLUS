#!/bin/bash -l

###############################################################################
# SLURM DIRECTIVES
###############################################################################
#SBATCH --job-name="SCENIC+_${CELL_TYPE}_${SAMPLE_NAME}"
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=32G

###############################################################################
# ENVIRONMENT SETUP
###############################################################################
set -euo pipefail
trap "echo 'An error occurred. Exiting...'; exit 1;" ERR

conda activate scenicplus
SCRIPT_DIR="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS"
cd "${SCRIPT_DIR}"

echo "Python executable: $(which python)"
echo ""

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

# Optional: Use precomputed cisTarget database
USE_PRECOMPUTED_CISTARGET_DB=true
# # Or create your own cisTarget motif database
# STEP_05_CREATE_CISTARGET_MOTIF_DATABASES=false

STEP_06_RUN_SNAKEMAKE_PIPELINE=true
STEP_07_FORMAT_INFERRED_GRN=true

###############################################################################
# PATH SETUP
###############################################################################


###############################################################################
# INPUT FILES & DIRECTORIES
###############################################################################
LOG_DIR="${SCRIPT_DIR}/LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs/"

INPUT_DIR="${SCRIPT_DIR}/input/${CELL_TYPE}/${SAMPLE_NAME}"
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
    ANNOTATION_VERSION="v9-nr"
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
echo ""

###############################################################################
# FUNCTION DEFINITIONS
###############################################################################
check_if_running() {
    echo ""
    echo "Checking for SCENIC+ jobs running for the current sample..."
    # Use the SLURM job name for comparison
    JOB_NAME="${SLURM_JOB_NAME}"  # Dynamically retrieve the job name from SLURM

    # Check for running jobs with the same name, excluding the current job
    RUNNING_COUNT=$(squeue --name="$JOB_NAME" --noheader | wc -l)

    # If other SCENIC+ jobs are running and there are lock files, exit
    if [ "$RUNNING_COUNT" -gt 1 ]; then
        echo "    A job with the name '"$JOB_NAME"' is already running:"
        echo "    Exiting to avoid conflicts."
        exit 1
    else
        echo "    No jobs already running"
    fi
    echo ""
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
    --chromsizes "${CHROMSIZES}" \
    --search_space "${OUTPUT_DIR}/search_space.tsv" \
    --tf_to_gene_adjacencies "${OUTPUT_DIR}/tf_to_gene_adj.tsv" \
    --region_to_gene_adjacencies "${OUTPUT_DIR}/region_to_gene_adj.tsv" \
    --eRegulons_direct "${OUTPUT_DIR}/eRegulons_direct.tsv" \
    --eRegulons_extended "${OUTPUT_DIR}/eRegulons_extended.tsv" \
    --AUCell_direct "${OUTPUT_DIR}/AUCell_direct.h5mu" \
    --AUCell_extended "${OUTPUT_DIR}/AUCell_extended.h5mu" \
    --scplus_mdata "${OUTPUT_DIR}/scplusmdata.h5mu" \
    --temp_dir "${SCRIPT_DIR}/tmp" \
    --n_cpu "${NUM_CPU}" \
    --seed 666 \
    --ensembl_species "${ENSEMBL_SPECIES}" \
    --motif_enrichment_species "${MOTIF_ENRICHMENT_SPECIES}" \
    --annotation_version "${ANNOTATION_VERSION}" \
    --output_config_path "${SCRIPT_DIR}/scplus_pipeline/Snakemake/config/${CELL_TYPE}_${SAMPLE_NAME}_config.yaml"
}


###############################################################################
# CHECK PATHS AND DEPENDENCIES
###############################################################################

check_if_running

echo "Checking for all required directories and files"
# Check required directories
check_dir_exists "$CISTARGET_SCRIPT_DIR"
check_dir_exists "$SCRIPT_DIR"

# Check required files
check_file_exists "$INPUT_DIR/$ATAC_FILE_NAME"
check_file_exists "$INPUT_DIR/$RNA_FILE_NAME"
check_file_exists "$BLACKLIST"
check_file_exists "$GENOME_FASTA"

# Check to see if SCENIC+ generated directories exist or create them
check_or_create_dir "$LOG_DIR"
check_or_create_dir "$OUTPUT_DIR"
check_or_create_dir "$TEMP_DIR"
check_or_create_dir "$QC_DIR"

# Generate the config file for the cell type and sample
generate_config

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

# # Update the n_cpu value in the config.yaml file
# CONFIG_FILE="${SCRIPT_DIR}/scplus_pipeline/Snakemake/config/config.yaml"
# sed -i "s/^\(\s*n_cpu:\s*\).*/\1${NUM_CPU}/" "${CONFIG_FILE}"
# echo "Updated config.yaml with n_cpu: ${NUM_CPU}"
# echo ""

# Download the precomputed cisTarget database
if [ "$USE_PRECOMPUTED_CISTARGET_DB" = true ]; then

    echo "Using pre-computed cisTarget database"

    # Ensure destination directory exists
    mkdir -p "$INPUT_DIR"

    # File: rankings.feather
    if [ -f "${ORGANISM_DIR}/${CISTARGET_RANKINGS_PRECOMP}" ]; then
        echo "    Precomputed cisTarget ${PYCISTOPIC_SPECIES_CODE} regions_vs_motifs.rankings.feather file exists"
    else
        echo "    Downloading rankings.feather file..."
        if [ "$ORGANISM" == "human" ]; then
            RANKINGS_FEATHER_LINK="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
        
        elif [ "$ORGANISM" == "mouse" ]; then
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
        if [ "$ORGANISM" == "human" ]; then
            SCORES_FEATHER_LINK="https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/screen/mc_v10_clust/region_based/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
        
        elif [ "$ORGANISM" == "mouse" ]; then
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
        --blacklist "${BLACKLIST}" \
        --chromsize_file_path "${CHROMSIZES}";
fi

###############################################################################
# STEP 3: GET TSS DATA
###############################################################################
if [ "$STEP_03_GET_TSS_DATA" = true ]; then
    echo "Step 3: Getting Transcription Start Site data"
    pycistopic tss get_tss \
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

    # Update the configfile line in the Snakefile
    sed -i "s|^\s*configfile:.*|configfile: \"${NEW_CONFIG_PATH}\"|" "${SNAKEFILE}"

    # Print confirmation
    echo "Updated Snakefile with configfile: ${NEW_CONFIG_PATH}"
    echo "    Variable value: '$(grep 'configfile' "${SNAKEFILE}")'"
    echo ""

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
            echo "        Exiting to avoid conflicts."
            exit 1
        
        # If no other SCENIC+ jobs are running and there are lock files, remove them
        else
            echo "        No other jobs with the name '"$JOB_NAME"', removing locks"
            rm -rf "$LOCK_FILE/*"
        fi
        
    else
        echo "    No lock files detected."
    fi

    echo "    Running snakemake"

    cd "${SCRIPT_DIR}/scplus_pipeline/Snakemake"
    snakemake --cores ${NUM_CPU} --latency-wait 600 > "${LOG_DIR}/Step 6: Snakemake.log" 2>&1;
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
                --inferred_grn_file "scplusmdata.h5mu" \
                --cell_type "$CELL_TYPE" \
                --sample_name "$SAMPLE_NAMe"
        echo "    DONE! Formatted GRN saved as 'scenic_plus_inferred_grn_${CELL_TYPE}.tsv'"
    else
        echo "    ERROR! formatting inferred GRN 'scplusmdata.h5mu': File not found in the output directory"
    fi
fi