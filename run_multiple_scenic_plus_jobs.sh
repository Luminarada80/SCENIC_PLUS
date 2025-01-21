#!/bin/bash -l

submit_run_scenic_plus_job() {
    local SAMPLE_NAME=$1
    local CELL_TYPE=$2
    local SPECIES=$3
    local RNA_FILE_NAME=$4
    local ATAC_FILE_NAME=$5

    # Ensure the log directory exists
    mkdir -p "LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs"

    # Submit the job
    sbatch \
        --export=ALL,SAMPLE_NAME="$SAMPLE_NAME",CELL_TYPE="$CELL_TYPE",SPECIES="$SPECIES",RNA_FILE_NAME="$RNA_FILE_NAME",ATAC_FILE_NAME="$ATAC_FILE_NAME" \
        --output="LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs/scenic_plus_${CELL_TYPE}_${SAMPLE_NAME}.out" \
        --error="LOGS/${CELL_TYPE}_logs/${SAMPLE_NAME}_logs/scenic_plus_${CELL_TYPE}_${SAMPLE_NAME}.err" \
        --job-name="SCENIC+_${CELL_TYPE}_${SAMPLE_NAME}" \
        /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/run_scenic_plus.sh
}

run_macrophage() {
    local CELL_TYPE="macrophage"
    local SAMPLE_NAMES=("buffer1")  # Add more samples here
    local SPECIES="human"

    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        local RNA_FILE_NAME="macrophage_${SAMPLE_NAME}_filtered_RNA.csv"
        local ATAC_FILE_NAME="macrophage_${SAMPLE_NAME}_filtered_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"
    done
}

run_mESC(){
    local CELL_TYPE="mESC"
    local SAMPLE_NAMES=(
        "1000_cells_E7.5_rep1"
        )
    local SPECIES="mouse"

    # Submit each SAMPLE_NAME as a separate job
    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        local RNA_FILE_NAME="mESC_${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="mESC_${SAMPLE_NAME}_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"
    done
}

run_macrophage
run_mESC