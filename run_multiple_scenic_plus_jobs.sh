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
    local SAMPLE_NAMES=(
        "macrophage_buffer1_filtered"
        # "macrophage_buffer2_filtered"
        # "macrophage_buffer3_filtered"
        # "macrophage_buffer4_filtered"
        # "macrophage_buffer1_stability1"
        # "macrophage_buffer1_stability2"
        # "macrophage_buffer1_stability3"
        # "macrophage_buffer1_stability4"
        # "macrophage_buffer1_stability5"
        # "macrophage_buffer1_stability6"
        # "macrophage_buffer1_stability7"
        # "macrophage_buffer1_stability8"
        # "macrophage_buffer1_stability9"
        # "macrophage_buffer1_stability10"
        # "macrophage_buffer2_stability1"
        # "macrophage_buffer2_stability2"
        # "macrophage_buffer2_stability3"
        # "macrophage_buffer2_stability4"
        # "macrophage_buffer2_stability5"
        # "macrophage_buffer2_stability6"
        # "macrophage_buffer2_stability7"
        # "macrophage_buffer2_stability8"
        # "macrophage_buffer2_stability9"
        # "macrophage_buffer2_stability10"
        # "macrophage_buffer3_stability1"
        # "macrophage_buffer3_stability2"
        # "macrophage_buffer3_stability3"
        # "macrophage_buffer3_stability4"
        # "macrophage_buffer3_stability5"
        # "macrophage_buffer3_stability6"
        # "macrophage_buffer3_stability7"
        # "macrophage_buffer3_stability8"
        # "macrophage_buffer3_stability9"
        # "macrophage_buffer3_stability10"
        # "macrophage_buffer4_stability1"
        # "macrophage_buffer4_stability2"
        # "macrophage_buffer4_stability3"
        # "macrophage_buffer4_stability4"
        # "macrophage_buffer4_stability5"
        # "macrophage_buffer4_stability6"
        # "macrophage_buffer4_stability7"
        # "macrophage_buffer4_stability8"
        # "macrophage_buffer4_stability9"
        # "macrophage_buffer4_stability10"
        )
    local SPECIES="human"

    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        local RNA_FILE_NAME="${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="${SAMPLE_NAME}_ATAC.csv"

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
        "1000_cells_E7.5_rep1" # failed, found no eRegulons
        # "1000_cells_E7.5_rep2"
        # "1000_cells_E7.75_rep1"
        # "1000_cells_E8.0_rep1"
        # "1000_cells_E8.0_rep2"
        # "1000_cells_E8.5_CRISPR_T_KO"
        # "1000_cells_E8.5_CRISPR_T_WT"
        "2000_cells_E7.5_rep1"
        # "2000_cells_E8.0_rep1"
        # "2000_cells_E8.0_rep2"
        # "2000_cells_E8.5_CRISPR_T_KO"
        # "2000_cells_E8.5_CRISPR_T_WT"
        "3000_cells_E7.5_rep1"
        # "3000_cells_E8.0_rep1"
        # "3000_cells_E8.0_rep2"
        # "3000_cells_E8.5_CRISPR_T_KO"
        # "3000_cells_E8.5_CRISPR_T_WT"
        "4000_cells_E7.5_rep1"
        # "4000_cells_E8.0_rep1"
        # "4000_cells_E8.0_rep2"
        # "4000_cells_E8.5_CRISPR_T_KO"
        # "4000_cells_E8.5_CRISPR_T_WT"
        "5000_cells_E7.5_rep1"
        # "5000_cells_E8.5_CRISPR_T_KO"
        # "70_percent_subsampled_1"
        # "70_percent_subsampled_2"
        # "70_percent_subsampled_3"
        # "70_percent_subsampled_4"
        # "70_percent_subsampled_5"
        # "70_percent_subsampled_6"
        # "70_percent_subsampled_7"
        # "70_percent_subsampled_8"
        # "70_percent_subsampled_9"
        # "70_percent_subsampled_10"
        "filtered_L2_E7.5_rep1"
        # "filtered_L2_E7.5_rep2"
        # "filtered_L2_E7.75_rep1"
        # "filtered_L2_E8.0_rep1"
        # "filtered_L2_E8.0_rep2"
        # "filtered_L2_E8.5_CRISPR_T_KO"
        # "filtered_L2_E8.5_CRISPR_T_WT"
        # "filtered_L2_E8.5_rep1"
        # "filtered_L2_E8.5_rep2"
        # "filtered_L2_E8.75_rep1"
        # "filtered_L2_E8.75_rep2"
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

run_K562(){
    local CELL_TYPE="K562"
    local SAMPLE_NAMES=(
        # "K562_human_filtered"
    )
    local SPECIES="human"

    # Submit each SAMPLE_NAME as a separate job
    for SAMPLE_NAME in "${SAMPLE_NAMES[@]}"; do
        local RNA_FILE_NAME="${SAMPLE_NAME}_RNA.csv"
        local ATAC_FILE_NAME="${SAMPLE_NAME}_ATAC.csv"

        # Submit the job for each sample
        submit_run_scenic_plus_job \
            "$SAMPLE_NAME" \
            "$CELL_TYPE" \
            "$SPECIES" \
            "$RNA_FILE_NAME" \
            "$ATAC_FILE_NAME"
    done
}


# run_macrophage
run_mESC