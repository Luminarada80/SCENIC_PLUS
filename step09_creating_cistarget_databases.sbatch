#!/bin/bash
#SBATCH --job-name=step_09_preparing_fasta_from_consensus_regions
#SBATCH --output=step_09_preprocessing_output.txt
#SBATCH --error=step_09_preprocessing_error.txt
#SBATCH --mem-per-cpu=20G
#SBATCH --cpus-per-task=64
#SBATCH -p memory
#SBATCH -A lab_uzun
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --mail-user=habdelhalim@pennstatehealth.psu.edu
#SBATCH --mail-type=ALL
#srun source /gpfs/Home/haa5704/miniconda3/envs/scenicplus

source /gpfs/Home/haa5704/.bashrc
conda init
conda activate scenicplus

OUT_DIR=""${PWD}""
CBDIR="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/create_cisTarget_databases/aertslab_motif_colleciton/v10nr_clust_public/singletons"
FASTA_FILE="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/hg38.10x_brain.with_1kb_bg_padding.fa"
MOTIF_LIST="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/motifs.txt"
SCRIPT_DIR="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/create_cisTarget_databases/"

"${SCRIPT_DIR}/create_cistarget_motif_databases.py" \
    -f ${FASTA_FILE} \
    -M ${CBDIR} \
    -m ${MOTIF_LIST} \
    -o ${OUT_DIR}/${DATABASE_PREFIX} \
    --bgpadding 1000 \
    -t 20
