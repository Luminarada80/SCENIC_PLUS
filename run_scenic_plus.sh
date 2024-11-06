#!/bin/bash -l

#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem-per-cpu=16G
#SBATCH -o scenic_plus.log
#SBATCH -e scenic_plus.err

conda activate scenicplus
cd /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS

echo $(which python)

# Add the local pycisTopic to the python path so it's recognized as a module
export PYTHONPATH=$PYTHONPATH:/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/pycisTopic/src

# echo Running step00_RNA_preprocessing.py
# python3 /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/step00_RNA_preprocessing.py

echo Running step01_preprocessing_pseudobulk_profiles.py
python3 /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/step01_ATAC_preprocessing_pseudobulk_profiles.py

echo Running step02_ATAC_inferring_consensus_peaks.py
python3 /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/step02_ATAC_inferring_consensus_peaks.py

echo Running step03_ATAC_qc.py
python3 /gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/step03_ATAC_qc.py
