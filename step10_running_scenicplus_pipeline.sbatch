#!/bin/bash
#SBATCH --job-name=step_10_runing_scenicplus_pipeline
#SBATCH --output=step_10_runing_scenicplus_pipeline_output.txt
#SBATCH --error=step_10_runing_scenicplus_pipeline_error.txt
#SBATCH --mem-per-cpu=1000
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

snakemake --cores 20 --latency-wait 600
