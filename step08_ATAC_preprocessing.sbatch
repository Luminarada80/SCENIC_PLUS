#!/bin/bash
#SBATCH --job-name=step8_scenicplus
#SBATCH --output=step_08_preprocessing_output.txt
#SBATCH --error=step_08_preprocessing_error.txt
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH -c 16
#SBATCH --mem-per-cpu=16G
#SBATCH -A lab_uzun
#SBATCH --nodes 1
#srun source /gpfs/Home/esm5360/miniconda3/envs/scenicplus

conda activate scenicplus

module load cluster
module load bedtools/2.31.0

REGION_BED="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/outs/consensus_peak_calling/consensus_regions.bed"
GENOME_FASTA="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/data/hg38.fa"
CHROMSIZES="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/data/hg38.chrom.sizes"
SCRIPT_DIR="/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/create_cisTarget_databases/"

${SCRIPT_DIR}/create_fasta_with_padded_bg_from_bed.sh \
        ${GENOME_FASTA} \
        ${CHROMSIZES} \
        ${REGION_BED} \
        hg38.10x_brain.with_1kb_bg_padding.fa \
        1000 \
        yes
head -n 2 hg38.10x_brain.with_1kb_bg_padding.fa


