import pycisTopic
import os
import pandas as pd

import argparse

# Define command-line arguments
parser = argparse.ArgumentParser(description="Process scRNA-seq and scATAC-seq data for pseudo-bulk analysis.")

# Add arguments for file paths and directories
parser.add_argument("--input_dir", required=True, help="Path to the data input directory")
parser.add_argument("--output_dir", required=True, help="path to the output directory")
parser.add_argument("--atac_file_name", required=True,  help="path to the temporary directory")
parser.add_argument("--mm10_blacklist", required=True,  help="path to the the pycisTopic mm10 blacklist bed file")


# Parse arguments
args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir
atac_file_name = args.atac_file_name
tmp_dir = args.tmp_dir
mm10_blacklist = args.mm10_blacklist

cell_data = pd.read_table(f"{input_dir}/{atac_file_name}", header=[0])
print(cell_data.head())

chromsizes = pd.read_table(
    'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes',
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

chromsizes.to_csv(f'{output_dir}/chromsizes.tsv', index=False, sep='\t')
print(chromsizes.head())

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(output_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "barcode",
    sample_id_col = "sample",
    chromsizes = chromsizes,
    bed_path = os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(output_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 10,
    normalize_bigwig = True,
    temp_dir = tmp_dir,
    split_pattern = "-"
)



with open(os.path.join(output_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")

with open(os.path.join(output_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")


