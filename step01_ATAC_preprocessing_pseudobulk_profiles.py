import pycisTopic
import os
import pandas as pd

import shared_variables

cell_data = pd.read_table(shared_variables.cell_data, header=[0])
print(cell_data.head())

out_dir = shared_variables.out_dir

fragments_dict = {
    "E7.5_rep1": shared_variables.fragments_dict
}

chromsizes = pd.read_table(
    shared_variables.mm10_chrom_size_link,
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
print(chromsizes.head())

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data,
    variable = "barcode",
    sample_id_col = "sample",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 10,
    normalize_bigwig = True,
    temp_dir = shared_variables.temp_dir,
    split_pattern = "-"
)

with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")

with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")

