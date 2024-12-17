import pycisTopic
import os
import pandas as pd

import argparse

# Define command-line arguments
parser = argparse.ArgumentParser(description="Process scRNA-seq and scATAC-seq data for pseudo-bulk analysis.")

# Add arguments for file paths and directories
parser.add_argument("--input_dir", required=True, help="Path to the data input directory")
parser.add_argument("--output_dir", required=True, help="path to the output directory")
parser.add_argument("--tmp_dir", required=True,  help="path to the temporary directory")
parser.add_argument("--mm10_blacklist", required=True,  help="path to the the pycisTopic mm10 blacklist bed file")


# Parse arguments
args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir
tmp_dir = args.tmp_dir
mm10_blacklist = args.mm10_blacklist

chromsizes = pd.read_table(
    'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes',
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

chromsizes.to_csv(f'{output_dir}/chromsizes.tsv', index=False, sep='\t')

bw_paths = {}
with open(os.path.join(output_dir, "consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})

bed_paths = {}
with open(os.path.join(output_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"

os.makedirs(os.path.join(output_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(output_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    n_cpu = 16,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    skip_empty_peaks = True,
    _temp_dir = tmp_dir
)

from pycisTopic.iterative_peak_calling import get_consensus_peaks
# Other param
peak_half_width=250

# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = mm10_blacklist)

consensus_peaks.to_bed(
    path = os.path.join(output_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)

