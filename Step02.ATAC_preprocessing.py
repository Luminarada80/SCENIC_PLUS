import pycisTopic
import pandas as pd
import numpy as np
import os
import pickle
import joblib
import subprocess
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
from pycisTopic.utils import fig2img
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.lda_models import run_cgs_models_mallet
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.lda_models import evaluate_models
from pycisTopic.topic_binarization import binarize_topics
from pycisTopic.utils import region_names_to_coordinates
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

import logging

import argparse

# Define command-line arguments
parser = argparse.ArgumentParser(description="Process scRNA-seq and scATAC-seq data for pseudo-bulk analysis.")

# Add arguments for file paths and directories
parser.add_argument("--input_dir", required=True, help="Path to the data input directory")
parser.add_argument("--output_dir", required=True, help="path to the output directory")
parser.add_argument("--tmp_dir", required=True,  help="path to the temporary directory")
parser.add_argument("--atac_file_name", required=True,  help="path to the ATAC data file")
parser.add_argument("--blacklist", required=True,  help="path to the the pycisTopic organism blacklist bed file")
parser.add_argument("--chromsize_file_path", required=True, help="Path to the chromosome size file for the organism")

# Parse arguments
args = parser.parse_args()

input_dir = args.input_dir
output_dir = args.output_dir
tmp_dir = args.tmp_dir
atac_file_name = args.atac_file_name
blacklist = args.blacklist
chromsize_file_path = args.chromsize_file_path

chromsizes = pd.read_table(
    chromsize_file_path,
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

chromsizes.to_csv(f'{output_dir}/chromsizes.tsv', index=False, sep='\t')

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

# Path to your peaks x barcodes CSV file
atac_csv_path = f"{input_dir}/{atac_file_name}"

logging.info(f'Preprocessing ATACseq data')

# Read in the ATACseq data as a pandas DataFrame
atac_data = pd.read_csv(atac_csv_path, sep=",", header=0, index_col=0)

# Ensure all data is numeric
atac_data = atac_data.apply(pd.to_numeric, errors='coerce')

# Fill missing values with 0
atac_data.fillna(0, inplace=True)

# Extract regions from row names (assuming "chr:start-end" format)
regions = atac_data.index.to_series().str.extract(r"(chr[^:]+):(\d+)-(\d+)")
regions.columns = ["Chrom", "Start", "End"]

# Debug: Check indices and shape
logging.info(f"Initial atac_data shape: {atac_data.shape}")
logging.info(f"Initial regions shape: {regions.shape}")

nonzero_cells = (atac_data.sum(axis=0) > 0)
logging.info(f"Number of non-zero coverage cells: {nonzero_cells.sum()} / {atac_data.shape[1]}")
atac_data = atac_data.loc[:, nonzero_cells]

# Handle missing values in extracted data
regions["Start"] = pd.to_numeric(regions["Start"], errors="coerce")
regions["End"] = pd.to_numeric(regions["End"], errors="coerce")

# Drop rows with NaN in Start or End
valid_indices = regions.dropna(subset=["Start", "End"]).index
# Debug: Log invalid indices
invalid_indices = set(atac_data.index) - set(valid_indices)
logging.info(f"Invalid indices removed: {len(invalid_indices)}")

# Debug: Verify alignment
logging.info(f"Post-alignment atac_data shape: {atac_data.shape}")
logging.info(f"Post-alignment regions shape: {regions.shape}")

# Convert Start and End to integers
regions["Start"] = regions["Start"].astype(int)
regions["End"] = regions["End"].astype(int)

# Calculate the mean signal for each region
regions["Score"] = atac_data.mean(axis=1).values

# Log a sample of the resulting regions DataFrame
logging.info(regions.head())

if not os.path.exists(f'{output_dir}/consensus_peak_calling/'):
    os.makedirs(f'{output_dir}/consensus_peak_calling/')
print(f'Regions (consensus_regions.bed):')
print(regions[["Chrom", "Start", "End", "Score"]].head())

# Save to BED format
regions[["Chrom", "Start", "End", "Score"]].to_csv(f'{output_dir}/consensus_peak_calling/consensus_regions.bed', sep="\t", header=False, index=False)
logging.info(f'Saved consensus_regions.bed to {output_dir}')

logging.info(f"Number of regions before filtering: {regions.shape[0]}")

# Create the cisTopic object
logging.info(f'\tCreating cistopic_obj from csv file')
cistopic_obj = create_cistopic_object(
    fragment_matrix=atac_data,
    path_to_blacklist=blacklist,
)

logging.info(f"Created cistopic_obj with n_cells × n_regions = {cistopic_obj.fragment_matrix.shape}")

# Inspect the cisTopic object
logging.info(f'\t{cistopic_obj}')

logging.info(f'\tSaving cistopic_obj.pkl')
pickle.dump(
    cistopic_obj,
    open(os.path.join(output_dir, "cistopic_obj.pkl"), "wb")
)

logging.info(f'\tStarting MALLET')
os.environ['MALLET_MEMORY'] = '200G'

# Configure path Mallet
mallet_path="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/scenicplus/Mallet-202108/bin/mallet"

# n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=64,
    n_iter=150,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=tmp_dir,
    save_path=tmp_dir,
    mallet_path=mallet_path,
)


logging.info(f'\t\tDone! Saving MALLET models as "models.pkl"')
pickle.dump(
    models,
    open(os.path.join(output_dir, "models.pkl"), "wb")
)

print(type(models))
print(len(models))

logging.info(f'\tEvaluating models...')
model = evaluate_models(
    models,
    select_model = 40,
    return_model = True
)

logging.info(f'\tAdding selected model to the cistopic_obj...')
cistopic_obj.add_LDA_model(model)

logging.info(f'\nSELECTED MODEL = {cistopic_obj.selected_model}')

pickle.dump(
    cistopic_obj,
    open(os.path.join(output_dir, "cistopic_obj.pkl"), "wb")
)

logging.info(f'\tBinarizing topics using "ntop" for top 3k topics')
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)    
plt.savefig(os.path.join(output_dir, 'Figure_4_topic_region_distribution_ntop.png'))
plt.close()

logging.info(f'\tBinarizing topics using otsu')
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)    
plt.savefig(os.path.join(output_dir, 'Figure_5_topic_region_distribution_otsu.png'))
plt.close()

logging.info(f'\tBinarizing the topics using the "li" method')
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100
)    
plt.savefig(os.path.join(output_dir, 'Figure_6_topic_region_distrbution_li.png'))
plt.close()

logging.info(f'\tComputing topic metrics on the cistopic_obj')
topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.savefig(os.path.join(output_dir, 'Figure_7_topic_qc.png'))
plt.close()

logging.info(f'\tImputing accessibility of the cistopic object')
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

logging.info(f'\tNormalizing imputed accessibility scores')
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

logging.info(f'\tFinding highly variable feature regions')
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)
plt.savefig(os.path.join(output_dir, 'Figure_7_DAR.png'))
plt.close()

logging.info(f'\tCreating "region_sets" directories in the output_dir')
os.makedirs(os.path.join(output_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(output_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(output_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(output_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

logging.info(f'\tCreating DataFrame with region IDs to coordinates mapping for the otsu binarized topics')
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(output_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
    
logging.info(f'\tCreating DataFrame with region IDs to coordinates mapping for the binarized top 3k topics')
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(output_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

