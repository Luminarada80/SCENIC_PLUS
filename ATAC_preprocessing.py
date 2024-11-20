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

import shared_variables

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(message)s')

out_dir = shared_variables.out_dir

# Path to your peaks x barcodes CSV file
atac_csv_path = "/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER/LINGER_MESC_SC_DATA/FULL_MESC_SAMPLES/multiomic_data_1000_cells_E7.5_rep1_ATAC.csv"

logging.info(f'Preprocessing ATACseq data')

# Read in the ATACseq data as a pandas DataFrame
atac_data = pd.read_csv(atac_csv_path, sep=",", header=0, index_col=0)

# Ensure all data is numeric
atac_data = atac_data.apply(pd.to_numeric, errors='coerce')

# Fill missing values with 0
atac_data.fillna(0, inplace=True)

# Extract regions from row names (assuming "chr:start-end" format)
regions = atac_data.index.str.extract(r"(chr[^:]+):(\d+)-(\d+)")
regions.columns = ["Chrom", "Start", "End"]
regions["Start"] = regions["Start"].astype(int)
regions["End"] = regions["End"].astype(int)
print(atac_data.sum(axis=1).head())  # Summation output for rows


# Calculate total signal for each region
regions["Score"] = atac_data.sum(axis=1).values
print(regions.head())  # Check final result

# Save to BED format
regions[["Chrom", "Start", "End", "Score"]].to_csv(f'{out_dir}/consensus_peak_calling/consensus_regions.bed', sep="\t", header=False, index=False)
print(f'Saved consensus_regions.bed to {out_dir}')

# chrom_sizes_path = "/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/mm10.chrom.sizes"
# bedgraph_dir = f"{out_dir}/consensus_peak_calling/"
# bigwig_dir = f"{out_dir}/consensus_peak_calling/"
# os.makedirs(bedgraph_dir, exist_ok=True)
# os.makedirs(bigwig_dir, exist_ok=True)

# # Loop through barcodes (columns of atac_data)
# # for barcode in atac_data.columns:
# # Create a BEDGraph DataFrame for the barcode
# barcode_bedgraph = regions.copy()
# # barcode_bedgraph["Score"] = regions["Score"] = atac_data.sum(axis=1).values

# # Save as a BEDGraph file
# bedgraph_path = os.path.join(bedgraph_dir, f"consensus_regions.bedgraph")
# barcode_bedgraph[["Chrom", "Start", "End", "Score"]].to_csv(
#     bedgraph_path, sep="\t", header=False, index=False
# )

# # Load the bedGraph file
# try:
#     bedgraph = pd.read_csv(bedgraph_path, sep="\t", header=None, names=["Chrom", "Start", "End", "Score"])
#     print(bedgraph.head())
# except Exception as e:
#     print("Error loading bedGraph:", e)

# # Check for missing values or incorrect rows
# if bedgraph.shape[1] != 4 or bedgraph.isnull().any().any():
#     print("Invalid rows found:")
#     print(bedgraph[bedgraph.isnull().any(axis=1)])
# else:
#     print("bedGraph file is valid.")

# print(f"Saved BEDGraph files for all barcodes in {bedgraph_dir}")

# bedgraph_path = os.path.join(bedgraph_dir, bedgraph_path)
# bigwig_path = os.path.join(bigwig_dir, f"consensus_regions.bw")

# # Sort the BEDGraph file
# sorted_bedgraph_path = os.path.join(bedgraph_dir, f"consensus_region_sorted.bedgraph")
# subprocess.run(
#     ["sort", "-k1,1", "-k2,2n", bedgraph_path, "-o", sorted_bedgraph_path],
#     check=True
# )

# # Convert to BigWig
# subprocess.run(
#     ["/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/bedGraphToBigWig", sorted_bedgraph_path, chrom_sizes_path, bigwig_path],
#     check=True
# )

# # Optionally clean up the sorted BEDGraph file
# os.remove(sorted_bedgraph_path)

# print(f"BigWig files saved in {bigwig_dir}")

# print(f"Saved BEDGraph files for all barcodes in {out_dir}")

# # Save the DataFrame as a BEDGraph file
# bedgraph_path = f'{out_dir}/consensus_peak_calling/consensus_regions.bedGraph'
# bigwig_path = f'{out_dir}/consensus_peak_calling/pseudobulk_bw_files'
# chrom_sizes_path = "/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/mm10.chrom.sizes"  # Path to a file containing chromosome sizes

# # Save regions as BEDGraph
# regions[["Chrom", "Start", "End", "Score"]].to_csv(
#     bedgraph_path, sep="\t", header=False, index=False
# )

# # Sort the BEDGraph file
# sorted_bedgraph_path = "regions_sorted.bedgraph"
# subprocess.run(
#     ["sort", "-k1,1", "-k2,2n", bedgraph_path, "-o", sorted_bedgraph_path],
#     check=True
# )

# # Convert the sorted BEDGraph to BigWig
# subprocess.run(
#     ["bedGraphToBigWig", sorted_bedgraph_path, chrom_sizes_path, bigwig_path],
#     check=True
# )

# # Cleanup temporary files if desired
# os.remove(bedgraph_path)
# os.remove(sorted_bedgraph_path)

# print(f"BigWig file created: {bigwig_path}")

# Create the cisTopic object
logging.info(f'\tCreating cistopic_obj from csv file')
cistopic_obj = create_cistopic_object(
    fragment_matrix=atac_data,
    path_to_blacklist=shared_variables.path_to_mm10_blacklist,
)

# Inspect the cisTopic object
logging.info(f'\t{cistopic_obj}')

logging.info(f'\tSaving cistopic_obj.pkl')
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)

logging.info(f'\tStarting MALLET')
os.environ['MALLET_MEMORY'] = '200G'

# Configure path Mallet
mallet_path="/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/scenicplus/Mallet-202108/bin/mallet"

# n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],

# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5],
    n_cpu=16,
    n_iter=150,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=shared_variables.temp_dir,
    save_path=shared_variables.temp_dir,
    mallet_path=mallet_path,
)


logging.info(f'\t\tDone! Saving MALLET models as "models.pkl"')
pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)

print(type(models))
print(len(models))

model = evaluate_models(
    models,
    select_model = 2,
    return_model = True
)

cistopic_obj.add_LDA_model(model)

logging.info(f'\tEvaluating models')

logging.info(f'\tBinarizing topics using "ntop" for top 3k topics')
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)    
plt.savefig(os.path.join(out_dir, 'Figure_4_topic_region_distribution_ntop.png'))
plt.close()

logging.info(f'\tBinarizing topics using otsu')
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)    
plt.savefig(os.path.join(out_dir, 'Figure_5_topic_region_distribution_otsu.png'))
plt.close()

logging.info(f'\tBinarizing the topics using the "li" method')
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100
)    
plt.savefig(os.path.join(out_dir, 'Figure_6_topic_region_distrbution_li.png'))
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
plt.savefig(os.path.join(out_dir, 'Figure_7_topic_qc.png'))
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
plt.savefig(os.path.join(out_dir, 'Figure_7_DAR.png'))
plt.close()

logging.info(f'\tCreating "region_sets" directories in the out_dir')
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

logging.info(f'\tCreating DataFrame with region IDs to coordinates mapping for the otsu binarized topics')
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
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
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

