import argparse
import os
from ruamel.yaml import YAML
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

# Initialize ruamel.yaml
yaml = YAML()
yaml.default_flow_style = False
yaml.preserve_quotes = True

# Define argparse arguments
parser = argparse.ArgumentParser(description="Generate a config.yaml file for SCENIC+ pipeline.")

parser.add_argument("--cisTopic_obj_fname", required=True, help="Path to cisTopic object.")
parser.add_argument("--GEX_anndata_fname", required=True, help="Path to gene expression AnnData file.")
parser.add_argument("--region_set_folder", required=True, help="Folder containing region set files.")
parser.add_argument("--ctx_db_fname", required=True, help="Path to motif ranking database.")
parser.add_argument("--dem_db_fname", required=True, help="Path to motif score database.")
parser.add_argument("--path_to_motif_annotations", required=True, help="Path to motif annotations.")
parser.add_argument("--combined_GEX_ACC_mudata", required=True, help="Path to combined GEX-ACC Mudata file.")
parser.add_argument("--dem_result_fname", required=True, help="Path to DEM results file.")
parser.add_argument("--ctx_result_fname", required=True, help="Path to CisTarget results file.")
parser.add_argument("--output_fname_dem_html", required=True, help="Path to DEM results HTML file.")
parser.add_argument("--output_fname_ctx_html", required=True, help="Path to CisTarget results HTML file.")
parser.add_argument("--cistromes_direct", required=True, help="Path to cistromes (direct) H5AD file.")
parser.add_argument("--cistromes_extended", required=True, help="Path to cistromes (extended) H5AD file.")
parser.add_argument("--tf_names", required=True, help="Path to TF names list.")
parser.add_argument("--genome_annotation", required=True, help="Path to genome annotation file.")
parser.add_argument("--chromsizes", required=True, help="Path to chromsizes file.")
parser.add_argument("--search_space", required=True, help="Path to search space file.")
parser.add_argument("--tf_to_gene_adjacencies", required=True, help="Path to TF-to-gene adjacencies file.")
parser.add_argument("--region_to_gene_adjacencies", required=True, help="Path to region-to-gene adjacencies file.")
parser.add_argument("--eRegulons_direct", required=True, help="Path to eRegulons (direct) file.")
parser.add_argument("--eRegulons_extended", required=True, help="Path to eRegulons (extended) file.")
parser.add_argument("--AUCell_direct", required=True, help="Path to AUCell (direct) file.")
parser.add_argument("--AUCell_extended", required=True, help="Path to AUCell (extended) file.")
parser.add_argument("--scplus_mdata", required=True, help="Path to SCENIC+ Mudata file.")
parser.add_argument("--ensembl_species", required=True, help="Species to use for downloading genome annotations from ensembl.")
parser.add_argument("--motif_enrichment_species", required=True, help="Species to use for 'aertslab_motif_colleciton'.")
parser.add_argument("--annotation_version", required=True, help="Annotation version for the 'aertslab_motif_colleciton'.")

# General parameters
parser.add_argument("--temp_dir", required=True, help="Path to temporary directory.")
parser.add_argument("--n_cpu", type=int, required=True, help="Number of CPUs.")
parser.add_argument("--seed", type=int, default=666, help="Random seed (default: 666).")

# Path to save config file
parser.add_argument("--output_config_path", required=True, help="Path to save the generated config.yaml file.")

# Parse arguments
args = parser.parse_args()

# Create the configuration dictionary
new_config = {
    "input_data": {
        "cisTopic_obj_fname": DoubleQuotedScalarString(args.cisTopic_obj_fname),
        "GEX_anndata_fname": DoubleQuotedScalarString(args.GEX_anndata_fname),
        "region_set_folder": DoubleQuotedScalarString(args.region_set_folder),
        "ctx_db_fname": DoubleQuotedScalarString(args.ctx_db_fname),
        "dem_db_fname": DoubleQuotedScalarString(args.dem_db_fname),
        "path_to_motif_annotations": DoubleQuotedScalarString(args.path_to_motif_annotations),
    },
    "output_data": {
        "combined_GEX_ACC_mudata": DoubleQuotedScalarString(args.combined_GEX_ACC_mudata),
        "dem_result_fname": DoubleQuotedScalarString(args.dem_result_fname),
        "ctx_result_fname": DoubleQuotedScalarString(args.ctx_result_fname),
        "output_fname_dem_html": DoubleQuotedScalarString(args.output_fname_dem_html),
        "output_fname_ctx_html": DoubleQuotedScalarString(args.output_fname_ctx_html),
        "cistromes_direct": DoubleQuotedScalarString(args.cistromes_direct),
        "cistromes_extended": DoubleQuotedScalarString(args.cistromes_extended),
        "tf_names": DoubleQuotedScalarString(args.tf_names),
        "genome_annotation": DoubleQuotedScalarString(args.genome_annotation),
        "chromsizes": DoubleQuotedScalarString(args.chromsizes),
        "search_space": DoubleQuotedScalarString(args.search_space),
        "tf_to_gene_adjacencies": DoubleQuotedScalarString(args.tf_to_gene_adjacencies),
        "region_to_gene_adjacencies": DoubleQuotedScalarString(args.region_to_gene_adjacencies),
        "eRegulons_direct": DoubleQuotedScalarString(args.eRegulons_direct),
        "eRegulons_extended": DoubleQuotedScalarString(args.eRegulons_extended),
        "AUCell_direct": DoubleQuotedScalarString(args.AUCell_direct),
        "AUCell_extended": DoubleQuotedScalarString(args.AUCell_extended),
        "scplus_mdata": DoubleQuotedScalarString(args.scplus_mdata),
    },
    "params_general": {
        "temp_dir": DoubleQuotedScalarString(args.temp_dir),
        "n_cpu": args.n_cpu,  # No need for quotes around integers
        "seed": args.seed,    # No need for quotes around integers
    },
    "params_data_preparation": {
        "bc_transform_func": DoubleQuotedScalarString("\"lambda x: f'{x}'\""),
        "is_multiome": True,  # No need for quotes around booleans
        "key_to_group_by": DoubleQuotedScalarString(""),
        "nr_cells_per_metacells": 10,
        "direct_annotation": DoubleQuotedScalarString("Direct_annot"),
        "extended_annotation": DoubleQuotedScalarString("Orthology_annot"),
        "species": DoubleQuotedScalarString(args.ensembl_species),
        "biomart_host": DoubleQuotedScalarString("http://www.ensembl.org"),
        "search_space_upstream": DoubleQuotedScalarString("1000 5000"),
        "search_space_downstream": DoubleQuotedScalarString("1000 5000"),
        "search_space_extend_tss": DoubleQuotedScalarString("10 10"),
    },
    "params_motif_enrichment": {
        "species": DoubleQuotedScalarString(args.motif_enrichment_species),
        "annotation_version": DoubleQuotedScalarString(args.annotation_version),
        "motif_similarity_fdr": 0.001,
        "orthologous_identity_threshold": 0.0,
        "annotations_to_use": DoubleQuotedScalarString("Direct_annot Orthology_annot"),
        # DEM parameters
        "fraction_overlap_w_dem_database": 0.4,
        "dem_max_bg_regions": 500,
        "dem_balance_number_of_promoters": True,
        "dem_promoter_space": 1000,
        "dem_adj_pval_thr": 0.50,
        "dem_log2fc_thr": 0.5,
        "dem_mean_fg_thr": 0.0,
        "dem_motif_hit_thr": 1.00,
        # CisTarget parameters
        "fraction_overlap_w_ctx_database": 0.4,
        "ctx_auc_threshold": 0.01,
        "ctx_nes_threshold": 0.01,
        "ctx_rank_threshold": 0.02,
    },
    "params_inference": {
        "tf_to_gene_importance_method": DoubleQuotedScalarString("GBM"),
        "region_to_gene_importance_method": DoubleQuotedScalarString("GBM"),
        "region_to_gene_correlation_method": DoubleQuotedScalarString("SR"),
        "order_regions_to_genes_by": DoubleQuotedScalarString("importance"),
        "order_TFs_to_genes_by": DoubleQuotedScalarString("importance"),
        "gsea_n_perm": 1000,
        "quantile_thresholds_region_to_gene": DoubleQuotedScalarString("0.85 0.90 0.95"),
        "top_n_regionTogenes_per_gene": DoubleQuotedScalarString("5 10 15"),
        "top_n_regionTogenes_per_region": DoubleQuotedScalarString(""),
        "min_regions_per_gene": 0,
        "rho_threshold": 0.05,
        "min_target_genes": 5,
    },
}

# # Load the existing config file if it exists
# if os.path.exists(args.output_config_path):
#     with open(args.output_config_path, 'r') as config_file:
#         existing_config = yaml.load(config_file)  # Use `yaml.load` for `ruamel.yaml`
# else:
existing_config = {}

# Merge new_config into existing_config
for key, value in new_config.items():
    if key in existing_config:
        if isinstance(existing_config[key], dict) and isinstance(value, dict):
            # Merge dictionaries
            existing_config[key].update(value)
        else:
            # Replace non-dictionary values
            existing_config[key] = value
    else:
        existing_config[key] = value

# Write the updated configuration back to the file
with open(args.output_config_path, 'w') as config_file:
    yaml.dump(existing_config, config_file)

print(f"Config file written to {args.output_config_path}")
