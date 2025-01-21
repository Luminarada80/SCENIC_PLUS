import os
import mudata
import pandas as pd
import argparse
import logging

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Path to the output directory"
    )

    parser.add_argument(
        "--inferred_grn_file",
        type=str,
        required=True,
        help="Path to the scplusmdata.h5mu inferred GRN file from the SCENIC+ snakemake pipeline"
    )
    parser.add_argument(
        "--cell_type",
        type=str,
        required=True,
        help="Cell type analyzed, used for naming the output file"
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        required=True,
        help="Name of the sample being processed"
    )

    args = parser.parse_args()
    
    return args

def main():
    args = parse_args()
    output_dir = args.output_dir
    inferred_grn_file = args.inferred_grn_file
    cell_type = args.cell_type
    sample_name = args.sample_name

    inferred_grn = mudata.read(f'{output_dir}/{inferred_grn_file}')

    inferred_grn_data = inferred_grn.uns["direct_e_regulon_metadata"]

    inferred_grn_data["Source"] = inferred_grn_data["TF"]
    inferred_grn_data["Target"] = inferred_grn_data["Gene"]
    inferred_grn_data["Score"] = inferred_grn_data["importance_x_abs_rho"]

    subset_inferred_grn = pd.DataFrame(inferred_grn_data[["Source", "Target", "Score"]])

    output_file_name = f"{output_dir}/scenic_plus_inferred_grn_{cell_type}_{sample_name}.tsv"

    subset_inferred_grn.to_csv(output_file_name, sep="\t", header=True, index=False)

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    main()