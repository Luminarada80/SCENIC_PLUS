import os
import mudata
import pandas as pd

inferred_grn_file = "/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/outs/scplusmdata.h5mu"

inferred_grn = mudata.read(inferred_grn_file)

inferred_grn_data = inferred_grn.uns["direct_e_regulon_metadata"]

inferred_grn_data["Source"] = inferred_grn_data["TF"]
inferred_grn_data["Target"] = inferred_grn_data["Gene"]
inferred_grn_data["Score"] = inferred_grn_data["importance_x_abs_rho"]

subset_inferred_grn = pd.DataFrame(inferred_grn_data[["Source", "Target", "Score"]])

print(subset_inferred_grn.head())

output_file_name = "/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/outs/scenic_plus_inferred_grn_K562.tsv"

subset_inferred_grn.to_csv(output_file_name, sep="\t", header=True, index=False)