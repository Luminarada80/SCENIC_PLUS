import scanpy as sc 
import pandas as pd

import shared_variables

adata = sc.read_10x_mtx(
    shared_variables.rna_data_dir,
    var_names = "gene_symbols"
)

adata.var_names_make_unique()

cell_data = pd.read_csv(shared_variables.cell_data, index_col=0)

cell_data.index = [cb.rsplit("-", 1)[0] for cb in cell_data.index]

#adata = adata[list(set(adata.obs_names) & set(cell_data.index))].copy()

#adata.obs = cell_data.loc[adata.obs_names]

adata.obs['cell_type'] = 'E7.5_rep1'

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)

adata.write("adata_final.h5ad")
