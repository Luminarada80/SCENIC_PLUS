import scanpy as sc
import pandas as pd
import numpy as np

# Paths to preprocessed RNA and ATAC data
rna_csv_path = "/gpfs/Labs/Uzun/DATA/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/LINGER/LINGER_MESC_SC_DATA/FULL_MESC_SAMPLES/multiomic_data_1000_cells_E7.5_rep1_RNA.csv"

# Read RNA-seq and ATAC-seq data
rna_data = pd.read_csv(rna_csv_path, index_col=0)  # Genes x Cells

# Check dimensions
print(f"RNA-seq shape: {rna_data.shape}")

# Create AnnData objects
adata_rna = sc.AnnData(rna_data.T)  # Transpose to Cells x Features
adata_rna.var_names_make_unique()

# Assign metadata if available
adata_rna.obs['modality'] = 'RNA'

adata_rna.raw = adata_rna.copy()

adata_rna.write("adata_final.h5ad")
