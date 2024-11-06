data_dir: str = '/gpfs/Labs/Uzun/DATA/PROJECTS/2024.SC_MO_TRN_BENCHMARKING.MIRA/SCENIC_PLUS.HABIBA/data'
results_dir: str = '/gpfs/Labs/Uzun/RESULTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS'

cell_data: str = f'{data_dir}/GSE205117_cell_metadata_filtered.tsv'

fragments_dict: str = f'{data_dir}/GSM6205427_E7.5_rep1_ATAC_fragments.tsv.gz'

rna_data_dir: str = f'{data_dir}/filtered_feature_bc_matrix/'

out_dir: str = f'{results_dir}/scenicplus/mESC_new_scenicplus/outs/'

temp_dir: str = f'{results_dir}/scenicplus/tmp_mESC/'

mm10_chrom_size_link: str = 'https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'

path_to_mm10_blacklist: str = '/gpfs/Labs/Uzun/SCRIPTS/PROJECTS/2024.GRN_BENCHMARKING.MOELLER/SCENIC_PLUS/pycisTopic/blacklist/mm10-blacklist.v2.bed'