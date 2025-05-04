import os
import sys
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models
import anndata as ad
import scipy.io
from pathlib import Path

if len(sys.argv) != 3:
    print("Usage: seurat_annotate_single.py <Path_to_humansolo_Solo.out> <Path_to_Output_CSV>")
    sys.exit(1)

# Input paths
sample_dir = sys.argv[1]
output_csv = sys.argv[2]
filtered_path = os.path.join(sample_dir, "Gene", "filtered")

# Check required files
required_files = ["matrix.mtx", "features.tsv", "barcodes.tsv"]
missing = [f for f in required_files if not os.path.exists(os.path.join(filtered_path, f))]
if missing:
    print(f"Missing files in {filtered_path}: {', '.join(missing)}")
    sys.exit(1)

# === Step 1: Manual load of 10X output ===
print(f"Loading data from {filtered_path}")
mtx = Path(filtered_path) / "matrix.mtx"
features = Path(filtered_path) / "features.tsv"
barcodes = Path(filtered_path) / "barcodes.tsv"

X = scipy.io.mmread(mtx).T.tocsr()  # Transpose to cells x genes
var_names = pd.read_csv(features, header=None, sep="\t")[1].astype(str).values
obs_names = pd.read_csv(barcodes, header=None)[0].astype(str).values

adata = ad.AnnData(X=X)
adata.var_names = var_names
adata.obs_names = obs_names

# === Step 2: QC filtering ===
print("Performing QC filtering")
adata.var['mt'] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]

# === Step 3: Normalize and dimensionality reduction ===
print("Normalizing and computing PCA")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
sc.tl.leiden(adata, resolution=1.0)

# === Step 4: Cell type annotation with CellTypist ===
print("Annotating cell types with CellTypist")
models.download_models(model = 'Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model='Immune_All_Low.pkl')
adata.obs['cell_type'] = predictions.predicted_labels.values

# === Step 5: Save results ===
print(f"Saving annotations to {output_csv}")
cell_info = pd.DataFrame({
    'Barcode': adata.obs_names,
    'CellType': adata.obs['cell_type'].values
})
cell_info.to_csv(output_csv, index=False)
