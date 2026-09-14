import sys
import os
import scanpy as sc
import torch
import numpy as np
import cell2location as c2l
import matplotlib.pyplot as plt

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT   = os.path.dirname(SCRIPT_DIR)
WORKING_DIR = SCRIPT_DIR
DATA_DIR    = os.path.join(SCRIPT_DIR, "data")
ADATA       = os.path.join(REPO_ROOT, "quality_control", "primary-cohort", "adata.h5ad")
OUT_DIR     = os.path.join(WORKING_DIR, "spot-deconvolution_zheng")

os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR

results_folder = OUT_DIR
ref_run_name   = f"{results_folder}/reference_signatures"
run_name       = f"{results_folder}/cell2location_map"

# ── Load data ──────────────────────────────────────────────────────────────────
adata_sc = sc.read_h5ad(f"{DATA_DIR}/OvC_adata_zheng_preprocessed.h5ad")
adata_sc.X = adata_sc.X.tocsr().astype("float32")

adata_st = sc.read_h5ad(ADATA)

# ── Filter to shared genes ─────────────────────────────────────────────────────
shared_features = [f for f in adata_st.var_names if f in adata_sc.var_names]

adata_sc = adata_sc[:, adata_sc.var.index.isin(shared_features)].copy()
adata_st = adata_st[:, shared_features].copy()
adata_sc = adata_sc[:, shared_features].copy()

print(f"Number of shared genes between spatial data and reference sc: {adata_st.shape[1]}")

# ── Gene filtering ─────────────────────────────────────────────────────────────
selected = c2l.utils.filtering.filter_genes(
    adata_sc, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
)

adata_sc = adata_sc[:, selected].copy()
adata_st = adata_st[:, selected].copy()

print(adata_sc.obs["cell_type"].value_counts())
print(f"The spatial AnnData object has {adata_st.shape[0]} obs and {adata_st.shape[1]} features.")
print(f"The single-cell AnnData object has {adata_sc.shape[0]} obs and {adata_sc.shape[1]} features.")

# ── Set up and train the regression model ──────────────────────────────────────
c2l.models.RegressionModel.setup_anndata(
    adata=adata_sc,
    batch_key="Patients",
    labels_key="cell_type",
)

model = c2l.models.RegressionModel(adata_sc)

model.train(max_epochs=300, batch_size=2500, train_size=1, lr=0.002)

fig = plt.figure(figsize=(4,3))
model.plot_history(20)
plt.savefig(f"{ref_run_name}/ELBO_loss.pdf", bbox_inches='tight')

fig = plt.figure(figsize=(4,3))
model.plot_QC()
plt.savefig(f"{ref_run_name}/reconstruction_accuracy.pdf", bbox_inches='tight')

# ── Save model ─────────────────────────────────────────────────────────────────
model.save(ref_run_name, overwrite=True)

# ── Export posterior & save results ────────────────────────────────────────────
adata_sc = model.export_posterior(
    adata_sc, sample_kwargs={"num_samples": 1000, "batch_size": model.adata.n_obs}
)

adata_sc.write(f"{ref_run_name}/adata_sc.h5ad")
adata_st.write(f"{ref_run_name}/adata_st.h5ad")

print("Done.")
