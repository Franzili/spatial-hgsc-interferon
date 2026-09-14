import os
import numpy as np
import scanpy as sc
import pandas as pd
import cell2location as c2l
from cell2location.utils import select_slide
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from aquarel import load_theme
import cmcrameri

theme = (
    load_theme("umbra_light").set_overrides({
    "figure.facecolor": 'white',
    "axes.facecolor": 'white'
})
    .set_grid(draw=False)
)

SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT   = os.path.dirname(SCRIPT_DIR)
WORKING_DIR = SCRIPT_DIR
DATA_DIR    = os.path.join(SCRIPT_DIR, "data")
OUT_DIR     = os.path.join(WORKING_DIR, "spot-deconvolution_zheng")

os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR

results_folder = OUT_DIR
ref_run_name   = f"{results_folder}/reference_signatures"
run_name       = f"{results_folder}/cell2location_map_reg20"

adata_sc = sc.read_h5ad(os.path.join(ref_run_name, 'adata_sc.h5ad'))
adata_st = sc.read_h5ad(os.path.join(ref_run_name, 'adata_st.h5ad'))

mod = c2l.models.RegressionModel.load(f"{ref_run_name}", adata_sc)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_sc.varm.keys():
    inf_aver = adata_sc.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_sc.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_sc.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_sc.uns['mod']['factor_names']

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_st.var_names, inf_aver.index)
adata_vis = adata_st[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
print(f"Shape intersect: {intersect.shape}")

# prepare anndata for cell2location model
c2l.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
model = c2l.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=15,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)

model.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
         )

# export the estimated cell abundance (summary of the posterior distribution).
adata_vis = model.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': model.adata.n_obs}
)

# Save model
model.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/spatial.h5ad"
adata_vis.write(adata_file)

# plot ELBO loss history during training, removing first 100 epochs from the plot
with theme:
    fig = plt.figure(figsize=(4,3))
    model.plot_history(100)
    plt.legend(labels=['full data training']);
    plt.title("ELBO loss full data training")
    plt.savefig(f"{run_name}/ELBO-loss_full-training.pdf", bbox_inches='tight')

    fig = plt.figure(figsize=(4,3))
    model.plot_QC()
    plt.savefig(f"{run_name}/model_qc.pdf", bbox_inches='tight')

    fig = plt.figure(figsize=(4,3))
    model.plot_spatial_QC_across_batches()
    plt.savefig(f"{run_name}/model_qc_across_batches.pdf", bbox_inches='tight')

