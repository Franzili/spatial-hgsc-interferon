library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(sceasy)
library(reticulate)
# Point RETICULATE_CONDA_ENV at a conda env that has anndata installed,
# e.g. Sys.setenv(RETICULATE_CONDA_ENV = "convert-seurat-anndata")
use_condaenv(Sys.getenv("RETICULATE_CONDA_ENV", unset = "convert-seurat-anndata"))

obj <- readRDS("data/raw_object.rds")
obj <- UpdateSlots(object = obj)
obj <- Seurat::UpdateSeuratObject(obj)

obj[["RNA"]] <- CreateAssayObject(counts = obj[["RNA"]]@counts)

Assays(obj)
DefaultAssay(obj) <- "RNA"
head(slotNames(obj))
head(obj@assays$RNA@counts)

sceasy::convertFormat(obj, from="seurat", to="anndata",
                      outFile='data/OvC_adata_zheng.h5ad')
