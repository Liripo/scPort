library(Matrix)
library(scPort)
library(Seurat)

# pbmc
counts <- Read10X(data.dir = "../../data/scRNA/filtered_gene_bc_matrices/hg19/")

# 似乎R代码更高效
tx <- bench::mark(
  clr_norm(counts,lang = "R"),
  clr_norm(counts,lang = "Rcpp"),
  as(Seurat::NormalizeData(counts,normalization.method = "CLR"),Class = "dgCMatrix"),
  iterations = 1
)



# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = counts, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

head(VariableFeatures(pbmc))

# 基因方差
bench::mark(
  Seurat:::SparseRowVar2(counts,hvf.info$mean,display_progress = F),
  sparseMatrixStats::rowVars(counts)
)


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# scale_data


pbmc <- ScaleData(pbmc)

bench::mark(
  scPort:::regress_out_matrix(
  LayerData(pbmc,layer = "data")[VariableFeatures(pbmc),],
  latent.data = pbmc[["percent.mt"]]
),

Seurat:::RegressOutMatrix(
  LayerData(pbmc,layer = "data")[VariableFeatures(pbmc),],
  latent.data = pbmc[["percent.mt"]],
  model.use = "linear"
),iterations = 1
)


bench::mark(
  {
    data <- Seurat:::FastSparseRowScale(LayerData(pbmc,layer = "data"))
    dimnames(data) <- dimnames(LayerData(pbmc,layer = "data"))
    data
  },
  Seurat:::FastRowScale(as.matrix(LayerData(pbmc,layer = "data"))),
  iterations = 1
)


bench::mark(
  Seurat:::FastSparseRowScale(LayerData(pbmc,layer = "data")),
  FastSparseRowScale_rcpp(LayerData(pbmc,layer = "data")),
  iterations = 1
)
