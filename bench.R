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


bench::mark(
  as.numeric(colSums(counts)),
  sparseMatrixStats::colSums2(counts)
)

sparseMatrixStats::rowSds(counts[1:10,])

Seurat:::SparseRowVarStd(
  mat = counts,
  mu = hvf.info$mean,
  sd = sqrt(hvf.info$variance.expected),
  vmax = sqrt(ncol(counts)),
  display_progress = F
)
