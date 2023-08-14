#' Seurat CLR 算法在R中与在Rcpp的高效实现
#' @return 返回标准化数据
#' @param mat 稀疏矩阵
#' @export
clr_norm <- function(counts,lang = c("R","Rcpp")) {
  if (!inherits(counts, 'dgCMatrix')){
    stop('`clr_norm` only supports for dgCMatrix.')
  }
  lang <- match.arg(lang)
  if (lang == "R") {
    log_mat <- counts
    log_mat@x <- log1p(counts@x)
    row_exp <- as.numeric(exp(rowMeans(log_mat)))
    counts <- counts/row_exp
    counts@x <- log1p(counts@x)
    counts
  } else {
    rn <- rownames(counts)
    cn <- colnames(counts)
    counts <- clr_norm_rcpp(counts)
    rownames(counts) <- rn
    colnames(counts) <- cn
  }
  counts
}

#' Seurat 选择高可变基因算法
#' vst 使用的是原始数据，不是标准化的数据
#' 本质上就是选取平均值高且
vst <- function(counts,
                clip.max = 'auto',
                loess.span = 0.3,
                nfeatures = 2000) {
  if (clip.max == 'auto') {
    # 在进行VST标准化后,将标准化后值大于clip.max的部分截断为clip.max
    # 这样可以减轻高表达基因过大值的影响
    # 默认的细胞数开方作为阈值通常可以适用大多数情况。
    clip.max <- sqrt(x = ncol(counts))
  }
  #
  hvf.info <- data.frame(mean = rowMeans(x = counts))
  # 计算基因方差,sparseMatrixStats更快
  hvf.info$variance <- sparseMatrixStats::rowVars(counts)
  # 这2列初始化为0
  hvf.info$variance.expected <- 0
  hvf.info$variance.standardized <- 0

  not.const <- hvf.info$variance > 0

  # 1.拟合平滑曲线模型
  fit <- loess(
    formula = log10(x = variance) ~ log10(x = mean),
    data = hvf.info[not.const, ],
    span = loess.span
  )
  # 2.使用模型计算的值
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted

  hvf.info$variance.standardized <- SparseRowVarStd_rcpp(
    mat = counts,
    mu = hvf.info$mean,
    sd = sqrt(hvf.info$variance.expected),
    vmax = clip.max
  )
  rownames(hvf.info) <- rownames(counts)
  hvf.info <- hvf.info[order(hvf.info$vst.variance.standardized, decreasing = TRUE), , drop = FALSE]
  head(rownames(hvf.info), n = nfeatures)
}
