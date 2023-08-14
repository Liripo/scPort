#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::sp_mat LogNorm(arma::sp_mat data, int scale_factor = 10000) {

  arma::SpRow<double> col_sums = arma::sum(data, 0);
  for (int i = 0; i < data.n_cols; ++i) {
    double sum = col_sums(i);
    // arma::sp_mat::iterator 是Armadillo库中稀疏矩阵arma::sp_mat的迭代器类型。
    // 允许你以类似迭代器的方式遍历稀疏矩阵的非零元素。

    for (arma::sp_mat::iterator it = data.begin_col(i); it != data.end_col(i); ++it) {
      *it = std::log1p(*it / sum * scale_factor);
    }
  }
  return data;
}

// Seurat CLR 算法在 Armadillo 中的实现
// [[Rcpp::export]]
arma::sp_mat clr_norm_rcpp(arma::sp_mat mat) {

  arma::sp_mat log_mat = mat;
  // 遍历修改稀疏矩阵中的非零元素
  for(arma::sp_mat::iterator it = log_mat.begin(); it != log_mat.end(); ++it) {
    *it = std::log1p(*it);
  }
  arma::mat row_means = arma::mat(arma::mean(log_mat, 1));
  row_means = arma::exp(row_means);
  // 转置后遍历列，行遍历效率很差
  mat = mat.t();
  for (int i = 0; i < mat.n_cols; ++i) {
    double exp_mean = row_means(i);
    for (arma::sp_mat::iterator it = mat.begin_col(i); it != mat.end_col(i); ++it) {
      *it = std::log1p(*it/exp_mean);
    }
  }
  return mat.t();
}


// Seurat vst 标准化
// [[Rcpp::export]]
NumericVector SparseRowVarStd_rcpp(
  arma::sp_mat mat,
  NumericVector mu,
  NumericVector sd,
  double vmax) {
  mat = mat.t();
  NumericVector allVars(mat.n_cols);
  for (int i = 0; i < mat.n_cols; i++) {
    if (sd[i] == 0) continue;
    double colSum = 0;
    // 由于迭代的是非零值，还需要去加上 0值的标准化，可以统一乘以同个数
    int nZero = mat.n_rows;
    for (arma::sp_mat::iterator it = mat.begin_col(i); it != mat.end_col(i); ++it) {
      nZero -= 1;
      colSum += std::pow(std::min(vmax, (*it - mu[i]) / sd[i]), 2);
    }
    colSum += std::pow((0 - mu[i]) / sd[i], 2) * nZero;
    allVars[i] = colSum / (mat.n_rows - 1);
  }
  return allVars;
}
