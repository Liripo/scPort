#include <RcppArmadillo.h>
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
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

