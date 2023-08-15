#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::sp_mat LogNorm(arma::sp_mat data, int scale_factor = 10000) {

  arma::SpRow<double> col_sums = arma::sum(data, 0);
  for (int i = 0; i < (int)data.n_cols; ++i) {
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
  for (int i = 0; i < (int)mat.n_cols; ++i) {
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
  for (int i = 0; i < (int)mat.n_cols; i++) {
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


// 回归残差的 Rcpp 实现
// [[Rcpp::export]]
NumericMatrix regress_out_matrix_rcpp(NumericMatrix mat, Rcpp::List qr) {
  // Rcpp 调用 R 函数
  Function qr_resid("qr.resid");
  // 行遍历效率不高? 转置之后列遍历似乎差别不大，快了一些
  for (int i = 0; i < mat.rows(); i++) {
    Rcpp::NumericVector y = mat.row(i);
    Rcpp::NumericVector regression_mat = qr_resid(Named("qr") = qr,Named("y") = y);
    mat.row(i) = regression_mat;
  }
  return mat;
}


// [[Rcpp::export(rng = false)]]
arma::mat FastSparseRowScale_rcpp(arma::sp_mat mat, bool scale = true, bool center = true,
                             double scale_max = 10){
  mat = mat.t();
  arma::SpRow<double> colMeans = arma::mean(mat, 0);
  arma::mat scaled_mat(mat.n_rows, mat.n_cols);
  for (int i=0; i < (int)mat.n_cols ; ++i){
    
    double colMean = colMeans(i);
    double colSdev = 0;

    if (scale == true){
      int nnZero = 0;
      if(center == true){
        for (arma::sp_mat::iterator it = mat.begin_col(i); it != mat.end_col(i); ++it ) {
          nnZero += 1;
          colSdev += pow((*it - colMean), 2);
        }
        colSdev += pow(colMean, 2) * (mat.n_rows - nnZero);
      }
      else{
        for (arma::sp_mat::iterator it = mat.begin_col(i); it != mat.end_col(i); ++it) {
          colSdev += pow(*it, 2);
        }
      }
      colSdev = sqrt(colSdev / (mat.n_rows - 1));
    }
    else{
      colSdev = 1;
    }
    if(center == false){
      colMean = 0;
    }
    scaled_mat.col(i) = (mat.col(i) - colMean) / colSdev;
    for(int s=0; s< (int)scaled_mat.n_rows; ++s){
      if (scaled_mat(s,i) > scale_max){
        scaled_mat(s,i) = scale_max;
      }
    }
  }
  return scaled_mat.t();
}