#ifndef DATA_MANIPULATION
#define DATA_MANIPULATION

#include <RcppArmadillo.h>


arma::sp_mat LogNorm(arma::sp_mat data, int scale_factor = 10000);

#endif//DATA_MANIPULATION
