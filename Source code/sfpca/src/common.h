#ifndef _COMMON_H_
#define _COMMON_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double GAMMA = 3.7;
arma::mat UCOEF;
arma::vec UCONST_VEC;
arma::mat VCOEF;
arma::vec VCONST_VEC;
using namespace std;
using namespace Rcpp;
using namespace arma;

#endif 
