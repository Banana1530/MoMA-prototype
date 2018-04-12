#ifndef _COMMON_H_
#define _COMMON_H_
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// modified in commmon.cpp
extern double GAMMA;
extern arma::mat UCOEF;
extern arma::vec UCONST_VEC;
extern arma::mat VCOEF;
extern arma::vec VCONST_VEC;
extern bool DEBUG;


using namespace std;
using namespace Rcpp;
using namespace arma;


int myassert(bool flag, std::string info);
double mat_norm(arma::vec u, arma::mat S_u);
#endif 
