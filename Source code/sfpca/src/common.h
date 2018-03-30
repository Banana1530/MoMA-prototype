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
const bool DEBUG = 0;

int myassert(bool flag, std::string info)
{
    if (flag == 0)
    {
        throw std::invalid_argument(info);
    }
    return 0;
};
double mat_norm(arma::vec u, arma::mat S_u)
{
    return sqrt(as_scalar(u.t() * S_u * u));
}
#endif 
