#include "common.h"
#include "model.h"
#include "solver.h"
#include "train.h"


// [[Rcpp::export]]
extern "C" SEXP sfpca(
    std::string model_type,
    arma::mat X ,
    arma::mat Y,
    arma::mat Omega_u,
    arma::mat Omega_v,
    double alpha_u = 0,
    double alpha_v = 0,
    double lambda_u = 0,
    double lambda_v = 0,
    std::string P_u = "l1",
    std::string P_v = "l1",
    double scad_a = 3.7,
    int non_neg = 0, // 1 means activate non-negativity constraint
    double EPS = 1e-4,
    long MAX_ITER = 1e+3,
    std::string solver = "ISTA",
    bool SVD = 1
)
{
  

  int nx = X.n_rows;
  int ny = Y.n_rows;
// Input checking
  myassert(nx==ny, "Sample size should match!\n");
  myassert(scad_a>2,"SCAD should have a>2!\n");
  myassert(solver.compare("ISTA")==0 || solver.compare("FISTA")==0,"Not supported solver!\n");

  if(DEBUG) cout<<"Before buildmod"<<endl;
  Model mod = build_model(X,Y,model_type);
  if (DEBUG)
	  cout << "After buildmod"<<endl;
	if(DEBUG) cout<<"Before build solver"<<endl;
		  Solver sol = build_solver(mod, Omega_u, Omega_v, alpha_u, alpha_v, lambda_u, lambda_v, P_u, P_v, scad_a, non_neg, EPS, MAX_ITER, solver, SVD);
	if(DEBUG) cout<<"After build solver";
  Rcpp::List res = train(mod,sol);
	return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be autoarma::matically
// run after the compilation.
//

/*** R

*/
