#include "common.h"
#include "model.h"
#include "solver.h"




typedef struct{
  Prox_op prox_u; // Model, sparse penalty, non_neg
  Prox_op prox_v;
  Grad grad_u; // alpha ==0, S_us?
  Grad grad_v;
  double prox_u_step; //lambda_u, L_u(S_u(Omeg,alpha_a)
  double prox_v_step;
  double grad_u_step; // L_u
  double grad_v_step;
  arma::vec u; // bool SVD, transformed matrix
  arma::vec v;
  arma::mat Su;
  arma::mat Sv;
  SOLVER_TYPE solver_type;
  long MAX_ITER;
  double EPS;
} Solver;

Model build_model(arma::mat X, arma::mat Y, std::string model_type){
  Model model={.model_type = string_to_MT(model_type)
  };
  
  switch(model.model_type){
    case PCA:
      if(DEBUG) cout << "Chosing case PCA: " << model.model_type << endl;

      model.X = X;
      break;
    case LDA:
      throw std::invalid_argument(model_type + " is not currently supported");
      break;
    case PLS:
      model.X = X.t() * Y;
      break;
    case CCA:
      model.X = X.t() * Y;
      break;
    default: 
      throw std::invalid_argument(model_type + " is not currently supported");
  }
  return model;
}

Solver build_solver(
                    Model mod,
                    arma::mat Omega_u = arma::mat(0),
                    arma::mat Omega_v = arma::mat(0),
                    double alpha_u = 0,
                    double alpha_v = 0,
                    double lambda_u = 0,
                    double lambda_v = 0,
                    std::string P_u = "l1", // assume for now they have same type of penalty
                    std::string P_v = "l1",
                    double scad_a = 3.7,
                    int non_neg = 0, // 1 means activate non-negativity constraint
                    double EPS = 1e-4,
                    long MAX_ITER = 1e+3,
                    std::string solver_type_string = "ISTA",
                    bool SVD = 1){

    int n = mod.X.n_rows, p = mod.X.n_cols;
   
    if(DEBUG) cout<<n << p<<endl;
    // Find S_u --> L_u
    arma::mat U;
    arma::vec s;
    arma::mat V;
    if(DEBUG) cout << "here Before  svd(U, s, V, mod.X);" << endl;
    svd(U, s, V, mod.X);
    if(DEBUG) cout << " here After svd(U,s,V,mod.X)" << endl;
    Solver sol;
    if(DEBUG) cout << "here Before sol.Su.eye(size(Omega_u))" <<endl;
    
    sol.Su.eye(size(Omega_u));
    sol.Su += n * alpha_u * Omega_u;
    if(DEBUG) myassert(Omega_v.n_cols == 200, "Wrong Omega_v input");
    sol.Sv.eye(size(Omega_v));
    sol.Sv += p * alpha_v * Omega_v;
    if (DEBUG){
      myassert(sol.Su.n_cols == sol.Su.n_cols && sol.Su.n_cols == 199, "Wrong Su");
      myassert(sol.Sv.n_cols == 200, "Wrong Sv");
    }
      

    double Lu = eig_sym(sol.Su).max() + 0.01;
    double Lv = eig_sym(sol.Sv).max() + 0.01;

    sol.MAX_ITER = MAX_ITER;
    sol.EPS = EPS;
    sol.grad_u_step = 1 / Lu;
    sol.grad_v_step = 1 / Lv;
    sol.prox_u_step = lambda_u / Lu;
    sol.prox_v_step = lambda_v / Lv;
    sol.solver_type = string_to_ST(solver_type_string);
    if (SVD)
    {
      sol.v = V.col(0);
      if(DEBUG) cout<< sol.v.n_elem<<endl;  
      sol.u = U.col(0);
      if(DEBUG) cout<< sol.u.n_elem<<endl;
      if(DEBUG) myassert(sol.u.n_elem == 199, "WRONG DIM OF U");
    }
    else
    {
      sol.v = randu<arma::vec>(p);
      sol.u = randu<arma::vec>(n);
    }

    
    // match gradient
    UCONST_VEC = mod.X * sol.v;
    VCONST_VEC = mod.X.t() * sol.u;
    if(abs(alpha_u) <1e-9 && abs(alpha_v)<1e-9){
      
      sol.grad_u = Grad_nosmoothu;
      sol.grad_v = Grad_nosmoothv;
    }else{
      if(DEBUG) cout<< "CHOOSE Grad_line"<<endl;
      UCOEF = sol.Su;
      VCOEF = sol.Sv;
      myassert(VCOEF.n_rows ==200, "VCOEF sucks");
      sol.grad_u = Grad_lineu;
      sol.grad_v = Grad_linev;
    }

    // match proximal operator
    if (P_u.compare("l1") == 0 || P_u.compare("L1") == 0)
    {

      if (non_neg == 1)
      {
        sol.prox_u = Lasso_pos;
        sol.prox_v = Lasso_pos;
      }
      else{
        if(DEBUG) cout<< "CHOOSING Lasso"<<endl;
        sol.prox_u = Lasso;
        sol.prox_v = Lasso;
      }
    }
    else if (P_u.compare("scad") == 0 || P_u.compare("SCAD") == 0)
    {
      GAMMA = scad_a;
      if(non_neg==0)
      {
        
        sol.prox_u = Scad;
        sol.prox_v = Scad;
      }
      else
        throw std::invalid_argument("NONNEG_SCAD Not currently supported");
    }
    else
      throw std::invalid_argument("Not currently supported");

    return sol;
  }     

// arma::vec lasso(arma::vec x, double l)
// {
//   cout<<"in lasso"<<endl;
//   return  sign(x)%max(abs(x) - l,zeros(size(x)));
// }

// arma::vec lasso_pos(arma::vec x, double l)
// {
//   cout<<"in lasso_pos"<<endl;
//   return  max(x - l, zeros(size(x)));
// }


// arma::vec scad(arma::vec x, double l)
// {
//  // cout<<"scadding";
//   int n = x.n_elem;
//   arma::vec z(n);
//   arma::vec abs = arma::abs(x);
//   arma::vec sgn = sign(x);
//   for(int i=0;i<n;i++){
//         z(i) = abs(i) > SCAD_A * l ? x(i)
//           : (abs(i) > 2 * l ? ((SCAD_A - 1) * x(i) - sgn(i) * SCAD_A * l) / (SCAD_A - 2) 
//           : sgn(i) * max((abs(i) - l), double(0)));
//   }
//   return z;
// }

void train(){}


// Prox_op match_prox_op(std::string type, int non_neg, double a){
//   if(type.compare("l1")==0 || type.compare("L1")==0){
//     if(non_neg==1)
//       return lasso_pos;
//     else
//       return lasso;
//   }
//   else if(type.compare("scad")==0 || type.compare("SCAD")==0){
//     if(non_neg==0)
//     {
//       SCAD_A = a;
//       return scad;
//     }
//   }
//   cout<<"here";
//   throw std::invalid_argument(type + " is not currently supported");
// }




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
  if(DEBUG) cout << "here Before mode = build_model" << endl;
  Model mod = build_model(X,Y,model_type);
  if(DEBUG) cout << "here After mod = build_mode" << endl;
  Solver sol = build_solver(mod,Omega_u,Omega_v,alpha_u,alpha_v,lambda_u,lambda_v,P_u,P_v,scad_a,non_neg,EPS,MAX_ITER,solver,SVD);
  if(DEBUG) cout<< "here after sol = build_solver"<<endl;
  // arma::mat U;  arma::vec s;  arma::mat V;  svd(U,s,V,X);
  // arma::mat Su; Su.eye(size(Omega_u)); Su += n*alpha_u*Omega_u;
  // arma::mat Sv; Sv.eye(size(Omega_v)); Sv += p*alpha_v*Omega_v;

  // double Lu = 1 / sol.grad_u_step;
  // double Lv = 1 / sol.grad_u_step;

  // Lasso prox_u;
  // Lasso prox_v;
  // Grad_line grad_u;
  // Grad_line grad_v;
  // grad_u.coef = Su;
  // grad_v.coef = Sv;

  // arma::vec v; //In mod now
  // arma::vec u;
  // // final result
  // if(SVD){
  //   v = V.col(0);
  //   u = U.col(0);
  // }
  // else{
  //   v = randu<arma::vec>(p);
  //   u = randu<arma::vec>(n);
  // }

  // value of the last round of update
  cout << size(sol.u) << endl;
  arma::vec oldu = zeros(size(sol.u));
  arma::vec oldv = zeros(size(sol.v));
  // last step
  arma::vec oldui = zeros(size(sol.u));
  arma::vec oldvi = zeros(size(sol.v));

  // stopping tolerance
  int iter = 0;
  int indu = 1;
  int indv = 1;
  int indo = 1;
  
  if(solver.compare("ISTA")==0)
  {
    while (indo > EPS && iter < MAX_ITER)
    {
      // ready for a new round of updates
      oldu = sol.u;
      oldv = sol.v;
      indu = 1;
      indv = 1;
      while (indu > EPS)
      {
        

        oldui = sol.u;
        if(DEBUG) cout<< "Before UCONST_VEC = X*sol.v; "<<endl;
        if(DEBUG) cout<< X.n_rows<<endl;
        UCONST_VEC = X*sol.v;
        if (DEBUG)
          cout << "AFTER UCONST_VEC = X*sol.v; " << endl;

        // change stepsize
        if(DEBUG)
          cout << "ol.u = (sol.prox_u(s" << endl;
        sol.u = (sol.prox_u(sol.u - sol.grad_u(sol.u)*sol.grad_u_step, sol.prox_u_step));
        if(DEBUG)
          cout << "AFTER ol.u = (sol.prox_u(s" << endl;
        norm(sol.u) > 0 ? sol.u /= mat_norm(sol.u, sol.Su) : sol.u.zeros();
        
        indu = norm(sol.u - oldui) / norm(oldui);
        if (DEBUG)
          cout << "AFTER iindu = norm(sol.u - oldui) / norm" << endl;
      }

      
      while (indv > EPS)
      {
        oldvi = sol.v;
        arma::mat Xt = X.t();
        if(DEBUG)
        cout<<Xt.n_cols <<'\t' << Xt.n_rows<<endl;
        VCONST_VEC = X.t() * sol.u;
        if(DEBUG)
          cout << "sol.v  = sol.prox_v(" << endl;
        sol.v = sol.prox_v(sol.v - sol.grad_v(sol.v)*sol.grad_v_step, 
                            sol.prox_v_step);
        norm(sol.v) > 0 ? sol.v /= mat_norm(sol.v, sol.Sv): sol.v.zeros();

        indv = norm(sol.v - oldvi) / norm(oldvi);
      }

      indo = norm(oldu - sol.u) / norm(oldu) + norm(oldv - sol.v) / norm(oldv);
      iter++;
    }
  }
  // else {
  //   // Doing FISTA

  //   arma::vec yu = u;
  //   arma::vec yv = v;
  //   double oldt = 1;
  //   double newt = 1;

  //   while (indo > EPS && iter < MAX_ITER)
  //   {
  //     // ready for a new round of updates
  //     oldu = u;
  //     oldv = v;
  //     indu = 1;
  //     indv = 1;
 
  //     oldt = 1;
  //     yu = u;
  //     while (indu > EPS)
  //     {
  //       oldui = u;
  //       oldt = newt;
  //       u = prox_u.f(yu + (X * v - Su * yu) / Lu, lambda_u / Lu);
  //       norm(u) > 0 ? u /= mat_norm(u, Su) : u.zeros();
  //       cout<<oldt<<endl;
  //       newt = 0.5 * (1 + sqrt(1 + 4 * oldt*oldt));
  //       yu = u + (oldt - 1) / newt * (u - oldui);
  //       indu = norm(u - oldui) / norm(oldui);
  //     }

  //     oldt = 1;
  //     yv = v;
  //     while (indv > EPS)
  //     {
       
  //       oldvi = v;
  //       oldt = newt;
  //       v = prox_v.f(yv + (X.t() * u - Sv * yv) / Lv, lambda_v / Lv);
  //       norm(v) > 0 ? v /= mat_norm(v, Sv) : v.zeros();
        
  //       newt = 0.5 * (1 + sqrt(1 + 4 * oldt*oldt));

  //       indv = norm(v - oldvi) / norm(oldvi);
  //     }

  //     indo = norm(oldu - u) / norm(oldu) + norm(oldv - v) / norm(oldv);
  //     iter++;
  //   }
  // }




  sol.u = sol.u/norm(sol.u);
  sol.v = sol.v/norm(sol.v);
  double d = as_scalar(sol.u.t() * mod.X * sol.v);
  return Rcpp::List::create(
    Rcpp::Named("u") = sol.u,
    Rcpp::Named("v") = sol.v,
    Rcpp::Named("d") = d,
    Rcpp::Named("DeflatedX") = mod.X - d * sol.u * sol.v.t());

}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be autoarma::matically
// run after the compilation.
//

/*** R

*/
