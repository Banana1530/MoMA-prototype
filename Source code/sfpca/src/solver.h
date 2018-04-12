#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "common.h"
#include "model.h"
typedef arma::vec (*Prox_op)(arma::vec, double);
typedef arma::vec (*Grad)(arma::vec);
typedef enum { ISTA, FISTA } SOLVER_TYPE;

arma::vec Grad_nosmoothv(arma::vec x);
arma::vec Grad_nosmoothu(arma::vec x);

arma::vec Grad_linev(arma::vec x);
arma::vec Grad_lineu(arma::vec x);

arma::vec Lasso(arma::vec x, double l);

arma::vec Lasso_pos(arma::vec x, double l);

arma::vec Scad(arma::vec x, double l);

SOLVER_TYPE string_to_ST(std::string solver_type_string);

typedef struct
{
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

Solver build_solver(
    Model mod,
    arma::mat Omega_u,
    arma::mat Omega_v,
    double alpha_u,
    double alpha_v,
    double lambda_u,
    double lambda_v,
    std::string P_u, // assume for now they have same type of penalty
    std::string P_v,
    double scad_a,
    int non_neg, // 1 means activate non-negativity constraint
    double EPS,
    long MAX_ITER,
    std::string solver_type_string,
    bool SVD);

#endif