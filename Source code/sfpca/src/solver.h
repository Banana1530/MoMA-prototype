#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "common.h"
#include "model.h"
typedef arma::vec (*Prox_op)(arma::vec, double);
typedef arma::vec (*Grad)(arma::vec);

typedef enum { ISTA,
               FISTA } SOLVER_TYPE;

arma::vec Grad_nosmoothv(arma::vec x)
{
    return -VCONST_VEC + x;
}
arma::vec Grad_nosmoothu(arma::vec x)
{
    return -UCONST_VEC + x;
}

arma::vec Grad_linev(arma::vec x)
{
    if (DEBUG)
        cout << "Grad_linv" << endl;
    return -VCONST_VEC + VCOEF * x;
}
arma::vec Grad_lineu(arma::vec x)
{
    return -UCONST_VEC + UCOEF * x;
}

arma::vec Lasso(arma::vec x, double l)
{
    if (DEBUG)
        cout << "lassoing" << endl;
    return sign(x) % max(abs(x) - l, zeros(size(x)));
};

arma::vec Lasso_pos(arma::vec x, double l)
{
    return max(x - l, zeros(size(x)));
};

arma::vec Scad(arma::vec x, double l)
{
    if (DEBUG)
        cout << "scadding";

    int n = x.n_elem;
    arma::vec z(n);
    arma::vec abs = arma::abs(x);
    arma::vec sgn = sign(x);
    for (int i = 0; i < n; i++)
    {
        z(i) = abs(i) > GAMMA * l ? x(i)
                                  : (abs(i) > 2 * l ? ((GAMMA - 1) * x(i) - sgn(i) * GAMMA * l) / (GAMMA - 2)
                                                    : sgn(i) * max((abs(i) - l), double(0)));
    }
    return z;
}

SOLVER_TYPE string_to_ST(std::string solver_type_string)
{
    if (solver_type_string.compare("ISTA") == 0)
        return ISTA;
    else if (solver_type_string.compare("FISTA") == 0)
        return FISTA;
    else
        throw std::invalid_argument(solver_type_string + " is not currently supported");
}

#endif