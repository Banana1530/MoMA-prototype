#include "solver.h"


double GAMMA;
arma::mat UCOEF;
arma::vec UCONST_VEC;
arma::mat VCOEF;
arma::vec VCONST_VEC;

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
    bool SVD)
{

    int n = mod.X.n_rows, p = mod.X.n_cols;

    if (DEBUG)
        cout << n << p << endl;
    // Find S_u --> L_u
    arma::mat U;
    arma::vec s;
    arma::mat V;
    if (DEBUG)
        cout << "here Before  svd(U, s, V, mod.X);" << endl;
    svd(U, s, V, mod.X);
    if (DEBUG)
        cout << " here After svd(U,s,V,mod.X)" << endl;
    Solver sol;
    if (DEBUG)
        cout << "here Before sol.Su.eye(size(Omega_u))" << endl;

    sol.Su.eye(size(Omega_u));
    sol.Su += n * alpha_u * Omega_u;
    if (DEBUG)
        myassert(Omega_v.n_cols == 200, "Wrong Omega_v input");
    sol.Sv.eye(size(Omega_v));
    sol.Sv += p * alpha_v * Omega_v;
    if (DEBUG)
    {
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
        if (DEBUG)
            cout << sol.v.n_elem << endl;
        sol.u = U.col(0);
        if (DEBUG)
            cout << sol.u.n_elem << endl;
        if (DEBUG)
            myassert(sol.u.n_elem == 199, "WRONG DIM OF U");
    }
    else
    {
        sol.v = randu<arma::vec>(p);
        sol.u = randu<arma::vec>(n);
    }

    // match gradient
    UCONST_VEC = mod.X * sol.v;
    VCONST_VEC = mod.X.t() * sol.u;
    if (abs(alpha_u) < 1e-9 && abs(alpha_v) < 1e-9)
    {

        sol.grad_u = Grad_nosmoothu;
        sol.grad_v = Grad_nosmoothv;
    }
    else
    {
        if (DEBUG)
            cout << "CHOOSE Grad_line" << endl;
        UCOEF = sol.Su;
        VCOEF = sol.Sv;
        myassert(VCOEF.n_rows == 200, "VCOEF sucks");
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
        else
        {
            if (DEBUG)
                cout << "CHOOSING Lasso" << endl;
            sol.prox_u = Lasso;
            sol.prox_v = Lasso;
        }
    }
    else if (P_u.compare("scad") == 0 || P_u.compare("SCAD") == 0)
    {
        GAMMA = scad_a;
        if (non_neg == 0)
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