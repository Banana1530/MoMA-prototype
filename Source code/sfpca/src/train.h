#ifndef _TRAIN_H_
#define _TRAIN_H_
#include "common.h"
#include "model.h"
#include "solver.h"

Rcpp::List train(Model mod, Solver sol)
{ // value of the last round cout

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

    if (sol.solver_type == ISTA)
    {
        while (indo > sol.EPS && iter < sol.MAX_ITER)
        {
            // ready for a new round of updates
            oldu = sol.u;
            oldv = sol.v;
            indu = 1;
            indv = 1;
            while (indu > sol.EPS)
            {

                oldui = sol.u;
                UCONST_VEC = mod.X * sol.v;
                // change stepsize
                sol.u = (sol.prox_u(sol.u - sol.grad_u(sol.u) * sol.grad_u_step, sol.prox_u_step));
                norm(sol.u) > 0 ? sol.u /= mat_norm(sol.u, sol.Su) : sol.u.zeros();

                indu = norm(sol.u - oldui) / norm(oldui);
                if (DEBUG)
                    cout << "AFTER iindu = norm(sol.u - oldui) / norm" << endl;
            }

            while (indv > sol.EPS)
            {
                oldvi = sol.v;
                VCONST_VEC = mod.X.t() * sol.u;
                sol.v = sol.prox_v(sol.v - sol.grad_v_step * sol.grad_v(sol.v),
                                   sol.prox_v_step);
                norm(sol.v) > 0 ? sol.v /= mat_norm(sol.v, sol.Sv) : sol.v.zeros();

                indv = norm(sol.v - oldvi) / norm(oldvi);
            }

            indo = norm(oldu - sol.u) / norm(oldu) + norm(oldv - sol.v) / norm(oldv);
            iter++;
        }
    }
    // else {
    // Doing FISTA

    //   arma::vec yu = u;
    //   arma::vec yv = v;
    //   double oldt = 1;
    //   double newt = 1;

    //   while (indo > sol.EPS && iter < sol.MAX_ITER)
    //   {
    //     // ready for a new round of updates
    //     oldu = u;
    //     oldv = v;
    //     indu = 1;
    //     indv = 1;

    //     oldt = 1;
    //     yu = u;
    //     while (indu > sol.EPS)
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
    //     while (indv > sol.EPS)
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

    // form final result
    sol.u = sol.u / norm(sol.u);
    sol.v = sol.v / norm(sol.v);
    double d = as_scalar(sol.u.t() * mod.X * sol.v);
    return Rcpp::List::create(
        Rcpp::Named("u") = sol.u,
        Rcpp::Named("v") = sol.v,
        Rcpp::Named("d") = d,
        Rcpp::Named("DeflatedX") = mod.X - d * sol.u * sol.v.t());
}

#endif 