#include "common.h"

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