#include "pmopt/optimization/screening.hpp"

using namespace PMOpt;

void Optimization::norm(double & norm_l1, double & norm_l2, const Variables & vars) noexcept
{
    norm_l1 = 0.0;
    norm_l2 = 0.0;
    for (auto j : vars._coef_indices)
    {
        norm_l1 += vars._reg_weight[j] * std::abs(vars._coef[j]);
        norm_l2 += vars._reg_weight[j] * vars._coef[j] * vars._coef[j];
    }
}