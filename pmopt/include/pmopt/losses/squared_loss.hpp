#pragma once

#include <iostream>
#include "pmopt/parameters.hpp"
#include "pmopt/variables.hpp"

namespace PMOpt
{

    struct SquaredLoss
    {
        static constexpr double SMOOTHNESS = 1.0;
        static const bool NUMERIC_OBJECTIVE = true;

        /**
         * @brief get dual variables
         * @param vars 
         * @param params
         * @param i_instance
         */
        inline static double dual(
            const Variables & vars,
            const Parameters & params,
            unsigned int i_instance
        )
        noexcept
        {
            return params._w[i_instance] * (params._y[i_instance] - vars._y_pred[i_instance]);
        }


        /**
         * @brief compute new intercept value
         * @param vars 
         * @param params 
         * @return double 
         */
        inline static double intercept(
            const Variables & vars,
            const Parameters & params
        )
        noexcept
        {
            double intercept = 0;
            for (auto i_instance : * params._instance_indices)
                intercept += params._mask[i_instance] * SquaredLoss::dual(vars, params, i_instance);
            return intercept + vars._intercept;
        }

        /**
         * @brief compute new intercept value
         * @param vars 
         * @param params 
         * @return double 
         */
        inline static double intercept_init(
            const Variables & vars,
            const Parameters & params
        )
        noexcept
        {
            return SquaredLoss::intercept(vars, params);
        }


        /**
         * @brief compute new coefficient
         * @param x_indices
         * @param x_values
         * @param vars 
         * @param params 
         * @param i_feature
         * @return double 
         */
        inline static double coefficient(
            const std::vector<unsigned int> & x_indices,
            const std::vector<double> & x_values,
            const Variables & vars,
            const Parameters & params,
            unsigned int i_feature
        )
        noexcept
        {
            double xTr = 0.0;
            double xTx = 0.0;
            for (unsigned int i_occ = 0; i_occ < x_indices.size(); ++i_occ)
            {

                auto i_instance = x_indices[i_occ];
                auto x = x_values[i_occ];
                auto w = params._w[i_instance];
                auto r = params._y[i_instance] - vars._y_pred[i_instance];
                auto m = params._mask[i_instance];
                
                xTr += x * r * w * m;
                xTx += x * x * w * m;
            }
            xTr += vars._coef[i_feature] * xTx;

            auto lam_l1 = params._lambda_l1 * vars._reg_weight[i_feature];
            auto lam_l2 = params._lambda_l2 * vars._reg_weight[i_feature];
            if (xTr > lam_l1) 
                return (xTr - lam_l1) / (xTx + lam_l2);
            else if (xTr < - lam_l1)
                return (xTr + lam_l1) / (xTx + lam_l2);
            else
                return 0.0;
        }


        /**
         * @brief evaluate loss function
         * @param[out] primal
         * @param[out] dual
         * @param[out] scale
         * @param vars
         * @param params
         */
        inline static void evaluate(
            double & primal,
            double & dual,
            double & scale,
            const Variables & vars,
            const Parameters & params
        )
        noexcept
        {
            double dual_first = 0.0;
            double dual_second = 0.0;
            for (auto i_instance : * params._instance_indices)
            {
                auto y = params._y[i_instance];
                auto w = params._w[i_instance];
                auto r = y - vars._y_pred[i_instance];
                auto m = params._mask[i_instance];
                dual_first += r * r * w * m;
                dual_second += r * w * y * m;
            }

            primal = dual_first * 0.5;

            if (params._l2_scale || params._lambda_l2 < 1e-16)
            {
                double max_scale = params._lambda_l1 / vars._max_score;
                scale = std::min(std::max(dual_second / dual_first, - max_scale), max_scale);
            }

            dual_first *= vars._scale * vars._scale * 0.5;
            dual_second *= vars._scale;
            dual = - dual_first + dual_second;
        }
    };
}