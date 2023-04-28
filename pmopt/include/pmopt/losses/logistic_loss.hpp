#pragma once

#include <iostream>
#include "pmopt/losses/squaredhinge_loss.hpp"
#include "pmopt/parameters.hpp"
#include "pmopt/variables.hpp"

namespace PMOpt
{

    struct LogisticLoss
    {
        static constexpr double SMOOTHNESS = 0.25;
        static constexpr double ZERO = 1e-14;
        static const bool NUMERIC_OBJECTIVE = false;

        static unsigned int _maxiter_middle;
        static unsigned int _maxiter_inner;
        static double _criterion_scale;
        static double _linsearch_update;

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
            auto w = params._w[i_instance];
            auto y = params._y[i_instance];
            auto r = std::exp( - y * vars._y_pred[i_instance]);
            return w * y * r / (1 + r);
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
            const auto & x_indices = * params._instance_indices;
            const auto & x_values = * params._ones;
            double total_update = 0.0;
            for (unsigned int it_mid = 0; it_mid < _maxiter_middle; ++it_mid)
            {
                // compute loss and its derivatives
                double loss, d_loss, dd_loss;
                evaluate_new(loss, d_loss, dd_loss, x_indices, x_values, vars, params, total_update);
                dd_loss = std::max(dd_loss, ZERO);

                // compute newton direction
                auto d = - d_loss / dd_loss;

                // if direction is small then break
                if (std::abs(d) < ZERO)
                    break;

                // find optimal update amount by line search
                auto criterion = d_loss * d;
                for (unsigned int it_in = 0; it_in < _maxiter_inner; ++it_in)
                {
                    // compute updated loss
                    double loss_new, dummy;
                    evaluate_new(loss_new, dummy, dummy, x_indices, x_values, vars, params, d + total_update);

                    // check optimality condition
                    if (loss_new - loss <= _criterion_scale * criterion)
                        break;
                    
                    // update direction and criterion
                    d *= _linsearch_update;
                    criterion *= _linsearch_update;
                }

                // update total_update
                total_update += d;
            }

            return vars._intercept + total_update;
        }


        /**
         * @brief compute new intercept value first time
         * @param vars (unused)
         * @param params 
         * @return double 
         */
        inline static double intercept_init(const Variables &, const Parameters & params) noexcept
        {
            double sum_w_positive = 0.0;
            double sum_w_negative = 0.0;
            for (auto i_instance : * params._instance_indices)
            {
                auto y = params._y[i_instance];
                auto w = params._w[i_instance];
                auto m = params._mask[i_instance];
                if (y > 0)
                    sum_w_positive += w * m;
                else
                    sum_w_negative += w * m;
            }
            return std::log(sum_w_positive / sum_w_negative);
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
            // compute loss and its derivatives 
            double loss, d_loss, dd_loss;
            evaluate_new(loss, d_loss, dd_loss, x_indices, x_values, vars, params, 0);

            // add regularization
            auto lam_l1 = params._lambda_l1 * vars._reg_weight[i_feature];
            auto lam_l2 = params._lambda_l2 * vars._reg_weight[i_feature];
            auto coef = vars._coef[i_feature];
            loss = (loss + lam_l2 * coef * coef * 0.5) / lam_l1;
            d_loss = (d_loss + lam_l2 * coef) / lam_l1;
            dd_loss = std::max((dd_loss + lam_l2) / lam_l1, ZERO);

            // compute newton direction            
            auto d = SquaredHingeLoss::newton_direction(d_loss, dd_loss, coef);

            // new coefficient
            auto coef_new = coef + d;

            // if direction is very small then return
            if (std::abs(d) < ZERO)
                return coef_new;

            // find optimal update amount by line search
            auto criterion = d_loss * d + std::abs(coef_new) - std::abs(coef);
            for (unsigned int it = 0; it < _maxiter_inner; ++it)
            {
                auto d_norm = std::abs(coef_new) - std::abs(coef);

                // compute loss after update
                double loss_new, dummy;
                evaluate_new(loss_new, dummy, dummy, x_indices, x_values, vars, params, d);

                // add regularization
                loss_new = (loss_new + lam_l2 * coef_new * coef_new * 0.5) / lam_l1;

                // check optimality
                if (d_norm + loss_new - loss <= _criterion_scale * criterion)
                    break;

                // update direction and criterion
                d *= _linsearch_update;
                criterion *= _linsearch_update;
                coef_new = coef + d;
            }

            return coef_new;
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
            primal = 0.0;
            dual = 0.0;
            scale = 1.0;

            if (params._l2_scale || params._lambda_l2 < 1e-16)
                scale = std::min(params._lambda_l1 / vars._max_score, 1.0);

            for (auto i_instance : * params._instance_indices)
            {
                auto w = params._w[i_instance];
                auto r = std::exp( - params._y[i_instance] * vars._y_pred[i_instance]);
                auto m = params._mask[i_instance];

                // primal loss function
                primal += w * std::log(1 + r) * m;

                // dual loss function 
                auto t = scale * r / (1 + r);

                if (t <= 0 || t >= 1) [[unlikely]]
                    continue;

                dual -= w * t * std::log(t) * m;
                dual -= w * (1 - t) * std::log(1 - t) * m;
            }
        }


        /**
         * @brief compute loss value and its derivatives (updated by difference)
         * @param[out] loss 
         * @param[out] d_loss 
         * @param[out] dd_loss 
         * @param x_indices 
         * @param x_values 
         * @param vars
         * @param params 
         * @param difference
         */
        inline static void evaluate_new(
            double & loss,
            double & d_loss,
            double & dd_loss,
            const std::vector<unsigned int> & x_indices,
            const std::vector<double> & x_values,
            const Variables & vars,
            const Parameters & params,
            double difference
        )
        noexcept
        {
            loss = 0.0;
            d_loss = 0.0;
            dd_loss = 0.0;

            for (unsigned int i_occ = 0; i_occ < x_indices.size(); ++i_occ)
            {
                auto i_instance = x_indices[i_occ];
                auto x = x_values[i_occ];
                auto w = params._w[i_instance];
                auto y = params._y[i_instance];
                auto r = std::exp( - y * (vars._y_pred[i_instance] + x * difference));
                auto m = params._mask[i_instance];

                loss += w * std::log(1 + r) * m;
                d_loss -= w * y * x * r / (1 + r) * m;
                dd_loss += w * x * x / (1 + r) / (1 + r) * r * m;
            }
        }
    };
}