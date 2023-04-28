#pragma once

#include "pmopt/losses/squared_loss.hpp"
#include "pmopt/utils/containers.hpp"
#include "pmopt/parameters.hpp"
#include "pmopt/variables.hpp"

namespace PMOpt
{

    struct SquaredHingeLoss
    {
        static constexpr double SMOOTHNESS = 1.0;
        static constexpr double INF = 1e+16;
        static constexpr double ZERO = 1e-8;
        static constexpr double ZERO_CHECK = 1e-14;
        static const bool NUMERIC_OBJECTIVE = false;

        static double _criterion_scale;
        static double _linsearch_update;
        static unsigned int _maxiter_inner;

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
            if (1 - params._y[i_instance] * vars._y_pred[i_instance] > 0)
                return params._w[i_instance] * (params._y[i_instance] - vars._y_pred[i_instance]);
            else
                return 0;
        }


        /**
         * @brief compute new intercept value
         * @param vars 
         * @param params 
         * @return double 
         */
        inline static double intercept(
            Variables & vars,
            const Parameters & params
        )
        noexcept
        {
            auto & sorted = vars._sorted_instances;
            sorted.clear();
            sorted.reserve(params._n_instances);

            for (auto i_instance : * params._instance_indices)
                if (params._mask[i_instance]) [[likely]]
                    sorted.push_back(i_instance);

            auto n = sorted.size();
            std::sort(sorted.begin(), sorted.end(), [&](unsigned int i1, unsigned int i2)
            {
                return params._y[i1] - vars._y_pred[i1] < params._y[i2] - vars._y_pred[i2];
            });

            // First, the misclassified instances at intercept = - infty will be computed.
            // In this case, all instances are classified as negative,
            // so all positive instances are misclassified
            double w_sum_lower = 0.0;
            double r_sum_lower = 0.0;
            const auto & y = params._y;
            const auto & w = params._w;
            for (auto i_instance : * params._instance_indices)
            {
                if (y[i_instance] < 0)
                    continue;

                // # of instances that 1 - y_i * (-infty) > 0
                // i.e. # of positive instances.
                auto r_i = y[i_instance] - vars._y_pred[i_instance];
                auto m = params._mask[i_instance];
                w_sum_lower += w[i_instance] * m;
                r_sum_lower += r_i * w[i_instance] * m;
            }

            double r_lower = - INF;
            double w_sum_in_upper = w_sum_lower;
            double r_sum_in_upper = r_sum_lower;
            double w_sum_in_lower = w_sum_lower;
            double r_sum_in_lower = r_sum_lower;
            double w_sum_upper = w_sum_lower;
            double r_sum_upper = r_sum_lower;
            for (unsigned int k = 0; k < n + 1; ++ k)
            {
                double r_k1 = - INF, r_k2 = - INF;
                if (k < n) [[likely]]
                    r_k1 = y[sorted[k]] - vars._y_pred[sorted[k]];
                if (k < n - 1) [[likely]]
                    r_k2 = y[sorted[k + 1]] - vars._y_pred[sorted[k + 1]];

                // Compute upper bound of the current residual
                double r_upper = (k < n) ? r_k1 : INF;

                // Update indices of upper bound
                // Newly added to the indices
                if (k < n && y[sorted[k]] < 0)
                {
                    w_sum_upper += w[sorted[k]];
                    r_sum_upper += r_k1 * w[sorted[k]];
                    w_sum_in_upper += w[sorted[k]];
                    r_sum_in_upper += r_k1 * w[sorted[k]];
                }

                // Update indices of upper bound
                // Newly remove from the indices
                if (k < n && y[sorted[k]] > 0)
                {
                    w_sum_in_upper -= w[sorted[k]];
                    r_sum_in_upper -= r_k1 * w[sorted[k]];
                }


                // If next two are the same, continue accumulating upper bound                
                if (k < n - 1 && std::abs(r_k2 - r_k1) < ZERO)
                    continue;

                // Compute the gradients at lower and upper bound
                double grad_lower = r_lower * w_sum_lower - r_sum_lower;
                double grad_upper = r_upper * w_sum_upper - r_sum_upper;

                // If the stationary is in this region
                if (grad_lower <= 0 && grad_upper > 0) [[unlikely]]
                {
                    assert(r_sum_in_lower / w_sum_in_lower >= r_lower);
                    assert(r_sum_in_lower / w_sum_in_lower < r_upper);
                    return r_sum_in_lower / w_sum_in_lower + vars._intercept;
                }

                // Next region
                r_lower = r_upper;
                r_sum_lower = r_sum_in_upper;
                w_sum_lower = w_sum_in_upper;
                r_sum_upper = r_sum_in_upper;
                w_sum_upper = w_sum_in_upper;
                r_sum_in_lower = r_sum_in_upper;
                w_sum_in_lower = w_sum_in_upper;

            }
            assert(false);
        }


        /**
         * @brief compute new intercept value first time
         * @param params 
         * @return double 
         */
        inline static double intercept_init(const Variables & vars, const Parameters & params) noexcept
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
            // compute loss and its derivatives
            double loss, d_loss, dd_loss, xTx;
            evaluate_new(loss, d_loss, dd_loss, xTx, x_indices, x_values, vars, params, 0);

            // compute objective by adding regularization
            auto lam_l1 = params._lambda_l1 * vars._reg_weight[i_feature];
            auto lam_l2 = params._lambda_l2 * vars._reg_weight[i_feature];
            auto coef = vars._coef[i_feature];
            loss = (loss + lam_l2 * coef * coef) / lam_l1 * 0.5;
            d_loss = (d_loss + lam_l2 * coef) / lam_l1;
            dd_loss = std::max((dd_loss + lam_l2) / lam_l1, ZERO);

            // compute newton direction
            auto d = newton_direction(d_loss, dd_loss, coef);

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
                auto optimality = d_norm + d_loss * d + xTx * d * d;

                // cheap optimality checking
                if (optimality <= _criterion_scale * criterion)
                    break;
                
                // expensive optimality checking
                double loss_new, dummy;
                evaluate_new(loss_new, dummy, dummy, dummy, x_indices, x_values, vars, params, d);

                // add regularization
                loss_new = (loss_new + lam_l2 * coef_new * coef_new) / lam_l1 * 0.5;

                // check optimality condition
                optimality = d_norm + loss_new - loss;
                if (optimality <= _criterion_scale * criterion)
                    break;
                
                // update direction
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
            double dual_first = 0.0;
            double dual_second = 0.0;
            for (auto i_instance : * params._instance_indices)
            {
                auto y = params._y[i_instance];
                auto w = params._w[i_instance];
                auto m = params._mask[i_instance];
                if (1 - y * vars._y_pred[i_instance] < 0)
                    continue;

                auto r = y - vars._y_pred[i_instance];
                dual_first += r * r * w * m;
                dual_second += r * y * w * m;
            }

            primal = dual_first * 0.5;

            if (params._l2_scale || params._lambda_l2 < 1e-16)
            {
                double max_scale = params._lambda_l1 / vars._max_score;
                scale = std::min(std::max(dual_second / dual_first, 0.0), max_scale);
            }

            dual_first *= vars._scale * vars._scale * 0.5;
            dual_second *= vars._scale;
            dual = - dual_first + dual_second;
        }


        /**
         * @brief compute loss value and its derivatives
         * @param[out] loss 
         * @param[out] d_loss 
         * @param[out] dd_loss 
         * @param[out] x_sum_square 
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
            double & x_sum_square,
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
            x_sum_square = 0.0;

            for (unsigned int i_occ = 0; i_occ < x_indices.size(); ++i_occ)
            {
                auto i_instance = x_indices[i_occ];
                auto x = x_values[i_occ];
                auto w = params._w[i_instance];
                auto y = params._y[i_instance];
                auto r = y - vars._y_pred[i_instance];
                auto m = params._mask[i_instance];

                x_sum_square += w * x * x * m;

                if (r * y < 0.0)
                    continue;

                r -= x * difference;

                loss += r * r * w * m;
                d_loss -= r * x * w * m;
                dd_loss += w * x * x * m;
            }
        }


        /**
         * @brief compute update direction of newton method
         * @param d_loss 
         * @param dd_loss 
         * @param coef 
         * @return double 
         */
        inline static double newton_direction(double d_loss, double dd_loss, double coef) noexcept
        {
            if (d_loss + 1 < coef * dd_loss)
                return - (d_loss + 1) / dd_loss;
            else if (d_loss - 1 > coef * dd_loss)
                return - (d_loss - 1) / dd_loss;
            else
                return - coef;
        }
    };
}