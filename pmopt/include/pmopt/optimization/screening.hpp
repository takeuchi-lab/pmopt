#pragma once

#include <cmath>
#include "pmopt/parameters.hpp"
#include "pmopt/variables.hpp"

namespace PMOpt
{
    namespace Optimization
    {

        /**
         * @brief compute upper / lower bound of dot product between X and dual
         * @tparam Loss: type of loss function
         * @param[out] upper upper bound of dot
         * @param[out] lower lower bound of dot
         * @param x_indices if x is pattern feature then this is occurrence else all training instances
         * @param x_values feature vector
         * @param vars
         * @param params
         */
        template < class Loss >
        inline void dot_bound(
            double & upper,
            double & lower,
            const std::vector<unsigned int> & x_indices,
            const std::vector<double> & x_values,
            const Variables & vars,
            const Parameters & params
        )
        noexcept
        {
            upper = 0.0;
            lower = 0.0;
            for (unsigned int i_occ = 0; i_occ < x_indices.size(); ++i_occ)
            {
                auto i_instance = x_indices[i_occ];
                auto m = params._mask[i_instance];
                auto dot = x_values[i_occ] * Loss::dual(vars, params, i_instance);
                (dot > 0.0) ? upper += dot * m : lower += dot * m;
            }
        }


        /**
         * @brief compute dot product between X and dual
         * @tparam Loss: type of loss function
         * @param[out] value dot product value
         * @param x_indices
         * @param x_values
         * @param vars
         * @param params
         */
        template < class Loss >
        inline void dot(
            double & value,
            const std::vector<unsigned int> & x_indices,
            const std::vector<double> & x_values,
            const Variables & vars,
            const Parameters & params
        )
        noexcept
        {
            value = 0.0;
            for (unsigned int i_occ = 0; i_occ < x_indices.size(); ++i_occ)
            {
                auto i_instance = x_indices[i_occ];
                auto m = params._mask[i_instance];
                value += x_values[i_occ] * Loss::dual(vars, params, i_instance) * m;
            }
        }


        /**
         * @brief compute pruning / screening criteria
         * @param[out] pruning_criterion 
         * @param[out] screening_criterion 
         * @param scale scale factor of dual variables
         * @param radius screening sphere radius 
         * @param dot_upper 
         * @param dot_lower 
         * @param x_sum 
         * @param x_sum_sq
         * @param n number of training instances 
         */
        inline void safe_pruning(
            double & pruning_criterion,
            double & screening_criterion,
            double scale,
            double radius,
            double dot_upper,
            double dot_lower,
            double x_sum,
            double x_sum_sq,
            unsigned int n
        )
        noexcept
        {
            pruning_criterion = std::max(dot_upper, - dot_lower) * scale;
            pruning_criterion += std::sqrt(x_sum_sq) * radius;
            screening_criterion = std::abs(dot_upper + dot_lower) * scale;
            screening_criterion += std::sqrt(x_sum_sq - x_sum * x_sum / n) * radius;
        }


        /**
         * @brief compute screening crierion
         * @param criterion 
         * @param scale
         * @param radius
         * @param dot_value 
         * @param x_sum
         * @param x_sum_sq
         * @param n 
         */
        inline void safe_screening(
            double & criterion,
            double scale,
            double radius,
            double dot_value,
            double x_sum,
            double x_sum_sq,
            unsigned int n
        )
        noexcept
        {
            criterion = std::abs(dot_value) * scale;
            criterion += std::sqrt(x_sum_sq - x_sum * x_sum / n) * radius;
        }


        /**
         * @brief compute l1-norm and l2-norm
         * @param[out] norm_l1
         * @param[out] norm_l2
         * @param vars
         */
        void norm(double &, double &, const Variables &) noexcept;

    }
}