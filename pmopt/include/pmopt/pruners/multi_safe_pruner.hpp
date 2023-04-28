#pragma once

#include "pmopt/optimization/multi_screening.hpp"
#include "pmopt/pruners/safe_pruner.hpp"
#include "pmopt/logger.hpp"

namespace PMOpt
{
    template < class Loss, class XGenerator, class WGenerator >
    struct MultiSafePruner : public SafePruner<Loss, XGenerator, WGenerator>
    {
        const std::vector<Variables> & _vars;
        std::vector<std::vector<double>> _distance_matrix;
        std::vector<unsigned int> _including;
        unsigned int _n_vars;

        /**
         * @brief Construct a new Safe Pruner object
         * @param logger
         * @param params 
         * @param n_vars
         * @param iter
         */
        MultiSafePruner(
            const std::vector<Variables> & vars,
            const Parameters & params,
            unsigned int n_vars,
            const std::string & loginfo
        ) : 
            SafePruner<Loss, XGenerator, WGenerator>(vars[0], params),
            _vars(vars),
            _n_vars(n_vars)
        {
            // compute distance matrix and inclusion relationship
            Optimization::check_inclusion<Loss>(_distance_matrix, _including, vars, params, n_vars);
            
            Logger logger(std::cout);
            if (Logger::_verbose <= Logger::INFO)
            {
                for (unsigned int i_vars_1 = 0; i_vars_1 < n_vars; ++i_vars_1)
                {
                    for (unsigned int i_vars_2 = i_vars_1 + 1; i_vars_2 < n_vars; ++i_vars_2)
                    {
                        const auto & vars_1 = _vars[i_vars_1];
                        const auto & vars_2 = _vars[i_vars_2];
                        logger(Logger::INFO, loginfo, " INTERSECTION (", i_vars_1, ", ", i_vars_2, ')',
                            " include: ", _including[i_vars_1], ' ', _including[i_vars_2],
                            " radius: ", vars_1._radius, ' ', vars_2._radius,
                            " distance: ", _distance_matrix[i_vars_1][i_vars_2]);
                    }
                }
            }
        }


        /**
         * @brief check screening pruning
         * @tparam PatternType 
         * @param[out] object
         * @param projected 
         * @return bool
         */
        template < class Symbol, class Positions >
        bool evaluate(
            Activeset<Symbol> & activeset,
            const ProjectedDatabase<Symbol, Positions> & projected
        )
        const noexcept
        {
            auto n = this->_params._n_training;
            auto x = XGenerator::generate(projected, this->_params);
            auto w = WGenerator::generate(projected, this->_params);
            const auto & x_i = projected._occurrence;
            const auto & x_v = * x;

            // compute sum of x for radial compoent of screening sphere
            double x_sum, x_sum_sq;
            Optimization::sum_square(x_sum, x_sum_sq, x_v, x_i, this->_params._mask);

            bool is_screened = false;

            // single reference pruning
            std::vector<double> dot_upper(_n_vars), dot_lower(_n_vars);
            for (unsigned int i_vars = 0; i_vars < _n_vars; ++i_vars)
            {
                if (_including[i_vars])
                    continue;

                const auto & vars = _vars[i_vars];

                // compute bound dot product
                Optimization::dot_bound<Loss>(dot_upper[i_vars], dot_lower[i_vars], x_i, x_v, vars, this->_params);

                // compute crieria
                double pruning_criterion, screening_criterion;
                Optimization::safe_pruning(pruning_criterion, screening_criterion,
                    vars._scale, vars._radius, dot_upper[i_vars], dot_lower[i_vars], x_sum, x_sum_sq, n);

                // the pattern is pruned if pruned by at least one variables
                if (pruning_criterion < this->_params._lambda_l1 * w)
                    return true;

                // ths pattern is screened if screening by at least one variables
                // update flag
                is_screened = is_screened || (screening_criterion < this->_params._lambda_l1 * w);
            }

            // two reference pruning
            for (unsigned int i_vars_1 = 0; i_vars_1 < _n_vars; ++i_vars_1)
            {
                if (_including[i_vars_1])
                    continue;

                for (unsigned int i_vars_2 = i_vars_1 + 1; i_vars_2 < _n_vars; ++i_vars_2)
                {
                    if (_including[i_vars_2])
                        continue;

                    const auto & vars_1 = _vars[i_vars_1];
                    const auto & vars_2 = _vars[i_vars_2];
                    double d = _distance_matrix[i_vars_1][i_vars_2];

                    // compute critera
                    double pruning_criterion, screening_criterion;
                    Optimization::safe_pruning_multi(pruning_criterion, screening_criterion,
                        d, vars_1._scale, vars_2._scale, vars_1._radius, vars_2._radius,
                        dot_upper[i_vars_1], dot_upper[i_vars_2], dot_lower[i_vars_1], dot_lower[i_vars_2], x_sum, x_sum_sq, n);

                    // the pattern is pruned if pruned by at least one variables
                    if (pruning_criterion < this->_params._lambda_l1 * w)
                        return true;

                    // ths pattern is screened if screening by at least one variables
                    // update flag
                    is_screened = is_screened || (screening_criterion < this->_params._lambda_l1 * w);
                }
            }


            if (!is_screened)
            {
                activeset._patterns.push_back(projected._pattern);
                activeset._x_indices.emplace_back(
                    std::make_shared<const std::vector<unsigned int>>(projected._occurrence)
                );
                activeset._x_values.emplace_back(x);
                activeset._reg_weight.push_back(w);
                activeset._x_sum_sq.push_back(x_sum_sq);
                activeset._x_sum.push_back(x_sum);
            }

            return false;
        }


        /**
         * @brief Evaluate pruning / screening score for numerical input
         * @param[out] activeset
         * @param i_feature
         */
        template < class Pattern >
        void evaluate_numeric(
            Activeset<Pattern> & activeset,
            unsigned int i_feature
        )
        const noexcept
        {
            auto n = this->_params._instance_indices->size();
            auto x = this->_params._X_num[i_feature];
            const auto & x_i = * this->_params._instance_indices;
            const auto & x_v = * x;

            // compute sum of x for radial compoent of screening sphere
            double x_sum, x_sum_sq;
            Optimization::sum_square(x_sum, x_sum_sq, x_v, x_i, this->_params._mask);

            // single reference pruning
            std::vector<double> dot_value(_n_vars);
            for (unsigned int i_vars = 0; i_vars < _n_vars; ++i_vars)
            {
                // skip inactive sphere
                if (_including[i_vars])
                    continue;

                const auto & vars = _vars[i_vars];

                // compute bound dot product
                Optimization::dot<Loss>(dot_value[i_vars], x_i, x_v, vars, this->_params);

                // compute screenign criterion
                double criterion;
                Optimization::safe_screening(criterion,
                    vars._scale, vars._radius, dot_value[i_vars], x_sum, x_sum_sq, n);

                // ths feature is screened if screening by at least one variables
                if (criterion < this->_params._lambda_l1)
                    return;
            }

            // two reference pruning
            for (unsigned int i_vars_1 = 0; i_vars_1 < _n_vars; ++i_vars_1)
            {
                // skip inactive sphere
                if (_including[i_vars_1])
                    continue;

                for (unsigned int i_vars_2 = i_vars_1 + 1; i_vars_2 < _n_vars; ++i_vars_2)
                {
                    // skip inactive sphere
                    if (_including[i_vars_2])
                        continue;

                    const auto & vars_1 = _vars[i_vars_1];
                    const auto & vars_2 = _vars[i_vars_2];
                    double d = _distance_matrix[i_vars_1][i_vars_2];

                    // compute critera
                    double criterion;
                    Optimization::safe_screening_multi(criterion,
                        d, vars_1._scale, vars_2._scale, vars_1._radius, vars_2._radius,
                        dot_value[i_vars_1], dot_value[i_vars_2], x_sum, x_sum_sq, n);

                    // ths feature is screened if screening by at least one variables
                    if (criterion < this->_params._lambda_l1)
                        return;
                }
            }

            // not screened
            activeset._keys.push_back(this->_params._X_num_names[i_feature]);
            activeset._x_indices.push_back(this->_params._instance_indices);
            activeset._x_values.push_back(x);
            activeset._x_sum_sq.push_back(x_sum_sq);
            activeset._x_sum.push_back(x_sum);
            activeset._reg_weight.push_back(1.0);
        }
    };
}