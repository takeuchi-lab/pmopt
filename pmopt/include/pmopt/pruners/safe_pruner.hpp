#pragma once

#include <iostream>
#include "pmopt/optimization/basic.hpp"
#include "pmopt/optimization/screening.hpp"
#include "pmopt/parameters.hpp"
#include "pmopt/pattern.hpp"
#include "pmopt/variables.hpp"

namespace PMOpt
{
    template < class Loss, class XGenerator, class WGenerator >
    struct SafePruner
    {

        const static unsigned int PRIORITY = 10;
        const std::string INDENT = "  -- ";

        const Variables & _vars;
        const Parameters & _params;


        /**
         * @brief Construct a new Safe Pruner object
         * @param vars
         * @param params 
         */
        SafePruner(
            const Variables & vars,
            const Parameters & params
        ) : _vars(vars), _params(params) {}


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
            auto n = _params._n_training;
            auto x = XGenerator::generate(projected, _params);
            const auto & x_i = projected._occurrence;
            const auto & x_v = * x;

            // compute sum of x for radial compoent of screening sphere
            double x_sum, x_sum_sq;
            Optimization::sum_square(x_sum, x_sum_sq, x_v, x_i, _params._mask);

            // compute bound dot product
            double dot_upper, dot_lower;
            Optimization::dot_bound<Loss>(dot_upper, dot_lower, x_i, x_v, _vars, _params);

            // compute crieria
            double pruning_criterion, screening_criterion;
            Optimization::safe_pruning(pruning_criterion, screening_criterion,
                _vars._scale, _vars._radius, dot_upper, dot_lower, x_sum, x_sum_sq, n);

            // evaluate pruning and screening condition
            auto w = WGenerator::generate(projected, this->_params);
            bool is_pruned = pruning_criterion < _params._lambda_l1 * w;
            bool is_screened = screening_criterion < _params._lambda_l1 * w;

            // save active pattern
            if (!is_screened)
            {
                // Here not adding string into activeset._keys,
                // because we cannot use enumerator which have symbol information.
                // After mining, keys will be add by using enumerator.pat_to_str.
                activeset._x_indices.push_back(
                    std::make_shared<const std::vector<unsigned int>>(projected._occurrence)
                );
                activeset._x_values.push_back(x);
                activeset._patterns.push_back(projected._pattern);
                activeset._x_sum_sq.push_back(x_sum_sq);
                activeset._x_sum.push_back(x_sum);
                activeset._reg_weight.push_back(w);
            }

            return is_pruned;
        }


        /**
         * @brief Evaluate screening score for numerical input
         * @param[out] activeset
         * @param index
         */
        template < class Symbol >
        void evaluate_numeric(
            Activeset<Symbol> & activeset,
            unsigned int i_feature
        )
        const noexcept
        {
            auto n = _params._instance_indices->size();
            auto x = _params._X_num[i_feature];
            const auto & x_i = * _params._instance_indices;
            const auto & x_v = * x;

            // compute sum of x for radial compoent of screening sphere
            double x_sum, x_sum_sq;
            Optimization::sum_square(x_sum, x_sum_sq, x_v, x_i, _params._mask);

            // compute bound dot product
            double dot_value;
            Optimization::dot<Loss>(dot_value, x_i, x_v, _vars, _params);

            // compute screening criterion
            double criterion;
            Optimization::safe_screening(criterion, _vars._scale, _vars._radius, dot_value, x_sum, x_sum_sq, n);

            // evaluate screening condition
            bool is_screened = criterion < _params._lambda_l1;


            if (!is_screened)
            {
                activeset._keys.push_back(_params._X_num_names[i_feature]);
                activeset._x_indices.push_back(_params._instance_indices);
                activeset._x_values.push_back(x);
                activeset._x_sum_sq.push_back(x_sum_sq);
                activeset._x_sum.push_back(x_sum);
                activeset._reg_weight.push_back(1.0);
            }
        }
    };


    // struct for save lammax result
    template < class Symbol >
    struct LammaxObject
    {
        std::shared_ptr<const Pattern<Symbol>> _pattern;
        std::string _name;
        double _max_score = 0.0;
    };


    template < class Loss, class XGenerator, class WGenerator >
    struct LammaxCalculator : public SafePruner< Loss, XGenerator, WGenerator >
    {
        using SafePruner<Loss, XGenerator, WGenerator>::SafePruner;

        /**
         * @brief Evaluate pruning / screening score for pattern
         * @tparam PatternType 
         * @param[out] mined
         * @param projected 
         * @return bool
         */
        template < class Symbol, class Positions >
        bool evaluate(
            LammaxObject<Symbol> & mined,
            const ProjectedDatabase<Symbol, Positions> & projected
        )
        const noexcept
        {
            const auto & x_i = projected._occurrence;
            const auto & x_v = * XGenerator::generate(projected, this->_params);

            // compute bound dot product
            double dot_upper, dot_lower;
            Optimization::dot_bound<Loss>(dot_upper, dot_lower, x_i, x_v, this->_vars, this->_params);

            // compute crieria
            // in lammax mining, no radial component required
            double pruning_criterion, screening_criterion;
            Optimization::safe_pruning(pruning_criterion, screening_criterion,
                1.0, 0.0, dot_upper, dot_lower, 0.0, 0.0, 1);

            // evaluate pruning and screening condition
            auto w = WGenerator::generate(projected, this->_params);
            bool is_pruned = pruning_criterion / w <= mined._max_score;
            bool is_screened = screening_criterion / w <= mined._max_score;

            // update max score
            if (!is_screened)
            {
                mined._max_score = screening_criterion;
                mined._pattern = projected._pattern;
                mined._name.clear();
            }

            return is_pruned;
        }


        /**
         * @brief Evaluate screening score for numerical input
         * @tparam PatternType 
         * @param[out] mined
         * @param i_feature
         */
        template < class Symbol >
        void evaluate_numeric(
            LammaxObject<Symbol> & mined,
            unsigned int i_feature
        )
        const noexcept
        {
            const auto & key = this->_params._X_num_names[i_feature];
            const auto & x_i = * this->_params._instance_indices;
            const auto & x_v = * this->_params._X_num[i_feature];

            // compute bound dot product
            double dot_value;
            Optimization::dot<Loss>(dot_value, x_i, x_v, this->_vars, this->_params);

            // compute screening criterion
            // no scaling and no radial component
            double criterion;
            Optimization::safe_screening(criterion, 1.0, 0.0, dot_value, 0.0, 0.0, 1);

            // evaluate screening condition
            bool is_screened = criterion <= mined._max_score;


            if (!is_screened)
            {
                mined._max_score = criterion;
                mined._pattern.reset();
                mined._name = key;
            }
        }
    };


    // struct for accumulating dual norm
    struct DualNormObject
    {
        double _sum_score = 0.0;
    };


    template < class Loss, class XGenerator, class WGenerator >
    struct DualNormCalculator : public SafePruner< Loss, XGenerator, WGenerator >
    {
        using SafePruner<Loss, XGenerator, WGenerator>::SafePruner;        

        /**
         * @brief compute dual norm
         * @tparam PatternType 
         * @param[out] mined
         * @param projected 
         * @return bool
         */
        template < class Symbol, class Positions >
        bool evaluate(
            DualNormObject & mined,
            const ProjectedDatabase<Symbol, Positions> & projected
        )
        const noexcept
        {

            const auto & x_i = projected._occurrence;
            const auto & x_v = * XGenerator::generate(projected, this->_params);

            // compute bound dot product
            double dot_upper, dot_lower;
            Optimization::dot_bound<Loss>(dot_upper, dot_lower, x_i, x_v, this->_vars, this->_params);

            // compute crieria
            // dual norm calculation, no scaling and no radial component
            double pruning_criterion, screening_criterion;
            Optimization::safe_pruning(pruning_criterion, screening_criterion,
                1.0, 0.0, dot_upper, dot_lower, 0.0, 0.0, 1);

            // evaluate pruning and screening condition
            auto w = WGenerator::generate(projected, this->_params);
            auto lam_l1 = this->_params._lambda_l1 * w;
            bool is_pruned = pruning_criterion <= lam_l1;
            bool is_screened = screening_criterion <= lam_l1;


            // add score
            if (!is_screened)
                mined._sum_score -= (screening_criterion - lam_l1) * (screening_criterion - lam_l1) / w;

            return is_pruned;
        }



        /**
         * @brief Evaluate screening score for numerical input
         * @param[out] mined
         * @param i_feature
         */
        void evaluate_numeric(
            DualNormObject & mined,
            unsigned int i_feature
        )
        const noexcept
        {
            const auto & x_i = * this->_params._instance_indices;
            const auto & x_v = * this->_params._X_num[i_feature];

            // compute bound dot product
            double dot_value;
            Optimization::dot<Loss>(dot_value, x_i, x_v, this->_vars, this->_params);

            // compute screening criterion
            // no scaling and no radial component
            double criterion;
            Optimization::safe_screening(criterion, 1.0, 0.0, dot_value, 0.0, 0.0, 1);
            
            // evaluate screening condition
            auto lam_l1 = this->_params._lambda_l1;
            bool is_screened = criterion <= lam_l1;

            // add score
            if (!is_screened)
                mined._sum_score -= (criterion - lam_l1) * (criterion - lam_l1);
        }
    };
}