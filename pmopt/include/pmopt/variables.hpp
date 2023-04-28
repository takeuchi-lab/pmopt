#pragma once

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include "pmopt/utils/containers.hpp"
#include "pmopt/parameters.hpp"
#include "pmopt/pattern.hpp"

namespace PMOpt
{
    // classes for storing variables to be optimized.
    // there are three classes for storing optimization and mining results.
    //
    // 1. Variables
    //      this is for management of updating primal and dual variables
    //      and also storing intermediate calculation like l2-norm.
    //      in optimization and mining, it is mainly used.
    // 
    // 2. Refernce
    //      this is for storing optimal solution.
    //      it will be used to screen features by being loaded into Variables class.
    //
    // 3. Activeset
    //      this is for mining algorithm and is only depend Pattern type.
    //      if we need to access Pattern, then it is used

    struct Reference;

    template < class Symbol >
    struct Activeset;

    struct Variables
    {

        constexpr static double ZERO_COEF = 1e-16;

        double _intercept = 0.0;
        double _scale = 1.0;
        double _radius = 0.0;
        double _max_score = 0.0;
        double _norm_primal = 0.0;
        double _norm_l1 = 0.0;
        double _norm_l2 = 0.0;
        double _norm_dual = 0.0;
        double _loss_primal = 0.0;
        double _loss_dual = 0.0;
        double _obj_primal = 0.0;
        double _obj_dual = 0.0;
        double _gap = 0.0;
        double _relative_gap = 0.0;

        // following vectors have the size : # of numerical and pattern feature
        // note 
        // - _coef, _score, _reg_weight may contain inactive (screened) feature
        // - _coef_indices contains only active features
        std::vector<double> _coef; // model coefficient for each feature
        std::vector<double> _score; // used in screening, mathematically: |X_j^T dual_raw|
        std::vector<double> _reg_weight; // regularization weight for each feature
        std::vector<unsigned int> _coef_indices; // indices of coef

        std::vector<unsigned int> _sorted_instances;

        // following vectors have the size : # of instances
        std::vector<double> _y_pred;

        // used in reloading coef from mined activeset
        std::vector<unsigned int> _loaded;


        /**
         * @brief load values from reference 
         * @param ref 
         */
        void load(const Reference & ref) noexcept;


        /**
         * @brief reload coef and weight according to mining results
         * @tparam Pattern
         * @param activeset 
         * @param ref 
         */
        template < class Pattern >
        void reload(const Activeset<Pattern> &, const Reference &) noexcept;


        /**
         * @brief store variables into reference
         * @param activeset
         * @return Reference 
         */
        template < class Pattern >
        Reference dump(const Activeset<Pattern> &) const noexcept;


        /**
         * @brief store variables into reference 
         * this is used in the case variables converged before activeset mining
         * @return Reference 
         */
        Reference dump(const Reference &) const noexcept;


        /**
         * @brief to string
         * @param int 
         * @return std::string 
         */
        std::string str(unsigned int) const noexcept;
    };


    struct Reference
    {
        double _max_score;
        double _intercept;
        std::vector<double> _coef;
        std::vector<double> _y_pred;
        std::vector<double> _reg_weight;
        std::map<std::string, unsigned int> _key2index;
        std::vector<std::shared_ptr<const std::vector<unsigned int>>> _x_indices;
        std::vector<std::shared_ptr<const std::vector<double>>> _x_values;

        /**
         * @brief find indices from vector of reference
         * @param[out] outer index
         * @param[out] inner index
         * @param references
         * @param key
         * @return find or not
         */
        static bool find(
            unsigned int &,
            unsigned int &,
            std::vector<std::reference_wrapper<const Reference>> &,
            const std::string &
        ) noexcept;
    };


    template < class Symbol >
    struct Activeset
    {
        std::vector<std::string> _keys;
        std::vector<std::shared_ptr<const Pattern<Symbol>>> _patterns;
        std::vector<std::shared_ptr<const std::vector<unsigned int>>> _x_indices;
        std::vector<std::shared_ptr<const std::vector<double>>> _x_values;
        std::vector<double> _x_sum;
        std::vector<double> _x_sum_sq;
        std::vector<double> _reg_weight;


        /**
         * @brief clear all components
         */
        void clear() noexcept
        {
            _keys.clear();
            _patterns.clear();
            _x_indices.clear();
            _x_values.clear();
            _x_sum.clear();
            _x_sum_sq.clear();
            _reg_weight.clear();
        }
    };


}



template < class Pattern >
void PMOpt::Variables::reload(
    const Activeset<Pattern> & activeset,
    const Reference & ref
)
noexcept
{
    auto size = activeset._keys.size();

    // save flag whether coef is loaded or not
    // they will be used in removing coef not loaded
    // to make _coef and _y_pred consistent  
    _loaded.clear();
    _loaded.resize(ref._key2index.size(), 0);


    // some initialization
    _coef_indices.resize(size);
    _coef.resize(size);
    _score.resize(size);
    _reg_weight.resize(size);


    // load coefficient of reference solution
    for (unsigned int j = 0; j < size; ++j)
    {
        // j : index of this variables class
        // index: index of ref from which this class loads

        const auto & key = activeset._keys[j];

        // if there exists reference coef
        // ref._key2index[key] cannot be used because ref is const
        if (ref._key2index.find(key) != ref._key2index.end())
        {
            auto index = ref._key2index.at(key);
            _coef[j] = ref._coef[index];
            _loaded[index] = true;
        }
        else
        {
            _coef[j] = 0.0;
        }

        // initialize some other parameters
        _score[j] = 0.0;
        _reg_weight[j] = activeset._reg_weight[j];
        _coef_indices[j] = j;
    }

    // remove effect of not loaded coef from dual_raw
    for (unsigned int index = 0; index < ref._key2index.size(); ++index)
    {
        if (_loaded[index])
            continue;

        // if not loaded, then remove effect of coef from dual_raw
        // this indices might contain the index not in training set
        // but it is no problem because _y_pred of test indices are never accessed
        // and we assume no new sample will be add from reference
        const auto & x_i = * ref._x_indices[index];
        const auto & x_v = * ref._x_values[index];
        for (unsigned int i = 0; i < x_i.size(); ++i)
            _y_pred[x_i[i]] -= x_v[i] * ref._coef[index];
    }
}


template < class Pattern >
PMOpt::Reference PMOpt::Variables::dump(const Activeset<Pattern> & activeset) const noexcept
{
    Reference ref;
    ref._max_score = _max_score;
    ref._intercept = _intercept;

    // copy
    ref._y_pred = _y_pred;

    // need for loop on _coef_indices
    // because inactive features are included in _coef and _reg_weight
    // ref._key2index.reserve(_coef_indices.size());
    ref._coef.reserve(_coef_indices.size());
    ref._reg_weight.reserve(_coef_indices.size());
    ref._x_indices.reserve(activeset._x_indices.size());
    ref._x_values.reserve(activeset._x_values.size());
    for (auto j : _coef_indices)
    {
        if (std::abs(_coef[j]) < ZERO_COEF)
            continue;

        ref._key2index[activeset._keys[j]] = ref._key2index.size();
        ref._coef.push_back(_coef[j]);
        ref._reg_weight.push_back(_reg_weight[j]);
        ref._x_indices.push_back(activeset._x_indices[j]);
        ref._x_values.push_back(activeset._x_values[j]);
    }

    return ref;
}
