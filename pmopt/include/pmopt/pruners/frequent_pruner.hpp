#pragma once

#include <memory>
#include <vector>
#include "pmopt/pattern.hpp"

namespace PMOpt
{

    template < class Symbol >
    struct FrequentObject
    {
        std::vector<std::shared_ptr<const Pattern<Symbol>>> _patterns;
        std::vector<unsigned int> _supports;
    };


    /**
     * @brief pruning condition used in frequent mining
     */
    struct FrequentPruner
    {
        unsigned int _minsup = 0;

        /**
         * @brief Construct a new Frequent Pruner object
         * @param minsup : minmimum support of pattern, which can be considered as frequent
         */
        FrequentPruner(unsigned int minsup) : _minsup(minsup) {}


        /**
         * @brief evaluate support of a pattern and check frequent or not
         * @tparam ProjectedDatabase
         * @param[out] mined support and pattern will be saved in it if frequent
         * @param projetced stores which instance the pattern is in
         * @return bool : true if pattern is not frequent
         */
        template < class Symbol, class Positions >
        bool evaluate(
            FrequentObject<Symbol> & mined,
            const ProjectedDatabase<Symbol, Positions> & projected
        )
        const noexcept
        {
            unsigned int support = projected._occurrence.size();
            bool is_pruned = support < _minsup;
            if (!is_pruned)
            {
                mined._supports.push_back(support);
                mined._patterns.push_back(projected._pattern);
            }
            return is_pruned;
        }
    };
}