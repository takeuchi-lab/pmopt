#pragma once

#include <chrono>
#include <stack>
#include <vector>
#include "pmopt/logger.hpp"
#include "pmopt/pattern.hpp"

namespace PMOpt
{

    template < class Enumerator >
    struct Miner
    {
        using Symbol = typename Enumerator::Symbol;
        using Positions = typename Enumerator::Positions;
        using PDB = ProjectedDatabase<Symbol, Positions>;

        const static unsigned int PRIORITY = 5;
        const std::string LOG_PREFIX_NEWLINE = "\n  ";
        const std::string LOG_PREFIX_OUTER = "  - ";
        const std::string LOG_PREFIX_INNER = "  -- ";

        Enumerator & _enumerator;
        unsigned int _maxpat = 3;
        unsigned int _n_children;
        double _time;
        std::vector<PDB> _children;
        std::vector<std::shared_ptr<const PDB>> _history;
        std::stack<std::shared_ptr<const PDB>> _next;

        /**
         * @brief Construct a new Miner object
         * @param enumerator 
         * @param maxpat maximum length of pattern to be mined
         */
        Miner(Enumerator & enumerator, unsigned int maxpat)
            : _enumerator(enumerator), _maxpat(maxpat) {}

        /**
         * @brief Mine patterns with pruner's pruning condition from the database
         * @tparam Pruner 
         * @param mined container for storing mined patterns and related properties
         * @param pruner 
         */
        template < class Pruner, class MinedObject >
        void mine(
            MinedObject & mined,
            const Pruner & pruner
        );
    };
}



template < class Enumerator >
template < class Pruner, class MinedObject >
void PMOpt::Miner<Enumerator>::mine(
    MinedObject & mined,
    const Pruner & pruner
)
{
    using namespace std::chrono;
    using seconds = std::chrono::duration<double>;
    auto start = system_clock::now();
    Logger logger(std::cout);

    // initialize dfs stack
    _history.resize(_maxpat + 1);
    _next.push(std::make_shared<const PDB>(Enumerator::ROOT_SYMBOL));


    // DFS
    while (!_next.empty())
    {
        // pop from dfs stack
        auto parent = _next.top();
        _history[parent->_length] = parent;
        _next.pop();

        logger(PRIORITY, LOG_PREFIX_NEWLINE, "depth:", parent->_length);
        logger(PRIORITY, LOG_PREFIX_OUTER, "parent:", * parent.get());

        // enumerating children
        _n_children = _enumerator.enumerate(_children, * parent.get());

        logger(PRIORITY+1, LOG_PREFIX_OUTER, "# of children:", _n_children);

        // Check pruning condition for each child
        for (unsigned int i_child = 0; i_child < _n_children; ++i_child)
        {

            if (!_children[i_child]._is_valid)
                continue;


            // evaluate prunign condition
            if (pruner.evaluate(mined, _children[i_child]))
                continue;

            logger(PRIORITY+2, LOG_PREFIX_INNER, _children[i_child], " is not pruned");


            // if under maxpat then add stack
            if (_children[i_child]._length < _maxpat)
                _next.push(std::make_shared<const PDB>(_children[i_child]));
        }
    }

    auto end = system_clock::now();
    _time = duration_cast<seconds>(end - start).count();
}