#pragma once

#include <memory>

namespace PMOpt
{

    /**
     * @brief Pattern itself
     * @tparam Symbol : the type of minimum element of pattern.
     * E.g.
     *  Itemset pattern and sequence pattern is implemented as Pattern<unsigned int>.
     *  Graph pattern is implemented as Pattern<Edge>
     */
    template < class Symbol >
    struct Pattern
    {
        Symbol _symbol;

        // pointer for parent pattern
        std::shared_ptr<const Pattern<Symbol>> _parent;
    };


    /**
     * @brief Information related to the pattern.
     *  E.g. occurrence of pattern or position that pattern occurred in.
     * 
     * @tparam Symbol : the type of minimum element of pattern
     * @tparam Position : the type represents how pattern occurs in samples
     */
    template < class Symbol, class Positions >
    struct ProjectedDatabase
    {
        std::shared_ptr<const Pattern<Symbol>> _pattern;

        unsigned int _length = 0;
        bool _is_valid = false;
        std::vector<unsigned int> _occurrence;
        std::vector<Positions> _positions;

        ProjectedDatabase() = default;

        ProjectedDatabase(Symbol root) : _pattern(new Pattern<Symbol>({
            root,
            std::shared_ptr<const Pattern<Symbol>>()
        })) {}
    };


    /**
     * @brief 
     * @note This function should be used only for debugging.
     *  To obtain exact patttern string representation, use Enumerator::pat_to_str instead.
     * @tparam Symbol 
     * @param os 
     * @param pattern 
     * @return std::ostream& 
     */
    template < class Symbol >
    inline std::ostream & operator <<(
        std::ostream & os,
        const Pattern<Symbol> & pattern
    )
    noexcept
    {
        const auto * p = & pattern;
        bool is_first = true;
        os << '{';
        while (p != nullptr)
        {
            if (!is_first) os << ',';
            is_first = false;

            os << p->_symbol;
            p = p->_parent.get();
        }
        os << '}';
        return os;
    }


    template < class Symbol, class Positions >
    inline std::ostream & operator <<(
        std::ostream & os,
        const ProjectedDatabase<Symbol, Positions> & pdb
    )
    {
        return os << * pdb._pattern.get();
    }

}