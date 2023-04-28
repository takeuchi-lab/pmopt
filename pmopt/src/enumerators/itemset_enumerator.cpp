#include "pmopt/enumerators/itemset_enumerator.hpp"

#include <cassert>
#include <sstream>
#include "pmopt/utils/containers.hpp"

using namespace PMOpt;

unsigned int ItemsetEnumerator::enumerate(
    std::vector<PDB> & children,
    const PDB & parent
)
const noexcept
{
    auto n_items = _database._items.size() - parent._pattern->_symbol - 1;
    auto n_instances = parent._occurrence.size();

    initialize_children(children, parent, n_items, n_instances);

    // If `parent` is ROOT then a parent of `parent` is nullptr
    if (!parent._pattern->_parent)
        search_from_root(children);
    else
        search(children, parent);

    for (auto & child : children)
        child._is_valid = !child._occurrence.empty();

    return n_items;
}


void ItemsetEnumerator::search(
    std::vector<PDB> & children,
    const PDB & parent
)
const noexcept
{
    for (unsigned int i_occ = 0; i_occ < parent._occurrence.size(); ++i_occ)
    {
        auto i_instance = parent._occurrence[i_occ];
        auto pos = parent._positions[i_occ];

        for (; pos < _database[i_instance].size(); ++pos)
        {
            auto item = _database[i_instance][pos] - parent._pattern->_symbol - 1;
            children[item]._occurrence.push_back(i_instance);
            children[item]._positions.push_back(pos + 1);
        }
    }
}


void ItemsetEnumerator::search_from_root(
    std::vector<PDB> & children
)
const noexcept
{
    for (unsigned int i_instance = 0; i_instance < _database.size(); ++i_instance)
    {
        for (unsigned int pos = 0; pos < _database[i_instance].size(); ++pos)
        {
            auto item = _database[i_instance][pos] - 1; // substract root item
            children[item]._occurrence.push_back(i_instance);
            children[item]._positions.push_back(pos + 1);
        }
    }
}

void ItemsetEnumerator::initialize_children(
    std::vector<PDB> & children,
    const PDB & parent,
    unsigned int n_items,
    unsigned int n_instances
)
const noexcept
{
    children.resize(std::max((unsigned int) children.size(), n_items));

    if (n_instances == 0)
        n_instances = _database.size();

    for (unsigned int i_child = 0; i_child < n_items; ++i_child)
    {
        auto & child = children[i_child];
        child._length = parent._length + 1;
        child._occurrence.clear();
        child._occurrence.reserve(n_instances);
        child._positions.clear();
        child._positions.reserve(n_instances);
        child._pattern = std::make_shared<const Pattern<unsigned int>>(
            Pattern<unsigned int>{i_child + parent._pattern->_symbol + 1, parent._pattern}
        );
    }
}


Itemset ItemsetEnumerator::pat_to_struct(const Pattern<unsigned int> & pattern) const noexcept
{
    if (pattern._symbol == _database.ROOT)
        return Itemset();

    Itemset itemset{pattern._symbol};

    auto p = pattern._parent;
    while (p->_parent)
    {
        itemset.push_back(p->_symbol);
        p = p->_parent;
    }
    std::sort(itemset.begin(), itemset.end());
    return itemset;
}


std::string ItemsetEnumerator::pat_to_str(const Pattern<unsigned int> & pattern) const noexcept
{
    return _database.struct_to_str(pat_to_struct(pattern));
}