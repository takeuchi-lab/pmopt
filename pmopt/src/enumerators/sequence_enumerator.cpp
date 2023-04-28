#include "pmopt/enumerators/sequence_enumerator.hpp"

using namespace PMOpt;

unsigned int SequenceEnumerator::enumerate(
    std::vector<PDB> & children,
    const PDB & parent
)
const noexcept
{
    auto n_items = _database._items.size() - 1;
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


void SequenceEnumerator::search(
    std::vector<PDB> & children,
    const PDB & parent
)
const noexcept
{
    (void) children;
    for (unsigned int i_occ = 0; i_occ < parent._occurrence.size(); ++i_occ)
    {
        auto i_instance = parent._occurrence[i_occ];
        const auto & pos_list = parent._positions[i_occ];

        for (auto pos : pos_list) // for each starting positions
        {
            for (unsigned int skip = 0; skip <= _max_skip; ++skip)
            {
                // if position reaches the end of the sequence
                if (pos + skip >= _database[i_instance].size())
                    continue;

                auto item = _database[i_instance][pos + skip] - 1;

                if (children[item]._occurrence.empty() || children[item]._occurrence.back() != i_instance)
                {
                    // item is found in the indexed instance first time
                    children[item]._occurrence.push_back(i_instance);
                    children[item]._positions.push_back(std::vector<unsigned int>{pos + skip + 1});
                }
                else
                {
                    // item is already found in the indexed instance
                    children[item]._positions.back().push_back(pos + skip + 1);
                }
            }
        }
    }
}


void SequenceEnumerator::search_from_root(
    std::vector<PDB> & children
)
const noexcept
{
    for (unsigned int i_instance = 0; i_instance < _database.size(); ++i_instance)
    {
        for (unsigned int pos = 0; pos < _database[i_instance].size(); ++pos)
        {
            auto item = _database[i_instance][pos] - 1;
            if (children[item]._occurrence.empty() || children[item]._occurrence.back() != i_instance)
            {
                children[item]._occurrence.push_back(i_instance);
                children[item]._positions.push_back(std::vector<unsigned int>{pos + 1});
            }
            else
            {
                children[item]._positions.back().push_back(pos + 1);
            }
        }
    }
}


void SequenceEnumerator::initialize_children(
    std::vector<PDB> & children,
    const PDB & parent,
    unsigned int n_items,
    unsigned int n_instances
)
const noexcept
{
    children.resize(n_items);

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
            Pattern<unsigned int>{i_child + 1, parent._pattern}
        );
    }
}


Sequence SequenceEnumerator::pat_to_struct(const Pattern<unsigned int> & pattern) const noexcept
{
    if (pattern._symbol == _database.ROOT)
        return Sequence();

    Sequence sequence{pattern._symbol};

    const auto * p = pattern._parent.get();
    while (p->_parent)
    {
        sequence.push_back(p->_symbol);
        p = p->_parent.get();
    }
    std::reverse(sequence.begin(), sequence.end());
    return sequence;
}


std::string SequenceEnumerator::pat_to_str(const Pattern<unsigned int> & pattern) const noexcept
{
    return _database.struct_to_str(pat_to_struct(pattern));
}