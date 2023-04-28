#include "pmopt/databases/itemset_database.hpp"

#include <iostream>
#include <cassert>
#include <sstream>
#include <algorithm>

using namespace PMOpt;


void ItemsetDatabase::push(std::string && buffer)
{
    Itemset itemset;

    std::stringstream ss(std::move(buffer));
    while (std::getline(ss, buffer, ' '))
    {
        if (_item2id.find(buffer) == _item2id.end())
        {
            // found new item
            _item2id[buffer] = _items.size();
            _items.push_back(buffer);
        }
        itemset.push_back(_item2id[buffer]);
    }

    // itemset must be sorted
    std::sort(itemset.begin(), itemset.end());

    emplace_back(std::move(itemset));
}


std::string ItemsetDatabase::struct_to_str(const Itemset & structure) const noexcept
{
    std::stringstream ss;

    bool is_first = true;
    for (auto item : structure)
    {
        assert(item != ROOT);

        if (!is_first)
            ss << " ";

        ss << _items[item];
        is_first = false;
    }

    return ss.str();
}