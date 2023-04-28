#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include "pmopt/structures/itemset.hpp"

namespace PMOpt
{
    
    struct ItemsetDatabase : public std::vector<Itemset>
    {

        static const unsigned int ROOT = 0;

        std::string _format;
        std::vector<std::string> _items;
        std::unordered_map<std::string, unsigned int> _item2id;


        /**
         * @brief Construct a new Itemset Database object
         * @param format name of column that stores structured instance in input csv
         */
        ItemsetDatabase(const std::string & format) : _format(format)
        {
            // item for root
            // In fact, root pattern has no item,
            // but in the implementation, it has an ROOT item
            _items.push_back(std::string());
            _item2id[std::string()] = ROOT;
        }


        /**
         * @brief Push database to new instance that represented by string
         * @param text new instance
         */
        void push(std::string &&);


        /**
         * @brief convert structure instance into string
         * @param structure
         * @return std::string 
         */
        std::string struct_to_str(const Itemset &) const noexcept;

    };
}
