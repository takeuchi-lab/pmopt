#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include "pmopt/structures/sequence.hpp"

namespace PMOpt
{

    struct SequenceDatabase : public std::vector<Sequence>
    {
        static const unsigned int ROOT = 0;

        const std::string FORMAT_STRING = "string";
        const std::string FORMAT_SEQUENCE = "sequence";

        std::string _format;
        std::vector<std::string> _items;
        std::unordered_map<std::string, unsigned int> _item2id;


        /**
         * @brief Construct a new Sequence Database object
         * @param format name of column that stores structured instance in input csv
         */
        SequenceDatabase(const std::string & format) : _format(format)
        {
            // In fact, root pattern has no item,
            // but in the implementation, it has an ROOT item, same as itemset mining
            _items.push_back(std::string());
            _item2id[std::string()] = ROOT;
        }


        /**
         * @brief Push database to new instance that represented by string
         * @param text new instance
         */
        void push(std::string &&);


        /**
         * @brief Push instance (string ver.)
         */
        void push_string(std::string &&);


        /**
         * @brief push instance (sequence ver.)
         */
        void push_sequence(std::string &&);


        /**
         * @brief convert structure instance into string
         * @param structure
         * @return std::string 
         */
        std::string struct_to_str(const Sequence &) const noexcept;


        /**
         * @brief convert structure instance into string (string ver.)
         * @param structure
         * @return std::string 
         */
        std::string struct_to_str_string(const Sequence &) const noexcept;


        /**
         * @brief convert structure instance into string (sequence ver.)
         * @param structure
         * @return std::string 
         */
        std::string struct_to_str_sequence(const Sequence &) const noexcept;

    };
}