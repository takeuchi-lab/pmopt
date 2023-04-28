#pragma once

#include <iostream>
#include <memory>
#include "pmopt/databases/sequence_database.hpp"
#include "pmopt/logger.hpp"
#include "pmopt/pattern.hpp"

namespace PMOpt
{

    struct SequenceEnumerator
    {
        using Database = SequenceDatabase;
        using Symbol = unsigned int;
        using Positions = std::vector<unsigned int>;
        using PDB = ProjectedDatabase<Symbol, Positions>;

        static constexpr unsigned int ROOT_SYMBOL = 0;

        static unsigned int _max_skip;

        const Database & _database;

        std::shared_ptr<const PDB> _root;

        /**
         * @brief Construct a new Sequence Enumerator object
         * @param database 
         * @param logger 
         * @param option
         */
        SequenceEnumerator(const Database & database) : _database(database) {}


        /**
         * @brief enumerate child patterns of the parent
         * @param[out] children std::vector<PDB> &
         * @param parent
         * @return unsigned int 
         */
        unsigned int enumerate(std::vector<PDB> &, const PDB &) const noexcept;


        /**
         * @brief search child pattern of parent
         * @param[out] children
         * @param parent
         */
        void search(std::vector<PDB> &, const PDB &) const noexcept;

        
        /**
         * @brief search child pattern of root
         * @param[out] children
         */
        void search_from_root(std::vector<PDB> &) const noexcept;


        /**
         * @brief initialize buffer for projected database
         * @param[out] children
         * @param parent
         * @param n_items
         * @param n_instances
         */
        void initialize_children(
            std::vector<PDB> &,
            const PDB &,
            unsigned int,
            unsigned int
        ) const noexcept;


        /**
         * @brief convert pattern into structure type
         * @return Structure 
         */
        Sequence pat_to_struct(const Pattern<unsigned int> &) const noexcept;


        /**
         * @brief convert pattern into string
         * @return std::string 
         */
        std::string pat_to_str(const Pattern<unsigned int> &) const noexcept;
    };
}