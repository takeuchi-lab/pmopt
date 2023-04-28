#pragma once

#include <iostream>
#include <memory>
#include "pmopt/databases/itemset_database.hpp"
#include "pmopt/pattern.hpp"

namespace PMOpt
{

    struct ItemsetEnumerator
    {
        using Database = ItemsetDatabase;
        using Symbol = unsigned int;
        using Positions = unsigned int;
        using PDB = ProjectedDatabase<Symbol, Positions>;

        static constexpr unsigned int ROOT_SYMBOL = 0;

        const ItemsetDatabase & _database;

        std::shared_ptr<const PDB> _root;


        /**
         * @brief Construct a new Itemset Enumerator object
         * @param database 
         */
        ItemsetEnumerator(const ItemsetDatabase & database) : _database(database) {}


        /**
         * @brief enumerate child patterns of the parent
         * @param[out] children std::vector<ProjectedDatabase> &
         * @param parent const ProjectedDatabase &
         * @return unsigned int number of children 
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
        Itemset pat_to_struct(const Pattern<unsigned int> &) const noexcept;


        /**
         * @brief convert pattern into string
         * @return std::string 
         */
        std::string pat_to_str(const Pattern<unsigned int> &) const noexcept;
    };
}