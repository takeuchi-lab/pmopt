#pragma once

#include <fstream>
#include "pmopt/pruners/frequent_pruner.hpp"
#include "pmopt/parameters.hpp"
#include "pmopt/variables.hpp"

namespace PMOpt
{

    namespace Utils
    {

        /**
         * @brief write mined frequent patterns with its support
         * @tparam Pattern 
         * @param mined
         * @param enumerator
         * @param filename 
         */
        template < class Enumerator >
        void write_frequent(
            FrequentObject<typename Enumerator::Symbol> & mined,
            const Enumerator & enumerator,
            const std::string & filename
        )
        {
            std::ofstream ofs(filename);
            if (!ofs)
            {
                std::cerr << "Error: File cannot open `" << filename << "`.\n";
                exit(1);
            }

            ofs << "pattern,support\n";
            for (unsigned int i = 0; i < mined._patterns.size(); ++i)
            {
                const auto & database = enumerator._database;
                const auto & pattern = * mined._patterns[i].get();
                ofs << database.struct_to_str(enumerator.pat_to_struct(pattern));
                ofs << ",";
                ofs << mined._supports[i];
                ofs << "\n";
            }
        }


        /**
         * @brief write model
         * @param reference
         * @param params
         * @param filename
         * @param trial
         * @param nol2
         */
        void write_model(
            const Reference &,
            const Parameters &,
            const std::string &,
            unsigned int,
            bool
        );
    }

}