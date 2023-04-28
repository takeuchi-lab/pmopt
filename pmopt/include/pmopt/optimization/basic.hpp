#pragma once

#include <vector>

namespace PMOpt
{
    namespace Optimization
    {

        /**
         * @brief standardize vector to make mean 0, stddev 1
         * @param[out] vec
         * @param[out] mean
         * @param[out] stddev
         * @param mask
         */
        void standardize(
            std::vector<double> &,
            double &,
            double &,
            const std::vector<unsigned int> &
        );
 

        /**
         * @brief scale vector to make sum equals to one
         * @param[out] vec
         * @param[out] sum
         * @param[out] max
         * @param mask
         */
        void scale(
            std::vector<double> &,
            double &,
            double &,
            const std::vector<unsigned int> &
        );


        /**
         * @brief compute sum and square of vec
         * expect vec to be sparse (i.e. vec have no zero entries)
         * @param[out] sum
         * @param[out] sum_square
         * @param vec
         * @param indices
         * @param mask
         */
        void sum_square(
            double &,
            double &,
            const std::vector<double> &,
            const std::vector<unsigned int> &,
            const std::vector<unsigned int> &
        ) noexcept;

    }
}