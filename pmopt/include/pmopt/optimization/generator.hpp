#pragma once

#include <cmath>
#include "pmopt/parameters.hpp"

namespace PMOpt
{

    /**
     * @brief class for generating binary vector as feature
     */
    struct BinaryXGenerator
    {
        static const bool IS_BINARY = true;

        /**
         * @brief compute feature vector of pattern
         * @tparam ProjectedDatabase 
         * @param projected 
         * @return Eigen vector type
         */
        template < class ProjectedDatabase >
        static auto generate(const ProjectedDatabase &, const Parameters & params)
        {
            return params._ones;
        }
    };


    /**
     * @brief class for generating weight that is power of pattern length
     */
    struct LengthWGenerator
    {
        static const bool IS_BINARY = false;

        /**
         * @brief generate regularization weight of pattern
         * @tparam ProjectedDatabase 
         * @param projected 
         * @param params
         * @return double 
         */
        template < class ProjectedDatabase >
        static double generate(const ProjectedDatabase & projected, const Parameters & params)
        {
            return std::pow(projected._length, params._reg_power);
        }
    };


}