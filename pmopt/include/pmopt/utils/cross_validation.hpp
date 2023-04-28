#pragma once

#include <cassert>
#include <vector>

namespace PMOpt
{
    
    namespace Utils
    {
        /**
         * @brief Make leave-group-out cv indices
         * @param[out] indices
         * @param group
         * @return unsigned int
         */
        unsigned int make_leave_group_out(
            std::vector<unsigned int> &,
            const std::vector<unsigned int> &
        ) noexcept;
        

        /**
         * @brief Make leave-one-out cv indices
         * @param[out] indices
         * @param size
         * @param seed
         * @param shuffle
         * @return unsigned int
         */
        unsigned int make_leave_one_out(
            std::vector<unsigned int> &,
            unsigned int,
            unsigned int,
            bool
        ) noexcept;


        /**
         * @brief Make K-fold cv indices
         * @param[out] indices
         * @param size 
         * @param n_folds 
         * @param seed 
         * @return unsigned int
         */
        unsigned int make_k_fold(
            std::vector<unsigned int> &,
            unsigned int, unsigned int, unsigned int) noexcept;



        /**
         * @brief Make indices
         * @param[out] mask
         * @param[out] n_training
         * @param fold_ids
         * @param fold_id
         */
        void make_test_mask(
            std::vector<unsigned int> &,
            unsigned int &,
            const std::vector<unsigned int> &,
            unsigned int
        ) noexcept;
    }
}