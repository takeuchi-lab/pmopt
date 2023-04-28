#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

namespace PMOpt
{

    struct Parameters
    {
        std::shared_ptr<const std::vector<unsigned int>> _instance_indices;
        std::shared_ptr<const std::vector<double>> _ones;
        std::vector<unsigned int> _mask;
        std::vector<double> _y;
        std::vector<double> _w;
        std::vector<std::shared_ptr<std::vector<double>>> _X_num;
        std::vector<std::string> _X_num_names;
        std::vector<unsigned int> _group;
        std::vector<double> _X_mean;
        std::vector<double> _X_std;

        unsigned int _n_instances = 0;
        unsigned int _n_training = 0;
        double _y_mean = 0.0;
        double _y_std = 1.0;
        double _w_sum = 1.0;
        double _w_max = 1.0;
        double _lambda_l1 = 0.0;
        double _lambda_l2 = 0.0;
        double _reg_power = 0.0;

        // whether dual variables made to be scaled in lambda_l2 != 0
        bool _l2_scale = false;

        /**
         * @brief set reference of loaded vectors into eigen map
         * must be called at once
         * @param n_instances
         * @param is_numeric
         */
        void initialize(unsigned int, bool);


        /**
         * @brief standardize _y, _w, ... according to _instance_indices
         * @param is_numeric if standardization of y is needed
         */
        void standardize(bool);


        /**
         * @brief check classes of objective vector for classification case.
         */
        void check_classes() const;
    };

}