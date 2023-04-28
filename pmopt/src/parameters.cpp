#include "pmopt/parameters.hpp"

#include <algorithm>
#include <iostream>
#include <unordered_set>
#include "pmopt/optimization/basic.hpp"
#include "pmopt/utils/containers.hpp"

using namespace PMOpt;

void Parameters::initialize(unsigned int n_instances, bool is_numeric)
{
    _n_instances = n_instances;
    assert(!_instance_indices);

    // Initialize following members of this class
    // _instance_indices
    // _ones
    // _w
    // _y

    // _instance_indices
    // This stores all instance indices: {0, 1, ..., n-1}
    std::vector<unsigned int> indices(n_instances);
    std::iota(indices.begin(), indices.end(), 0);
    _instance_indices = std::make_shared<std::vector<unsigned int>>(std::move(indices));

    // _ones
    // This is a vector of ones: {1, ..., 1}
    _ones = std::make_shared<const std::vector<double>>(n_instances, 1.0);

    // _w
    // This is a vector of weights for each instance
    // It is loaded from the input file
    // If not loaded, we set ones: {1, ..., 1}
    if (_w.empty())
        _w.resize(n_instances, 1.0);

    assert(_w.size() == n_instances);

    // _y
    // This is a vector of objective variables
    // Here we check variety of _y
    std::unordered_set<double> y_set;
    for (auto y : _y)
        y_set.insert(y);
    if (y_set.size() == 1)
    {
        std::cerr << "Error: All y values are the same.\n";
        exit(1);
    }

    assert(_y.empty() || _y.size() == n_instances);

    // If the loss function expects non numeric objective
    // _y must be +1 or -1
    if (!is_numeric)
        for (auto y : _y)
            if (y != +1 && y != -1)
            {
                std::cerr << "Error: y must be +1/-1 if the loss function is for classification.\n";
                exit(1);
            }

    if (_X_num.empty())
        return;

    // _X_mean and _X_std
    // These are vectors of mean and standard deviation
    _X_mean.resize(_X_num.size(), 0.0);
    _X_std.resize(_X_num.size(), 1.0);
    assert(std::all_of(_X_num.begin(), _X_num.end(), [&](auto _X_j){
        return _X_j->size() == n_instances;
    }));
}


void Parameters::standardize(bool is_numeric)
{
    // standardize y
    if (is_numeric)
        try
        {
            Optimization::standardize(_y, _y_mean, _y_std, _mask);
        }
        catch(const std::runtime_error & e)
        {
            std::cerr << e.what() << " in standardization of y."
                << " Please make sure that objective `y` have large enough variance." << '\n';
            exit(1);
        }
    else
    {
        // Check both classes are in the training set
        check_classes();   
    }

    // standardize w
    try
    {
        Optimization::scale(_w, _w_sum, _w_max, _mask);
    }
    catch(const std::runtime_error & e)
    {
        std::cerr << e.what() << " in scaling of w."
            << " Please make sure that weight `w` have large enough variance." << '\n';
        exit(1);
    }

    // standardize x
    for (unsigned int j = 0; j < _X_num.size(); ++j)
    {
        try
        {
            Optimization::standardize(* _X_num[j], _X_mean[j], _X_std[j], _mask);
        }
        catch(const std::runtime_error & e)
        {
            std::cerr << e.what() << " in standardization of X_" << j+1
                << ". Please make sure that j-th numeric feature have large enough variance." << '\n';
            exit(1);
        }
    }   
        
}


void Parameters::check_classes() const
{
    // Error if y has a value neither +1 nor -1
    // Error if y does not have both +1 and -1
    bool has_negative = false;
    bool has_positive = false;
    for (unsigned int i_instance = 0; i_instance < _n_instances; ++i_instance)   
    {
        // exclude testing instances
        if (!_mask[i_instance])
            continue;
        
        if (_y[i_instance] == -1)
            has_negative = true;
        else if (_y[i_instance] == +1)
            has_positive = true;
        else
        {
            std::cerr << "Error: y must be +1/-1 for classification loss function.\n";
            exit(1);
        }

        if (has_negative && has_positive)
            return;
    }
    
    std::cerr << "Error: y must have both +1 and -1 for classification loss function.\n";
    exit(1);
}