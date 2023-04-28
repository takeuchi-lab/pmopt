#include "pmopt/optimization/basic.hpp"

#include <cmath>
#include <iostream>

using namespace PMOpt;

namespace
{
    const double DIV_TOL = 1e-8;
}

void Optimization::standardize(
    std::vector<double> & vec,
    double & mean,
    double & stddev,
    const std::vector<unsigned int> & mask
)
{
    // Standardization multiple times
    // Let v_1, v_2, ..., v_n be n-dimensional vector
    // Here we consider some standardized vector v' = v'_1, v'_2, ..., v'_n such that:
    //      v'_i = (v_i - m) / s
    // where
    //      m: the mean of v'_i
    //      s: the standard deviation of v'_i
    //
    // and then, we can restore original v_i by calculating
    //      v_i = v'_i s + m
    //
    // Let m', s' be the mean and standard deviation of v'_i respectively.
    // We can standardize v'_i to obtain v''_i
    //      v''_i = (v'_i - m') / s'
    //            = ((v_i - m) / s - m')/s'
    //            = (v_i - m - sm')/(ss')
    // In this time, we can restore original v_i by calculating
    //      v_i = v' ss' + m + sm'
    // and we can consider m + sm' as new mean and ss' as new standard deviation,
    // so we update the mean and the standard deviation:
    //      m' <- m + sm'
    //      s' <- ss'

    unsigned int n = 0;
    double current_mean = 0.0;
    double current_stddev = 0.0;
    for (unsigned int i = 0; i < mask.size(); ++i)
    {
        current_mean += vec[i] * mask[i];
        current_stddev += vec[i] * vec[i] * mask[i];
        n += mask[i];
    }
    current_mean /= n;
    current_stddev /= n;

    // Check variance of y
    current_stddev = current_stddev - current_mean * current_mean;
    if (std::abs(current_stddev) < DIV_TOL)
        std::cerr << "Warning: Variance of y vector is too small.\n";
    if (current_stddev != current_stddev)
        throw std::runtime_error("NaN detected");

    current_stddev = std::max(current_stddev, DIV_TOL);
    current_stddev = std::sqrt(current_stddev);
    for (auto & v : vec)
        v = (v - current_mean) / current_stddev;

    mean += stddev * current_mean;
    stddev *= current_stddev;
}


void Optimization::scale(
    std::vector<double> & vec,
    double & sum,
    double & max,
    const std::vector<unsigned int> & mask
)
{
    // Scale multiple times
    // Let v_i be n-dimensional vector and v'_i be scaled by m.
    //      v'_i = v_i / m
    // Then we can restore original v_i by calculating
    //      v_i = v'_i m
    // 
    // Here, we consider further scaling of v'_i:
    //      v''_i = v'_i / m'
    //            = v_i / (m m')
    // Then we can restore original v_i by calculating
    //      v_i = v''_i m m'
    // Thus mm' can be considered as new scaling factor

    double current_sum = 0.0;
    max = 0.0;
    for (auto & v : vec)
        v = std::max(v, DIV_TOL);

    for (unsigned int i = 0; i < mask.size(); ++i)
    {
        current_sum += vec[i] * mask[i];
        max = std::max(vec[i] * mask[i], max);
    }

    for (auto & v : vec)
        v /= current_sum;

    if (current_sum != current_sum)
        throw std::runtime_error("NaN detected");

    max /= current_sum;
    sum *= current_sum;
}



void Optimization::sum_square(
    double & sum,
    double & sum_square,
    const std::vector<double> & vec,
    const std::vector<unsigned int> & indices,
    const std::vector<unsigned int> & mask
)
noexcept
{
    sum = 0.0;
    sum_square = 0.0;
    for (unsigned int i_occ = 0; i_occ < indices.size(); ++i_occ)
    {
        auto m = mask[indices[i_occ]];
        sum += vec[i_occ] * m;
        sum_square += vec[i_occ] * vec[i_occ] * m;
    }
}