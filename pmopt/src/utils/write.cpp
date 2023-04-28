#include "pmopt/utils/write.hpp"

#include <iomanip>
#include <unordered_map>

using namespace PMOpt;

const unsigned int PRECISION = 16;

void Utils::write_model(
    const Reference & reference,
    const Parameters & params,
    const std::string & filename,
    unsigned int trial,
    bool nol2
)
{
    std::ofstream ofs(filename, (trial == 0) ? std::ios_base::out : std::ios_base::app);
    if (!ofs)
    {
        std::cerr << "Error: File cannot open: " << filename << "\n";
        exit(1);
    }

    if (trial == 0)
    {
        ofs << "id,lambda_l1";
        if (!nol2)
            ofs << ",lambda_l2";
        ofs << ",coefficient,pattern\n";
    }

    // Here coefficients and an intercept of the model will be written.
    // They are abtained from scaled dataset,
    // so it must be unscaled to be consistent with the original dataset.
    // Let m_y, s_y be the mean / standard deviation of y and y' be standardized one such that:
    //      y'_i = (y_i - m_y) / s_y
    // And let m_k, s_k be the mean / standard deviation of x_k, where k-th numeric input feature
    // and x'_k be standardized one like
    //      x'_ik = (x_ik - m_k) / s_k
    // Also we denote x_j as j-th pattern input feature.
    // Then we can write the model as:
    //      y'_i = \sum_k x'_ik w_k + \sum_j x_ij w_j + b
    // So,
    //      (y_i - m_y) / s_y = \sum_j x_ij w_j +  \sum_k (x_ik - m_k) / s_k w_k + b
    //      y_i = \sum_j x_ij s_y w_j + \sum_k (x_ik - m_k) / s_k s_y w_k + b s_y + m_y
    //      y_i = \sum_j x_ij s_y w_j + \sum_k x_ik / s_k s_y w_k + b s_y + m_y - \sum_k m_k / s_k s_y w_k
    // Here, let w'_j, w'_k, b' be
    //      w'_j = s_y w_j
    //      w'_k = s_y / s_k w_j
    //      b' = b s_y + m_y - \sum_k m_k / s_k s_y w_k
    // Now the model for the original dataset can be written as
    //      y_i = \sum_j x_ij w'_j + \sum_k x_ik w'_k + b'

    std::unordered_map<std::string, unsigned int> name2index;
    for (unsigned int i = 0; i < params._X_num_names.size(); ++i)
        name2index[params._X_num_names[i]] = i;

    ofs << std::fixed;
    ofs << std::setprecision(PRECISION);

    // write coefficients
    double intercept_of_intercept = 0.0;
    for (const auto & pair : reference._key2index)
    {
        ofs << trial << ',' << params._lambda_l1;
        if (!nol2)
            ofs << ',' << params._lambda_l2;

        const auto & key = pair.first;
        auto coef = reference._coef[pair.second];
        if (name2index.find(key) != name2index.end()) // if numeric feature
        {
            auto index = name2index[key];
            ofs << ',' << coef * params._y_std / params._X_std[index];

            // accumulate standardization effect of intercept
            intercept_of_intercept += coef * params._X_mean[index] / params._X_std[index];
        }
        else // if pattern feature
            ofs << ',' << coef * params._y_std;
        ofs << ',' << key << '\n';
    }

    // write an intercept
    ofs << trial << ',' << params._lambda_l1;
    if (!nol2)
        ofs << ',' << params._lambda_l2;
    auto intercept = reference._intercept;
    ofs << ',' << (intercept - intercept_of_intercept) * params._y_std + params._y_mean << ",\n";

}