#include "pmopt/variables.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>

using namespace PMOpt;

void Variables::load(const Reference & ref) noexcept
{
    _intercept = ref._intercept;
    _scale = 1.0;
    _radius = 0.0;
    _max_score = 0.0;
    _norm_primal = 0.0;
    _norm_l1 = 0.0;
    _norm_l2 = 0.0;
    _norm_dual = 0.0;
    _loss_primal = 0.0;
    _loss_dual = 0.0;
    _obj_primal = 0.0;
    _obj_dual = 0.0;
    _gap = 0.0;
    _relative_gap = 0.0;

    _score.clear();
    _coef.resize(ref._coef.size());
    _y_pred.resize(ref._y_pred.size());
    _reg_weight.resize(ref._reg_weight.size());
    _coef_indices.resize(ref._coef.size());
    std::iota(_coef_indices.begin(), _coef_indices.end(), 0); 
    std::copy(ref._coef.begin(), ref._coef.end(), _coef.begin());
    std::copy(ref._y_pred.begin(), ref._y_pred.end(), _y_pred.begin());
    std::copy(ref._reg_weight.begin(), ref._reg_weight.end(), _reg_weight.begin());
}


Reference Variables::dump(const Reference & ref) const noexcept
{
    Reference ref_new;
    ref_new._max_score = _max_score;
    ref_new._intercept = _intercept;
    
    // copy
    ref_new._y_pred =  _y_pred;
    ref_new._key2index = ref._key2index;
    ref_new._coef = ref._coef;
    ref_new._reg_weight = ref._reg_weight;
    ref_new._x_indices = ref._x_indices;
    ref_new._x_values = ref._x_values;
    return ref_new;
}


std::string Variables::str(unsigned int precision) const noexcept
{
    const unsigned int WIDTH = 13;
    std::stringstream ss;
    ss << std::scientific << std::setprecision((int) precision)
        << "\n\tP : " << std::setw(WIDTH) << _obj_primal
        << "  PL: " << std::setw(WIDTH) << _loss_primal
        << "  PR: " << std::setw(WIDTH) << _norm_primal
        << "\n\tD : " << std::setw(WIDTH) << _obj_dual
        << "  DL: " << std::setw(WIDTH) << _loss_dual
        << "  DR: " << std::setw(WIDTH) << _norm_dual
        << "\n\tG : " << std::setw(WIDTH) << _gap
        << "  L1: " << std::setw(WIDTH) << _norm_l1
        << "  L2: " << std::setw(WIDTH) << _norm_l2
        << "\n\tR : " << std::setw(WIDTH) << _radius
        << "  S : " << std::setw(WIDTH) << _scale
        << "  M : " << std::setw(WIDTH) << _max_score
        << "\n";
    return ss.str();
}


bool Reference::find(
    unsigned int & outer,
    unsigned int & inner,
    std::vector<std::reference_wrapper<const Reference>> & references,
    const std::string & key
)
noexcept
{
    for (outer = 0; outer < references.size(); ++outer)
    {
        const Reference & ref = references[outer];
        if (ref._key2index.find(key) != ref._key2index.end())
        {
            inner = ref._key2index.at(key);
            return true;
        }
    }
    return false;
}
