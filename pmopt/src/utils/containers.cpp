#include "pmopt/utils/containers.hpp"

using namespace PMOpt;

bool Utils::contains(const std::string & key, const std::vector<std::string> & values)
{
    bool is_contained = false;
    for (const auto & val : values)
        is_contained = is_contained || key == val;
    return is_contained;
}