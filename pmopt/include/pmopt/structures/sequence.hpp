#pragma once

#include <vector>
#include "pmopt/utils/containers.hpp"

namespace PMOpt
{

    struct Sequence : public std::vector<unsigned int>
    {
        using std::vector<unsigned int>::vector;
    };

    inline std::ostream & operator<<(std::ostream & os, const Sequence & sequence) noexcept
    {
        return os << Utils::to_string(sequence);
    }
}