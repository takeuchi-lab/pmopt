#pragma once

#include <exception>
#include <fstream>
#include <sstream>
#include <string>

namespace PMOpt
{
    namespace Utils
    {

        std::string get_structure_type(
            const std::vector<std::string> & candidates,
            const std::string & filename
        )
        {
            std::ifstream ifs(filename);
            
            std::string buffer;
            std::getline(ifs, buffer);

            std::string structure_type;
            std::stringstream ss(std::move(buffer));
            while (std::getline(ss, buffer, ','))
            {
                for (const auto & candidate : candidates)
                {
                    if (buffer != candidate)
                        continue;

                    if (structure_type == std::string())
                        structure_type = buffer;
                    else
                        throw std::runtime_error("multiple structure column found in" + filename);
                }
            }
            return structure_type;
        }


        unsigned int stopi(const std::string & text)
        {
            auto value = std::stoi(text);
            if (value <= 0)
                throw std::invalid_argument("not strictly positive");
            return (unsigned int) value;
        }

        unsigned int stoui(const std::string & text)
        {
            auto value = std::stoi(text);
            if (value < 0)
                throw std::invalid_argument("smaller that zero");
            return (unsigned int) value;
        }
    }
}