#pragma once

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <exception>

#include "pmopt/parameters.hpp"

namespace PMOpt
{

    namespace Loader
    {
        const static char SEP = ',';
        const static char COL_NUMERICAL = '#';
        const static std::string COL_OBJECTIVE = "objective";
        const static std::string COL_WEIGHT = "weight";
        const static std::string COL_GROUP = "group";


        // Load database from given file
        template < class Database >
        void load(Database & database, Parameters & params, const std::string & filename);
    }

}


template < class Database >
void PMOpt::Loader::load(Database & database, Parameters & params, const std::string & filename)
{

    std::ifstream ifs(filename);
    if (!ifs)
    {
        std::cerr << "Error: File cannot open `" << filename << "`.\n";
        exit(1);
    }


    std::string buffer;

    // Read header of csv
    std::vector<std::string> columns;
    std::getline(ifs, buffer);
    std::stringstream ss(std::move(buffer));
    std::unordered_map<std::string, unsigned int> num2id;
    while (std::getline(ss, buffer, SEP))
    {
        if (buffer[0] == COL_NUMERICAL)
            params._X_num_names.push_back(buffer);
        columns.emplace_back(std::move(buffer));
    }

    // initialization for numerical features
    for (const auto & name : params._X_num_names)
        num2id[name] = num2id.size();
    if (params._X_num.empty())
    {
        params._X_num.reserve(num2id.size());
        for (unsigned int j = 0; j < num2id.size(); ++j)
            params._X_num.push_back(std::make_shared<std::vector<double>>());
    }

    // Read data
    unsigned int line = 2;
    while (std::getline(ifs, buffer))
    {

        // Read one line
        unsigned int index = 0;
        ss = std::stringstream(std::move(buffer));
        while (std::getline(ss, buffer, SEP))
        {
            try
            {
                if (columns[index] == database._format)
                    database.push(std::move(buffer));
                else if (columns[index] == COL_OBJECTIVE)
                    params._y.push_back(std::stod(buffer));
                else if (columns[index] == COL_WEIGHT)
                    params._w.push_back(std::stod(buffer));
                else if (columns[index] == COL_GROUP)
                    params._group.push_back(std::stod(buffer));
                else if (num2id.find(columns[index]) != num2id.end())
                    params._X_num[num2id[columns[index]]]->push_back(std::stod(buffer));
            }
            catch (std::invalid_argument &)
            {
                std::cerr << "Error: Invalid type of column `" << columns[index] 
                    << "` in line: " << line << " of " << filename << '\n';
                exit(1);
            }
            catch (std::runtime_error & error)
            {
                std::cerr << "Error: Invalid input of column `" << columns[index]
                    << "` in line: " << line << " of " << filename << '\n';
                std::cerr << "\twhat: " << error.what() << '\n';
                exit(1);
            }

            
            ++index;
        }

        if (index != columns.size() && columns.back() == database._format)
        {
            database.push(std::string());
            ++index;
        }

        if (index != columns.size())
        {
            std::cerr << "Error: Missing value in line: " << line << " of " << filename << '\n';
            exit(1);
        }

        ++line;
    }
}