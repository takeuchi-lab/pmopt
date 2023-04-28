#pragma once

#include <iostream>
#include <boost/version.hpp>

namespace PMOpt
{
    namespace Utils
    {

        static const std::string VERSION = "0.1.0";

        // print version and exit
        inline void exit_with_version()
        {
            std::cout << "PMOptz : " << VERSION << std::endl;
            std::cout << "Boost  : "     
                << BOOST_VERSION / 100000     << "."  // major version
                << BOOST_VERSION / 100 % 1000 << "."  // minor version
                << BOOST_VERSION % 100                // patch level
                << std::endl;
            exit(0);
        }


    }
   
}
