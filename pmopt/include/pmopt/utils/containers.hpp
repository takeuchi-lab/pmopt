#pragma once

#include <sstream>
#include <string>
#include <vector>
#include <concepts>

namespace PMOpt
{
    namespace Utils
    {

        /**
         * @brief Concept for printable class (i.e. has overridden operator<<)
         * @tparam T 
         */
        template < class T >
        concept Printable = requires ( const T & x )
        {
            { std::declval< std::ostream & >() << x };
        };


        /**
         * @brief convert iterable container into string
         * @tparam T type of element of container
         * @tparam U type of container
         * @param iterable 
         * @return std::string 
         */
        template < Utils::Printable T, template < class... > class U >
        std::string to_string(const U<T> & iterable, std::string delim = " ")
        {
            std::stringstream ss;
            bool is_first = true;
            for (const T & it : iterable)
            {
                if (!is_first)
                    ss << delim;
                ss << it;
                is_first = false;
            }
            return ss.str();
        }


        /**
         * @brief Check if `key` is contained in `values`
         * @param key 
         * @param values 
         * @return bool 
         */
        bool contains(const std::string & key, const std::vector<std::string> & values);
    }


    // template < Utils::Printable T, template < class... > class U >
    // inline std::ostream & operator<<(std::ostream & os, const U<T> & iterable)
    // {
    //     bool is_first = true;
    //     for (const T & it : iterable)
    //     {
    //         if (!is_first)
    //             os << std::string(" ");
    //         os << it;
    //         is_first = false;
    //     }

    //     return os;
    // }
}