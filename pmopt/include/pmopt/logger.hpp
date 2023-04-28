#pragma once

#include <chrono>
#include <iostream>
#include <iomanip>
#include "pmopt/utils/containers.hpp"

namespace PMOpt
{
    struct Logger
    {

        const static unsigned int ESSENTIAL = 1;
        const static unsigned int INFO = 1;
        const static unsigned int DEBUG = 50;
        const static unsigned int WIDTH = 10;

        static unsigned int _verbose;

        std::ostream & _os;
        unsigned int _indent;

        Logger(std::ostream & os, unsigned int indent) : _os(os), _indent(indent) {} 

        Logger(std::ostream & os) : Logger(os, 0) {}

        
        /**
         * @brief Logging function
         * @tparam T : printable
         * @param priority Level of log, smaller is more important 
         * @param args Array of objects to be printed
         */
        template < Utils::Printable ... Args >
        inline void log(unsigned int priority, const Args & ... args) const noexcept
        {
            if (priority > _verbose)
                return;

            // _os << std::string(_indent, ' ');

            // _os << "[ " << std::setw(WIDTH) << key << " ]";
            __log(std::forward<const Args &>(args)...);
            _os << std::endl;
        }


        template < Utils::Printable ... Args >
        inline void operator()(unsigned int priority, const Args & ... args) const noexcept
        {
            log(priority, std::forward<const Args &>(args)...);
        }


        // logging for multiple arguments
        template < Utils::Printable T, Utils::Printable ... Args >
        inline void __log(const T & head, const Args & ... args) const noexcept
        {
            _os << head;
            __log(std::forward<const Args &>(args)...);
        }

        // logging for last one
        template < Utils::Printable T >
        inline void __log(const T & last) const noexcept
        {
            _os << last;
        }


        // logging for empty arguments
        // inline void __log() const noexcept {}

        // print current time
        inline std::string now() const noexcept
        {
            auto now_clock = std::chrono::system_clock::now();
            auto now_time = std::chrono::system_clock::to_time_t(now_clock);
            std::string now_string = std::ctime(& now_time);
            return now_string.substr(0, now_string.size() - 1);
        }
    };
}