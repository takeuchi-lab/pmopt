include(CheckCXXCompilerFlag)

# CMake version
cmake_minimum_required(VERSION 3.21)
enable_testing()

# Project name
project(pmopt CXX)

# Language setting
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF) # without gnu extensions

# Set compile flags for release
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -Wsign-compare -Wsign-conversion -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Check sanitizer availability
set(FLAGS_SANITIZER "-fsanitize=address,undefined")
set(CMAKE_REQUIRED_FLAGS ${FLAGS_SANITIZER})
check_cxx_compiler_flag(${FLAGS_SANITIZER} SANITIZER_SUPPORTED)
unset(CMAKE_REQUIRED_FLAGS)

# Set compile flags for debug
message(STATUS "Sanitizer Availability: ${SANITIZER_SUPPORTED}")
if(SANITIZER_SUPPORTED)
    set(CMAKE_CXX_FLAGS_DEBUG "-g ${FLAGS_SANITIZER}")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif()


# Add include directories globally
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/sysinclude)

# Add project sub-directories
add_subdirectory(src)