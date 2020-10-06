CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
project(cartdg VERSION 0.0
        DESCRIPTION "Discontinuous Galerkin numerics for CFD on a Cartesian grid"
        LANGUAGES CXX)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_library(tecio REQUIRED)
add_library(cartdg SHARED)

add_subdirectory(src)
target_include_directories(cartdg PUBLIC include)
add_subdirectory(test)
add_subdirectory(demo)
