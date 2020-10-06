CMAKE_MINIMUM_REQUIRED(VERSION 3.17)
project(cartdg VERSION 0.0
        DESCRIPTION "Discontinuous Galerkin numerics for CFD on a Cartesian grid"
        LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(cartdg SHARED)

add_subdirectory(src)
target_include_directories(cartdg PUBLIC include)
