CMAKE_MINIMUM_REQUIRED(VERSION 3.16)
project(cartdg VERSION 0.0
        DESCRIPTION "Discontinuous Galerkin numerics for CFD on a Cartesian grid"
        LANGUAGES CXX)
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Select Release for high performance and Debug for easier debugging" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

if(NOT MAX_BASIS_RANK)
   set(MAX_BASIS_RANK 8 CACHE STRING
       "Maximum degree of basis polynomials + 1. Must be an integer greater than 1." FORCE)
endif()
execute_process(COMMAND python3 ${CMAKE_SOURCE_DIR}/script/auto_generate.py ${CMAKE_BINARY_DIR} ${MAX_BASIS_RANK})

find_library(tecio REQUIRED)
add_library(cartdg SHARED)
target_sources(cartdg PUBLIC ${CMAKE_BINARY_DIR}/Gauss_lobatto.cpp
                             ${CMAKE_BINARY_DIR}/Solution__get_kernel.cpp)

add_subdirectory(src)
target_include_directories(cartdg PUBLIC include)
add_subdirectory(test)
add_subdirectory(demo)
