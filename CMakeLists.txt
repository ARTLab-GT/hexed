CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

project(cartdg VERSION 0.0
        DESCRIPTION "Discontinuous Galerkin numerics for CFD on a Cartesian grid"
        LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# determine build type
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Select Release for high performance and Debug for easier debugging" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS "-Wall")

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "script/auto_generate.py;script/benchmark.cpp.in;script/benchmark.py")

# determine whether to compile with sanitizers
option(sanitize "If ON, compile with sanitizer flags." OFF)
if (sanitize)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=leak -fsanitize=undefined")
endif()

# determine whether to parallelize
option(threaded "If ON, allow multiple threads. Substantially reduces execution time." ON)
if(threaded)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# determine rank of polynomial basis
if(NOT DEFINED MAX_BASIS_RANK)
   set(MAX_BASIS_RANK 8 CACHE STRING
       "Maximum degree of basis polynomials + 1. Must be an integer greater than 1." FORCE)
endif()
execute_process(COMMAND python3 ${CMAKE_SOURCE_DIR}/script/auto_generate.py ${CMAKE_BINARY_DIR} ${MAX_BASIS_RANK})

# determine rank of polynomial basis
if(NOT DEFINED TECIO_DIR)
   set(TECIO_DIR /usr/local/tecplot/360ex_2018r2 CACHE STRING
       "Where to find TECIO (probably the root directory of your tecplot installation)" FORCE)
endif()

# add a configure file to pass CMake variables to the code
configure_file(cartdgConfig.hpp.in cartdgConfig.hpp)

# define main target, which is a library
add_library(cartdg SHARED)

# add tecio as external library
add_library(tecio SHARED IMPORTED)
set_target_properties(tecio PROPERTIES IMPORTED_LOCATION ${TECIO_DIR}/bin/libtecio.so)
target_include_directories(cartdg PRIVATE ${TECIO_DIR}/include)

# add automatically generated source files
target_sources(cartdg PRIVATE ${CMAKE_BINARY_DIR}/Gauss_lobatto.cpp)

target_include_directories(cartdg PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
                           $<INSTALL_INTERFACE:include>
)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(demo)

# object library to contain the kernels
add_library(kernels OBJECT)
set_property(TARGET kernels PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(kernels PRIVATE ${CMAKE_SOURCE_DIR}/kernels ${CMAKE_SOURCE_DIR}/include)
target_include_directories(kernels PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/kernels/include>)
add_subdirectory(${CMAKE_BINARY_DIR}/kernels/src)
target_link_libraries(cartdg PRIVATE kernels)

# add executable for benchmarking
add_executable(benchmark benchmark/main.cpp)
set_property(TARGET benchmark PROPERTY RUNTIME_OUTPUT_DIRECTORY benchmark)
target_link_libraries(benchmark PRIVATE cartdg kernels tecio)
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/script/benchmark.py benchmark.py SYMBOLIC)

install(TARGETS cartdg DESTINATION lib EXPORT cartdgTargets)
install(DIRECTORY include/ DESTINATION include/cartdg)
install(FILES ${CMAKE_BINARY_DIR}/cartdgConfig.hpp DESTINATION include/cartdg)
install(EXPORT cartdgTargets
        NAMESPACE cartdg::
        FILE cartdgTargets.cmake
        DESTINATION lib/cmake/cartdg
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
"${CMAKE_CURRENT_BINARY_DIR}/cartdgConfig.cmake"
INSTALL_DESTINATION "lib/cmake"
NO_SET_AND_CHECK_MACRO
NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# generate the version file for the config file
write_basic_package_version_file(
"${CMAKE_CURRENT_BINARY_DIR}/cartdgConfigVersion.cmake"
VERSION "${cartdg_VERSION_MAJOR}.${cartdg_VERSION_MINOR}"
COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES
${CMAKE_CURRENT_BINARY_DIR}/cartdgConfig.cmake
DESTINATION lib/cmake/cartdg
)
