CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

project(cartdg VERSION 0.0
        DESCRIPTION "Discontinuous Galerkin numerics for CFD on a Cartesian grid"
        LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

# determine build type
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Select Release for high performance and Debug for easier debugging" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS "-Wall")

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "script/auto_generate.py;script/Basis.py")

# determine whether to compile with sanitizers
option(sanitize "If ON, compile with sanitizer flags." OFF)
if (sanitize)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=bounds-strict -fsanitize=undefined -fsanitize=address -fsanitize=leak")
    if (CMAKE_CXX_COMPILER_VERSION AND CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 9.0.0)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=pointer-compare -fsanitize=pointer-subtract")
    endif()
endif()

# determine whether to parallelize
include(ProcessorCount)
option(threaded "If ON, allow multiple threads. Substantially reduces execution time." ON)
if(threaded)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    ProcessorCount(n_threads)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    set(n_threads, 1)
endif()

# determine row_size of polynomial basis
if(NOT DEFINED MAX_BASIS_ROW_SIZE)
   set(MAX_BASIS_ROW_SIZE 8 CACHE STRING
       "Maximum degree of basis polynomials + 1. Must be an integer greater than 1." FORCE)
endif()
execute_process(COMMAND python3 ${CMAKE_SOURCE_DIR}/script/auto_generate.py ${CMAKE_BINARY_DIR} ${MAX_BASIS_ROW_SIZE})

option(use_otter "include capability to use otter for plotting" ON)
option(use_tecplot "include capability to use tecplot for plotting" ON)

# add a configure file to pass CMake variables to the code
configure_file(config_build.hpp.in config_build.hpp)
configure_file(config_build.hpp.in config.hpp)

# define main target, which is a library
add_library(cartdg SHARED)
# add external libraries
if (use_otter)
    find_package(otter REQUIRED)
    target_link_libraries(cartdg PUBLIC otter)
endif()
if (use_tecplot)
    target_link_libraries(cartdg PRIVATE tecio)
endif()
# add automatically generated source files
target_sources(cartdg PRIVATE ${CMAKE_BINARY_DIR}/Gauss_lobatto.cpp ${CMAKE_BINARY_DIR}/Gauss_legendre.cpp)
target_include_directories(cartdg PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
                           $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
                           $<INSTALL_INTERFACE:include>
)

add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(demo)

install(TARGETS cartdg DESTINATION lib EXPORT cartdgTargets)
install(DIRECTORY include/ DESTINATION include/cartdg)
install(FILES ${CMAKE_BINARY_DIR}/config.hpp DESTINATION include/cartdg/)
install(PROGRAMS script/cartdg-post-process script/cartdg-wrapped-post-process TYPE BIN)
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
