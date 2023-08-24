CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

project(hexed VERSION 0.0.0
        DESCRIPTION "Discontinuous Galerkin engine for CFD with automated, unstructured quad/hex meshing."
        LANGUAGES CXX)
# get info about current build for documentation purposes
string(TIMESTAMP DATE "%Y-%m-%d") # get the current date for reference purposes
execute_process(
    COMMAND git log -1 --format=%H
    OUTPUT_VARIABLE COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS False)

# determine build type
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
       "Select Release for high performance and Debug for easier debugging" FORCE)
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug")
endif()
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# default installation directory is one that we have access to on lab machines
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "~/.local" CACHE PATH "where to install hexed" FORCE)
endif()
if (NOT DEFINED PYTHON_INSTALL_PREFIX)
    set(PYTHON_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib CACHE PATH "where to install Python module" FORCE)
endif()

option(build_tests "would you like to build the unit tests? (requires Catch2)" OFF)
option(build_docs "build documentation? (requires doxygen)" OFF)
option(sanitize "If ON, compile with sanitizer flags." OFF)
option(profile "if ON, compile with gprof profiling flags (in which case `threaded` should be OFF)")
option(threaded "If ON, allow multiple threads. Substantially reduces execution time." ON)
option(use_tecplot "include capability to use visualize in Tecplot's native format" OFF)
option(use_xdmf "include capability to use visualize in XDMF format" ON)
option(use_cppcheck "use the cppcheck static analysis tool (requires cppcheck)" OFF)
option(use_nlopt "enable features that require root finding with NLopt" OFF)
option(use_occt "enable CAD geometry processing with OpenCASCADE Technology" ON)
option(obsessive_timing "turn on additional timing measurements, potentially to the point of slowing down the code" OFF)

# determine whether to compile with sanitizers
if (sanitize)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=bounds-strict -fsanitize=undefined -fsanitize=address -fsanitize=leak")
    if (CMAKE_CXX_COMPILER_VERSION AND CMAKE_CXX_COMPILER_VERSION GREATER_EQUAL 9.0.0)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=pointer-compare -fsanitize=pointer-subtract")
    endif()
endif()

# determine whether to parallelize
include(ProcessorCount)
if(threaded)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    ProcessorCount(n_threads)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    set(n_threads, 1)
endif()

if (profile)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
endif()

# determine row_size of polynomial basis
if(NOT DEFINED MAX_BASIS_ROW_SIZE)
   set(MAX_BASIS_ROW_SIZE 8 CACHE STRING
       "Maximum degree of basis polynomials + 1. Must be an integer greater than 1." FORCE)
endif()

# configure files
configure_file(hexed-config.cmake.in hexed-config.cmake)
configure_file(hexed_config.hpp.in include/auto/hexed_config.hpp)
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include/hexed SYMBOLIC)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipsoid.igs ${CMAKE_BINARY_DIR}/ellipsoid.igs)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipsoid.stp ${CMAKE_BINARY_DIR}/ellipsoid.stp)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipse.IGS ${CMAKE_BINARY_DIR}/ellipse.IGS)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipse.STEP ${CMAKE_BINARY_DIR}/ellipse.STEP)
add_custom_command(
    OUTPUT Gauss_lobatto.cpp Gauss_legendre.cpp
    COMMAND python3 ${CMAKE_SOURCE_DIR}/script/auto_generate.py ${CMAKE_BINARY_DIR} ${MAX_BASIS_ROW_SIZE}
    DEPENDS
        ${CMAKE_SOURCE_DIR}/script/auto_generate.py
        ${CMAKE_SOURCE_DIR}/script/Basis.py
)
set(colormap_files
    colormaps/cividis.map
    colormaps/inferno.map
    colormaps/magma.map
    colormaps/plasma.map
    colormaps/turbo.map
    colormaps/twilight.map
    colormaps/viridis.map
)
configure_file(script/hexed-wrapped-post-process.in hexed-wrapped-post-process)
configure_file(script/format.mcr.in format.mcr)

# generate colormaps to use in tecplot
add_custom_command(
    OUTPUT ${colormap_files}
    COMMAND python3 ${CMAKE_SOURCE_DIR}/script/colormaps.py
    DEPENDS ${CMAKE_SOURCE_DIR}/script/colormaps.py
)
add_custom_target(
    colormaps ALL
    DEPENDS ${colormap_files}
)

# define main target, which is a library
add_library(hexed SHARED)

# add external libraries
if (use_tecplot)
    target_link_libraries(hexed PRIVATE tecio)
endif()
if (use_nlopt)
    find_package(NLopt REQUIRED)
    target_link_libraries(hexed PRIVATE nlopt)
endif()
if (use_occt)
    find_package(OpenCASCADE REQUIRED)
    target_link_libraries(hexed PUBLIC TKBRep TKOpenGl TKV3d TKIGES TKSTEP TKSTL)
endif()
if (use_xdmf)
    find_package(Xdmf REQUIRED)
    find_package(libxml2 REQUIRED)
    target_link_libraries(hexed PRIVATE Xdmf)
    target_include_directories(hexed PRIVATE ${LIBXML2_INCLUDE_DIRS})
endif()

# add include/source files
target_sources(hexed PRIVATE ${CMAKE_BINARY_DIR}/Gauss_lobatto.cpp ${CMAKE_BINARY_DIR}/Gauss_legendre.cpp)
target_include_directories(hexed PUBLIC ${CMAKE_BINARY_DIR}/include/auto)
target_include_directories(hexed PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(hexed INTERFACE ${CMAKE_BINARY_DIR}/include)
add_subdirectory(src)
if (build_tests)
    add_subdirectory(test)
endif()
if(build_docs)
    add_subdirectory(doc)
endif()

# run cppcheck
if (use_cppcheck)
    add_custom_command(OUTPUT cppcheck_output.txt COMMAND cppcheck --project=compile_commands.json --inline-suppr > cppcheck_output.txt DEPENDS hexed)
    add_custom_target(cppcheck ALL DEPENDS cppcheck_output.txt)
endif()

add_custom_command(
    OUTPUT constants.hil
    COMMAND python3 ${CMAKE_SOURCE_DIR}/script/translate.py ${CMAKE_SOURCE_DIR}/include/constants.hpp
    DEPENDS         ${CMAKE_SOURCE_DIR}/script/translate.py ${CMAKE_SOURCE_DIR}/include/constants.hpp
)
add_custom_target(constants ALL DEPENDS constants.hil)

# install
install(PROGRAMS script/hexed-post-process ${CMAKE_BINARY_DIR}/hexed-wrapped-post-process script/hexed-fit-geom TYPE BIN)
install(TARGETS hexed hexecute hil)
export(TARGETS hexed FILE hexed-targets.cmake)
install(FILES ${CMAKE_BINARY_DIR}/hexed-config.cmake DESTINATION lib/cmake/hexed)
install(FILES ${CMAKE_SOURCE_DIR}/python/hexed_utils.py DESTINATION ${PYTHON_INSTALL_PREFIX})
