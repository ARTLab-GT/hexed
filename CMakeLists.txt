CMAKE_MINIMUM_REQUIRED(VERSION 3.16)

project(hexed VERSION 0.2.1
        DESCRIPTION "discontinuous Galerkin CFD solver with automated, unstructured quad/hex meshing"
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

option(autoinstall_python "automatically use pip to install python dependencies and hexedpy in your current environment" OFF)
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
if (use_tecplot)
    set(TECPLOT_DIR "/opt/tecplot/360ex_2020r2" CACHE STRING "where to find tecplot")
    set(CMAKE_INSTALL_RPATH ${TECPLOT_DIR}/bin)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "-g3 -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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

# install Python dependencies
if (autoinstall_python)
    add_custom_command(OUTPUT py_deps_installed
        COMMAND ${CMAKE_SOURCE_DIR}/script/install_py_deps.sh
        COMMAND touch py_deps_installed
        DEPENDS ${CMAKE_SOURCE_DIR}/script/install_py_deps.sh
    )
    add_custom_target(install_py_deps ALL DEPENDS py_deps_installed)
else()
    add_custom_target(install_py_deps ALL)
endif()

# configure files
configure_file(hexed-config.cmake.in hexed-config.cmake)
configure_file(hexed_config.hpp.in include/auto/hexed_config.hpp)
configure_file(config.cpp.in config.cpp)
configure_file(script/benchmark.py.in benchmark.py)
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR}/include/hexed SYMBOLIC)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipsoid.igs ${CMAKE_BINARY_DIR}/ellipsoid.igs)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipsoid.stp ${CMAKE_BINARY_DIR}/ellipsoid.stp)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipse.IGS ${CMAKE_BINARY_DIR}/ellipse.IGS)
file(COPY_FILE ${CMAKE_SOURCE_DIR}/test/ellipse.STEP ${CMAKE_BINARY_DIR}/ellipse.STEP)
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
add_custom_command(
    OUTPUT constants.hil
    COMMAND python3 ${CMAKE_SOURCE_DIR}/script/translate.py ${CMAKE_SOURCE_DIR}/include/constants.hpp
    DEPENDS         ${CMAKE_SOURCE_DIR}/script/translate.py ${CMAKE_SOURCE_DIR}/include/constants.hpp
)
add_custom_target(constants ALL DEPENDS constants.hil)

# build Python package
file(COPY ${CMAKE_SOURCE_DIR}/python DESTINATION ${CMAKE_BINARY_DIR})
file(COPY ${CMAKE_SOURCE_DIR}/LICENSE.txt DESTINATION ${CMAKE_BINARY_DIR}/python/)
configure_file(python/pyproject.toml.in python/pyproject.toml)
set(wheel python/dist/hexedpy-${CMAKE_PROJECT_VERSION}-py3-none-any.whl)
set(pysrc python/dist/hexedpy-${CMAKE_PROJECT_VERSION}.tar.gz)
add_custom_command(
    OUTPUT ${wheel} ${pysrc}
    COMMAND ${CMAKE_SOURCE_DIR}/script/build_python.sh
    DEPENDS constants.hil python ${CMAKE_SOURCE_DIR}/LICENSE.txt ${CMAKE_SOURCE_DIR}/script/build_python.sh py_deps_installed
)
add_custom_target(hexedpy_package ALL DEPENDS ${wheel} ${pysrc})
if (autoinstall_python)
    # note: Both install commands are necessary to ensure that `hexedpy` is reinstalled with each new build.
    # The first one installs `hexedpy` with all its dependencies (but won't reinstall anything if it has already been installed)
    # and the second forces the reinstallation of `hexedpy` itself, but not its dependencies.
    # Thus `hexedpy` will always be reinstalled, but the dependencies will only be installed if they aren't already.
    add_custom_command(
        OUTPUT hexedpy_installed
        COMMAND pip3 install --upgrade ${wheel}
        COMMAND pip3 install --force-reinstall --no-deps ${wheel}
        COMMAND touch hexedpy_installed
        DEPENDS ${wheel}
    )
    add_custom_target(install_hexedpy ALL DEPENDS hexedpy_installed)
endif()

# autogenerate C++ files
configure_file(script/autogen.sh.in autogen.sh)
file(CHMOD autogen.sh PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
add_custom_command(
    OUTPUT Gauss_lobatto.cpp Gauss_legendre.cpp
    COMMAND ./autogen.sh
    DEPENDS
        autogen.sh
        ${CMAKE_SOURCE_DIR}/script/auto_generate.py
        python
        py_deps_installed
)

# generate colormaps to use in tecplot
add_custom_command(
    OUTPUT ${colormap_files}
    COMMAND python3 ${CMAKE_SOURCE_DIR}/script/colormaps.py
    DEPENDS ${CMAKE_SOURCE_DIR}/script/colormaps.py
)
add_custom_target(
    colormaps ALL
    DEPENDS ${colormap_files} py_deps_installed
)

# define main target, which is a library
add_library(hexed SHARED)

# add external libraries
find_package(HDF5 REQUIRED)
target_link_libraries(hexed PRIVATE hdf5_cpp)

if (use_tecplot)
    target_include_directories(hexed PRIVATE ${TECPLOT_DIR}/include)
    target_link_directories(hexed PRIVATE ${TECPLOT_DIR}/bin)
    target_link_libraries(hexed PRIVATE tecio)
endif()
if (use_nlopt)
    find_package(NLopt REQUIRED)
    target_link_libraries(hexed PRIVATE nlopt)
endif()
if (use_occt)
    find_package(OpenCASCADE REQUIRED)
    target_link_libraries(hexed PUBLIC TKBRep TKOpenGl TKV3d)
    if (OpenCASCADE_VERSION_MINOR GREATER 7)
        target_link_libraries(hexed PUBLIC TKDEIGES TKDESTEP TKDESTL)
    else()
        target_link_libraries(hexed PUBLIC TKIGES TKSTEP TKSTL)
    endif()
endif()
if (use_xdmf)
    find_package(Xdmf REQUIRED)
    find_package(libxml2 REQUIRED)
    target_link_libraries(hexed PRIVATE Xdmf)
    target_include_directories(hexed PRIVATE ${LIBXML2_INCLUDE_DIRS})
endif()

# add include/source files
target_sources(hexed PRIVATE ${CMAKE_BINARY_DIR}/Gauss_lobatto.cpp ${CMAKE_BINARY_DIR}/Gauss_legendre.cpp ${CMAKE_BINARY_DIR}/config.cpp)
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

# install
install(PROGRAMS script/hexed-post-process ${CMAKE_BINARY_DIR}/hexed-wrapped-post-process script/hexed-fit-geom TYPE BIN)
install(TARGETS hexed hexecute hil)
export(TARGETS hexed FILE hexed-targets.cmake)
install(FILES
    ${CMAKE_SOURCE_DIR}/hil/builtin.hil
    ${CMAKE_SOURCE_DIR}/hil/interactive.hil
    ${CMAKE_SOURCE_DIR}/hil/hexed.hil
    ${CMAKE_BINARY_DIR}/constants.hil
    ${CMAKE_SOURCE_DIR}/LICENSE.txt
    DESTINATION lib/hexed
)
install(FILES ${CMAKE_BINARY_DIR}/hexed-config.cmake DESTINATION lib/cmake/hexed)

# uninstall
configure_file(${CMAKE_SOURCE_DIR}/script/uninstall.sh.in uninstall.sh)
add_custom_target(uninstall
    COMMAND ./uninstall.sh
)
