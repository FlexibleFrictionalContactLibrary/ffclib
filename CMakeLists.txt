# ffclib: fclib's global frictional contact problem with compliant bilateral rows (FORMAT.md). The library
# (src/ffclib/, included as <ffclib/...>) reads and writes the files, evaluates the error measures and solves
# the problems (src/ffclib/solver/); examples/solve.cpp shows its use.
#
#   cmake -S . -B build [-DFFCLIB_WITH_CLARABEL=ON -DFFCLIB_WITH_SCS=ON -DFFCLIB_WITH_QOCO=ON -DFFCLIB_WITH_CONICXX=ON]
#         [-DFFCLIB_WITH_CUOPT=ON -DCUOPT_ROOT=<venv>]
#   cmake --build build -j && ctest --test-dir build && cmake --install build --prefix <prefix>
#
# Eigen and HDF5 are taken from an installation or fetched (cmake/Dependencies.cmake). The conic solvers Clarabel
# (needs Rust/cargo), SCS, QOCO and ConicXX are fetched, cuOpt is found in an installation. At top level ffclib is
# a shared library that carries the solvers and a fetched HDF5 (SCS's shared library is installed next to it); as a
# subproject (add_subdirectory, FetchContent) it is static and builds neither tests nor examples. Link
# ffclib::ffclib; installed: find_package(ffclib).
cmake_minimum_required(VERSION 3.21)
project(ffclib VERSION 0.1.0 LANGUAGES C CXX)

option(FFCLIB_BUILD_TESTS "Build the tests (fetches fclib for the compatibility test)" ${PROJECT_IS_TOP_LEVEL})
option(FFCLIB_BUILD_EXAMPLES "Build the examples" ${PROJECT_IS_TOP_LEVEL})
option(FFCLIB_SHARED "Build ffclib as a shared library (needed to install it with conic solvers)" ${PROJECT_IS_TOP_LEVEL})
option(FFCLIB_WITH_CLARABEL "Conic solver Clarabel (fetched; needs Rust/cargo)" ON)
option(FFCLIB_WITH_SCS "Conic solver SCS (fetched)" ON)
option(FFCLIB_WITH_QOCO "Conic solver QOCO (fetched)" ON)
option(FFCLIB_WITH_CONICXX "Conic solver ConicXX (fetched)" ON)
option(FFCLIB_WITH_CUOPT "Conic solver NVIDIA cuOpt (GPU; an installed cuOpt >= 26.10.1, e.g. pip libcuopt-mathopt-cu12)" OFF)
option(FFCLIB_INSTALL "Install the library and its CMake package" ${PROJECT_IS_TOP_LEVEL})
set(_ffclib_conic OFF)
if(FFCLIB_WITH_CLARABEL OR FFCLIB_WITH_SCS OR FFCLIB_WITH_QOCO OR FFCLIB_WITH_CONICXX OR FFCLIB_WITH_CUOPT)
    set(_ffclib_conic ON)
endif()
# A static ffclib would hand the solvers' libraries on to its users, which the CMake package cannot carry.
if(FFCLIB_INSTALL AND _ffclib_conic AND NOT FFCLIB_SHARED)
    message(FATAL_ERROR "ffclib: installing it with conic solvers needs FFCLIB_SHARED=ON")
endif()

include(FetchContent)
include(GNUInstallDirs)

if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ffclib's own modules by their path: a parent project may have modules of the same name in its CMAKE_MODULE_PATH.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")  # FindCUOPT.cmake
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/Dependencies.cmake")
# The CMake package of a static ffclib hands HDF5 on to its users and finds it with find_dependency().
if(FFCLIB_INSTALL AND FFCLIB_FETCHED_HDF5 AND NOT FFCLIB_SHARED)
    message(FATAL_ERROR "ffclib: installing a static ffclib needs an installed HDF5 (or FFCLIB_SHARED=ON, which carries the fetched one)")
endif()

set(FFCLIB_WARNINGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(FFCLIB_WARNINGS -Wall -Wextra)
endif()

if(FFCLIB_SHARED)
    set(_ffclib_type SHARED)
else()
    set(_ffclib_type STATIC)
endif()
add_library(ffclib ${_ffclib_type} src/ffclib/problem.cpp src/ffclib/io.cpp src/ffclib/error.cpp src/ffclib/solver/gauss_seidel.cpp src/ffclib/solver/conic/cone_program.cpp)
add_library(ffclib::ffclib ALIAS ffclib)
if(FFCLIB_SHARED)
    # The installed library finds SCS's next to itself and cuOpt's where it was linked.
    set_target_properties(ffclib PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR} INSTALL_RPATH "$ORIGIN" INSTALL_RPATH_USE_LINK_PATH TRUE)
    # Only ffclib's API is exported; the statically linked solvers (QOCO's AMD and QDLDL, Clarabel's Rust
    # runtime, ConicXX) stay local and cannot clash with other copies in the process.
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(ffclib PRIVATE "LINKER:--version-script=${CMAKE_CURRENT_SOURCE_DIR}/cmake/ffclib_exports.map")
    endif()
endif()
target_compile_features(ffclib PUBLIC cxx_std_17)
set_target_properties(ffclib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(ffclib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(ffclib PUBLIC Eigen3::Eigen PRIVATE hdf5::hdf5)
target_compile_options(ffclib PRIVATE ${FFCLIB_WARNINGS})
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/ConicSolvers.cmake")

if(FFCLIB_BUILD_TESTS)
    enable_testing()
    # fclib itself (C, needs HDF5's high-level API) checks that both read each other's files. Only its source is
    # compiled: its CMake project builds its own tests and needs SuiteSparse.
    FetchContent_Declare(fclib GIT_REPOSITORY https://github.com/FrictionalContactLibrary/fclib.git
        GIT_TAG f894f34f02093910ceb49943d69f4901142bb880 SOURCE_SUBDIR do-not-configure)
    FetchContent_MakeAvailable(fclib)
    add_library(ffclib_test_fclib STATIC "${fclib_SOURCE_DIR}/src/fclib.c")
    set_target_properties(ffclib_test_fclib PROPERTIES C_STANDARD 99)
    target_include_directories(ffclib_test_fclib SYSTEM PUBLIC "${fclib_SOURCE_DIR}/src")
    target_link_libraries(ffclib_test_fclib PUBLIC hdf5::hdf5 hdf5::hdf5_hl)

    foreach(_test test_io test_error test_gauss_seidel test_conic)
        add_executable(${_test} tests/${_test}.cpp)
        target_link_libraries(${_test} PRIVATE ffclib::ffclib)
        target_compile_options(${_test} PRIVATE ${FFCLIB_WARNINGS})
        add_test(NAME ${_test} COMMAND ${_test})
    endforeach()
    target_link_libraries(test_io PRIVATE ffclib_test_fclib)
endif()

if(FFCLIB_BUILD_EXAMPLES)
    add_executable(ffclib_solve examples/solve.cpp)
    target_link_libraries(ffclib_solve PRIVATE ffclib::ffclib)
    target_compile_options(ffclib_solve PRIVATE ${FFCLIB_WARNINGS})
endif()

if(FFCLIB_INSTALL)
    include(CMakePackageConfigHelpers)
    install(TARGETS ffclib EXPORT ffclibTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    if(FFCLIB_WITH_SCS)
        install(TARGETS scsdir LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})  # found through ffclib's RPATH $ORIGIN
    endif()
    install(FILES LICENSE NOTICE DESTINATION ${CMAKE_INSTALL_DATADIR}/doc/ffclib)
    install(DIRECTORY src/ffclib DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.hpp")
    install(EXPORT ffclibTargets NAMESPACE ffclib:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ffclib)
    configure_package_config_file(cmake/ffclibConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/ffclibConfig.cmake" INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ffclib)
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/ffclibConfigVersion.cmake" COMPATIBILITY SameMinorVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ffclibConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/ffclibConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ffclib)
endif()
