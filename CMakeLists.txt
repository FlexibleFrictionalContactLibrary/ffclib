# ffclib: fclib's global frictional contact problem with compliant bilateral rows (FORMAT.md). The library
# (src/ffclib/, included as <ffclib/...>) reads and writes the files, evaluates the error measures and solves
# the problems (src/ffclib/solver/); examples/solve.cpp shows its use.
#
#   cmake -S . -B build [-DFFCLIB_WITH_CUOPT=ON -DCUOPT_ROOT=<venv>]
#   cmake --build build -j && ctest --test-dir build
#
# The CPU conic solvers (Clarabel: needs Rust/cargo; SCS, QOCO, ConicXX) are fetched and built by default at top
# level; -DFFCLIB_WITH_<NAME>=OFF leaves one out.
#
# As a subproject: add_subdirectory(<this repository>) and link ffclib::ffclib. Installed: find_package(ffclib).
cmake_minimum_required(VERSION 3.21)
project(ffclib VERSION 0.1.0 LANGUAGES C CXX)

option(FFCLIB_BUILD_TESTS "Build the tests (fetches fclib for the compatibility test)" ON)
option(FFCLIB_BUILD_EXAMPLES "Build the examples" ON)
option(FFCLIB_WITH_CLARABEL "Conic solver Clarabel (fetched; needs Rust/cargo)" OFF)
option(FFCLIB_WITH_SCS "Conic solver SCS (fetched)" OFF)
option(FFCLIB_WITH_QOCO "Conic solver QOCO (fetched)" OFF)
option(FFCLIB_WITH_CONICXX "Conic solver ConicXX (fetched)" OFF)
option(FFCLIB_WITH_CUOPT "Conic solver NVIDIA cuOpt (GPU; an installed cuOpt >= 26.10.1, e.g. pip libcuopt-mathopt-cu12)" OFF)
# The CMake package does not carry the conic solvers' libraries yet: no installation with them.
set(_ffclib_conic OFF)
if(FFCLIB_WITH_CLARABEL OR FFCLIB_WITH_SCS OR FFCLIB_WITH_QOCO OR FFCLIB_WITH_CONICXX OR FFCLIB_WITH_CUOPT)
    set(_ffclib_conic ON)
endif()
if(PROJECT_IS_TOP_LEVEL AND NOT _ffclib_conic)
    set(_ffclib_install ON)
else()
    set(_ffclib_install OFF)
endif()
option(FFCLIB_INSTALL "Install the library and its CMake package (not with the conic solvers)" ${_ffclib_install})

include(FetchContent)
include(GNUInstallDirs)

if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.4 QUIET NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(STATUS "ffclib: Eigen3 not found, fetching Eigen 3.4.0")
        set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(Eigen3 GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0 GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(Eigen3)
    endif()
endif()
find_package(HDF5 REQUIRED COMPONENTS C)

set(FFCLIB_WARNINGS "")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(FFCLIB_WARNINGS -Wall -Wextra)
endif()

add_library(ffclib src/ffclib/problem.cpp src/ffclib/io.cpp src/ffclib/error.cpp src/ffclib/solver/gauss_seidel.cpp src/ffclib/solver/conic/cone_program.cpp)
add_library(ffclib::ffclib ALIAS ffclib)
target_compile_features(ffclib PUBLIC cxx_std_17)
set_target_properties(ffclib PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(ffclib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(ffclib PUBLIC Eigen3::Eigen PRIVATE hdf5::hdf5)
target_compile_options(ffclib PRIVATE ${FFCLIB_WARNINGS})
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(ConicSolvers)

if(FFCLIB_BUILD_TESTS)
    enable_testing()
    # fclib itself (C, needs HDF5's high-level API) checks that both read each other's files. Only its source is
    # compiled: its CMake project builds its own tests and needs SuiteSparse.
    find_package(HDF5 REQUIRED COMPONENTS C HL)
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

if(FFCLIB_INSTALL AND _ffclib_conic)
    message(STATUS "ffclib: not installed, the CMake package does not carry the conic solvers' libraries yet")
elseif(FFCLIB_INSTALL)
    include(CMakePackageConfigHelpers)
    install(TARGETS ffclib EXPORT ffclibTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(DIRECTORY src/ffclib DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} FILES_MATCHING PATTERN "*.hpp")
    install(EXPORT ffclibTargets NAMESPACE ffclib:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ffclib)
    configure_package_config_file(cmake/ffclibConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/ffclibConfig.cmake" INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ffclib)
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/ffclibConfigVersion.cmake" COMPATIBILITY SameMinorVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/ffclibConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/ffclibConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/ffclib)
endif()
