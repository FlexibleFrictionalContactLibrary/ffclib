# cfc: a lightweight library for fclib's global frictional contact problems with compliant bilateral
# rows (FORMAT.md): reader and writer of the HDF5 files, the common residuals and the cone programs of the
# convexified problem. Dependencies: Eigen (>= 3.4) and the HDF5 C library.
#
#   cmake -S . -B build && cmake --build build && ctest --test-dir build
#
# As a subproject: add_subdirectory(cfc) and link cfc::cfc; installed: find_package(cfc) and cfc::cfc.
cmake_minimum_required(VERSION 3.20)
project(cfc VERSION 0.2.0 LANGUAGES C CXX)

option(CFC_BUILD_TOOLS "Build cfc_info" ON)
option(CFC_BUILD_TESTS "Build the tests (fetches fclib for the interoperability test)" ${PROJECT_IS_TOP_LEVEL})
option(CFC_INSTALL "Install the library and its CMake package" ${PROJECT_IS_TOP_LEVEL})
set(CFC_FCLIB_FILES "" CACHE STRING "Files of the FCLIB collection that test_cfc reads in addition (optional)")

if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)
include(GNUInstallDirs)

if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.4 QUIET NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(STATUS "cfc: Eigen3 not found, fetching Eigen 3.4.0")
        set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_TESTING OFF CACHE BOOL "" FORCE)
        set(EIGEN_BUILD_PKGCONFIG OFF CACHE BOOL "" FORCE)
        FetchContent_Declare(Eigen3 GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0 GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(Eigen3)
    endif()
endif()
find_package(HDF5 REQUIRED COMPONENTS C)

add_library(cfc src/problem.cpp src/file.cpp src/residuals.cpp src/cone_program.cpp)
add_library(cfc::cfc ALIAS cfc)
target_compile_features(cfc PUBLIC cxx_std_17)
set_target_properties(cfc PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(cfc PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(cfc PUBLIC Eigen3::Eigen PRIVATE hdf5::hdf5)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(cfc PRIVATE -Wall -Wextra)
endif()

if(CFC_BUILD_TOOLS)
    add_executable(cfc_info tools/cfc_info.cpp)
    target_link_libraries(cfc_info PRIVATE cfc::cfc)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(cfc_info PRIVATE -Wall -Wextra)
    endif()
endif()

if(CFC_BUILD_TESTS)
    enable_testing()
    # fclib itself (C, needs HDF5's high-level API) checks that both read each other's files. Only its
    # source is compiled: its CMake project builds tests and needs SuiteSparse.
    find_package(HDF5 REQUIRED COMPONENTS C HL)
    FetchContent_Declare(fclib GIT_REPOSITORY https://github.com/FrictionalContactLibrary/fclib.git
        GIT_TAG f894f34f02093910ceb49943d69f4901142bb880 SOURCE_SUBDIR do-not-configure)
    FetchContent_MakeAvailable(fclib)
    add_library(cfc_test_fclib STATIC "${fclib_SOURCE_DIR}/src/fclib.c")
    set_target_properties(cfc_test_fclib PROPERTIES C_STANDARD 99)
    target_include_directories(cfc_test_fclib SYSTEM PUBLIC "${fclib_SOURCE_DIR}/src")
    target_link_libraries(cfc_test_fclib PUBLIC hdf5::hdf5 hdf5::hdf5_hl)

    add_executable(test_cfc tests/test_cfc.cpp)
    target_link_libraries(test_cfc PRIVATE cfc::cfc cfc_test_fclib)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(test_cfc PRIVATE -Wall -Wextra)
    endif()
    add_test(NAME test_cfc COMMAND test_cfc ${CFC_FCLIB_FILES} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(CFC_INSTALL)
    include(CMakePackageConfigHelpers)
    install(TARGETS cfc EXPORT cfcTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
    install(DIRECTORY include/cfc DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    if(CFC_BUILD_TOOLS)
        install(TARGETS cfc_info RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()
    install(EXPORT cfcTargets NAMESPACE cfc:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cfc)
    configure_package_config_file(cmake/cfcConfig.cmake.in "${CMAKE_CURRENT_BINARY_DIR}/cfcConfig.cmake" INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cfc)
    write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/cfcConfigVersion.cmake" COMPATIBILITY SameMinorVersion)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/cfcConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/cfcConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cfc)
endif()
