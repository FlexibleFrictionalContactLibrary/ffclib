# Eigen and HDF5: a parent project's targets, an installation, or fetched (Eigen 3.4.0 header-only; HDF5 1.14.6
# with zlib, built as static libraries). Settings for the fetched dependencies are normal variables of this
# directory scope, so they do not leak into the cache of a parent project.
include(FetchContent)

# Build only what ffclib links and skip the dependencies' own install rules (CMake >= 3.28; older versions also
# build and install their other libraries, headers and packages). Used for HDF5 and the conic solvers.
set(_ffclib_exclude "")
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.28)
    set(_ffclib_exclude EXCLUDE_FROM_ALL)
endif()
# option() and set(CACHE) in the dependencies honor the normal variables below (their own cmake_minimum_required
# predates these policies).
set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(CMAKE_POLICY_DEFAULT_CMP0126 NEW)
set(BUILD_TESTING OFF)
# Extracted files get the time of extraction (rebuilds after a new download).
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

# --- Eigen ---
# Eigen's CMake project registers its build tree in the user package registry (~/.cmake/packages/Eigen3), where
# every later find_package(Eigen3) on the machine finds it, even when the configure step failed before
# Eigen3Targets.cmake was generated there. A fetched Eigen does not register itself (CMP0090), and the search
# ignores the registry, so stale entries of other build trees do no harm.
set(CMAKE_POLICY_DEFAULT_CMP0090 NEW)
set(CMAKE_EXPORT_PACKAGE_REGISTRY OFF)
if(NOT TARGET Eigen3::Eigen)
    find_package(Eigen3 3.4 QUIET NO_MODULE NO_CMAKE_PACKAGE_REGISTRY NO_CMAKE_SYSTEM_PACKAGE_REGISTRY)
    if(NOT Eigen3_FOUND)
        message(STATUS "ffclib: Eigen3 not found, fetching Eigen 3.4.0")
        set(EIGEN_BUILD_DOC OFF)
        set(EIGEN_BUILD_PKGCONFIG OFF)
        FetchContent_Declare(Eigen3 GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0 GIT_SHALLOW TRUE)
        FetchContent_MakeAvailable(Eigen3)
        # A fetched Eigen's headers are system headers (an installed one is already).
        get_target_property(_ffclib_eigen_include_dirs eigen INTERFACE_INCLUDE_DIRECTORIES)
        set_target_properties(eigen PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_ffclib_eigen_include_dirs}")
    endif()
endif()

# --- HDF5 ---
# The C library; the tests also need the high-level library (fclib in test_io).
set(_ffclib_hdf5_components C)
if(FFCLIB_BUILD_TESTS)
    list(APPEND _ffclib_hdf5_components HL)
endif()
set(FFCLIB_FETCHED_HDF5 OFF)
if(NOT TARGET hdf5::hdf5)
    find_package(HDF5 QUIET COMPONENTS ${_ffclib_hdf5_components})
    if(NOT HDF5_FOUND)
        set(FFCLIB_FETCHED_HDF5 ON)
    endif()
endif()
if(FFCLIB_FETCHED_HDF5)
    message(STATUS "ffclib: HDF5 (${_ffclib_hdf5_components}) not found, fetching HDF5 1.14.6 and zlib 1.3.1")
    # Static and position independent: they end up in the shared ffclib.
    set(BUILD_SHARED_LIBS OFF)
    set(BUILD_STATIC_LIBS ON)
    set(_ffclib_pic "${CMAKE_POSITION_INDEPENDENT_CODE}")
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(HDF_PACKAGE_NAMESPACE "hdf5::")  # the namespace of HDF5's own targets (its cache defaults, not read here)
    set(HDF5_BUILD_HL_LIB ${FFCLIB_BUILD_TESTS})
    set(HDF5_BUILD_TOOLS OFF)
    set(HDF5_BUILD_UTILS OFF)
    set(HDF5_BUILD_EXAMPLES OFF)
    set(HDF5_BUILD_CPP_LIB OFF)
    set(HDF5_BUILD_FORTRAN OFF)
    set(HDF5_BUILD_JAVA OFF)
    set(HDF5_ENABLE_SZIP_SUPPORT OFF)
    # ffclib's arrays are deflate compressed. HDF5 fetches zlib and links it statically (no system zlib needed).
    set(HDF5_ENABLE_Z_LIB_SUPPORT ON)
    set(HDF5_ALLOW_EXTERNAL_SUPPORT TGZ)
    set(ZLIB_USE_LOCALCONTENT OFF)
    set(ZLIB_TGZ_ORIGPATH https://github.com/madler/zlib/releases/download/v1.3.1)
    set(ZLIB_TGZ_NAME zlib-1.3.1.tar.gz)
    # zlib as part of HDF5, like HDF5 itself as a subproject: no global output directories in the cache (the lower
    # case prefix, HDF_DIR_PATHS(zlib)) and no install, export or package registry entry (upper case).
    set(zlib_EXTERNALLY_CONFIGURED 1)
    set(ZLIB_EXTERNALLY_CONFIGURED 1)
    FetchContent_Declare(hdf5 URL https://github.com/HDFGroup/hdf5/releases/download/hdf5_1.14.6/hdf5-1.14.6.tar.gz
        URL_HASH SHA256=e4defbac30f50d64e1556374aa49e574417c9e72c6b1de7a4ff88c4b1bea6e9b ${_ffclib_exclude})
    FetchContent_MakeAvailable(hdf5)
    set(CMAKE_POSITION_INDEPENDENT_CODE "${_ffclib_pic}")
    unset(BUILD_SHARED_LIBS)
    unset(BUILD_STATIC_LIBS)
    # The names of FindHDF5's imported targets.
    add_library(hdf5::hdf5 ALIAS hdf5-static)
    if(FFCLIB_BUILD_TESTS)
        add_library(hdf5::hdf5_hl ALIAS hdf5_hl-static)
    endif()
endif()
